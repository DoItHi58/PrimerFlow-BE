import os
import time
from datetime import datetime, timezone

import pysam
from fastapi import APIRouter, HTTPException, status

from app.algorithms.PrimerDesigner import PrimerDesigner
from app.schemas.request import PrimerDesignRequest
from app.schemas.response import PrimerDesignResponse

router = APIRouter()


def _project_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))


def _resolve_paths() -> tuple[str, str]:
    root = _project_root()
    db_path = os.getenv("DB_PATH") or os.path.join(root, "database", "annotations.db")
    genome_path = os.getenv("GENOME_PATH") or os.path.join(
        root, "database", "raw_data", "GRCh38.primary_assembly.genome.fa.gz"
    )
    return db_path, genome_path


def _validate_genome_fasta(genome_path: str) -> None:
    if not os.path.exists(genome_path):
        raise HTTPException(status_code=503, detail="genome FASTA 파일을 찾을 수 없습니다.")

    if genome_path.endswith(".gz"):
        fai_path = f"{genome_path}.fai"
        gzi_path = f"{genome_path}.gzi"
        missing = [path for path in (fai_path, gzi_path) if not os.path.exists(path)]
        if missing:
            raise HTTPException(
                status_code=503,
                detail=(
                    "GENOME_PATH가 .gz인 경우 bgzip FASTA와 인덱스(.fai, .gzi)가 모두 필요합니다. "
                    f"누락 파일: {', '.join(missing)}"
                ),
            )

    try:
        fasta = pysam.FastaFile(genome_path)
        fasta.close()
    except OSError as exc:
        msg = str(exc)
        if "Cannot index files compressed with gzip" in msg or "File truncated" in msg:
            raise HTTPException(
                status_code=503,
                detail=(
                    "현재 genome 파일은 일반 gzip 형식으로 보입니다. "
                    "pysam 사용을 위해서는 bgzip으로 압축된 FASTA(.fa.gz)와 "
                    "인덱스(.fai, .gzi)가 필요합니다."
                ),
            ) from exc
        raise HTTPException(status_code=503, detail=f"genome FASTA 파일을 열 수 없습니다: {exc}") from exc


def _normalize_gc_range(gc_min: float, gc_max: float) -> tuple[float, float]:
    if gc_min > 1 or gc_max > 1:
        return gc_min / 100.0, gc_max / 100.0
    return gc_min, gc_max


def _validate_db_path(db_path: str) -> None:
    if not os.path.exists(db_path):
        raise HTTPException(status_code=503, detail="annotations.db 파일을 찾을 수 없습니다.")


def _intron_size_range(request: PrimerDesignRequest) -> tuple[int, int] | None:
    if not request.position.intronSize:
        return None
    return int(request.position.intronSize.min), int(request.position.intronSize.max)


def _to_candidate(candidate: dict, index: int) -> dict:
    seq = candidate["seq"]
    gc_percent = ((seq.count("G") + seq.count("C")) / len(seq)) * 100 if seq else 0.0
    return {
        "id": f"cand_{index}",
        "sequence": seq,
        "start_bp": int(candidate["start"]),
        "end_bp": int(candidate["end"]),
        "strand": "forward" if candidate["strand"] == "+" else "reverse",
        "metrics": {
            "tm_c": float(candidate.get("tm", 0.0)),
            "gc_percent": float(gc_percent),
            "penalties": candidate.get("penalty"),
        },
    }


def _filter_candidates_by_template(
    designer: PrimerDesigner,
    request: PrimerDesignRequest,
    candidates: list[dict],
    template_info: dict,
) -> list[dict]:
    mapped_candidates = [
        designer.map_to_genomic_coords(candidate, template_info) for candidate in candidates
    ]

    if request.position.searchRange.from_ <= request.position.searchRange.to:
        mapped_candidates = [
            candidate
            for candidate in mapped_candidates
            if request.position.searchRange.from_
            <= candidate["genomic_start"]
            <= request.position.searchRange.to
        ]

    filtered_candidates = [
        candidate
        for candidate in mapped_candidates
        if designer.local_db_filter(
            chrom=template_info["chrom"],
            primer=candidate,
            junction_mode=request.position.exonJunctionSpan,
            restriction_enzymes=request.position.restrictionEnzymes,
            intron_inclusion=request.position.intronInclusion,
            intron_size_range=_intron_size_range(request),
        )
    ]

    if request.specificity.checkEnabled and filtered_candidates:
        end_strict = request.specificity.endMismatchStrictness
        mismatch_cutoff = end_strict.minMismatch if end_strict else 2
        filtered_candidates = designer.filter_specific_primers(
            primers=filtered_candidates,
            target_chrom=template_info["chrom"],
            target_start=template_info["genomic_start"],
            target_end=template_info["genomic_start"] + template_info["template_length"] - 1,
            mispriming_library=request.specificity.misprimingLibrary,
            snp_exclusion=request.specificity.snpExclusion,
            splice_variant_handling=request.specificity.spliceVariantHandling,
            mismatch_cutoff=mismatch_cutoff,
        )

    return filtered_candidates


def _build_response(
    request: PrimerDesignRequest,
    template_info: dict | None,
    candidates: list[dict],
    started: float,
) -> PrimerDesignResponse:
    response = {
        "genome": {
            "id": (template_info["chrom"] if template_info else request.basic.targetOrganism),
            "name": request.basic.targetOrganism,
            "sequence": request.basic.templateSequence,
            "length_bp": len(request.basic.templateSequence),
        },
        "candidates": [_to_candidate(candidate, i) for i, candidate in enumerate(candidates, start=1)],
        "meta": {
            "params": request,
            "timestamp": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            "execution_time_ms": int((time.perf_counter() - started) * 1000),
        },
    }
    return PrimerDesignResponse(**response)


def _close_designer(designer: PrimerDesigner | None) -> None:
    if designer is None:
        return
    try:
        designer.cur.close()
    except Exception:
        pass
    try:
        designer.db.close()
    except Exception:
        pass
    try:
        designer.genome.close()
    except Exception:
        pass


@router.post("/design", response_model=PrimerDesignResponse, status_code=status.HTTP_200_OK)
async def design(request: PrimerDesignRequest) -> PrimerDesignResponse:
    """프라이머 설계"""
    started = time.perf_counter()
    db_path, genome_path = _resolve_paths()

    _validate_db_path(db_path)
    _validate_genome_fasta(genome_path)

    designer: PrimerDesigner | None = None

    try:
        designer = PrimerDesigner(genome_fasta=genome_path, annotation_db=db_path)

        tm_range = (request.basic.primerTm.min, request.basic.primerTm.max)
        gc_range = _normalize_gc_range(request.properties.gcContent.min, request.properties.gcContent.max)

        candidates = designer.generate_candidates(
            template=request.basic.templateSequence,
            tm_range=tm_range,
            gc_range=gc_range,
            max_poly_x=request.properties.maxPolyX,
            gc_clamp=request.properties.gcClamp,
        )

        template_info = designer.locate_template_in_genome(request.basic.templateSequence)

        if template_info:
            candidates = _filter_candidates_by_template(designer, request, candidates, template_info)

        return _build_response(request, template_info, candidates, started)
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"프라이머 설계 중 오류가 발생했습니다: {exc}") from exc
    finally:
        _close_designer(designer)
