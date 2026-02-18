import time
from datetime import datetime, timezone

from fastapi import APIRouter, status

from app.schemas.request import PrimerDesignRequest
from app.schemas.response import PrimerDesignResponse

router = APIRouter()

@router.post("/design", response_model=PrimerDesignResponse, status_code=status.HTTP_200_OK)
async def design(request: PrimerDesignRequest) -> PrimerDesignResponse:
    """프라이머 설계"""
    started = time.perf_counter()

    response = {
        "genome": {
            "id": "dummy_genome_001",
            "name": "dummy_genome",
            "sequence": "ACGT" * 30,
            "length_bp": 120,
        },
        "candidates": [
            {
                "id": "cand_1",
                "sequence": "ACGTACGTACGTACGTACGT",
                "start_bp": 1,
                "end_bp": 20,
                "strand": "forward",
                "metrics": {"tm_c": 60.0, "gc_percent": 50.0, "penalties": None},
            },
            {
                "id": "cand_2",
                "sequence": "TGCATGCATGCATGCATGCA",
                "start_bp": 101,
                "end_bp": 120,
                "strand": "reverse",
                "metrics": {"tm_c": 60.0, "gc_percent": 50.0, "penalties": None},
            },
        ],
        "meta": {
            "params": request,
            "timestamp": "",
            "execution_time_ms": 0,
        },
    }

    response["meta"]["timestamp"] = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    response["meta"]["execution_time_ms"] = int((time.perf_counter() - started) * 1000)

    return PrimerDesignResponse(**response)
