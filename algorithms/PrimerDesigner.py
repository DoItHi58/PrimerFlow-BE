import sqlite3
import numpy as np
import pysam
from typing import List, Dict, Tuple, Optional, Literal
############################################
# Thermodynamics & Alignment Utilities
############################################
RC_MAP = str.maketrans("ATCG", "TAGC")
def reverse_complement(seq: str) -> str:
    return seq.translate(RC_MAP)[::-1]
NN_PARAMS = {
    'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2), 'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
    'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7), 'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
    'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0), 'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
    'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4), 'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9)
}
def calc_tm_nn(seq: str, dna_nM=50.0, salt_mM=50.0) -> float:
    dh, ds = 0.0, 0.0
    for i in range(len(seq) - 1):
        h, s = NN_PARAMS.get(seq[i:i+2], (0, 0))
        dh += h
        ds += s
    ds_corrected = ds + 0.368 * (len(seq) - 1) * np.log(salt_mM / 1000.0)
    R = 1.987
    c = dna_nM * 1e-9
    denominator = ds_corrected + R * np.log(c / 4)
    if denominator == 0: return 0.0
    return ((dh * 1000) / denominator) - 273.15
def needleman_wunsch_mismatch(seq1: str, seq2: str) -> int:
    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((n + 1, m + 1))
    for i in range(n + 1): score_matrix[i][0] = -i
    for j in range(m + 1): score_matrix[0][j] = -j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = 1 if seq1[i-1] == seq2[j-1] else -1
            score_matrix[i][j] = max(score_matrix[i-1][j-1] + match,
                                    score_matrix[i-1][j] - 1,
                                    score_matrix[i][j-1] - 1)
    return max(len(seq1), len(seq2)) - int(score_matrix[n][m])
############################################
# Main Designer
############################################
class PrimerDesigner:
    def __init__(self, genome_fasta: str, annotation_db: str):
        self.genome = pysam.FastaFile(genome_fasta)
        self.db = sqlite3.connect(annotation_db)
        self.cur = self.db.cursor()
    def generate_candidates(
        self,
        template: str,
        k_min=18, k_max=25,
        tm_range=(57, 63),
        gc_range=(0.4, 0.6),
        max_poly_x=4,
        gc_clamp=True
    ) -> List[Dict]:
        """Stage 2.1: 물성 기반 후보군 생성"""
        candidates = []
        for k in range(k_min, k_max + 1):
            for i in range(len(template) - k + 1):
                seq = template[i:i+k]
                for strand, s in [('+', seq), ('-', reverse_complement(seq))]:
                    tm = calc_tm_nn(s)
                    if not (tm_range[0] <= tm <= tm_range[1]): continue
                    # 1. Poly-X 필터 반영
                    if any(base * max_poly_x in s for base in "ATCG"): continue
                   
                    # 2. GC Content 필터 반영
                    gc_content = (s.count('G') + s.count('C')) / len(s)
                    if not (gc_range[0] <= gc_content <= gc_range[1]): continue
                    # 3. 3' 말단 안정성 및 GC Clamp
                    dg3 = sum(NN_PARAMS.get(s[-5:][j:j+2], (0, 0))[0] for j in range(4))
                    if dg3 <= -10.0: continue
                    if gc_clamp and s[-1] not in "GC": continue
                    if s[-5:].count('G') + s[-5:].count('C') > 4: continue
                    candidates.append({
                        'seq': s, 'start': i, 'end': i + k - 1,
                        'strand': strand, 'tm': tm, 'dg3': dg3
                    })
        return candidates
    def local_db_filter(
        self,
        chrom: str,
        primer: Dict,
        junction_mode: Literal["none", "flanking", "spanning"] = "none",
        restriction_enzymes: List[str] = [],
        intron_inclusion: bool = True,
        intron_size_range: Optional[Tuple[int, int]] = None
    ) -> bool:
        """Stage 2.2: 위치 및 구조 기반 필터링"""
       
        # 1. SNP 필터링 (3' end strictness)
        s_pos, e_pos = (primer['end'] - 4, primer['end']) if primer['strand'] == '+' else (primer['start'], primer['start'] + 4)
        self.cur.execute("SELECT COUNT(*) FROM snp WHERE chrom=? AND pos BETWEEN ? AND ?", (chrom, s_pos, e_pos))
        if self.cur.fetchone()[0] > 0: return False
        # 2. 제한효소 필터링 (사용자 지정 리스트 반영)
        if restriction_enzymes:
            placeholders = ','.join(['?'] * len(restriction_enzymes))
            query = f"SELECT COUNT(*) FROM restriction_site WHERE chrom=? AND name IN ({placeholders}) AND NOT (end < ? OR start > ?)"
            self.cur.execute(query, (chrom, *restriction_enzymes, primer['start'], primer['end']))
            if self.cur.fetchone()[0] > 0: return False
        # 3. Exon/Intron 구조 필터링
        self.cur.execute("SELECT start, end FROM exon WHERE chrom=? ORDER BY start", (chrom,))
        exons = self.cur.fetchall()
        # Intron Inclusion 로직: 프라이머가 인트론 구간에 걸쳐 있는지 확인
        is_in_intron = any(primer['start'] > exons[i][1] and primer['end'] < exons[i+1][0] for i in range(len(exons)-1))
        if not intron_inclusion and is_in_intron: return False
       
        # Intron Size 제한 확인
        if is_in_intron and intron_size_range:
            for i in range(len(exons)-1):
                if primer['start'] > exons[i][1] and primer['end'] < exons[i+1][0]:
                    i_size = exons[i+1][0] - exons[i][1]
                    if not (intron_size_range[0] <= i_size <= intron_size_range[1]): return False
        # Exon Junction Spanning 로직
        if junction_mode == "spanning":
            is_on_junction = any(primer['start'] < e[1] and primer['end'] > exons[idx+1][0]
                                 for idx, e in enumerate(exons[:-1]))
            if not is_on_junction: return False
           
        return True
    def specificity_check(
        self,
        chrom: str,
        primer: Dict,
        target_start: int,
        target_end: int,
        mispriming_library: bool = False,
        snp_exclusion: bool = False,
        splice_variant_handling: bool = False,
        max_hits=50,
        mismatch_cutoff=2
    ) -> bool:
        """Stage 2.3: 게놈 전체 특이성 및 변이체 처리"""
       
        # 1. Mispriming Library (RepeatMasker 등 반복 서열 DB) 필터링
        if mispriming_library:
            self.cur.execute("SELECT COUNT(*) FROM repeats WHERE chrom=? AND NOT (end < ? OR start > ?)",
                             (chrom, primer['start'], primer['end']))
            if self.cur.fetchone()[0] > 0: return False
        hits = 0
        for ref in self.genome.references:
            full_seq = self.genome.fetch(ref)
            for search_seq in [primer['seq'], reverse_complement(primer['seq'])]:
                pos = full_seq.find(search_seq)
                while pos != -1:
                    # 의도된 타겟 구간 제외
                    if ref == chrom and target_start <= pos <= target_end:
                        pos = full_seq.find(search_seq, pos + 1)
                        continue
                   
                    # 2. Splice Variant Handling
                    if splice_variant_handling:
                        self.cur.execute("SELECT transcript_id FROM exon WHERE chrom=? AND start <= ? AND end >= ?", (ref, pos, pos+len(search_seq)))
                        if self.cur.fetchone(): # 동일 유전자의 다른 전사체라면 off-target에서 제외
                            pos = full_seq.find(search_seq, pos + 1)
                            continue
                    off_target = full_seq[pos:pos+len(primer['seq'])]
                   
                    # 3. SNP Exclusion: Off-target 위치에 SNP가 있으면 결합력이 약해지므로 hits 카운트에서 제외 가능성 판단
                    if snp_exclusion:
                        self.cur.execute("SELECT COUNT(*) FROM snp WHERE chrom=? AND pos BETWEEN ? AND ?", (ref, pos, pos+len(search_seq)))
                        if self.cur.fetchone()[0] > 0:
                            # SNP가 있으면 결합이 불안정해져서 실제 PCR 시 문제가 안 될 수 있음 (전략적 허용)
                            pos = full_seq.find(search_seq, pos + 1)
                            continue
                    # 3' 말단 미스매치 정밀 검사
                    mm = needleman_wunsch_mismatch(primer['seq'][-10:], off_target[-10:])
                    if mm < mismatch_cutoff: return False
                    hits += 1
                    if hits > max_hits: return False
                    pos = full_seq.find(search_seq, pos + 1)
        return True
    def pair_primers(
        self,
        primers: List[Dict],
        product_range=(100, 300),
        max_tm_diff=3.0,
        opt_tm=60.0
    ) -> List[Dict]:
        """Stage 2.4: 최종 페어링 및 Penalty 스코어링"""
        fwd = [p for p in primers if p['strand'] == '+']
        rev = [p for p in primers if p['strand'] == '-']
        pairs = []
        for f in fwd:
            for r in rev:
                size = r['start'] - f['end']
                if not (product_range[0] <= size <= product_range[1]): continue
               
                # Tm 차이 및 페널티 계산 (UI 설정값 opt_tm 반영)
                tm_diff = abs(f['tm'] - r['tm'])
                if tm_diff > max_tm_diff: continue
                penalty = (abs(f['tm'] - opt_tm) + abs(r['tm'] - opt_tm) +
                           abs(f['dg3'] + 8.0) + abs(r['dg3'] + 8.0) + (tm_diff * 2))
                pairs.append({
                    'fwd': f, 'rev': r, 'product_size': size, 'penalty': penalty
                })
        return sorted(pairs, key=lambda x: x['penalty'])