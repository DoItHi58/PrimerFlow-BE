import sqlite3
import numpy as np
import pysam
from typing import List, Dict, Tuple, Literal

############################################
# Basic utilities
############################################

RC_MAP = str.maketrans("ATCG", "TAGC")

def reverse_complement(seq: str) -> str:
    return seq.translate(RC_MAP)[::-1]

############################################
# Thermodynamics & NW Alignment
############################################

NN_PARAMS = {
    'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
    'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
    'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
    'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
    'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
    'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
    'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4),
    'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9)
}

def calc_tm_nn(seq: str, dna_nM=50.0, salt_mM=50.0) -> float:
    dh, ds = 0.0, 0.0
    for i in range(len(seq) - 1):
        h, s = NN_PARAMS.get(seq[i:i+2], (0, 0))
        dh += h
        ds += s

    # ds는 cal/(K·mol) 단위, dh는 kcal/mol 단위이므로 dh에 1000을 곱함
    # 염 농도 보정 (SantaLucia 1998)
    ds_corrected = ds + 0.368 * (len(seq) - 1) * np.log(salt_mM / 1000.0)
    R = 1.987 # Gas constant
    c = dna_nM * 1e-9
    
    # 분모 0 체크 및 계산
    denominator = ds_corrected + R * np.log(c / 4)
    if denominator == 0: return 0.0
    
    tm = (dh * 1000) / denominator
    return tm - 273.15

def needleman_wunsch_mismatch(seq1: str, seq2: str) -> int:
    """실제 NW 알고리즘을 이용한 3' 말단 부위 미스매치 계산 (단순화된 버전)"""
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
    # 전체 길이에서 최적 점수를 뺀 값을 미스매치 비용으로 환산
    return max(len(seq1), len(seq2)) - int(score_matrix[n][m])

############################################
# Main Designer
############################################

class PrimerDesigner:
    def __init__(self, genome_fasta: str, annotation_db: str):
        self.genome = pysam.FastaFile(genome_fasta)
        self.db = sqlite3.connect(annotation_db)
        self.cur = self.db.cursor()

    def generate_candidates(self, template: str, k_min=18, k_max=25, tm_range=(57, 63)) -> List[Dict]:
        candidates = []
        for k in range(k_min, k_max + 1):
            for i in range(len(template) - k + 1):
                seq = template[i:i+k]
                for strand, s in [('+', seq), ('-', reverse_complement(seq))]:
                    tm = calc_tm_nn(s)
                    if not (tm_range[0] <= tm <= tm_range[1]): continue

                    # Poly-X 필터 (예: AAAA, GGGG 등 4개 이상 반복 제거)
                    if any(base * 4 in s for base in "ATCG"): continue
                    
                    # GC Content 필터 (40-60% 권장)
                    gc_content = (s.count('G') + s.count('C')) / len(s)
                    if not (0.4 <= gc_content <= 0.6): continue

                    # 3' 말단 안정성 및 GC Clamp
                    dg3 = sum(NN_PARAMS.get(s[-5:][j:j+2], (0, 0))[0] for j in range(4))
                    if dg3 <= -10.0 or s[-1] not in "GC": continue
                    if s[-5:].count('G') + s[-5:].count('C') > 4: continue

                    candidates.append({
                        'seq': s, 'start': i, 'end': i + k - 1,
                        'strand': strand, 'tm': tm, 'dg3': dg3
                    })
        return candidates

    def local_db_filter(self, chrom: str, primer: Dict, 
                        junction_mode: Literal["none", "flanking", "spanning"] = "none") -> bool:
        # 1. SNP 필터링 (Strand 방향에 따른 3' 말단 정의)
        if primer['strand'] == '+':
            s_pos, e_pos = primer['end'] - 4, primer['end']
        else:
            s_pos, e_pos = primer['start'], primer['start'] + 4
            
        self.cur.execute("SELECT COUNT(*) FROM snp WHERE chrom=? AND pos BETWEEN ? AND ?", (chrom, s_pos, e_pos))
        if self.cur.fetchone()[0] > 0: return False

        # 2. 제한효소 필터링
        self.cur.execute("SELECT COUNT(*) FROM restriction_site WHERE chrom=? AND NOT (end < ? OR start > ?)",
                         (chrom, primer['start'], primer['end']))
        if self.cur.fetchone()[0] > 0: return False

        # 3. Exon Junction 필터링
        
        self.cur.execute("SELECT start, end FROM exon WHERE chrom=? ORDER BY start", (chrom,))
        exons = self.cur.fetchall()
        
        if junction_mode == "flanking":
            # 인트론을 사이에 두고 서로 다른 엑손에 위치해야 함 (Pairing 단계에서 주로 체크)
            pass 
        elif junction_mode == "spanning":
            # 프라이머 자체가 엑손 경계에 걸쳐 있어야 함
            is_on_junction = any(primer['start'] < end and primer['end'] > start 
                                 for start, end in zip([e[1] for e in exons[:-1]], [e[0] for e in exons[1:]]))
            if not is_on_junction: return False
        
        return True

    def specificity_check(self, chrom: str, primer: Dict, target_start: int, target_end: int,
                          max_hits=50, mismatch_cutoff=2) -> bool:
        hits = 0
        for ref in self.genome.references:
            full_seq = self.genome.fetch(ref)
            # 정방향 및 역상보 서열 모두 탐색
            for search_seq in [primer['seq'], reverse_complement(primer['seq'])]:
                pos = full_seq.find(search_seq)
                while pos != -1:
                    # 타겟 구간 내부의 결합은 허용 (의도된 결합)
                    if ref == chrom and target_start <= pos <= target_end:
                        pos = full_seq.find(search_seq, pos + 1)
                        continue
                    
                    off_target = full_seq[pos:pos+len(primer['seq'])]
                    # Needleman-Wunsch를 통한 3' 말단 미스매치 정밀 검사
                    mm = needleman_wunsch_mismatch(primer['seq'][-10:], off_target[-10:])
                    if mm < mismatch_cutoff: return False

                    hits += 1
                    if hits > max_hits: return False
                    pos = full_seq.find(search_seq, pos + 1)
        return True

    def pair_primers(self, primers: List[Dict], product_range=(100, 300), 
                     max_tm_diff=3.0, opt_tm=60.0) -> List[Dict]:
        fwd = [p for p in primers if p['strand'] == '+']
        rev = [p for p in primers if p['strand'] == '-']
        pairs = []

        for f in fwd:
            for r in rev:
                size = r['start'] - f['end']
                if not (product_range[0] <= size <= product_range[1]): continue
                if abs(f['tm'] - r['tm']) > max_tm_diff: continue

                # 최적 Tm과의 차이를 기반으로 페널티 계산
                penalty = (abs(f['tm'] - opt_tm) + abs(r['tm'] - opt_tm) +
                           abs(f['dg3'] + 8.0) + abs(r['dg3'] + 8.0))

                pairs.append({
                    'fwd': f, 'rev': r, 'product_size': size, 'penalty': penalty
                })
        return sorted(pairs, key=lambda x: x['penalty'])