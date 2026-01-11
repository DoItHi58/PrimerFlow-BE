from pydantic import BaseModel
from typing import Optional, Literal, Any

class Metrics(BaseModel):
    tm_c: Optional[float] = None            # 녹는점
    gc_percent: Optional[float] = None      # GC 함량
    penalties: Optional[Any] = None         # 패널티 점수 등

# 1) GenomeSequence
class GenomeSequence(BaseModel):
    id: str                                 # 고유 ID
    name: str                               # FASTA Header
    sequence: str                           # 정규화된 유전자 서열 (A/C/G/T/N...)
    length_bp: int                          # 서열 길이

# 2) PrimerCandidate
class PrimerCandidate(BaseModel):
    id: str                                 # 후보 ID
    sequence: str                           # 프라이머 서열
    start_bp: int                           # 프라이머 시작 위치 (1-based)
    end_bp: int                             # 프라이머 끝 위치
    strand: Literal["forward", "reverse"]   # 방향
    metrics: Metrics                    
