from pydantic import BaseModel
from typing import Optional
from schemas.primer import GenomeSequence, PrimerCandidate
from schemas.request import PrimerDesignRequest

class Meta(BaseModel):
    params: PrimerDesignRequest             # 요청 시 사용된 파라미터 (검증용)
    timestamp: str                          # 생성 시간
    execution_time_ms: Optional[int] = None # 실행 시간

class PrimerDesignResponse(BaseModel):
    genome: GenomeSequence                  # 분석된 게놈 정보
    candidates: list[PrimerCandidate]       # 후보 프라이머 목록
    meta: Meta                              # 메타 데이터