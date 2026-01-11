from pydantic import BaseModel
from typing import Optional, Literal
from schemas.common import Range, PrimerTm, EndMismatchStrictness, SearchRange

# 1
class BasicInput(BaseModel):
    templateSequence: str
    targetOrganism: str
    productSize: Range
    primerTm: PrimerTm

# 2
class PrimerProperties(BaseModel):
    gcContent: Range
    maxTmDifference: float
    gcClamp: bool
    maxPolyX: int
    concentration: float

# 3
class PrimerSpecificity(BaseModel):
    checkEnabled: bool
    spliceVariantHandling: bool
    snpExclusion: bool
    endMismatchStrictness: Optional[EndMismatchStrictness] = None
    misprimingLibrary: bool

# 4
class PrimerPosition(BaseModel):
    searchRange: SearchRange
    exonJunctionSpan: Literal["none", "flanking", "spanning"]
    intronInclusion: bool
    intronSize: Optional[Range] = None
    restrictionEnzymes: list[str]

class PrimerDesignRequest(BaseModel):
    basic: BasicInput
    properties: PrimerProperties
    specificity: PrimerSpecificity
    position: PrimerPosition