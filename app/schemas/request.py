from typing import Literal, Optional

from pydantic import BaseModel

from app.schemas.common import EndMismatchStrictness, PrimerTm, Range, SearchRange


class BasicInput(BaseModel):
    templateSequence: str
    targetOrganism: str
    productSize: Range
    primerTm: PrimerTm


class PrimerProperties(BaseModel):
    gcContent: Range
    maxTmDifference: float
    gcClamp: bool
    maxPolyX: int
    concentration: float


class PrimerSpecificity(BaseModel):
    checkEnabled: bool
    spliceVariantHandling: bool
    snpExclusion: bool
    endMismatchStrictness: Optional[EndMismatchStrictness] = None
    misprimingLibrary: bool


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
