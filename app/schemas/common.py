from pydantic import BaseModel, Field


class Range(BaseModel):
    min: float
    max: float


class PrimerTm(BaseModel):
    min: float
    opt: float
    max: float


class EndMismatchStrictness(BaseModel):
    regionSize: int
    minMismatch: int


class SearchRange(BaseModel):
    from_: int = Field(..., alias="from")
    to: int
