from fastapi import APIRouter
from fastapi import status, HTTPException
from schemas.request import PrimerDesignRequest
from schemas.response import PrimerDesignResponse

router = APIRouter()

@router.post("/design", response_model=PrimerDesignResponse, status_code=status.HTTP_200_OK)
async def design(request: PrimerDesignRequest) -> PrimerDesignResponse:
    """프라이머 설계"""
    raise HTTPException(status_code=501, detail="Not implemented")
    


    
