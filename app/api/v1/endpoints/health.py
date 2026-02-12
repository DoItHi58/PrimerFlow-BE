from fastapi import APIRouter

router = APIRouter()


@router.get("/health", summary="Health check")
async def health_check():
    """Simple health check endpoint."""
    return {"status": "ok"}
