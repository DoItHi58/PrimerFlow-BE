from fastapi import FastAPI
from api.v1.endpoints.health import router as health_router
from api.v1.endpoints.design import router as design_router

app = FastAPI(title="PrimerFlow API", version="0.1.0")

app.include_router(health_router)
app.include_router(design_router)

