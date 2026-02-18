<<<<<<< HEAD
"""호환용 엔트리포인트.
=======
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from api.v1.endpoints.health import router as health_router
from api.v1.endpoints.design import router as design_router
>>>>>>> origin/main

`app.main:app`을 재노출합니다.
"""

<<<<<<< HEAD
from app.main import app

__all__ = ["app"]
=======
# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://127.0.0.1:3000",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(health_router)
app.include_router(design_router)
>>>>>>> origin/main

