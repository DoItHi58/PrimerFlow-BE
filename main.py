"""호환용 엔트리포인트.

`app.main:app`을 재노출합니다.
"""

from app.main import app

__all__ = ["app"]

