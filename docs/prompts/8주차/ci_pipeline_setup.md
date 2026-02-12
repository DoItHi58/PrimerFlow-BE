# CI 파이프라인 구축 (GitHub Actions)

## 1. 배경 및 목적

- PR/푸시 시 코드 품질을 자동으로 검증(테스트/린트/타입체크/문법 체크/스모크 테스트)하여, “CI 통과 = 머지 가능”한 개발 흐름을 만들기 위함
- 배포(CD)는 아직 대상/환경이 확정되지 않아 범위에서 제외하고, 먼저 CI만 구축
- 로컬 개발 환경과 CI 환경에서 동일한 기준으로 검증되도록(커맨드/엔트리포인트/프로젝트 구조) 표준화를 진행

## 2. 프롬프트 (User Input)

```text
- 배포할 건 아직 안정했어 (CD는 빼고 CI만)
- ruff 추가하고 설명해줘 / F, I 규칙으로 하자
- health 체크 테스트 만들어봐
- 파일들에 따른 .gitignore 수정해줘
- main.py 관련해서 패키징 해볼래? (app/ 구조로 정리 도움)
```

## 3. AI 응답 요약 (AI Output)

- GitHub Actions 기반의 최소 CI 구성을 제안하고, “테스트/린트/타입체크”를 머지 게이트로 사용하는 흐름을 설명
- Ruff를 도입하되 초기 범위는 `F`(기본 오류)와 `I`(import 정렬)로 제한하여 잡음을 줄이고 점진적으로 강화하는 방식 제안
- Pyright 타입체크를 추가하되, 3rd-party stub 미제공으로 생기는 경고는 초기엔 완화한 뒤 점진적으로 강화하는 전략 제안
- pytest는 `python -m pytest`로 실행해 환경별 엔트리포인트 차이로 인한 import 문제를 줄이는 방향을 제안
- FastAPI 앱이 최소한 import 가능한지 확인하는 스모크 테스트(예: `from app.main import app`)를 CI에 포함하도록 제안

## 4. 결과 및 적용 (Result)

- CI 워크플로우 구성
	- 파일: `.github/workflows/ci.yml`
	- 트리거: `develop` 브랜치 push, `develop/main` 대상 PR, 수동 실행(`workflow_dispatch`)
	- 주요 단계
		- 의존성 설치 및 검증: `pip install -r requirements.txt` / `pip check`
		- 테스트: `python -m pytest -q`
		- 린트: `ruff check .`
		- 타입체크: `pyright`
		- 문법 체크: `python -m compileall -q .`
		- 스모크 테스트: `python -c "from app.main import app; assert app is not None"`

- 품질 도구 설정 고정
	- 파일: `pyproject.toml`
	- Ruff: `select = ["F", "I"]`, `known-first-party = ["app"]`
	- Pyright: `typeCheckingMode = "basic"` + 초기 잡음 완화(`reportMissingImports/reportMissingTypeStubs = none`)

- 애플리케이션 엔트리포인트/구조 정리
	- 파일: `main.py`
		- 호환(shim) 용도로 `from app.main import app`를 재노출

- 테스트 추가
	- 파일: `tests/test_health.py`
	- `/health`가 `200` 및 `{ "status": "ok" }`를 반환하는지 검증

- `.gitignore` 보강
	- 파일: `.gitignore`
	- ruff/pyright/pytest 캐시, 빌드 산출물, coverage 파일, 로그 파일 등을 제외 대상으로 추가


- 검증(로컬 재현)
	- CI와 동일한 커맨드로 로컬에서 재현 실행하여 통과 확인
		- `ruff check .`
		- `pyright`
		- `python -m pytest -q`
		- `python -m compileall -q .`
		- `python -c "from app.main import app; assert app is not None"`

