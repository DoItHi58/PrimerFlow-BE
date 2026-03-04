# 🧬 PrimerFlow

> **High-Performance PCR Primer Design & Visualization Platform**

## 프로젝트 개요
**PrimerFlow**는 생명과학 연구원들이 PCR 프라이머를 설계할 때 겪는 비효율을 해결하기 위한 웹 솔루션입니다.


## 프로젝트 구조

```text
PrimerFlow-BE/
├─ .github/               # GitHub 워크플로우 설정
├─ .husky/                
├─ app/                   # FastAPI 애플리케이션 코드
│  ├─ main.py
│  ├─ api/
│  │  └─ v1/
│  │     └─ endpoints/
│  │        ├─ design.py
│  │        └─ health.py
│  ├─ algorithms/
│  └─ schemas/
├─ database/              # DB 파일 및 원천 데이터
├─ docs/                  # 명세/가이드 문서
│  └─ prompts/
│  └─ strategy/
├─ scripts/               # DB 구축/점검 스크립트
├─ tests/                 # pytest 테스트 코드
├─ main.py
├─ requirements.txt
├─ README.md
└─ .gitignore
```


## 개발 환경 설정

### 1. 가상환경 생성 및 활성화

- Windows
    ```powershell
    python -m venv .venv
    .\.venv\Scripts\Activate.ps1    # PowerShell
    # 또는
    .\.venv\Scripts\activate.bat    # cmd
    ```
- macOS / Linux
    ```bash
    # (Ubuntu/WSL) venv 모듈이 없으면 먼저 설치
    sudo apt update
    sudo apt install -y python3-venv
    python3 -m venv .venv
    source .venv/bin/activate
    ```

### 2. 의존성 설치

```powershell
pip install -r requirements.txt
```

### 3. 개발 서버 실행

```powershell
uvicorn app.main:app --reload
```

### 4. 테스트 실행

```powershell
python -m pytest -q
```

- 기본 엔드포인트: http://localhost:8000/
- OpenAPI 문서: http://localhost:8000/docs
- ReDoc 문서: http://localhost:8000/redoc


## 배포 정보

- 서비스 URL: https://primerflow-be.onrender.com
- Health Check: https://primerflow-be.onrender.com/health
- OpenAPI 문서: https://primerflow-be.onrender.com/docs


### 5. Commit convention & commitlint

- 이 레포는 commitlint/husky를 사용합니다. 클론 후 한 번만 실행:
  - `npm install`
  - `git config core.hooksPath .husky` (로컬 기기당 1회)


## 기술 스택
### Backend
- **Framework** : FastAPI
- **Language** : Python
- **Validation** : Pydantic
- **API docs**: Swagger (OpenAPI)
- **Server**: Uvicorn
- **Database**: SQLite
- **Bioinformatics**: pysam

### Quality & Testing
- **Linter/Formatter**: Ruff
- **Type Checking**: Pyright
- **Test**: pytest

### Collaboration
- **Commit Rule**: commitlint
- **Git Hooks**: Husky

### AI Tools
- OpenAI Codex, Google Gemini, GitHub Copilot

## 주간 진행 상황

### Week 1 (25.12.22 - 12.28)
- **작업 내역** : [1주차 Commit](https://github.com/Seq-Lab/PrimerFlow-BE/commit/9b9bf9882e8c376c14fa8daf3cecbde0a3b4d911)
    - 백엔드 기본 구조 세팅
    - [협업 가이드 추가](https://github.com/Seq-Lab/PrimerFlow-BE/commit/9c5e5de6a69456014aa58f39036ce55c5d420dcc)


### Week 2 (25.12.29 - 26.1.4)
- **작업 내역** : [2주차 Commit](https://github.com/Seq-Lab/PrimerFlow-BE/commit/4d11e13ac10cfcb9f8f09c41035bbc8fe3148adf)
    - 명세서 작성 및 프롬프트 추가
    - 협업 편의성을 위해 commitlint 추가
- **AI 활용**
    - main 브랜치 PR 차단 워크플로우 추가 : `.github/workflows/check-main-pr.yml` 생성
     - Spec문서 작성 : GPT와 Gemini에 동일 프롬프트를 입력하고 결과를 통합해 정리


### Week 3 (26.1.5 - 1.11)
- **작업 내역** : [3주차 Commit](https://github.com/Seq-Lab/PrimerFlow-BE/commit/7289717ca93a8d654a0fe5c4d7c3a685a06dc616)
    - `schemas/` 폴더 내 Pydantic 모델 정의
    - `/design` 엔드포인트 정의 (구현 미완료)
    - 알고리즘 명세 및 아키텍처 다이어그램 문서 추가
- **AI 활용**
    - Copilot 리뷰 한국어 지침 추가 : `.github/copilot-instructions.md` 생성

### Week 4 (26.1.12 - 1.18)
- **작업 내역** : [4주차 Commit](https://github.com/Seq-Lab/PrimerFlow-BE/commit/0f66dc4ef55129c4381c34d9d081c29c09b1e388)
    - `PrimerDesigner.py` 구현: 핵심 PCR 프라이머 설계 로직 첫 번째 버전 추가
    - 블랙박스 테스트 케이스 설계(EP, BVA, CE) 및 파일 추가
    - 테스트 케이스 시나리오 초안 완성 (지속적 보완 예정)
- **AI 활용**
    - 블랙박스 테스트 케이스 설계 : `docs/prompts/4주차/test_design_prompt.md` 참고

### Week 8 (26.2.9 - 2.15)
- **작업 내역** : [8주차 Commit](https://github.com/Seq-Lab/PrimerFlow-BE/commit/235451f6e47103919f97dce966bac7b985a1955e)
    - 기존 루트(`api/, schemas/, algorithms/`)를 `app/` 하위로 이동하고 import 경로를 `app.*`로 정리
    - GitHub Actions CI 파이프라인 및 Ruff/Pyright 설정(`pyproject.toml`) 추가
    - `/health` 엔드포인트에 대한 기본 테스트 추가 및 실행 가이드 업데이트
- **AI 활용**
    - CI 파이프라인 구축 도움 : `docs/prompts/8주차/ci_pipeline_setup.md` 참고

### Week 9 (26.2.16 - 2.22)
- **작업 내역** : [9주차 Commit](https://github.com/Seq-Lab/PrimerFlow-BE/commit/d1fce29f1c9f054c69f6b34ecf6c6d8299618e2f)
    - 대용량 원천 데이터(GFF3/VCF/rmsk/FASTA)로부터 annotations.db를 생성하는 `scripts/build_db.py` 추가
    - DB 무결성/건수/샘플 미리보기를 위한 `scripts/check_db_detail.py` 및 기존 알고리즘 연동 확인용 `scripts/test_db_integration.py` 추가
    - DB 스키마/데이터 소스/설치 절차 문서(`docs/spec_database.md`) 추가
- **AI 활용**
    - DB 구축 도움 : `docs/prompts/9주차/DB_setup.md` 참고

### Week 10 (26.2.23 - 3.1)
- **작업 내역** : [10주차 Commit](https://github.com/Seq-Lab/PrimerFlow-BE/pull/40)
    - `health/db` 엔드포인트 신설: SQLite DB(`database/annotations.db`) 존재 여부, 읽기 권한 및 테이블 샘플 데이터 확인 로직 추가.
    - API 탐색 개선: 루트(/) 경로 접속 시 문서 및 헬스 체크 엔드포인트 링크를 반환하도록 수정.
