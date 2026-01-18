# 🧬 PrimerFlow

> **High-Performance PCR Primer Design & Visualization Platform**

## 프로젝트 개요
**PrimeFlow**는 생명과학 연구원들이 PCR 프라이머를 설계할 때 겪는 비효율을 해결하기 위한 웹 솔루션입니다.


## 프로젝트 구조

```text
PrimerFlow-BE/
├─ api/
│  ├─ deps.py
│  └─ v1/
│     └─ endpoints/        # 엔드포인트 모음
├─ schemas/                # Pydantic 모델 모음
├─ algorithms/             # 알고리즘 모음
├─ docs/                   # 협업 가이드, Spec 문서 모음
│  └─ prompts/
│  └─ strategy/
├─ main.py                 # FastAPI 앱 엔트리포인트 
├─ requirements.txt        # Python 패키지 목록
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
    python3 -m venv .venv
    source .venv/bin/activate
    ```

### 2. 의존성 설치

```powershell
pip install -r requirements.txt
```

### 3. 개발 서버 실행

```powershell
uvicorn main:app --reload
```

- 기본 엔드포인트: http://localhost:8000/
- OpenAPI 문서: http://localhost:8000/docs
- ReDoc 문서: http://localhost:8000/redoc


### 4. Commit convention & commitlint

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