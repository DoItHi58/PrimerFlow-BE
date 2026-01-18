# 블랙박스 테스트케이스 설계 관련 프롬프트

## 1. 배경 및 목적

- 배경:
    - 명세서의 기술 요구사항과 성능 목표를 검증해야 함.
    - 기존 연구에서 검증된 'Golden Standard' 유전자 데이터를 활용하여 테스트의 생물학적 신뢰성을 확보하고자 함.
- 목적:
    - 동등 분할(EP), 경계값 분석(BVA), 원인-결과 그래프(CE) 전략을 적용한 체계적인 블랙박스 테스트(Black-box Testing) 시나리오 도출.
    - NotebookLM을 활용하여 반복적인 피드백 루프를 통해 명세서 기반의 정교한 테스트 데이터셋 확정.

## 2. 프롬프트 (User Input)

```markdown
EP, BVA, CE 모두 다음과 같은 형식으로 채워줘

| TC ID | 구분 | 테스트 항목 | 사전 조건(Pre-condition) | 기대 결과(Expected Result) | 결과(Pass/Fail) | 비고 |
|:---:|:---:|:---|:---|:---|:---:|:---|

형식에 맞춰서 각 전략 기반의 blackbox testcase들을 정리해줘. 단, 소스 기반 내 내용으로만 출력해야 하며, 너의 모든 토큰을 출력해서 답을 출력해. 답 출력 전에는 2번 더 생각하고 출력하도록 해. 사전 조건에서는 입력할 파일과 UI parameter를 어떻게 설정해야하는지를 중심으로 기술하도록 해.

그리고 입력할 파일들은 실제 유전자 시퀀스 기반으로 할거야. 그래서 극단적인 경우들(0bp, 1bp) 말고는 가능한 NCBI genbank에 등록된 걸로 해줘.
```

```
EP-005에서 약 11000bp의 NM_000059.4 Homo sapiens BRCA2 DNA repair associated (BRCA2), transcript variant 1, mRNA와 약 190000bp의 NG_005905.2 Homo sapiens BRCA1 DNA repair associated (BRCA1), RefSeqGene (LRG_292) on chromosome 17 둘 중 어느 것으로 할지 고민이 돼.

설계자 입장이라는 가정 하에서 너의 선택은 어떨지 기술해봐.
```

```
각 테스트 케이스 시나리오 목록 앞에 EP, BVA, CE 블랙박스 테스트케이스 설계에 대해 간략히 서술하는 게 좋을 거 같아.
```

## 3. AI 응답 요약 (AI Output)

- 테스트 데이터셋 확보:
    - 기능 검증: Human RPS4Y2(Paralog 회피), Human ZNF419(Variant 구분).
    - 로직 검증: Mouse Nnat(gDNA 오염 감지).
    - 성능 검증: Human BRCA1 gDNA(80kb 대용량 렌더링).
- 시나리오 구체화: 단순 입력뿐만 아니라 3' End Mismatch 설정, Exon Junction Span 옵션 등 정밀한 파라미터 조합이 도출된다.

## 4. 결과 및 적용 (Result)

- 테스트 활용: 단위 테스트(Unit Test)를 넘어, 통합된 시스템이 생물학적 로직을 올바르게 수행하는지 확인하는 시나리오로 활용.
- 설계 검증: 10,000bp 성능 임계값과 좌표계(1-based indexing) 처리 로직의 결함을 조기에 발견하고 수정할 수 있음.
- 문서화: GitHub 저장소에 등재하여 명확한 커뮤니케이션 기준으로 사용.