**작성자** : 강두하 <br>
**작성일** : 2026-01-11~17 <br>

# Equivalence Partitioning(EP) 토대의 테스트 시나리오 목록
Equivalence Partitioning(EP)은 입력 데이터를 시스템이 동일하게 처리할 것으로 예상되는 그룹으로 나누어 테스트하는 기법입니다. 모든 값을 다 테스트할 수 없으므로, 각 그룹의 대표값을 선정하여 효율적으로 검증합니다.
- 유효 입력 그룹 (Valid): 표준 FASTA 포맷, 정상적인 생물학적 서열(cDNA, gDNA).
-  무효 입력 그룹 (Invalid): 비허용 문자 포함 등 명세서 위반 케이스.
- 생물학적 그룹: 유사 유전자(Paralog)가 많은 유전자군, 스플라이싱 변이(Variant)가 복잡한 유전자군. Variant의 경우에는 뒤에 서술할 Cause-Effect graphing(CE)과 겹쳐서 CE에서만 다룹니다.

| TC ID | 구분 | 테스트 항목 | 입력 | 기대 결과(Expected Result) | 결과(Pass/Fail) | 비고 |
|:---:|:---:|:---|:---|:---|:---:|:---|
| EP-001 | 표준 FASTA 포맷 파싱 및 유효 문자 검증 | DNA 서열 파일 유효성 확인 | 1. 입력 파일: `EP-001_01.fna` 또는 `EP-001_02.fna`<br>2. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``Off`` | 1. 파싱 성공 및 서열 길이 표시.<br>2. '설계' 버튼 활성화.<br>3. 캔버스에 시각화됨. | - | A/C/G/T/N 허용 및 헤더 파싱 로직 |
| EP-002 | 입력 유효성 | 비허용 문자 포함 파일 처리 (Invalid Class) | 1. 입력 파일: `EP-002.fna`<br>2. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``Off`` | 1. 즉시 차단 (백엔드 요청 없음).<br>2. "허용되지 않는 문자가 포함되어 있습니다" 에러 메시지 출력. | - | 유효하지 않은 문자 검증 |
| EP-003 | 입력 유효성 | 파일 업로드 포맷 검증 | 1. 입력 파일:`EP-003.png`<br>2. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``Off`` | 에러 메시지 | - | fasta나 txt 외 다른 확장자 파일을 입력으로 받음. |
| EP-004 | 생물학적 로직 | 유사 유전자(Paralog) 회피 능력 검증 | 1. 입력 파일: `EP-001_1.fna` 또는 `EP-001_2.fna`<br>2. UI 설정:<br>- Target Organism: **``Homo sapiens(RefSeqmRNA)``**<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: **`On`**<br>- **Mispriming library: `Human`** | 1. 결과 리스트에서 RPS4Y1을 Off-target으로 감지.<br>2. 유사성이 높은 구간을 피해서 프라이머 설계 또는 경고 메시지 출력. | - | Paralog 간의 단일 염기 차이 구분 |
| EP-005 | 성능 검증 | 대용량 유전자 데이터 로딩 (Valid) | 1. 입력 파일: EP-005.fna <br>2. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``Off``<br>- Search Range: `1-10000` | 1. 60fps 유지하며 렌더링.<br>2. 스크롤/줌 동작 시 캔버스 끊김 현상 없음. | - | 10,000bp 이상 대용량 시각화 요구사항 |

# Boundary Value Analysis(BVA) 토대의 테스트 시나리오 목록
오류가 발생하기 쉬운 입력 조건의 경계를 집중적으로 테스트하는 기법입니다 최소값, 최대값, 그리고 임계값(Threshold) 직전/직후의 값을 테스트하여 시스템의 안정성을 검증합니다.
- 데이터 길이 경계: 0bp(오류), 1bp(최소), 10,000bp(성능 최적화 기준).
- 생물학적 임계값: 3' 말단 Mismatch 허용 개수(1개 vs 2개) 

| TC ID | 구분 | 테스트 항목 | 사전 조건(Pre-condition) | 기대 결과(Expected Result) | 결과(Pass/Fail) | 비고 |
|:---:|:---:|:---|:---|:---|:---:|:---|
| BVA-001 | 입력 길이 | 최소 무효 길이 (0bp) 검증 | 1. 입력 파일: BVA-001_1.fna 또는 BVA-001_2.fna (공백 또는 헤더만 존재)<br>2. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``Off`` | 1. "서열이 비어 있습니다" 에러 메시지 출력.<br>2. 프로세스 중단. | - | 서열 길이 0 예외 처리 |
| BVA-002 | 입력 길이 | 최소 유효 길이 (1bp) 검증 | 1. 입력 파일: BVA-002.fna <br>2. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``Off`` | 1. 정상 로딩.<br>2. 서열 길이 1bp로 표시. | - | 1bp 이상 유효 |
| BVA-003 | 좌표계 | Search Range 시작점 경계 (0 vs 1) | 1. 입력 파일: BVA-003.fna<br>2. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``Off``<br>- **Search Range: `From 0`** | 1. Validation Error: "좌표는 1부터 시작해야 합니다" 출력.<br>2. 또는 자동으로 1로 보정됨. | - | 1-based Indexing 준수 |
| BVA-004 | 특이성 임계값 | 3' 말단 Mismatch 경계 (1bp vs 2bp) | 1. 입력 파일: BVA-004.fna<br>2. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``On``- **3' end strictness: `말단 5 bp 내에 2 mismatches`** | 1. Mismatch가 1개인 타겟(Variant 1)은 '증폭 위험(Non-specific)'으로 간주하여 필터링되거나 경고 표시.<br>2. Mismatch가 2개 이상이어야 안전함. | - | 3' 말단 2bp 불일치 시 증폭 방지 |
| BVA-005 | 성능 경계 | View Culling 동작 임계값 (10,000bp) | 1. 입력 파일: BVA-005.fna<br>2. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``Off`` | 1. View Culling 알고리즘이 동작하여 화면 밖 데이터 렌더링 제외 확인.<br>2. 성능 저하 없이 부드러운 이동. | - | 10,000bp 기준 성능 최적화 |

# Cause-Effect Graphing(CE) 토대의 테스트 시나리오 목록
입력 조건들의 논리적 조합이 시스템의 출력에 미치는 영향을 검증합니다. 단순한 입력값 테스트를 넘어, 복합적인 비즈니스 로직이나 옵션 간의 상호작용을 테스트하는 데 유용합니다.
- DB 옵션 조합: mRNA DB vs Genome DB 선택에 따른 결과 차이 (gDNA 오염 감지).
- 시각화 로직: 프라이머 구간 겹침(원인) → Y축 레이어 분리(결과).
- 상태 전이: 네트워크 오류(원인) → 재시도 UI 표시(결과).

| TC ID | 구분 | 테스트 항목 | 사전 조건(Pre-condition) | 기대 결과(Expected Result) | 결과(Pass/Fail) | 비고 |
|:---:|:---:|:---|:---|:---|:---:|:---|
| CE-001 | 생물학적 로직 검증 | gDNA 오염 감지 (mRNA + Genome DB) | 1. 입력 파일: CE-001.fna<br>2. UI 설정:<br>- Target Organism: `Mus musculus`<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``On``<br>- Intron Inclusion: `On`<br>- Intron Size: `1000-10000`<br>- Exon Junction span: `flanking` | 결과 리스트에 Target(432bp) 외에 1632bp (gDNA Hit)가 함께 표시됨 | - | 인트론 포함 시 증폭 산물 크기 변화 예측 |
| CE-002 | 생물학적 로직 검증 | Exon Junction Spanning + Variant | 1. 입력 파일: CE-002.fna<br>2. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``On``<br>- Exon Junction Span: `Spanning`<br>| 1. Variant 5에만 특이적인(Exon 2-3 등) 프라이머 생성.<br>2. 다른 변이체(Variant 1, 4 등)는 증폭되지 않음(Not amplified)으로 표시. | - | Spanning 옵션의 특이성 보장 효과. |
| CE-003 | UI/UX | 설계 요청 실패 시 상태 처리 | 1. 입력 파일: 유효한 FASTA 파일<br>2. 상태 조작: 네트워크 연결 끊기 또는 서버 500 에러 Mocking<br>3. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `70-200`(Default)<br>- Primer Tm(Min/Opt/Max): `57/60/63`(Default)<br>- Specificity Check Enable: ``Off`` | 1. 로딩 스피너 종료.<br>2. "서버와 연결할 수 없습니다" 오류 메시지 및 재시도 버튼 표시. | - | 예외 상황 UX 처리 |
| CE-004 | 시각화 | Auto Layout 겹침 방지 로직 | 1. 입력 파일: CE-004.fna<br>2. UI 설정:<br>- Target Organism: ``Homo sapiens``<br>- PCR Product Size: `50-60`(Default)<br>- Primer Tm(Min/Opt/Max): `57/55/60`(Default)<br>- Specificity Check Enable: ``Off`` | 1. 캔버스 상에서 겹치는 프라이머들이 서로 다른 Y축(Layer)에 배치됨.<br>2. 모든 프라이머가 시각적으로 구분 가능함. | - | Auto Layout 알고리즘 검증 |