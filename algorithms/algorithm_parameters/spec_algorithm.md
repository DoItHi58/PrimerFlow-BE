# Primer design 시 고려해야 하는 것들
##### 본 문서에서는 생물학적 단어와 내용은 최대한 배제하여 개발 시 숙지하여야 한다고 판단되는 단어만 첨언하여 작성하였습니다.

# Introduction
Primer 설계는 주형 DNA(template DNA)라는 거대한 문자열 텍스트 안에서, 우리가 원하는 특정 구간만을 복사하기 위해 "정확한 위치"를 찾아내는 "검색 패턴(Search Pattern)"을 만드는 과정과 유사합니다. Primer가 주형 DNA에 얼마나 "정확하고(Specificity)", "빠르게(Efficiency)" 붙느냐가 핵심입니다. 아래는 Primer 설계 시 고려해야 하는 변수들을 컴퓨터 공학적 관점(알고리즘 및 Pseudocode)에서 정리한 내용입니다.

## 1. UI 입력 파라미터 (Input Parameters)

사용자 인터페이스(UI)에서 입력받아야 할 항목들은 알고리즘 내 역할에 따라 다음과 같이 분류된다.

### 1.1 기본 입력 (Basic Input)
알고리즘 구동을 위한 필수 데이터이다.

| 항목명 (Label) | 필수 | 입력 타입 | 설명 및 로직 |
| :--- | :---: | :--- | :--- |
| **PCR Template Sequence** | ✅ | String / File | 프라이머를 설계할 표적 DNA 서열 (FASTA 포맷). |
| **Target Organism (DB)** | ✅ | Dropdown | 특이성 검증(Specificity Check)을 수행할 대상 유전체 DB 선택 (예: hg38, mm10). |
| **PCR Product Size** | ✅ | Range (Min-Max) | 원하는 증폭 산물의 크기 범위 (예: $100 \sim 300 \text{bp}$). |
| **Primer Tm** | ✅ | 3 Inputs | 프라이머의 목표 융해 온도 범위 (Min, Opt, Max).<br>*(일반적으로 $57 \sim 63^\circ\text{C}$)* |

### 1.2 프라이머 물성 설정 (Primer Property)
프라이머 자체의 물리화학적 성질 및 열역학적 안정성을 제어한다.

| 항목명 (Label) | 필수 | 입력 타입 | 알고리즘 처리 로직 |
| :--- | :---: | :--- | :--- |
| **Primer GC Content (%)** | | Range (Min-Max) | 프라이머 서열 내 G/C 비율 검사 (보통 40~60%). |
| **Max Tm Difference** | | Number | Forward/Reverse 프라이머 간 $T_m$ 차이 허용치 (권장 $\le 3^\circ\text{C}$). |
| **GC Clamp Requirement** | | Checkbox | 3' 말단 1~2bp가 G 또는 C로 끝나는지 확인. |
| **Max Poly-X Run** | | Number | 단일 염기 반복(예: AAAAA) 허용 횟수 제한. |
| **Concentration (nM)** | | Number | $T_m$ 및 $\Delta G$ 계산(Nearest-Neighbor 모델) 시 농도 보정 상수로 사용. |

### 1.3 결합 위치 및 구조 제어 (Binding Position)
유전자 구조(Exon/Intron), 변이(SNP), 제한효소 정보를 기반으로 프라이머 위치를 필터링한다.

| 항목명 (Label) | 필수 | 입력 타입 | 알고리즘 처리 로직 |
| :--- | :---: | :--- | :--- |
| **Search Range** | | 2 Inputs (From-To) | 템플릿 내에서 프라이머 탐색 구간 제한. |
| **Exon Junction Span** | | Dropdown | Annotation DB(SQLite) 조회. 프라이머가 Exon 경계에 걸치는지 확인 (mRNA 타겟 시 gDNA 증폭 방지). |
| **Intron Inclusion** | | Checkbox | Annotation DB 조회. Fwd/Rev 프라이머 사이에 인트론이 포함되는지 확인. |
| **Intron Size Range** | | Range (Min-Max) | 포함될 인트론의 크기 범위 지정. |
| **Restriction Enzymes** | | List / Search | **[사용자 정의]** 템플릿 내 해당 효소 절단 부위와 프라이머가 겹치면 제거(Overlap check). |

### 1.4 특이성 검증 (Specificity)
Genome DB(pysam)를 활용하여 비특이적 결합을 차단한다.

| 항목명 (Label) | 필수 | 입력 타입 | 알고리즘 처리 로직 |
| :--- | :---: | :--- | :--- |
| **Specificity Check Enable** | ✅ | Checkbox | Genome DB 전수 검색 수행 여부 결정 (OFF 시 속도 향상, 위험 증가). |
| **Splice Variant Handling** | | Checkbox | Annotation DB 조회. 동일 유전자의 다른 전사체(Variant)에 결합하는 경우 허용(Pass). |
| **SNP Exclusion** | | Checkbox | SNP DB 조회. 프라이머 위치(특히 3' 말단)에 SNP 존재 시 제거. |
| **3' End Mismatch Strictness** | | Dropdown / Number | 비특이적 타겟 발견 시, 3' 말단 $N$ bp 내에 미스매치가 $M$개 미만이면 '위험'으로 간주하고 탈락시킴. |
| **Mispriming Library** | | Checkbox | 반복 서열(Repeat DB) 유사도 검사. |

---

## 2. 알고리즘 설계 가이드 (Algorithm Logic)

입력된 파라미터들은 다음의 4단계 파이프라인을 통해 처리되어야 한다.

### 2.1 1단계: 후보군 생성 및 물성 필터링 (In-Memory Processing)
1.  **Sliding Window:** Template 서열을 지정된 길이(k-mer)로 자르고 Reverse Complement를 생성한다.
2.  **열역학 계산 (Thermodynamics):**
    * **$T_m$ 계산:** 입력된 Concentration과 Salt 농도를 반영하여 Nearest-Neighbor 모델로 계산한다 (SantaLucia, 1998).
    * **$\Delta G$ 평가:** 3' 말단의 자유 에너지($\Delta G$)를 계산하여 너무 안정적인($\le -10 \text{ kcal/mol}$) 경우 비특이적 결합 위험으로 간주하고 제거한다 (Qu et al., 2012).
3.  **GC Clamp/Run:** 3' 말단의 마지막 염기가 G/C인지 확인하고, 3' 말단 5bp 내의 G/C 개수를 제한하여 GC Run을 방지한다.

### 2.2 2단계: 위치 및 구조 기반 필터링 (Local DB Query)
1.  **Restriction Enzyme Check:**
    * 사용자가 선택한 제한효소의 인식 서열(Site)을 Template에서 찾아 좌표화한다.
    * 프라이머 후보의 좌표 범위($Start \sim End$)가 제한효소 좌표와 겹치면(Intersection) 제거한다.
2.  **SNP Filtering:**
    * 프라이머 좌표를 SNP DB에 쿼리하여 겹치는 변이가 있는지 확인한다.
    * 특히 **3' 말단 5bp 이내의 SNP**는 치명적이므로 제거한다 (Ye et al., 2012).
3.  **Gene Structure (Exon/Intron):**
    * SQLite Annotation DB에서 해당 유전자의 Exon 경계 좌표를 로드한다.
    * **Exon Junction Span** 옵션 활성화 시, 프라이머 서열이 Exon $i$의 끝과 Exon $i+1$의 시작을 동시에 포함하는지 검사한다.

### 2.3 3단계: 특이성 검증 (Genome-wide Specificity Check)
> **Note:** `Specificity Check Enable`이 활성화된 경우에만 수행하며, 가장 연산 비용이 높다.

1.  **Genome Search (pysam):**
    * 프라이머 서열을 Query로 하여 Genome DB를 검색한다.
    * 짧은 서열의 민감도를 높이기 위해 **Word Size는 작게(예: 7)**, **E-value는 높게(예: 30,000)** 설정한다 (Ye et al., 2012).
2.  **Global Alignment & 3' Check:**
    * 단순 BLAST 히트(Hit)만으로는 부족하므로, 검색된 유사 서열(Potential Off-target)과 프라이머를 **Needleman-Wunsch 알고리즘**으로 전역 정렬한다.
    * **판정 로직:** 비특이적 타겟의 3' 말단 5bp 이내에 미스매치가 설정된 Strictness 개수(예: 2개) 미만으로 존재하면, 해당 타겟에서 증폭이 일어날 수 있다고 판단하여 프라이머를 탈락시킨다 (Ye et al., 2012).
3.  **Variant 예외 처리:**
    * 발견된 타겟이 의도한 유전자의 다른 Splice Variant라면 비특이적 결합으로 카운트하지 않는다.

### 2.4 4단계: 페어링 및 결과 출력 (Pairing & Output)
1.  **Pairing:** 통과된 Forward, Reverse 후보들을 조합하여 `PCR Product Size`와 `Max Tm Difference` 조건을 만족하는 쌍을 구성한다.
2.  **Scoring:** 프라이머 쌍의 페널티 점수(목표 $T_m$과의 차이, 3' 안정성 등)를 매겨 최적의 순서대로 정렬하여 반환한다.

---

## References
* Altschul, S. F., et al. (1990). Basic local alignment search tool. *Journal of Molecular Biology*.
* Altschul, S. F., et al. (1997). Gapped BLAST and PSI-BLAST. *Nucleic Acids Research*.
* Qu, W., et al. (2012). MFEprimer-2.0: a fast thermodynamics-based program for checking PCR primer specificity. *Nucleic Acids Research*.
* SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS*.
* Ye, J., et al. (2012). Primer-BLAST: A tool to design target-specific primers for polymerase chain reaction. *BMC Bioinformatics*.