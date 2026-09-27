# m0ip3h_amplitude_observation — 진폭과 형광 관측식

**계산 목적**: 이완진동 참조 영역(m0ip3h_waveform_duty)은 주기·파형은 맞지만 C 피크/휴지 비가 수백~천 배이고
피크가 자극에 따라 4배 감소한다. 데이터의 ΔF/F₀ 진폭은 0.9–1.5로 용량에 거의 무관하다.
형광 관측식 F(C) = 1 + R·C/(K_d + C), ΔF/F₀ = (F(C) − F(C_rest))/F(C_rest)로
(i) 어떤 (K_d, R)이 모델 C를 데이터 진폭 범위로 보내는지, (ii) 그때 포화가 F-공간 파형을 데이터 듀티 밖으로 넓히지 않는지,
(iii) 진폭–용량 경향이 데이터와 맞는지 확인.

## 파일
| 파일 | 역할 |
|---|---|
| `m0ip3h_amplitude_sweep.c` | m0ip3h_waveform_sweep와 동일한 모델·인터페이스 (자체 포함용 사본) |
| `run_m0ip3h_amplitude_sweeps.sh` | 참조 영역 S̄ 스윕(C_peak, C_rest) + α·β·κ_S가 진폭에 주는 영향 → `results/` |
| `plot_m0ip3h_amplitude_observation.py` | 용량 매핑 적합, 관측식 (K_d, R) 격자 → `results/m0ip3h_observation_grid.csv`, 그림 → `figures/` |

## 실행
    make && ./run_m0ip3h_amplitude_sweeps.sh && python3 plot_m0ip3h_amplitude_observation.py

## 2026-09-21 결과
- 데이터: ΔF/F₀ 진폭 중앙값 0.87/1.23/1.50/0.90/1.41 (3/5/10/15/30 µM) — 용량 의존 경향 없음. 세포 내 corr(ISI, 진폭) 중앙값 0.09.
  원자료 15 µM 파일의 피크 ΔF/F₀ 중앙값 0.90 (처리 CSV와 일치).
- 모델(참조 영역, 용량 매핑 c=11.6, K=300 µM, T₀=0.72 s): C_peak 31→12, C_rest ≈ 0.02 (비 1500→580).
- 관측식 격자(15 µM): 진폭 IQR [0.42, 1.25]와 F-듀티 [0.05, 0.15]를 동시에 만족하는 (K_d, R)은 6점, 모두 **비포화 영역**
  (K_d ≳ C_peak, R ≈ 1–3). 포화 영역(K_d ≪ C_peak, R ≥ 10)은 진폭이 너무 크거나 F-공간 스파이크가 넓어짐.
  → 데이터의 좁은 F-스파이크는 지시약이 포화되지 않았음을 뜻하고, 따라서 C-공간 듀티를 데이터 듀티와 직접 비교한 앞 계산의 전제가 정당화됨.
- 선택 (K_d=32, R=3): 15 µM 진폭 0.94(데이터 0.90). 그러나 **진폭–용량 경향이 반대**: 모델 1.45→0.79 감소, 데이터는 평탄/약간 증가.
  원인은 모델에서 높은 자극일수록 H 회복이 빨라 더 낮은 H에서 스파이크가 점화되어 방출량이 줄기 때문. 실제 세포에서는 진폭이 자극에 둔감(전부-아니면-무 스파이크).
  → 다음 계산의 표적: 진폭 불변성(예: 방출 문턱 구조, κ_A·n_A, 또는 ER 저장량 변수)과 세포 내 ISI–진폭 무상관.

## Figures
- `figures/m0ip3h_amplitude_observation_main.{png,pdf,eps}` (paper): amplitude vs dose data vs model; ΔF/F0(t) model vs data
- `figures/m0ip3h_amplitude_observation_supp.{png,pdf,eps}` (supplement): data amplitude; C_peak/C_rest; amplitude and F-duty vs K_d
