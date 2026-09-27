# m0ip3h_stimulus_encoding — IP3→H 결합을 넣으면 자극이 주기를 부호화하는가

**계산 목적**: m0_period_waveform에서 원식(AIMS 반응항)은 자극 S̄가 주기를 바꾸지 못함이 확인됨
(P→H 경로 부재). H 식에 IP3 의존항을 최소로 추가했을 때 (i) 주기가 자극과 함께 **감소**하는지,
(ii) 데이터의 6.7배 범위를 낼 수 있는지, (iii) 파형(듀티)이 개선되는지 확인.
개정 계획 v4의 관문 G0 "모델 계약"(IP3→H 항 추가 여부) 결정 근거.

## 두 후보 (gH가 중립값이면 원식과 정확히 일치)
| mech | 식 | 중립값 | 의미 |
|---|---|---|---|
| 1 | κ_H(P) = κ_H (d_H + g_H P)/(d_H + P) | g_H = 1 | IP3가 Ca 불활성화 임계값을 올림 (Li–Rinzel형) |
| 2 | dH/dt = (1 + g_H P/(d_H+P)) (Q(C) − H)/ε_H | g_H = 0 | IP3가 수용체 회복·불활성화 속도를 가속 |

## 파일
| 파일 | 역할 |
|---|---|
| `m0ip3h_stimulus_sweep.c` | M0 + 결합항(mech, gH, dH) 적분기·스윕·특징 계산. 인터페이스는 m0_period_waveform_sweep와 동일 |
| `run_m0ip3h_sweeps.sh` | 고정 설정 S̄ 스윕 8개 → `results/m0ip3h_period_vs_stimulus_*.csv`, `..._vs_gH_*.csv` |
| `search_m0ip3h_dose_range.py` | mech 2의 (ε_H, g_H, d_H, ρ) 격자 192개 × S̄ 9점 탐색 (약 1.5분) → `results/m0ip3h_dose_range_search.csv`, `..._best_curves.csv` |
| `plot_m0ip3h_stimulus_encoding.py` | 요약 그림 → `figures/m0ip3h_stimulus_encoding_{main,supp}.{png,pdf,eps}` |

## 실행
    make && ./run_m0ip3h_sweeps.sh && python3 search_m0ip3h_dose_range.py && python3 plot_m0ip3h_stimulus_encoding.py

## 2026-09-21 결과
- 임계값 결합(mech 1)은 주기를 **증가**시킴(진폭↑ → 불활성화 깊어짐). 데이터와 반대.
- 속도 결합(mech 2)은 주기를 **감소**시킴. Table 1 ρ에서는 P의 자극 의존 범위가 좁아(바닥 ~2,000) 1.8배에 그침.
  ρ를 줄여 P 범위를 넓히면 3배 이상이나 높은 S̄에서 유효 ε_H가 Hopf 문턱(≈16) 아래로 떨어져 진동 소멸.
  기본 ε_H를 300–600으로 올리면 전 범위 진동 유지.
- 격자 탐색: 전 범위 진동·감소 조합 120/192. 최대 주기 비 **8.3배**(ε_H=600, g_H=100, d_H=10000, ρ=0.024), 데이터 6.7배 달성 가능.
- 용량 매핑 S̄ = cU/(K+U), t = T₀s 적합: T₀ ≈ 0.44 s, K ≈ 240 µM, log-RMS 0.22. 3 µM 쪽(206 s)을 과소 예측 — 데이터의 3→5 µM 급감이 모델 곡선보다 가파름.
- **듀티는 모든 조합에서 0.47–0.56**: 자극 부호화와 무관하게 파형 문제는 남음 → 다음 계산(파형: C 동역학의 시간척도 분리).

## Figures
- `figures/m0ip3h_stimulus_encoding_main.{png,pdf,eps}` (paper): coupling direction; best config mapped to the dose data
- `figures/m0ip3h_stimulus_encoding_supp.{png,pdf,eps}` (supplement): rate-coupling range; grid search; gain sweep; waveform

## Update 2026-09-22 (data source)
ISI medians and their 95 % CIs come from `../data_dose_units/results/dose_table.csv` (train medians; 3, 5, 15 and
30 µM from constant-dose records, 10 µM step-only and shown hollow, not fitted).  The dose map of the best config is
fitted to the four constant-dose levels (K grid up to 3000 µM): T0 = 0.56 s, log-RMS 0.30.
