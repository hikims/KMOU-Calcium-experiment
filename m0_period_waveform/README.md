# m0_period_waveform — M0(결정론 AIMS 반응계)의 주기·파형 스윕

**계산 목적**: AIMS(2026) 식 2.6–2.8의 반응항(확산 제거, Table 1 기본값)이
(i) 자극 의존 주기(데이터: 3→30 µM에서 ISI 206→31 s, 6.6배 감소)와
(ii) 스파이크형 파형(데이터: 반치폭/주기 ≈ 0.16)을 낼 수 있는가.
개정 계획 v4의 관문 G0·G1("기본 재현", 모델 계약 결정 근거)에 해당.

## 파일 (이름 = 계산 목적)
| 파일 | 역할 |
|---|---|
| `m0_period_waveform_sweep.c` | M0 적분기(RK4) + 한 파라미터 스윕 + 주기·반치폭·감쇠 특징 계산 |
| `run_m0_sweeps.sh` | 스윕 7개 실행 → `results/` |
| `plot_m0_period_waveform.py` | 요약 그림 → `figures/` (png·pdf·eps). scipy 필요. 원자료 TXT·처리 CSV를 찾으면 데이터 패널을 실제값으로 채움 |
| `results/m0_period_vs_stimulus_lambda0.001.csv` | 자극 S̄ 스윕, Table 1 λ |
| `results/m0_period_vs_stimulus_lambda1.csv` | 자극 S̄ 스윕, λ=1 (게이트 덜 포화) |
| `results/m0_period_vs_lambda_S1.csv` | λ 스윕, S̄=1 |
| `results/m0_period_vs_epsH.csv` | ε_H 스윕 |
| `results/m0_period_vs_beta.csv` | β 스윕 |
| `results/m0_waveform_trace_baseline_epsH70.csv` | 기본값 파형 t,C,P,H |
| `results/m0_waveform_trace_epsH20.csv` | ε_H=20 파형 |
| `figures/m0_period_waveform_{main,supp}.{png,pdf,eps}` | 요약 그림 6패널 (a 데이터 ISI·IQR·n, b–d M0 스윕, e 듀티, f 한 주기 비교) |

## 위치·실행
    code/m0_period_waveform/
    make && ./run_m0_sweeps.sh && python3 plot_m0_period_waveform.py

plot 스크립트는 `../../Fluorescence*Spike Trains_data/2/*.csv` 가 있으면
용량별 중앙 ISI를 다시 계산하고, 없으면 하드코딩 값(2026-09-21 계산)을 씁니다.

## 실행 파일 사용법
    ./m0_period_waveform_sweep                       # 기본값 1회, 특징 한 줄
    ./m0_period_waveform_sweep S=3 beta=0.2 epsH=50  # 파라미터 지정
    ./m0_period_waveform_sweep sweep=epsH lo=10 hi=400 n=25 log=1
    ./m0_period_waveform_sweep trace=1 Tend=4000 Ttrans=3000
키: S alpha beta gamma rho mu lambda kS kC kH epsH C0 P0 H0 dt Tend Ttrans
출력 열: value, osc, Tosc, Cmax, Cmin, W50, rise50, decay50, duty(=W50/Tosc), Pmean, Hmin, Hmax
(과도기 Ttrans 이후 마지막 완전한 주기에서 계산; 기저선 = 직전 스파이크 사이 국소 최소)

## 2026-09-21 결과
- 기본값 Tosc = 202.2 → AIMS ≈200 재현(코드 검증).
- S̄ 0→10에서 주기 불변(202.25): λ=0.001이라 P가 2,000–5,500까지 쌓여 P/(1+P)≈1 포화.
  λ=1이면 S̄≥2에서만 진동, 주기 166→191로 증가(데이터와 반대), 변화폭 15%.
- 주기는 ε_H가 결정: Tosc ≈ 2.8 ε_H (log 기울기 1.11). β는 271→194 (1.4배). 진폭도 ε_H와 함께 증가.
- 파형: 모든 스윕에서 duty 0.38–0.60(넓은 진동). 데이터 0.16.
- 결론: 원식에는 P→H 경로가 없어 자극이 주기를 부호화하지 못함. G0에서
  IP3 의존 불활성화(κ_H(P) 또는 회복률의 P 의존) 추가 여부를 결정해야 함 → 02 계산.

## Figures
- `figures/m0_period_waveform_main.{png,pdf,eps}` (paper): data ISI vs dose; M0 period vs stimulus; one cycle vs measured spike
- `figures/m0_period_waveform_supp.{png,pdf,eps}` (supplement): period vs ε_H; period vs β; duty cycle
