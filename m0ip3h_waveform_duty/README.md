# m0ip3h_waveform_duty — 스파이크를 좁게(듀티 ≈ 0.1) 만들면서 자극 부호화를 유지할 수 있는가

**계산 목적**: 앞의 두 계산에서 (i) 원식은 자극이 주기를 못 바꾸고, (ii) IP3→H 속도 결합으로 주기 비 6.7배는 되지만
(iii) 파형은 어떤 파라미터에서도 듀티 0.47–0.56(넓은 진동)임이 확인됨. 데이터는 0.05–0.15.
어떤 최소 성분이 스파이크형 파형을 만드는지, 그 상태에서 자극 부호화가 유지되는지 확인.

## 추가한 성분 (기본값에서는 모두 꺼져 원식과 정확히 일치)
| 성분 | 식 | 기본값 | 근거 |
|---|---|---|---|
| 협동 Ca 활성화 | C^nA/(kA^nA + C^nA) (원식은 C/(1+C), nA=1) | nA=1, kA=1 | De Young–Keizer/Li–Rinzel의 m∞³: CICR 문턱 → 흥분성 |
| 기저 유입 | +Jin in dC/dt | Jin=0 | 협동 활성화에서는 C=0이 흡수 상태가 되므로 leak/influx 필요 (모든 표준 모델에 있음) |
| Ca 의존 H 가속 | rate × (1 + rC C²/(κH²+C²)) | rC=0 | Li–Rinzel τ_h(C): 불활성화는 회복보다 빠름 |

## 파일
| 파일 | 역할 |
|---|---|
| `m0ip3h_waveform_sweep.c` | 앞 계산의 실행 파일 + nA, kA, Jin, rC 키 (인터페이스 동일) |
| `run_m0ip3h_waveform_sweeps.sh` | 성분별 1D 스윕(rC·β·α·Jin), 참조 영역의 S̄·gH 스윕, 파형 trace → `results/` |
| `search_m0ip3h_waveform_duty.py` | 9개 파라미터 LHS 무작위 탐색(기본 600표본, 전 코어 사용, 2코어 ~4분) → `results/m0ip3h_waveform_duty_search.csv`, 로그 `results/search_log.txt` |
| `plot_m0ip3h_waveform_duty.py` | 요약 그림 → `figures/m0ip3h_waveform_duty_{main,supp}.{png,pdf,eps}` |

## 실행
    make && ./run_m0ip3h_waveform_sweeps.sh && python3 search_m0ip3h_waveform_duty.py 600 2 && python3 plot_m0ip3h_waveform_duty.py

## 2026-09-21 결과
- 원식 활성화(nA=1)에서는 펌프(β, κS)를 올리면 어디서나 진동이 사라지고, 진동이 있는 곳은 듀티 0.46–0.48로 고정.
  이유: 방출항이 Hill 1이라 CICR 문턱이 없어 Hopf형 완만 진동만 가능(구조적).
- 협동 활성화(nA=2)+기저 유입+Ca 의존 H 가속을 넣으면 이완진동이 되고, **β가 듀티를 결정**(β 0.7→4.1에서 듀티 0.57→0.08).
  원식의 Hopf 경계(β*≈0.27)는 더 이상 진동을 막지 않음.
- 참조 영역 (mech=2, ρ=0.024, gH=30, dH=10⁴, εH=1000, nA=2, kA=0.4, Jin=0.02, kS=0.25, rC=30, β=4, α=20):
  듀티 0.13이 S̄ 0.1–30 전 구간에서 일정, 주기 비 4.3배, 한 주기 모양이 실측 스파이크(15 µM, ISI 58 s)와 거의 일치.
- 소거 실험(S̄=1): nA=1, β=0.26, α=3.5로 되돌리면 진동 소멸; rC=0이면 듀티 0.52로 복귀; Jin=0은 듀티 0.05로 유지.
- 무작위 탐색 600표본: 108개 양쪽 자극 진동, **14개가 듀티 ≤ 0.15 & 주기 비 ≥ 4** (최고 비 7.3, 듀티 0.10–0.12).
- 남은 불일치: 진폭. 모델 C 피크/휴지 비가 수백~수천인데 ΔF/F₀는 1–3. 관측식(형광 포화·완충) 또는 완충 Ca 스케일 문제 → 다음 계산.

## Figures
- `figures/m0ip3h_waveform_duty_main.{png,pdf,eps}` (paper): duty vs β; reference regime T & duty vs S̄; one cycle vs measured spike
- `figures/m0ip3h_waveform_duty_supp.{png,pdf,eps}` (supplement): ablation; random search; raw traces
