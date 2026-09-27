# 칼슘 논문 계산 코드 전체 묶음 (2026-09-22 기준 최신)

이 zip을 풀면 `code/` 폴더가 생깁니다. `Calcium_experiment/` 아래에 두면(데이터 폴더들과 나란히) 상대경로
`../../*Spike Trains_data/` 등이 그대로 맞습니다.

    Calcium_experiment/
      Fluorescence Intensity Data from HEK293 Ca²⁺ Spike Trains_data/   (1/, 2/)
      Omnitemporal choreographies Ca2+ signaling_data/                 (Source Data_Emrich.xlsx)
      code/                                                             <- 이 묶음

모든 폴더에 소스, 실행 스크립트, README, `results/`(CSV), `figures/`(png/pdf/eps)가 들어 있습니다.
`results/`는 제 쪽에서 계산한 값이라 그대로 그림을 만들 수 있고, 무거운 계산은 원하실 때 다시 돌리면 됩니다
(m0_joint_scan은 1 500점 결과; 4 000점은 `./run_m0_joint_scan.sh 4000`).

## 폴더 (계산 순서)
| 폴더 | 내용 | 처음부터 재계산 시간(8코어) |
|---|---|---|
| data_dose_units | 사건 원장, train/segment 표, 용량 지도(끝점 규약) — **다른 폴더가 읽는 기준 표** | 6 s |
| data_qc_splits | 원자료 QC, 15 µM dev/lock 분할 | 1 min |
| m0_period_waveform | 원식 M0 주기·파형 스윕 | 1 min |
| m0_joint_scan | 원식 M0 공동 파라미터 스캔(OpenMP) | 4 min (4000점) |
| m0ip3h_stimulus_encoding | IP3→H 결합 방식과 용량 부호화 | 10 min(탐색 포함) |
| m0ip3h_waveform_duty | 협동 활성화·펌프·duty | 10 min(탐색 포함) |
| m0ip3h_amplitude_observation | 형광 관측함수, 진폭 | 2 min |
| m0ip3h_amplitude_invariance | qoff와 정형성 | 1 min |
| m0ip3h_bifurcation_structure | 고정점·고유값·진동 영역 격자(OpenMP) | 10 min |
| m1_first_passage_theory | 환원 과정 Dynkin 모멘트, 전체 SDE 검증 | 10 min |
| m1_small_nu_asymptotics | 작은 ν 전개(sympy) + 고해상 구적·MC(OpenMP) | 3 min |
| m1_isi_variability | H/P/Ca 잡음 vs 데이터(segment 단위, 같은 전처리) | 10 min |
| m1_sigma_tav_slope | σ–T_av 기울기(세포 축·용량 축) | 15 min |
| m1_particle_filter_synthetic | 입자 필터 합성 시험(보충) | 30 min |
| step_fold_change | step 기록의 fold change | 10 s |
| emrich_replication | Emrich 자료 재분석(보충) | 30 s |
| emrich_jin_prediction | J_in 시나리오(보충, OpenMP) | 3 min |
| assemble_paper_figures | 본문 Fig 1–7, 보충 S1–S25 조립 | 1 min |
| step1_ode | 초기 단계 코드(보관용) | — |

## 실행
- 그림만 다시 만들기(저장된 results 사용): `./RUN_ALL.sh`  (~20분; 컴파일은 자동)
- 전부 처음부터: `./RUN_FULL.sh` (~2시간)
둘 다 `code/` 안에서 실행합니다. 필요: C 컴파일러(macOS: `xcode-select --install`), OpenMP(`brew install libomp`),
python3 + numpy pandas scipy matplotlib openpyxl sympy.
