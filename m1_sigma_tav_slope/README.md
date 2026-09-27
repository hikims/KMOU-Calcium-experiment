# m1_sigma_tav_slope — the variability law σ_ISI ∝ T_av across cells (Proposition 1 / 3 check)

## Question
Across trains the data give log σ_ISI – log T_av slope 1.09 (bootstrap 95 % CI 1.01–1.15) and sd = 0.335 (T_av − 24 s).
If cells differ in the recovery time scale ε_H, which noise convention reproduces slope 1?

## Model
Same SDE kernel as `m1_isi_variability` (`m1_slope_sde.c`), S̄ = 1, ε_H swept 200–8000 (period 18–640 model units =
10–340 s with T0 = 0.533 s).  Noise amplitudes as tuned in `m1_isi_variability` at ε_H = 1000:
- M1, ν common: σ_H = ν/√ε_H with ν = 0.026·√1000 = 0.822   (Proposition 1: constant CV, slope 1)
- M1, σ_H common: σ_H = 0.026                                (sd ∝ ε_H^{3/2}: slope 1.5)
- M2: σ_P = 0.0167;  Ca-noise: σ_C = 0.39

## Files
`m1_slope_sde.c`, `Makefile`, `run_m1_slope_sweeps.sh` (4 sweeps × 9 ε_H, ~70 s), `plot_m1_sigma_tav_slope.py`
→ `figures/m1_sigma_tav_slope_{main,supp}.{png,pdf,eps}` (1×3), `results/m1_slope_summary.csv`, `results/m1_slope_data_trains.csv`

## Run
    make && ./run_m1_slope_sweeps.sh && python3 plot_m1_sigma_tav_slope.py

## Findings
- slopes: data 1.09 [1.01, 1.15]; M1 ν-common 1.02; M1 σ_H-common 1.49 (excluded); M2 0.95; Ca-noise 1.27 (excluded).
- The plan's Proposition 3 (P-noise gives slope 0.5) does NOT hold for multiplicative noise in the slow P (λ = 0.001):
  P fluctuations modulate the recovery rate slowly, so sd ∝ T_av as well.  M2 is instead excluded by the serial ISI
  correlation (0.93 vs data 0.03 after detrending) and corr(ISI, amplitude) (−0.93 vs 0.10), see m1_isi_variability.
  -> The theory section should state the discriminator as {slope, serial correlation, amplitude coupling}, not slope alone.
- Intercept: data T_ref ≈ 24 s, M1 (ν common) ≈ 0 s.  The data contain a deterministic ~20 s component per ISI that the
  model lacks (refractory / spike-and-early-recovery time is only ~3 s here).  Candidate: the H-noise should vanish
  during the spike and the early recovery (e.g. noise ∝ √(H(1−H)) is already small at H ≈ 0.06; a stronger
  suppression such as noise ∝ (H−H_min) or a longer deterministic refractory phase).  To be tested in the calibration.

## Figures
- `figures/m1_sigma_tav_slope_main.{png,pdf,eps}` (paper): σ vs T_av law; slopes
- `figures/m1_sigma_tav_slope_supp.{png,pdf,eps}` (supplement): CV vs T_av

## Update 2026-09-22 (data source, dose axis)
Data from `../data_dose_units/results/train_table.csv` (union, detrended sd).  The pooled log σ – log T_av slope is
**1.31 [1.18, 1.42]** (raw 1.16; within dose with dose fixed effects 1.18; per dose 0.9–1.5 with wide CIs) —
no longer 1.09.  The data mix two axes: cells at one dose (ε_H) and the same cells at different doses (S̄).
New sweeps along S̄ at ε_H = 1000 (`m1_slope_M1_nu_common_dose.csv`, `m1_slope_M2_dose.csv`): M1 with ν common
gives slope 1.02 along the cell axis and 1.61 along the dose axis (the first-passage CV rises as S̄ falls because
H* rises), M2 gives 0.95 / 0.84, σ_H-common 1.49, Ca-noise 1.27.  The data value lies between the two M1 axes and
above every M2 value; the statement "slope 1" becomes "slope between the cell-axis value 1 and the dose-axis
value 1.6", which only M1 with ν common satisfies.  T0 = 0.527 s (from data_dose_units).

## Update 2026-09-22 (unit = authors' segments)
Units are the authors' segments (sd after detrending inside the segment): slope 1.33 [1.23, 1.42]
(union trains detrended 1.31, raw 1.20, within dose 1.26).  Model axes unchanged (M1 ν common: cells 1.02,
dose 1.61; M2: 0.95 / 0.84).
