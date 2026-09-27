# m1_isi_variability — M1: bounded H-noise and the within-cell ISI statistics

## Question
Does bounded noise in the receptor-recovery variable H (M1) reproduce the within-cell ISI statistics of the data
(CV ≈ 0.22, right-skewed ISI distribution, amplitude CV ≈ 0.15, corr(ISI, amplitude) ≈ 0.1, serial ISI correlation
≈ 0 after detrending), and do the same statistics rule out noise in P (M2) or in Ca?

## Model
Deterministic part = stereotyped reference of `m0ip3h_amplitude_invariance` (qoff=0, kR=1, rC=300), S̄ = 1.
    dH = (deterministic) dt + σ_H √(H(1−H)) dW      (M1, Jacobi-type, keeps 0 ≤ H ≤ 1)
    dP = (deterministic) dt + σ_P P dW              (M2)
    dC = (deterministic) dt + σ_C √C dW             (Ca-triggered comparison)
Euler–Maruyama, dt = 0.01, reflection at the boundaries.  `nuH` sets σ_H = ν/√ε_H (common rescaled intensity across
cells that differ in ε_H; used in `m1_sigma_tav_slope`).

## Files
- `m1_isi_sde.c` — kernel; modes: statistics line / `sweep=` / `isi=1` (spike list: t, ISI, H at onset, amplitude) / `trace=1`
- `run_m1_isi_sweeps.sh` — CV etc. vs σ_H, σ_C, σ_P (13 points each, ~1400 ISIs per point; ~45 s) → `results/`
- `plot_m1_isi_variability.py` — tunes each noise amplitude to the data CV, writes the ISI lists, the trace and
  `results/m1_noise_location_table.csv`; figure `figures/m1_isi_variability_{main,supp}.{png,pdf,eps}` (2×2)

## Run
    make && ./run_m1_isi_sweeps.sh && python3 plot_m1_isi_variability.py

## Data handling
Base csv files only (`*_2.csv`, `*_3.csv`… are segments of the same trains); a train = (file, ST, stimulus) with ≥ 10 spikes
(172 trains).  Serial correlation is reported raw and after removing a linear trend of ISI vs time within the train.

## Findings
- σ_H ≈ 0.026 gives CV 0.21, skewness 0.8 (data 0.52, IQR 0.26–0.95), amplitude CV 0.09 (data 0.15),
  corr(ISI, amp) −0.09 (data 0.10), serial corr 0.02 (data 0.30 raw, 0.03 detrended).  Spikes start at H ≈ H_on
  (sd of H at onset 0.02), which is why the amplitude is decoupled from the preceding ISI.
- Ca-noise tuned to the same CV: amplitude CV 0.26 and corr(ISI, amp) +0.37, skewness 1.8 — too variable and coupled.
- P-noise tuned to the same CV: amplitude CV 0.04, corr(ISI, amp) −0.93, serial corr +0.93 (P is slow, λ = 0.001,
  so its fluctuations modulate whole runs of ISIs) — the raw data serial correlation (0.30) is mostly linear drift;
  after detrending it is 0.03.  A residual slow component (M12 or drift, M1d) is left for the calibration step.
- ISI/mean distribution of the model matches the pooled data histogram and a gamma with the same CV.

## Figures
- `figures/m1_isi_variability_main.{png,pdf,eps}` (paper): trace; ISI distribution; noise-location statistics
- `figures/m1_isi_variability_supp.{png,pdf,eps}` (supplement): CV/skewness vs σ_H

## Update 2026-09-22 (data source)
The data statistics now come from `../data_dose_units/results/train_table.csv` (union of all processed files,
172 trains) and are **detrended** (linear trend in time removed per train) because the union makes trains up to
2 h long and slow drift inflates CV, skewness and the serial correlation.  Targets: CV 0.213 [0.14, 0.34],
skew 0.55, amp CV 0.145, corr(ISI, amp) 0.05, serial −0.02 (raw: 0.253 / 0.77 / — / 0.04 / 0.23).
Tuned σ_H* = 0.0253 (was 0.026); the noise-location table and the conclusions are unchanged
(H-noise: CV 0.21, skew 0.76, amp CV 0.09, corr −0.08, serial 0.03).

## Update 2026-09-22 (unit = authors' segments; same-pipeline comparison)
Data unit is now one of the authors' stationary segments (`segment_table.csv`, 196 segments, median 22 spikes),
detrended inside the segment: CV 0.196 [0.14, 0.31], skew 0.49, amp CV 0.135, corr 0.05, serial −0.03
(union trains, detrended: 0.213 / 0.55 / 0.145 / 0.05 / −0.02 — `unit_comparison.csv`).  The synthetic ISI lists
are cut into segments with the same length distribution and passed through the same detrending and formulas
("segmented" rows of `m1_noise_location_table.csv`).  With that, H-noise gives CV 0.185, skew 0.47, amp CV 0.08,
corr −0.13, serial −0.06 — the earlier skewness gap (0.76 vs 0.5) was a finite-sample effect.  P-noise: CV 0.09,
serial 0.66; Ca-noise: skew 0.74, amp CV 0.24, corr 0.50.  σ_H* = 0.0232.  This is a comparison of the three
specified noise models, not a general identification.
