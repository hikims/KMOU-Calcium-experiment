# m0_joint_scan — joint parameter scan of the ORIGINAL AIMS equations M0 (audit item 6-3)

## Question
Is there any parameter set of the unmodified reaction system (AIMS Math 2026, eqs. 2.6–2.8, no diffusion)
that gives, over some stimulus window, (C1) a period ratio equal to the data (train-level medians,
3 → 30 µM: 6.95; ratio of the CI ends 4.6), (C2) a narrow spike (duty W50/T ≤ 0.15) and (C3) an amplitude
that stays within a factor 2?  A finite scan cannot prove impossibility; it shows how close the best points
come and which criteria fail together.  This is the negative control for the modified terms used later.

## Method
Latin-hypercube sample in log10 space over the 10 parameters (one decade either side of Table 1; λ over
[1e-4, 10] because the P/(1+P) gate saturates at Table 1's λ = 0.001), 9 log-spaced stimulus values
S̄ ∈ [0.03, 30] per point.  The dose → S̄ map is monotone but unknown, so the period ratio is the maximum
T(S_i)/T(S_j) over any sub-window in which every point oscillates and T is non-increasing; duty and the
amplitude ratio are taken inside that window.  RK4, dt = 0.02, Tend = 4000 + 40 ε_H (period ≈ 2.8 ε_H);
half-width crossings interpolated; state clipping counted (`nclip`).  OpenMP over sample points.

## Result
The numbers depend on the number of sample points; `plot_m0_joint_scan.py` writes the current ones to
`results/m0_joint_scan_summary.md` (quote that file, not this README).  With 4 000 points (seed 1): 642 points
oscillate at >= 3 stimulus values, the largest period ratio over any monotone window is 4.45, one point is narrow
(duty <= 0.15), and **no point satisfies the period-ratio and duty criteria together** (with 1 500 points the
picture was the same).  Reading: in the original equations a stimulus-dependent period requires λ large enough
that the IP3 gate is not saturated, and then the period changes through the approach to a Hopf point (amplitude
→ 0), while a narrow spike needs a strong, fast Ca-inactivation that the term (Q(C) − H)/ε_H does not provide.
This is a finite scan of the stated box, not a proof of impossibility.

## Files
- `m0_joint_scan.c` — kernel: `npts=N` scan (OpenMP), `curve=1` stimulus curve of one point, `trace=1`
- `run_m0_joint_scan.sh [NPTS] [NS]` — scan + baseline curve → `results/m0_joint_scan.csv`, `results/m0_joint_baseline_curve.csv`
- `plot_m0_joint_scan.py` — criteria, pass rates, best candidates (re-runs curves/traces in parallel) →
  `results/m0_joint_scan_summary.csv`, `m0_joint_top20_by_T_ratio.csv`, `m0_joint_best_curves.csv`, `m0_joint_best_trace.csv`;
  `figures/m0_joint_scan_main.*` (a: ratio vs duty with the target box, b: period curves, c: pass rates),
  `figures/m0_joint_scan_supp.*` (a: ratio vs λ, b: min duty vs ε_H, c: oscillating points, d: best trace)
- Targets are read from `../data_dose_units/results/` (fallback: hard-coded 2026-09-22 values)

## Run
    make && ./run_m0_joint_scan.sh 4000 && python3 plot_m0_joint_scan.py
(~0.5 s per point per thread; 4 000 points ≈ 4 min on 8 threads)
