# m1_particle_filter_synthetic — gate 2: synthetic recovery of (ε_H, σ_H, β) with a bootstrap particle filter

## What it does
- `m1_pf_sde.c` — bootstrap PF for the full M1 SDE (stereotyped reference, noise in H) with observation
  F = b + a C + e, e ~ N(0, σ_obs²); particles start at random phases of the deterministic cycle (hidden initial state);
  systematic resampling at every observation.  Modes: `simulate=1` (synthetic record t,C,P,H,F), `obs=FILE`
  (marginal log-likelihood), `obs=FILE filt=1` (filtered means of C, P, H).
- `run_pf_synthetic.py [NPART NPOINTS NSEED]` — simulates one record (T_obs = 2000 model units ≈ 18 min, Δt_obs = 2 ≈ 1 s,
  ~24 spikes; truth ε_H = 1000, σ_H = 0.026, β = 4, a = 0.11, σ_obs = 0.05) and computes profile log-likelihoods over
  each parameter (×/÷ 2, others at truth), NSEED PF seeds per point, all cores.  ~10 s per likelihood at N = 1000.
- `plot_m1_particle_filter_synthetic.py` → `figures/m1_particle_filter_synthetic_summary.{png,pdf,eps}` (2×2)
- `pmmh_synthetic.py [NITER NPART STEP]` — PMMH driver (random walk on log parameters), optional; see cost note inside.

## Run
    make && python3 run_pf_synthetic.py 1000 9 2 && python3 plot_m1_particle_filter_synthetic.py
(the shipped results used NPOINTS = 7)

## Findings (this record, N = 1000)
- ε_H and β are sharply identified: the profile falls by > 100 log-units within ±25 % of the truth; the filter tracks
  the spikes and the hidden H (filtered H follows the true H to within ~0.05).
- σ_H is weakly identified from a single 24-spike record: log-likelihood flat (within the PF noise, sd ≈ 5) between
  0.8× and 2× the truth; only the lower side (< 0.6×) is excluded.  Expected — σ_H is only visible through the ISI
  scatter, and 24 ISIs give CV ± 30 %.  Consequences for the plan: (i) use the full-length records (5000–8000 s,
  60–150 spikes) for σ_H; (ii) treat σ_H (or ν) as a population parameter shared across cells (hierarchical), which
  is exactly the test "σ_H posteriors overlap" in §6.2; (iii) profile over ε_H and β per cell is cheap and sharp.
- PF likelihood noise: sd ≈ 4–5 log-units at N = 1000, ≈ 2 at N = 4000 (40 s per evaluation).  PMMH therefore needs
  N ≥ 4000 or a smoother observation model (larger σ_obs) — a numerical setting to fix on the 15 µM development set.
- ESS collapses to ~1 at spike times (sharp spikes vs σ_obs = 0.05): the filter survives through resampling, but an
  auxiliary/guided proposal would be the standard improvement if real records prove harder.

Gate 2 verdict: passed for (ε_H, β); σ_H requires pooling across cells or longer records — not a failure of the
model but a design constraint for the calibration step.

## Figures
- `figures/m1_particle_filter_synthetic_main.{png,pdf,eps}` (paper): hidden H recovered; profiles of ε_H and β
- `figures/m1_particle_filter_synthetic_supp.{png,pdf,eps}` (supplement): synthetic record with the filter; σ_H profile
