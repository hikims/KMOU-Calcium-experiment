# m1_first_passage_theory — first-passage reduction of the ISI (Proposition 1) and its check against the full SDE

## Reduction
Between spikes C ≈ C_rest and Q ≈ 1, so with s' = s/ε_eff, ε_eff = ε_H/[fP·fC·Q]_rest,
    dH = (1 − H) ds' + ν √(H(1−H)) dW,   ν = σ_H √ε_eff,   H0 = H_min → absorb at H* = H_on
    ISI = T_sp + ε_eff τ_FP,   T_sp = T_det − ε_eff ln((1−H0)/(1−H*))   (spike + early recovery)
Moments of τ_FP by the Dynkin/Siegert integrals (m = E τ, s = SD τ; third moment for the skewness) with the
substitution z = y w^{1/k}, k = 2/ν², which makes the quadrature stable for every ν (m1_fpt_jacobi.c, header).
Consequences: CV(ISI) = ε_eff s/(T_sp + ε_eff m) depends on (ν, H0, H*) and on T_sp/ε_eff only — cells that differ in
ε_H share the CV if ν and H* are common (Proposition 1, slope 1 in m1_sigma_tav_slope); the stimulus enters through
fP (ε_eff) and H*, so period ratios are ε_H-independent (Proposition 2).

## Files
- `m1_fpt_jacobi.c` — quadrature of mean/SD/skewness; `mc=N` Monte-Carlo check of the 1-D SDE; sweeps over nu, Hs, H0
- `m1_full_sde.c` — copy of the full 3-variable SDE kernel (m1_isi_variability) for the reduction check
- `run_m1_fpt.sh` — ν sweep, H* sweeps for ν = 0.1/0.3/1, four Monte-Carlo checks, two full-model sweeps (~3 min)
- `plot_m1_first_passage_theory.py` → `figures/m1_first_passage_theory_{main,supp}.{png,pdf,eps}` (1×3),
  `results/m1_fpt_reduction_check.csv`

## Run
    make && ./run_m1_fpt.sh && python3 plot_m1_first_passage_theory.py

## Findings
- Quadrature = Monte-Carlo (mean, SD, CV) to < 1 % for ν = 0.1–2.
- Full 3-variable SDE vs theory (y = x): T_av and σ_ISI over the σ_H sweep (ν 0.04–3.6) and the ε_H sweep (200–8000)
  agree to a few %, with T_sp taken from each cell's deterministic period.
- Data CV 0.22 (IQR 0.17–0.26) ↔ ν ≈ 0.3 at H* = 0.41; CV is almost flat in H* for ν ≤ 0.3, so the stimulus
  (which moves H* and ε_eff) leaves the CV unchanged — consistent with the constant CV across doses in the data.
- The mean τ_FP falls below the deterministic value for ν ≳ 0.5 (noise-induced early crossings); at ν ≈ 0.3 the
  shift is 2 %, so the deterministic dose→period calibration of M0 carries over to M1.

## Figures
- `figures/m1_first_passage_theory_main.{png,pdf,eps}` (paper): first-passage moments; reduction check
- `figures/m1_first_passage_theory_supp.{png,pdf,eps}` (supplement): CV vs H*

## Update 2026-09-22 (error panel, ν conventions)
`results/m1_fpt_reduction_error_summary.csv` and panel (c): relative error of the reduced-process prediction vs the
full SDE by ν range — ν ≤ 0.5 (16 points, the data range): T_av median 2.7 % (max 4.9 %), σ median 4.4 % (max 7.5 %);
0.5 < ν ≤ 1: 10–17 %; ν > 1: 21–29 %.  The statement "within a few per cent" applies to ν ≤ 0.5 only.  The check
uses each cell's deterministic period for T_sp, so it tests the stochastic correction, not an independent period
prediction.  Conventions: ν₀ = σ_H √ε_H (what the run scripts set, 0.822 at the reference) vs ν = σ_H √ε_eff
(what enters the reduced process, 0.315); 0 is an entrance boundary of the reduced process only for ν² ≤ 2.
