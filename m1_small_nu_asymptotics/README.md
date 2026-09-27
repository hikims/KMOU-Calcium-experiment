# m1_small_nu_asymptotics — small-ν expansion of the first-passage moments (Corollary of Theorem 2)

## Result (closed forms; a = ν²/2, z0 = 1/(1−H0), z* = 1/(1−H*), T0 = ln(z*/z0))
    I1 = z* − z0 − T0,          K3 = z*² − z0² − 4(z* − z0) + 2 T0
    E[τ]  = T0 − a I1 + a² (z*² − z0² − 3(z* − z0) + T0) + O(a³)
    Var τ = 2a I1 + a² (−4(z*² − z0²) + 14(z* − z0) − 6 T0) + O(a³)
    κ3    = 6 a² K3 + O(a³)
Hence to leading order  sd = ν √I1,  E[τ] = T0 − Var/2  (the mean shift is minus half the variance),
skew = (3/2) ν K3 / I1^{3/2},  and the shape ratio  skew/CV → (3/2) K3 T0 / I1²  is a number fixed by (H0, H*).
For the stereotyped regime (H0 = 0.064, H* = 0.408): I1 = 0.1627, K3 = 0.1449, sd/ν = 0.4034 (quadrature 0.4027),
shift/ν² = 0.0814 (0.0812), skew/ν = 3.31 (quadrature 3.5 at ν ≤ 0.1, converging with the grid), skew/CV → 3.76.

Derivation: the inner integral ∫₀¹ F(y w^a) dw = Σ_m (−a)^m (y∂_y)^m F(y) (∫₀¹ (ln w)^m dw = (−1)^m m!); in z = 1/(1−y)
every term is a polynomial in z times a power of ln z, so the nested Dynkin integrals are elementary
(`derive_small_nu.py`, sympy, < 2 s).

## Checks and consequences
- Quadrature (Simpson × Gauss–Legendre, OpenMP over the y grid) and Monte-Carlo (200 000 paths per ν, per-thread RNG)
  agree with the expansion up to ν ≈ 0.3 (CV ≈ 0.25, the data range); beyond, the exact moments saturate.
- The FPT is more skewed than a gamma law (2 CV) and close to an inverse Gaussian (3 CV) in the data range;
  the shape ratio of the model (3.76) is above the data's 2.4 (train medians): with H-noise alone the model
  predicts ISI distributions somewhat more skewed than observed.  `supp (c)` shows which (H0, H*) would give 2.4.
- The third cumulant is a cubic cancellation; at ν ≤ 0.02 the quadrature needs n_y ≥ 16 000 to be within 1 %
  (`supp (a)`), which is why m1_first_passage_theory's small-ν skewness was slightly high.

## Files
- `derive_small_nu.py` → `results/asymptotic_coefficients.{csv,txt}`
- `m1_fpt_asymptotics.c` — `sweep=nu` (quadrature + closed forms), `mc=N [hist=1]`, `ratio=1` (shape ratio vs H*)
- `run_m1_small_nu.sh` → `results/m1_asym_vs_nu.csv`, `m1_asym_convergence.csv`, `m1_asym_mc.csv`,
  `m1_asym_hist_nu{0.3,1.0}.csv`, `m1_asym_ratio_vs_Hs_H0*.csv`
- `plot_m1_small_nu_asymptotics.py` → `results/m1_asym_summary.csv`; `figures/m1_small_nu_asymptotics_main.*`
  (a: mean shift and sd vs ν, b: skewness vs ν with gamma / inverse-Gaussian references),
  `figures/m1_small_nu_asymptotics_supp.*` (a: grid convergence, b: FPT density at ν = 0.3, c: shape ratio vs H*)

## Run
    make && ./run_m1_small_nu.sh && python3 plot_m1_small_nu_asymptotics.py       (~6 min on 2 threads)
