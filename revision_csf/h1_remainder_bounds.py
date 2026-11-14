"""
h1_remainder_bounds.py — non-asymptotic enclosures for Theorem 5.2.

Jacobi form, a = nu^2/2, theta = y d/dy, g = 1/(1-y).  The bounded solutions of the cumulant equations are
    u = A g,   v = A[2 a y u^2],   w = A[6 a y u v],     A f(y) = int_0^inf e^{-s} f(y e^{-a s}) ds  = (1 + a theta)^{-1} f,
kappa_n = int_{H0}^{H*} (u, v, w).  A is a positive averaging operator over [0, y] that commutes with theta, and
A f = f - a A theta f exactly.  All theta^m g are non-negative and increasing, hence so are u, v, w and their theta-images,
and A f <= f for increasing f.  This gives, for EVERY a > 0 (with reflection at 0 when nu^2 > 2):
    T0 - a I1                 <=  E tau   <=  T0 - a I1 + a^2 m2
    2 a I1 - a^2 (2 I1 + 4 K3) <=  Var tau <=  2 a I1
    6 a^2 K3 - a^3 c3         <=  kappa_3 <=  6 a^2 K3
with m2 = dz2 - 3 dz + T0, c3 = 18 K3 + 96 J, J = int y^3/(1-y)^4 dy; the lower/upper correction terms are exactly the
next coefficients of the expansion, so every truncation of the series is an enclosure.
This script verifies the enclosures against the exact cumulants for random geometries and nu in [0.02, 3.5].
Output: results/h1_remainder_bounds.csv
"""
import os, sys, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import a2_exact_vs_expansion as A2
def coeff(h0, hs):
    z0, zs = 1 / (1 - h0), 1 / (1 - hs); T0 = np.log(zs / z0); dz, dz2 = zs - z0, zs**2 - z0**2
    I1 = dz - T0; K3 = dz2 - 4 * dz + 2 * T0; m2 = dz2 - 3 * dz + T0; Jf = lambda z: z**3 / 3 - 1.5 * z**2 + 3 * z - np.log(z); J = Jf(zs) - Jf(z0)
    return dict(T0=T0, I1=I1, K3=K3, m2=m2, cv=2 * I1 + 4 * K3, c3=18 * K3 + 96 * J)
rng = np.random.default_rng(3); geos = [(0.064, 0.408)] + [(h0, h0 + (0.97 - h0) * rng.random() ** 1.0) for h0 in rng.uniform(0.005, 0.5, 14)]
rows = []
for h0, hs in geos:
    c = coeff(h0, hs)
    for nu in [0.02, 0.05, 0.1, 0.2, 0.281, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 3.5]:
        a = nu * nu / 2; e = A2.exact("jacobi", nu, h0=h0, hs=hs)
        b = dict(mean_lo=c["T0"] - a * c["I1"], mean_hi=c["T0"] - a * c["I1"] + a * a * c["m2"], var_lo=2 * a * c["I1"] - a * a * c["cv"], var_hi=2 * a * c["I1"], k3_lo=6 * a * a * c["K3"] - a**3 * c["c3"], k3_hi=6 * a * a * c["K3"])
        tol = 1e-7
        ok = (b["mean_lo"] <= e["mean"] * (1 + tol)) and (e["mean"] <= b["mean_hi"] * (1 + tol)) and (b["var_lo"] <= e["var"] * (1 + tol)) and (e["var"] <= b["var_hi"] * (1 + tol)) and (b["k3_lo"] <= e["k3"] * (1 + 1e-5)) and (e["k3"] <= b["k3_hi"] * (1 + 1e-5))
        rows.append(dict(H0=h0, Hs=hs, nu=nu, mean=e["mean"], var=e["var"], k3=e["k3"], **b, ok=ok,
                         rel_width_mean=(b["mean_hi"] - b["mean_lo"]) / e["mean"], rel_width_var=(b["var_hi"] - b["var_lo"]) / e["var"], rel_width_k3=(b["k3_hi"] - b["k3_lo"]) / e["k3"]))
R = pd.DataFrame(rows); R.to_csv(os.path.join(HERE, "results", "h1_remainder_bounds.csv"), index=False)
print("enclosures hold in %d of %d cases (%d geometries x 15 noise levels)" % (R.ok.sum(), len(R), len(geos)))
if not R.ok.all(): print(R[~R.ok].round(5).to_string(index=False))
pd.set_option("display.width", 220); print("\nreference geometry (0.064, 0.408): coefficients", {k: round(v, 5) for k, v in coeff(0.064, 0.408).items()})
print(R[(R.H0 == 0.064)][["nu", "mean_lo", "mean", "mean_hi", "var_lo", "var", "var_hi", "k3_lo", "k3", "k3_hi", "rel_width_mean", "rel_width_var", "rel_width_k3"]].to_string(index=False, float_format=lambda v: f"{v:.5g}"))
