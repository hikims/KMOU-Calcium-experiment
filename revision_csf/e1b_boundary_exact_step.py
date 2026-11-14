"""
e1b_boundary_exact_step.py — second check of formula (11) for nu^2 > 2 with a scheme that treats the boundary 0 exactly.

Euler-Maruyama with reflection converges slowly when 0 is a regular boundary with k = 2/nu^2 < 1/2 (e1: the simulated
mean still drifts with the time step at nu >= 2.5).  Here each step freezes the factor (1-H) of the diffusion
coefficient and samples the resulting CIR transition exactly (non-central chi-square), which is the instantaneously
reflecting process at 0:  dH = (1-H) ds + s_n sqrt(H) dW,  s_n = nu sqrt(1-H_n).
Output: results/e1b_boundary_exact_step.csv
"""
import os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import a2_exact_vs_expansion as A
H0, HS = 0.064, 0.408
def mc(a):
    nu, dt, seed, npath = a; rng = np.random.default_rng(seed); h = np.full(npath, H0); alive = np.arange(npath); tau = np.empty(npath); n = 0; e = np.exp(-dt)
    while len(alive):
        x = h[alive]; s2 = nu * nu * (1 - x); c = s2 * (1 - e) / 4; x = c * rng.noncentral_chisquare(4 / s2, np.maximum(x * e / c, 1e-300)); n += 1
        hit = x >= HS; tau[alive[hit]] = n * dt; h[alive] = x; alive = alive[~hit]
    m, s = tau.mean(), tau.std(); return dict(nu=nu, dt=dt, mean=m, sd=s, skew=float(np.mean(((tau - m) / s) ** 3)), se=s / np.sqrt(npath))
if __name__ == "__main__":
    nus = [1.0, 2.0, 2.5, 3.0, 3.6]; dts = (1.6e-3, 4e-4, 1e-4)
    with Pool(15) as pool: R = pd.DataFrame(pool.map(mc, [(nu, dt, 11, 200000) for nu in nus for dt in dts]))
    rows = []
    for nu in nus:
        e = A.exact("jacobi", nu, h0=H0, hs=HS); g = R[R.nu == nu].sort_values("dt"); ex = {q: np.polyfit(np.sqrt(g.dt.values), g[q].values, 1)[1] for q in ("mean", "sd")}
        rows.append(dict(nu=nu, k=2 / nu**2, mean_formula=e["mean"], sd_formula=np.sqrt(e["var"]), **{f"mean_dt{d:g}": g[g.dt == d]["mean"].iloc[0] for d in dts}, mean_extrap=ex["mean"], sd_extrap=ex["sd"],
                         err_mean=ex["mean"] / e["mean"] - 1, err_sd=ex["sd"] / np.sqrt(e["var"]) - 1, se_rel=g.se.iloc[0] / e["mean"]))
    T = pd.DataFrame(rows); T.to_csv(os.path.join(HERE, "results", "e1b_boundary_exact_step.csv"), index=False); pd.set_option("display.width", 220); print(T.round(4).to_string(index=False))
