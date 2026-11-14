"""
e1_reflecting_boundary.py — the moment formula (11) beyond the entrance regime.

For nu^2 > 2 the boundary 0 of the Jacobi form is regular.  With instantaneous reflection the moment functions satisfy
L M_n = -n M_{n-1}, M_n(H*) = 0 and zero flux in the natural scale, lim_{x->0} M_n'(x)/s'(x) = 0; because the speed
measure of (0, y) is finite for every nu, two integrations give the same formula (11).  Check: formula (linear ODE
system of a2 = quadrature of (11)) against Euler-Maruyama simulation with reflection at 0, three time steps.
Output: results/e1_reflecting_boundary.csv
"""
import os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import a2_exact_vs_expansion as A
H0, HS = 0.064, 0.408
def mc(a):
    nu, dt, seed, npath = a; rng = np.random.default_rng(seed); h = np.full(npath, H0); alive = np.arange(npath); tau = np.empty(npath); sq = np.sqrt(dt); n = 0; hit0 = np.zeros(npath, bool)
    while len(alive):
        x = h[alive]; x = x + (1 - x) * dt + nu * np.sqrt(np.clip(x * (1 - x), 0, None)) * sq * rng.standard_normal(len(alive)); neg = x < 0; hit0[alive[neg]] = True; x = np.abs(x); n += 1
        hit = x >= HS; tau[alive[hit]] = n * dt; h[alive] = x; alive = alive[~hit]
    m, s = tau.mean(), tau.std(); return dict(nu=nu, dt=dt, mean=m, sd=s, skew=float(np.mean(((tau - m) / s) ** 3)), frac_touching_0=float(hit0.mean()), se_mean=s / np.sqrt(npath))
if __name__ == "__main__":
    nus = [1.0, 1.3, 1.5, 2.0, 2.5, 3.0, 3.6]; jobs = [(nu, dt, 7, 200000) for nu in nus for dt in (4e-4, 1e-4, 2.5e-5)]
    with Pool(21) as pool: R = pd.DataFrame(pool.map(mc, jobs))
    rows = []
    for nu in nus:
        e = A.exact("jacobi", nu, h0=H0, hs=HS); g = R[R.nu == nu].sort_values("dt")
        # Richardson-type extrapolation of the weak error (order 1/2 at the absorbing threshold): fit m(dt) = m0 + c sqrt(dt)
        ex = {}
        for q in ("mean", "sd"): c = np.polyfit(np.sqrt(g.dt.values), g[q].values, 1); ex[q] = c[1]
        rows.append(dict(nu=nu, nu2=nu * nu, k=2 / nu**2, boundary="entrance" if nu * nu <= 2 else "regular", mean_formula=e["mean"], sd_formula=np.sqrt(e["var"]), skew_formula=e["skew"],
                         mean_mc_fine=g["mean"].iloc[0], sd_mc_fine=g.sd.iloc[0], skew_mc_fine=g["skew"].iloc[0], mean_mc_extrap=ex["mean"], sd_mc_extrap=ex["sd"],
                         err_mean=ex["mean"] / e["mean"] - 1, err_sd=ex["sd"] / np.sqrt(e["var"]) - 1, frac_paths_touching_0=g.frac_touching_0.iloc[0], se_rel=g.se_mean.iloc[0] / e["mean"]))
    T = pd.DataFrame(rows); T.to_csv(os.path.join(HERE, "results", "e1_reflecting_boundary.csv"), index=False); pd.set_option("display.width", 220); print(T.round(4).to_string(index=False))
