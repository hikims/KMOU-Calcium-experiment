"""
b1_fpt_density.py — density of the first-passage time tau (H0 -> H*) of the reduced recovery process, tabulated in nu.

Backward equation for the survival function S(x, s) = P_x(tau > s):  S_s = b S_x + a D S_xx,  S(H*, s) = 0,  S(x, 0) = 1,
a = nu^2/2.  Boundary x = 0: entrance (Jacobi, D(0) = 0) -> the equation itself, S_s = b S_x; regular (D(0) > 0) -> reflection.
Non-uniform grid refined towards H* (boundary layer of width a D/b), central differences (upwind where the cell Peclet
number exceeds 1.9, only near x = 0 for the Jacobi form), Crank-Nicolson with four backward-Euler start-up half steps.
Density f(s) = -(A S)(H0).  Checks: normalisation and the first three cumulants against the exact values of a2.
Output: results/b1_fpt_density_<form>.pkl  (list over nu of dict(nu, mean, sd, skew, ds, f)), results/b1_density_check.csv
"""
import os, sys, pickle, numpy as np, pandas as pd
from scipy.linalg import solve_banded
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "results")
sys.path.insert(0, HERE)
from a2_exact_vs_expansion import FORMS, exact, H0, HS

def make_grid(a, b, D, growth=1.03):
    """grid on [xlo, H*]; xlo = H0/2 with reflection when the path cannot reach it (exponent > 40), else 0"""
    xlo = H0 / 2 if 2 * b(H0) * (H0 / 2) / (2 * a * max(D(H0), D(H0 / 2))) > 40 else 0.0
    hmax = min(2.0e-4, 1.5 * a * min(D(H0 / 2), D(HS), D(H0)) / b(H0 / 2))       # cell Peclet number < 1.9 above H0/2
    l = a * D(HS) / b(HS); h0 = min(hmax, max(l / 8, 1e-7)); xs = [HS]; h = h0
    while xs[-1] - h > xlo: xs.append(xs[-1] - h); h = min(hmax, h * growth)
    x = np.array(xs[::-1]); x[0] = xlo
    j0 = int(np.argmin(np.abs(x - H0))); x[j0] = H0
    return x, j0, xlo

def density(form, nu, nper=160, zmax=16.0):
    b, D = FORMS[form]; a = nu**2 / 2; x, j0, xlo = make_grid(a, b, D); N = len(x) - 1          # unknowns 0..N-1, S_N = 0
    e = exact(form, nu); m, sd = e["mean"], np.sqrt(e["var"]); ds = sd / nper; nstep = int(np.ceil((m + zmax * sd) / ds))
    hm = np.diff(x)[:-1]; hp = np.diff(x)[1:]; xi = x[1:-1]; bb, dd = b(xi), a * D(xi)
    lo = 2 * dd / (hm * (hm + hp)); di = -2 * dd / (hm * hp); up = 2 * dd / (hp * (hm + hp))
    pe = bb * np.maximum(hm, hp) / np.maximum(dd, 1e-300); cen = pe < 1.9
    lo += np.where(cen, -bb * hp / (hm * (hm + hp)), 0.0); di += np.where(cen, bb * (hp - hm) / (hm * hp), -bb / hp); up += np.where(cen, bb * hm / (hp * (hm + hp)), bb / hp)
    A = np.zeros((3, N))                                   # banded storage: row 0 upper, 1 diagonal, 2 lower
    A[1, 1:] = di; A[0, 2:] = up[:-1]; A[2, 0:N - 1] = lo  # sub-diagonal entry of row j sits at A[2, j-1]
    h1 = x[1] - x[0]
    if xlo > 0 or D(0.0) > 0: A[1, 0] = -2 * a * D(xlo) / h1**2; A[0, 1] = 2 * a * D(xlo) / h1**2      # reflection
    else: A[1, 0] = -b(0.0) / h1; A[0, 1] = b(0.0) / h1                                       # entrance: S_s = b S_x
    def matvec(S):
        r = A[1] * S; r[:-1] += A[0, 1:] * S[1:]; r[1:] += A[2, :-1] * S[:-1]; return r
    def lhs(theta_ds):
        M = -theta_ds * A.copy(); M[1] += 1.0; return M
    S = np.ones(N); f = np.zeros(nstep + 1)
    Mbe = lhs(ds / 2); Mcn = lhs(ds / 2)                   # backward Euler with step ds/2 has the same matrix as CN with step ds
    for _ in range(4): S = solve_banded((1, 1), Mbe, S)     # two full steps as four BE half steps
    f[2] = -matvec(S)[j0]
    for n in range(3, nstep + 1):
        S = solve_banded((1, 1), Mcn, S + (ds / 2) * matvec(S)); f[n] = -matvec(S)[j0]
    f = np.maximum(f, 0.0); s = ds * np.arange(nstep + 1)
    tz = lambda g: float(np.sum(g) * ds)                    # trapezoid (end values vanish)
    n0 = tz(f); m1 = tz(f * s) / n0; v = tz(f * (s - m1) ** 2) / n0; k3 = tz(f * (s - m1) ** 3) / n0
    return dict(form=form, nu=nu, mean=m, sd=sd, skew=e["skew"], ds=ds, f=f / n0, norm=n0, mean_num=m1, sd_num=np.sqrt(v), skew_num=k3 / v**1.5, npts=N)

NU = {"jacobi": np.geomspace(0.022, 1.40, 211), "one_minus_H": np.geomspace(0.011, 0.75, 211), "cle": np.geomspace(0.011, 0.75, 211)}
def _job(a): return density(*a)
if __name__ == "__main__":
    forms = sys.argv[1:] or ["jacobi", "cle"]
    chk = []
    with Pool(26) as pool:
        for form in forms:
            res = pool.map(_job, [(form, float(nu)) for nu in NU[form]], chunksize=2)
            pickle.dump(res, open(os.path.join(OUT, f"b1_fpt_density_{form}.pkl"), "wb"))
            for r in res: chk.append({k: r[k] for k in ("form", "nu", "npts", "norm", "mean", "mean_num", "sd", "sd_num", "skew", "skew_num")})
    C = pd.DataFrame(chk); C["CV"] = C.sd / C["mean"]; C["e_norm"] = C.norm - 1; C["e_mean"] = C.mean_num / C["mean"] - 1; C["e_sd"] = C.sd_num / C.sd - 1; C["e_skew"] = C.skew_num / C["skew"] - 1
    C.to_csv(os.path.join(OUT, "b1_density_check.csv"), index=False)
    pd.set_option("display.width", 200); print(C.iloc[::15][["form", "nu", "npts", "CV", "skew", "e_norm", "e_mean", "e_sd", "e_skew"]].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("max abs errors:", C[["e_norm", "e_mean", "e_sd", "e_skew"]].abs().max().to_dict())
