"""
d2_m0_directed_search.py — directed search for an M0 parameter set that joins the frequency range with a narrow spike
(a Latin-hypercube scan is sparse in ten dimensions).

Objective for one parameter set (9 stimulus values, as in the scan), minimised over every stimulus window in which all
points oscillate and the period does not increase:
    J = max(0, ln(4.61 / period ratio)) + max(0, ln(max duty / 0.15))           (J = 0: lenient C1 and C2 both met)
Nelder-Mead in log10 parameters inside the scan box, 300 evaluations, from 48 starting points of the 40 000-point scan
(12 with the largest period ratio, 12 with the smallest duty, 24 with the smallest J of the ratio-maximising window).
End points with J < 0.3 are traced at the two ends of their window and screened for damped transients.
Output: results/d2_m0_directed_search.csv
"""
import os, io, subprocess, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "results"); BIN = os.path.join(HERE, "..", "m0_joint_scan", "m0_joint_scan")
PN = ["alpha", "beta", "gamma", "rho", "mu", "lambda", "kS", "kC", "kH", "epsH"]
LO = np.log10([0.35, 0.026, 0.42, 0.24, 1e-4, 1e-4, 0.02, 0.008, 0.015, 7]); HI = np.log10([35.0, 2.6, 42.0, 24.0, 1e-2, 10.0, 2.0, 0.8, 1.5, 700])
ENV = {**os.environ, "OMP_NUM_THREADS": "1"}
def curve(x):
    p = subprocess.run([BIN, "curve=1", "ns=9"] + [f"{k}={10**v:.8g}" for k, v in zip(PN, x)], capture_output=True, text=True, env=ENV)
    txt = p.stdout[p.stdout.find("S,osc"):]
    return pd.read_csv(io.StringIO(txt)) if "S,osc" in p.stdout else None
def best_window(c):
    best = (9.0, 1.0, 1.0, -1, -1)
    if c is None: return best
    T, osc, duty = c["T"].values, c.osc.values, c.duty.values
    for i in range(9):
        if not osc[i]: continue
        for j in range(i + 1, 9):
            if not osc[j] or T[j] > T[j - 1] * 1.0001: break
            r = T[i] / T[j]; d = duty[i:j + 1].max(); J = max(0.0, np.log(4.61 / r)) + max(0.0, np.log(d / 0.15))
            if J < best[0]: best = (J, r, d, i, j)
    return best
def obj(u):
    x = LO + np.clip(u, 0, 1) * (HI - LO); pen = float(np.sum((u - np.clip(u, 0, 1)) ** 2)) * 10
    return best_window(curve(x))[0] + pen
def search(a):
    tag, x0 = a; u0 = (x0 - LO) / (HI - LO); r = minimize(obj, u0, method="Nelder-Mead", options=dict(maxfev=300, xatol=1e-3, fatol=1e-4, initial_simplex=None))
    x = LO + np.clip(r.x, 0, 1) * (HI - LO); c = curve(x); J, ratio, duty, i, j = best_window(c); sus = np.nan
    if J < 0.3 and i >= 0:
        rr = []
        for k in (i, j):
            p = subprocess.run([BIN, "trace=1"] + [f"{kk}={10**v:.8g}" for kk, v in zip(PN, x)] + [f"S={c.S.iloc[k]}"], capture_output=True, text=True, env=ENV)
            if "\n" in p.stdout.strip():
                txt = p.stdout[p.stdout.find("t,C"):]; y = pd.read_csv(io.StringIO(txt)).C.values; q = len(y) // 4; amp = [y[m * q:(m + 1) * q].max() - y[m * q:(m + 1) * q].min() for m in range(4)]; rr.append(amp[3] / amp[0] if amp[0] > 0 else 0)
        sus = min(rr) if rr else np.nan
    J0 = best_window(curve(x0))
    return dict(start=tag, J_start=J0[0], ratio_start=J0[1], duty_start=J0[2], J=J, ratio=ratio, duty=duty, amp_last_over_first_quarter=sus, nfev=r.nfev, **{k: 10**v for k, v in zip(PN, x)})
if __name__ == "__main__":
    S = pd.read_csv(os.path.join(OUT, "d1_m0_scan_40000.csv")); o = S[(S.n_osc >= 3) & (S.i_lo >= 0)].copy()
    o["J"] = np.maximum(0, np.log(4.61 / o.T_ratio)) + np.maximum(0, np.log(o.duty_max / 0.15))
    starts = [("ratio", r) for r in o.sort_values("T_ratio", ascending=False).head(12).itertuples()] + [("duty", r) for r in o.sort_values("duty_max").head(12).itertuples()] + [("J", r) for r in o.sort_values("J").head(24).itertuples()]
    jobs = [(f"{t}:{int(r.idx)}", np.log10([float(o.loc[r.Index, k]) for k in PN])) for t, r in starts]
    with Pool(28) as pool: R = pd.DataFrame(pool.map(search, jobs, chunksize=1))
    R = R.sort_values("J"); R.to_csv(os.path.join(OUT, "d2_m0_directed_search.csv"), index=False); pd.set_option("display.width", 230)
    print(R[["start", "J_start", "ratio_start", "duty_start", "J", "ratio", "duty", "amp_last_over_first_quarter", "nfev", "lambda", "epsH", "beta"]].round(3).head(20).to_string(index=False))
    print("searches reaching J = 0:", int((R.J <= 1e-9).sum()), "; of these sustained (ratio of last to first quarter amplitude >= 0.5):", int(((R.J <= 1e-9) & (R.amp_last_over_first_quarter >= 0.5)).sum()))
    print("smallest J among sustained end points: %.3f" % R[(R.amp_last_over_first_quarter >= 0.5) | R.amp_last_over_first_quarter.isna()].J.min())
