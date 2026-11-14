"""
d1_m0_scan_analysis.py — the M0 negative control with 40 000 Latin-hypercube points (manuscript: 1 500).

Same kernel, box and criteria as m0_joint_scan: C1 period ratio >= 6.95 (lenient: >= 4.61, ratio of the CI ends),
C2 duty <= 0.15 inside the window that attains the ratio, C3 amplitude ratio in [0.5, 2].  Every point that passes
C1(lenient) or C2 is traced at the low, middle and high stimulus of its window and screened for damped transients
(last-quarter amplitude below half of the first quarter), as in round2_checks/check4b.
Output: results/d1_m0_scan_summary.csv, results/d1_m0_candidates.csv
"""
import os, io, sys, subprocess, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "results"); BIN = os.path.join(HERE, "..", "m0_joint_scan", "m0_joint_scan")
PN = ["alpha", "beta", "gamma", "rho", "mu", "lambda", "kS", "kC", "kH", "epsH"]; LEVELS = np.geomspace(0.03, 30, 9); R_DATA, R_CI, DUTY = 6.95, 4.61, 0.15
def trace(a):
    idx, kv, Sv, tag = a; p = subprocess.run([BIN, "trace=1"] + [f"{k}={v}" for k, v in {**kv, "S": Sv}.items()], capture_output=True, text=True)
    if p.returncode or "\n" not in p.stdout.strip(): return dict(idx=idx, level=tag, S=Sv, ok=0, ratio_q4_q1=np.nan)
    x = pd.read_csv(io.StringIO(p.stdout)).C.values; q = len(x) // 4; amp = [x[j * q:(j + 1) * q].max() - x[j * q:(j + 1) * q].min() for j in range(4)]
    return dict(idx=idx, level=tag, S=Sv, ok=1, ratio_q4_q1=amp[3] / amp[0] if amp[0] > 0 else np.nan)
if __name__ == "__main__":
    S = pd.read_csv(os.path.join(OUT, "d1_m0_scan_40000.csv")); osc = S[S.n_osc >= 3].copy()
    osc["C1"] = osc.T_ratio >= R_DATA; osc["C1_ci"] = osc.T_ratio >= R_CI; osc["C2"] = (osc.duty_max <= DUTY) & (osc.i_lo >= 0); osc["C3"] = osc.amp_ratio.between(0.5, 2.0) & (osc.i_lo >= 0)
    cand = osc[(osc.C1_ci | osc.C2) & (osc.i_lo >= 0)]; jobs = []
    for r in cand.itertuples():
        kv = {k: float(getattr(r, k)) if k != "lambda" else float(cand.loc[r.Index, "lambda"]) for k in PN}
        for tag, i in (("lo", int(r.i_lo)), ("mid", int((r.i_lo + r.i_hi) // 2)), ("hi", int(r.i_hi))): jobs.append((int(r.idx), kv, float(LEVELS[i]), tag))
    with Pool(24) as pool: T = pd.DataFrame(pool.map(trace, jobs))
    worst = T.groupby("idx").ratio_q4_q1.min(); osc["sustained"] = osc.idx.map(worst).fillna(np.nan) >= 0.5
    rows = []
    for lab, m in (("points", pd.Series(True, index=S.index)),):
        rows.append(dict(criterion="points scanned", n=len(S)))
    for lab, m in (("oscillating at >= 3 stimulus values", pd.Series(True, index=osc.index)), ("C1: period ratio >= 6.95", osc.C1), ("C1 lenient: >= 4.61", osc.C1_ci), ("C2: duty <= 0.15", osc.C2), ("C3: amplitude ratio in [0.5, 2]", osc.C3),
                   ("C1 & C2", osc.C1 & osc.C2), ("C1 lenient & C2", osc.C1_ci & osc.C2), ("C1 lenient & C2, sustained at lo/mid/hi", osc.C1_ci & osc.C2 & osc.sustained),
                   ("C1 & C2 & C3", osc.C1 & osc.C2 & osc.C3), ("C1 lenient & C2 & C3, sustained", osc.C1_ci & osc.C2 & osc.C3 & osc.sustained), ("C2, sustained", osc.C2 & osc.sustained), ("C1 lenient, sustained", osc.C1_ci & osc.sustained)):
        rows.append(dict(criterion=lab, n=int(m.sum())))
    Sm = pd.DataFrame(rows); Sm.to_csv(os.path.join(OUT, "d1_m0_scan_summary.csv"), index=False); print(Sm.to_string(index=False))
    c = osc[(osc.C1_ci | osc.C2)].copy(); c["min_ratio_q4_q1"] = c.idx.map(worst); c.sort_values("T_ratio", ascending=False).to_csv(os.path.join(OUT, "d1_m0_candidates.csv"), index=False)
    pd.set_option("display.width", 230)
    print("max period ratio among points with duty <= 0.15 and sustained: %s" % (osc[osc.C2 & osc.sustained].T_ratio.max() if (osc.C2 & osc.sustained).any() else "none"))
    print("min duty among sustained points with ratio >= 4.61: %s" % (osc[osc.C1_ci & osc.sustained].duty_max.min() if (osc.C1_ci & osc.sustained).any() else "none"))
    print(c.sort_values("T_ratio", ascending=False).head(12)[["idx", "T_ratio", "duty_max", "amp_ratio", "S_lo", "S_hi", "lambda", "epsH", "min_ratio_q4_q1", "nclip"]].round(3).to_string(index=False))
