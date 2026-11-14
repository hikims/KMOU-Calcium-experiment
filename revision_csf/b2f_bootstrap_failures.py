"""
b2f_bootstrap_failures.py — why some replicates of the parametric bootstrap admit no fit of the families with offset.

Regenerates the synthetic segments of every replicate of b2_family_comparison.py (same seeds; the fits do not use the
random generator) and tests the criterion: a segment contains an interval at or below  offset - 0.8 q  (q the frame step;
0.8 q is the reach of the quadrature nodes of the frame-step kernel in isifam.py), so that every family with offset has
zero likelihood.  Compares the criterion with the stored failure flag (a total log-likelihood beyond 1e9 in
results/b2_E_parametric_bootstrap_raw.csv) and writes results/b2_F_bootstrap_failures.csv.
Needs results/b1_fpt_density_*.pkl (b1_fpt_density.py) and results/b2_A_segments.csv.
"""
import os, numpy as np, pandas as pd
import b2_family_comparison as B, isifam as F
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
if __name__ == "__main__":
    segs = B.load_segments(); P = pd.read_csv(os.path.join(OUT, "b2_A_segments.csv")); assert list(P.id) == [s["id"] for s in segs]
    meta = [(s["n"] + s["ndrop"], s["q"]) for s in segs]; raw = pd.read_csv(os.path.join(OUT, "b2_E_parametric_bootstrap_raw.csv"))
    reach = float(np.max(np.abs(F.NODES))); rows = []
    for fi, (fam, d) in enumerate(B.FAMS):
        nm = B.NAMES[fi]; pars = list(zip(P["m_" + nm], P["cv_" + nm])); g = raw[raw.truth == nm].reset_index(drop=True)
        for r in range(len(g)):
            rng = np.random.default_rng(1000 * fi + r); nbelow = 0
            for (m, cv), (n, q) in zip(pars, meta):
                x = B.synth_segment(rng, fam, d, m, cv, n, q); nbelow += int((x <= B.DELTA - (reach * q if q > 0 else 0.0)).sum())
            rows.append(dict(truth=nm, replicate=r, intervals_below_support=nbelow, criterion=nbelow > 0, stored_failure=bool(abs(g.loc[r, "invgauss+off"]) > 1e9)))
    R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "b2_F_bootstrap_failures.csv"), index=False)
    print("replicates:", len(R), " criterion equals stored failure flag in:", int((R.criterion == R.stored_failure).sum()))
    print(R.groupby("truth", sort=False).stored_failure.sum().to_string())
