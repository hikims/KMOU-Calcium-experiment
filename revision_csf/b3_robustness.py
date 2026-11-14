"""
b3_robustness.py — additions to b2: (1) doublets dropped plus a uniform contamination component (eps = 1 %, 2 %) that absorbs
isolated long intervals; (2) offset proportional to the segment mean (delta_i = phi * mean_i) against the common offset;
(3) fitted noise parameters by dose (nu of the Jacobi form; N = 1/nu^2 of the chemical-Langevin form).
Output: results/b3_*.csv
"""
import os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import isifam as F, b2_family_comparison as B
OUT = B.OUT
def fit_prop(a):
    x, q, phi, fams = a; d = phi * x.mean()
    return [F.fit(f, x, q, d)[0] for f in fams]
if __name__ == "__main__":
    pd.set_option("display.width", 250); fmt = lambda v: f"{v:.2f}"; segs = B.load_segments()
    with Pool(28) as pool:
        out = []
        for eps in (0.01, 0.02):
            L, _, _ = B.fit_all(pool, segs, eps=eps); c = B.summary(L, segs, B.NAMES, "fpt_jacobi+off"); c.insert(0, "eps", eps); out.append(c)
        C = pd.concat(out); C.to_csv(os.path.join(OUT, "b3_contamination.csv"), index=False); print("1. doublets dropped + uniform contamination"); print(C[["eps", "family", "logL", "dlogL_vs_best", "dlogL_vs_ref", "seg_lo", "seg_hi", "rec_lo", "rec_hi", "n_best"]].to_string(index=False, float_format=fmt))
        fams = ["gamma", "invgauss", "lognorm", "fpt_jacobi", "fpt_cle"]; rows = []
        for phi in (0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3):
            L = np.array(pool.map(fit_prop, [(s["x"], s["q"], phi, fams) for s in segs], chunksize=4)); rows.append(dict(phi=phi, **dict(zip(fams, L.sum(0)))))
        P = pd.DataFrame(rows); P.to_csv(os.path.join(OUT, "b3_proportional_offset.csv"), index=False); print("\n2. total logL with an offset proportional to the segment mean (compare with b2_B: common offset)"); print(P.to_string(index=False, float_format=fmt))
    A = pd.read_csv(os.path.join(OUT, "b2_A_segments.csv")); A["N_cle"] = 1 / A.nu_cle**2; A["cv_isi_jac"] = A["cv_fpt_jacobi+off"] * A["m_fpt_jacobi+off"] / (A["m_fpt_jacobi+off"] + B.DELTA)
    D = A.groupby("dose").agg(segments=("n", "size"), nu_jacobi_median=("nu_jacobi", "median"), nu_q25=("nu_jacobi", lambda v: v.quantile(.25)), nu_q75=("nu_jacobi", lambda v: v.quantile(.75)),
                              N_cle_median=("N_cle", "median"), N_q25=("N_cle", lambda v: v.quantile(.25)), N_q75=("N_cle", lambda v: v.quantile(.75)), mean_recovery_s=("m_fpt_jacobi+off", "median"))
    D.to_csv(os.path.join(OUT, "b3_noise_by_dose.csv")); print("\n3. fitted noise parameters by dose"); print(D.round(3).to_string())
    print("Spearman correlation of N with dose: %.3f; of nu_jacobi with dose: %.3f" % (A[["N_cle", "dose"]].corr(method="spearman").iloc[0, 1], A[["nu_jacobi", "dose"]].corr(method="spearman").iloc[0, 1]))
