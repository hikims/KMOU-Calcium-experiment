"""
b2c_sensitivity.py — part C of b2_family_comparison.py alone (sensitivity variants of the likelihood comparison), so that
it can be regenerated without the model-as-data runs and the parametric bootstrap of that script.
Output: results/b2_C_sensitivity.csv
"""
import os, sys, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b2_family_comparison as B
if __name__ == "__main__":
    pd.set_option("display.width", 250); segs = B.load_segments(); out = []
    with Pool(int(os.environ.get("NPROC", 12))) as pool:
        for lab, sg, eps in (("all intervals + 1% contamination", B.load_segments(trim=0), 0.01), ("no detrending, doublets dropped", B.load_segments(detr=False), 0.0),
                             ("two burst segments removed", [s for s in segs if not s["burst"]], 0.0)):
            Lc, _, _ = B.fit_all(pool, sg, eps=eps); c = B.summary(Lc, sg, B.NAMES, "fpt_jacobi+off"); c.insert(0, "variant", lab); out.append(c)
    Cc = pd.concat(out); Cc.to_csv(os.path.join(B.OUT, "b2_C_sensitivity.csv"), index=False)
    print(Cc[["variant", "family", "logL", "dlogL_vs_best", "dlogL_vs_ref", "seg_lo", "seg_hi", "rec_lo", "rec_hi"]].to_string(index=False, float_format=lambda v: f"{v:.2f}"))
