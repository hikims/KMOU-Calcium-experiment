"""
a4_slopes_noise_forms.py — sigma-T_av slopes of the full model for the three H-noise forms (cell axis: eps_H varied;
dose axis: S varied), amplitudes as matched in a3 at the reference point.

In the chemical-Langevin form the rescaled noise is nu = 1/sqrt(N) whatever the recovery rate, so a fixed number of units N
means a common nu on both axes; for the sigma_H forms the cell axis is run at common nu (sigma_H ∝ eps_H^-1/2) and at
common sigma_H, the dose axis at common sigma_H (as in the manuscript).
Output: results/a4_slope_sweeps.csv, results/a4_slopes.csv
"""
import os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "round2_checks"))
import common as cm
OUT = os.path.join(HERE, "results"); BIN = os.path.join(HERE, "m1_noiseform_sde")
star = pd.read_csv(os.path.join(OUT, "a3_matched_segmented.csv")).set_index("form").amp
def job(a):
    lab, axis, val, extra = a; par = {**cm.REF, ("epsH" if axis == "cell" else "S"): val, **extra}
    if "nu0" in par: par["sigH"] = par.pop("nu0") / np.sqrt(par["epsH"])
    r = pd.concat([cm.run(BIN, {**par, "seed": s, "Tend": 3000 + 700 * (84 * (par["epsH"] / 1000) * (3.5 if par["S"] < 0.1 else 1.0))}, None) for s in (1, 2, 3)]).mean()
    return dict(case=lab, axis=axis, value=val, Tav=r.Tav, sd=r.sd, CV=r.CV, n=r.n)
if __name__ == "__main__":
    eps = np.geomspace(200, 8000, 9); Sv = np.geomspace(0.05, 30, 9); sj, s1, N = float(star["jacobi"]), float(star["one_minus_H"]), float(star["cle"]); jobs = []
    for lab, extra_cell, extra_dose in (("jacobi, common nu / fixed sigma_H", dict(hnoise=0, nu0=sj * np.sqrt(1000)), dict(hnoise=0, sigH=sj)),
                                         ("jacobi, common sigma_H", dict(hnoise=0, sigH=sj), None),
                                         ("1-H, common nu / fixed sigma_H", dict(hnoise=1, nu0=s1 * np.sqrt(1000)), dict(hnoise=1, sigH=s1)),
                                         ("1-H, common sigma_H", dict(hnoise=1, sigH=s1), None),
                                         ("chemical Langevin, fixed N", dict(hnoise=2, Nrec=N), dict(hnoise=2, Nrec=N))):
        jobs += [(lab, "cell", float(e), extra_cell) for e in eps]
        if extra_dose: jobs += [(lab, "dose", float(s), extra_dose) for s in Sv]
    with Pool(8) as pool: R = pd.DataFrame(pool.map(job, jobs))
    R.to_csv(os.path.join(OUT, "a4_slope_sweeps.csv"), index=False)
    Sl = R.groupby(["case", "axis"], sort=False).apply(lambda g: pd.Series(dict(slope=np.polyfit(np.log(g.Tav), np.log(g.sd), 1)[0], CV_min=g.CV.min(), CV_max=g.CV.max(), Tav_min=g.Tav.min(), Tav_max=g.Tav.max())), include_groups=False).reset_index()
    Sl.to_csv(os.path.join(OUT, "a4_slopes.csv"), index=False); pd.set_option("display.width", 200); print(Sl.round(3).to_string(index=False))
