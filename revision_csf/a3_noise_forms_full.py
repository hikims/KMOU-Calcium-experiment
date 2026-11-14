"""
a3_noise_forms_full.py — full three-variable model with three diffusion coefficients for the H-noise, same pipeline as Table 5/S8.

forms:  jacobi       sigH sqrt(H(1-H))                      (manuscript)
        one_minus_H  sigH sqrt(1-H)
        cle          sqrt((aH(1-H)+bH H)/N), state-dependent rates of the model (chemical Langevin, N two-state units)
Step 1  amplitude grid -> segmented, detrended CV; interpolate the amplitude that gives the data median CV.
Step 2  at the matched amplitude: 6 seeds, segmented statistics (median over segments; seed sd; bootstrap CI).
Step 3  long unsegmented runs (Tend = 2e6, 4 seeds): CV, skew, skew/CV, amplitude CV, correlations, onset and reset values of H.
Output: results/a3_grid.csv, a3_matched_segmented.csv, a3_matched_long.csv
"""
import os, sys, io, subprocess, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "round2_checks"))
import common as cm
OUT = os.path.join(HERE, "results"); BIN = os.path.join(HERE, "m1_noiseform_sde")
FORM = {"jacobi": ("sigH", dict(hnoise=0), np.geomspace(0.014, 0.04, 9)),
        "one_minus_H": ("sigH", dict(hnoise=1), np.geomspace(0.007, 0.022, 9)),
        "cle": ("Nrec", dict(hnoise=2), np.geomspace(25, 160, 9))}
def job(a):
    form, val, seed, Tend, mode = a; key, extra, _ = FORM[form]
    return a, cm.run(BIN, {key: val, "seed": seed, "Tend": Tend, **extra}, mode)
def long_stats(lst):
    x, am = lst.ISI.values, lst.amp.values; m = x.mean(); s = x.std(ddof=1); sk = float(np.mean(((x - m) / x.std()) ** 3))
    return dict(n=len(x), Tav=m, CV=s / m, skew=sk, R=sk / (s / m), amp_cv=am.std(ddof=1) / am.mean(), corr_isi_amp=np.corrcoef(x, am)[0, 1], serial=np.corrcoef(x[:-1], x[1:])[0, 1],
                Hon_mean=lst.Hon.mean(), Hon_sd=lst.Hon.std(), Hoff_mean=lst.Hoff_prev.mean(), Hoff_sd=lst.Hoff_prev.std(), corr_Hoff_Hon=np.corrcoef(lst.Hoff_prev, lst.Hon)[0, 1],
                corr_ISI_Hon=np.corrcoef(x, lst.Hon)[0, 1], corr_ISI_Hoff=np.corrcoef(x, lst.Hoff_prev)[0, 1])
if __name__ == "__main__":
    target = cm.DMED["CV"]; print("data medians:", {k: round(v, 3) for k, v in cm.DMED.items()})
    with Pool(24) as pool:
        # step 1
        jobs = [(f, float(v), 1, 120000, "isi=1") for f, (_, _, g) in FORM.items() for v in g]
        rows = []
        for (f, v, seed, _, _), lst in pool.map(job, jobs):
            seg = cm.seg_medians(cm.segmented(lst)); x = lst.ISI.values; rows.append(dict(form=f, amp=v, n=len(lst), CV_unseg=x.std(ddof=1) / x.mean(), **{f"{s}_seg": seg[s] for s in cm.STATS}))
        G = pd.DataFrame(rows); G.to_csv(os.path.join(OUT, "a3_grid.csv"), index=False); print(G.round(4).to_string(index=False))
        star = {}
        for f in FORM:
            g = G[G.form == f].sort_values("CV_seg"); star[f] = float(np.exp(np.interp(target, g.CV_seg.values, np.log(g.amp.values))))
        print("matched amplitudes:", star)
        # step 2
        jobs = [(f, star[f], seed, 120000, "isi=1") for f in FORM for seed in range(1, 7)]
        segs = {f: [] for f in FORM}
        for (f, v, seed, _, _), lst in pool.map(job, jobs):
            d = cm.segmented(lst, seed=seed); d["seed"] = seed; segs[f].append(d)
        out = []
        for f in FORM:
            S = pd.concat(segs[f]); rec = dict(form=f, amp=star[f])
            for s in cm.STATS:
                lo, hi = cm.boot_median_ci(S[s].values); rec[s] = S[s].median(); rec[s + "_lo"] = lo; rec[s + "_hi"] = hi; rec[s + "_seed_sd"] = S.groupby("seed")[s].median().std(ddof=1)
            rec["skew_over_CV_seg_median"] = float((S["skew"] / S["CV"]).median()); out.append(rec)
        out.append(dict(form="data", amp=np.nan, **cm.DMED, skew_over_CV_seg_median=float((cm.SEGT.skew_detr / cm.SEGT.CV_detr).median())))
        M = pd.DataFrame(out); M.to_csv(os.path.join(OUT, "a3_matched_segmented.csv"), index=False)
        pd.set_option("display.width", 250); print(M[["form", "amp"] + cm.STATS + [s + "_seed_sd" for s in cm.STATS] + ["skew_over_CV_seg_median"]].round(4).to_string(index=False))
        # step 3
        jobs = [(f, star[f], seed, 2000000, "isi=1") for f in FORM for seed in range(11, 15)]
        rows = []
        for (f, v, seed, _, _), lst in pool.map(job, jobs): rows.append(dict(form=f, amp=v, seed=seed, **long_stats(lst)))
        L = pd.DataFrame(rows); L.to_csv(os.path.join(OUT, "a3_matched_long_seeds.csv"), index=False)
        A = L.groupby("form").mean(numeric_only=True).drop(columns="seed"); SD = L.groupby("form").std(numeric_only=True)
        A.to_csv(os.path.join(OUT, "a3_matched_long.csv")); print(A.round(4).T.to_string()); print("seed sd of R:", SD.R.round(3).to_dict(), " of skew:", SD["skew"].round(3).to_dict())
