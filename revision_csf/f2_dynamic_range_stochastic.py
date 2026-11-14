"""
f2_dynamic_range_stochastic.py — interval statistics by dose for the manuscript's rate coupling (n = 1, d_H = 10^4) and the
cooperative one of f1 (n = 2, d_H = 2000); dose maps.

Dose map: S = c U (linear) with c and the time unit T0 either from the two end points 3 and 30 µM (manuscript) or by
least squares over the five train-level medians; saturating map S = c U/(K + U) as a third option.
Stochastic part (H-noise, Jacobi form): one sigma_H for all doses, chosen so that the segmented, detrended CV at the
15 µM stimulus equals the data median of the 15 µM segments; then CV, skewness, amplitude CV by dose (same pipeline)
and the log-log slope of sd against mean interval across the five doses.
Output: results/f2_dose_map.csv, results/f2_by_dose.csv
"""
import os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "round2_checks"))
import common as cm
OUT = os.path.join(HERE, "results"); BIN = os.path.join(HERE, "m1_noiseform_sde")
VAR = {"n=1, dH=1e4 (manuscript)": dict(nP=1, dH=10000), "n=2, dH=2000": dict(nP=2, dH=2000)}
DOSE = np.array([3, 5, 10, 15, 30.0])
DT = pd.read_csv(os.path.join(HERE, "..", "data_dose_units", "results", "dose_table.csv"))
def data_T():
    out = []
    for d in DOSE:
        g = DT[DT.dose == d]; sub = g[g.subset.isin(["constant+pre", "all", "constant", "step"])]
        out.append(float(g.isi_train_median.iloc[0]))
    return np.array(out)
def run_pt(a):
    v, S, sig, seed = a; T = float(cm.run(BIN, {**cm.REF, **VAR[v], "S": S, "sigH": 1e-9, "Tend": 40000, "Ttrans": 20000}, None).iloc[0].Tdet)
    lst = cm.run(BIN, {**cm.REF, **VAR[v], "S": S, "sigH": sig, "seed": seed, "Tend": 3000 + 1450 * T}, "isi=1"); st = cm.run(BIN, {**cm.REF, **VAR[v], "S": S, "sigH": sig, "seed": seed, "Tend": 3000 + 1450 * T}, None).iloc[0]
    seg = cm.segmented(lst, seed=seed); return a, T, st, cm.seg_medians(seg)
if __name__ == "__main__":
    pd.set_option("display.width", 230); C = pd.read_csv(os.path.join(OUT, "f1_curves.csv")); Td = data_T(); print("data train-level median ISI by dose (s):", dict(zip(DOSE, Td.round(1))))
    Sg = cm.SEGT; dstat = Sg.groupby("dose").agg(CV=("CV_detr", "median"), skew=("skew_detr", "median"), amp_cv=("amp_cv", "median"), Tav=("Tav", "median"), sd=("sd_detr", "median"), nseg=("n", "size"))
    maps, rows = {}, []
    for v, par in VAR.items():
        c = C[(C.nP == par["nP"]) & (C.dH == par["dH"]) & (C.gH == 30) & (C.Tdet > 0)].sort_values("S"); lT = lambda S: np.interp(np.log(S), np.log(c.S.values), np.log(c.Tdet.values))
        def rms(p, kind):
            S = np.exp(p[0]) * DOSE if kind == "lin" else np.exp(p[0]) * DOSE / (np.exp(p[2]) + DOSE)
            if S.min() < c.S.min() or S.max() > c.S.max(): return 10.0
            return float(np.sqrt(np.mean((lT(S) + p[1] - np.log(Td)) ** 2)))
        # end points
        best = None
        for cc in np.geomspace(c.S.min() / 3, c.S.max() / 30, 400):
            r = (lT(cc * 3) - lT(cc * 30)); 
            if best is None or r > best[0]: best = (r, cc)
        cc = best[1]; t0 = np.exp(0.5 * ((np.log(Td[0]) - lT(cc * 3)) + (np.log(Td[-1]) - lT(cc * 30)))); pred = np.exp(lT(cc * DOSE)) * t0
        maps[v] = cc * DOSE; rows.append(dict(variant=v, map="end points, linear", c=cc, K=np.inf, T0_s=t0, log_rms=float(np.sqrt(np.mean(np.log(pred / Td) ** 2))), ratio_3_30=pred[0] / pred[-1], **{f"ISI_{int(d)}": p for d, p in zip(DOSE, pred)}))
        for kind, x0s in (("lin", [[np.log(cc), np.log(t0)]]), ("sat", [[np.log(cc * 100), np.log(t0), np.log(100.0)], [np.log(cc * 20), np.log(t0), np.log(20.0)], [np.log(cc * 1000), np.log(t0), np.log(1000.0)]])):
            r = min((minimize(rms, x0, args=(kind,), method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9, maxiter=4000)) for x0 in x0s), key=lambda q: q.fun)
            S = np.exp(r.x[0]) * DOSE if kind == "lin" else np.exp(r.x[0]) * DOSE / (np.exp(r.x[2]) + DOSE); pred = np.exp(lT(S) + r.x[1])
            rows.append(dict(variant=v, map="least squares, " + ("linear" if kind == "lin" else "saturating"), c=np.exp(r.x[0]), K=np.inf if kind == "lin" else np.exp(r.x[2]), T0_s=np.exp(r.x[1]), log_rms=r.fun, ratio_3_30=pred[0] / pred[-1], **{f"ISI_{int(d)}": p for d, p in zip(DOSE, pred)}))
    M = pd.DataFrame(rows); M.to_csv(os.path.join(OUT, "f2_dose_map.csv"), index=False); print(M.round(3).to_string(index=False))
    with Pool(20) as pool:
        out = []
        for v in VAR:
            S15 = maps[v][3]; grid = np.geomspace(0.012, 0.05, 8); res = pool.map(run_pt, [(v, S15, float(s), 1) for s in grid]); cv = np.array([r[3]["CV"] for r in res])
            sig = float(np.exp(np.interp(dstat.loc[15, "CV"], cv, np.log(grid)))); res = pool.map(run_pt, [(v, float(S), sig, seed) for S in maps[v] for seed in (1, 2, 3, 4)])
            for d, S in zip(DOSE, maps[v]):
                rr = [r for r in res if r[0][1] == float(S)]; med = pd.DataFrame([r[3] for r in rr]).mean(); st = pd.DataFrame([r[2] for r in rr]).mean()
                out.append(dict(variant=v, dose=d, S=S, sigma_H=sig, Tdet=rr[0][1], Tav=st.Tav, sd=st.sd, CV_unseg=st.CV, amp=st.amp_mean, CV=med.CV, skew=med["skew"], amp_cv=med.amp_cv, corr=med.corr_isi_amp, serial=med.serial,
                                data_CV=dstat.loc[d, "CV"], data_skew=dstat.loc[d, "skew"], data_amp_cv=dstat.loc[d, "amp_cv"], data_Tav=dstat.loc[d, "Tav"], nseg=dstat.loc[d, "nseg"]))
    B = pd.DataFrame(out); B.to_csv(os.path.join(OUT, "f2_by_dose.csv"), index=False); print(B.round(3).to_string(index=False))
    for v in VAR:
        g = B[B.variant == v]; print(v, ": dose-axis log-log slope of sd vs T_av (model) %.2f; CV ratio 3/30 µM model %.2f, data %.2f; amplitude ratio 30/3 µM model %.2f" % (np.polyfit(np.log(g.Tav), np.log(g.sd), 1)[0], g.CV.iloc[0] / g.CV.iloc[-1], g.data_CV.iloc[0] / g.data_CV.iloc[-1], g.amp.iloc[-1] / g.amp.iloc[0]))
    print("data: log-log slope of the dose medians (sd vs T_av): %.2f" % np.polyfit(np.log(dstat.Tav), np.log(dstat.sd), 1)[0])
