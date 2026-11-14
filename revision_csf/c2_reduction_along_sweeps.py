"""
c2_reduction_along_sweeps.py — reduction constants along the sweeps of the variability law (how T_sp
varies along the sweeps; range of nu along the dose axis; reduced vs full moments along that sweep).

For each point: deterministic cycle -> C_rest (minimum of C), P_bar (cycle mean of P), H_off (minimum of H), H_on (H at
onset), kappa = f_C f_P Q at (C_rest, P_bar), eps_eff = eps_H/kappa, T0 = ln[(1-H_off)/(1-H_on)], T_sp = T_det - eps_eff T0,
nu = sigma_H sqrt(eps_eff); reduced mean/sd from the exact Jacobi first-passage cumulants with that geometry; full-model
values from new runs (Tend = 4e5, two seeds) with the manuscript's sweep amplitude sigma_H = 0.026.
Output: results/c2_dose_axis.csv, results/c2_cell_axis.csv
"""
import os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "round2_checks"))
import common as cm, a2_exact_vs_expansion as A2
OUT = os.path.join(HERE, "results"); BIN = os.path.join(HERE, "m1_noiseform_sde")
P0 = dict(cm.REF); SIG = 0.026
def point(a):
    S, epsH, sig, nu0 = a; par = {**P0, "S": S, "epsH": epsH}
    if nu0: sig = nu0 / np.sqrt(epsH)
    st = cm.run(BIN, {**par, "sigH": 1e-9, "Tend": 20000 + 6 * 400 * epsH / 1000}, None).iloc[0]
    if st.Tdet == 0: return None
    tr = cm.run(BIN, {**par, "sigH": 1e-9, "tspan": 3 * st.Tdet}, "trace=1"); lst = cm.run(BIN, {**par, "sigH": 1e-9, "Tend": 3000 + 12 * st.Tdet}, "isi=1")
    Cr, Pb, Hoff, Hon = tr.C.min(), tr.P.mean(), lst.Hoff_prev.iloc[-1], st.Hon_det
    fP = 1 + P0["gH"] * Pb / (P0["dH"] + Pb); fC = 1 + P0["rC"] * Cr**2 / (P0["kR"]**2 + Cr**2); Q = 0.15**2 / (0.15**2 + Cr**2); kap = fC * fP * Q
    eeff = epsH / kap; T0 = np.log((1 - Hoff) / (1 - Hon)); Tsp = st.Tdet - eeff * T0; nu = sig * np.sqrt(eeff)
    A2.FORMS["_tmp"] = A2.FORMS["jacobi"]; e = A2.exact("jacobi", nu, h0=Hoff, hs=Hon)
    full = pd.concat([cm.run(BIN, {**par, "sigH": sig, "seed": s, "Tend": max(400000, 600 * st.Tdet)}, None) for s in (1, 2)]).mean()
    return dict(S=S, epsH=epsH, sigH=sig, Tdet=st.Tdet, C_rest=Cr, P_bar=Pb, fP=fP, fC=fC, Q=Q, kappa=kap, eps_eff=eeff, H_off=Hoff, H_on=Hon, T0=T0, T_sp=Tsp, Tsp_over_Tdet=Tsp / st.Tdet, nu=nu,
                Tav_red=Tsp + eeff * e["mean"], sd_red=eeff * np.sqrt(e["var"]), Tav_full=full.Tav, sd_full=full.sd, CV_full=full.CV, n_full=full.n,
                err_mean=(Tsp + eeff * e["mean"]) / full.Tav - 1, err_sd=eeff * np.sqrt(e["var"]) / full.sd - 1, local_slope_nu_common=st.Tdet / (st.Tdet - Tsp))
if __name__ == "__main__":
    pd.set_option("display.width", 250)
    with Pool(24) as pool:
        Sv = sorted(set(np.round(np.geomspace(0.05, 30, 9), 6)) | {0.064, 0.11, 0.21, 0.32, 0.64, 1.0})
        D = pd.DataFrame([r for r in pool.map(point, [(S, 1000.0, SIG, 0) for S in Sv]) if r]); D.to_csv(os.path.join(OUT, "c2_dose_axis.csv"), index=False)
        Cc = pd.DataFrame([r for r in pool.map(point, [(1.0, e, 0, 0.822) for e in np.geomspace(200, 8000, 9)]) if r]); Cc.to_csv(os.path.join(OUT, "c2_cell_axis.csv"), index=False)
    cols = ["S", "epsH", "Tdet", "kappa", "eps_eff", "H_off", "H_on", "T0", "T_sp", "Tsp_over_Tdet", "nu", "Tav_full", "CV_full", "err_mean", "err_sd", "local_slope_nu_common"]
    print("dose axis (eps_H = 1000, sigma_H = 0.026)"); print(D[cols].round(4).to_string(index=False))
    print("\ncell axis (S = 1, common nu0 = 0.822)"); print(Cc[cols].round(4).to_string(index=False))
    for nm, T in (("dose", D), ("cell", Cc)):
        print(nm, "axis: log-log slope full %.3f, reduced %.3f" % (np.polyfit(np.log(T.Tav_full), np.log(T.sd_full), 1)[0], np.polyfit(np.log(T.Tav_red), np.log(T.sd_red), 1)[0]))
