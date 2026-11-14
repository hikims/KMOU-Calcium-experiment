"""
j1_linear_law.py — slope of the linear variability law sigma = alpha (T_av - offset), for the comparison with the model
of Azzoni et al. (Table 10 of the manuscript).

Data: median over the 196 segments of sd_detr/(T_av - 8.64 s), in total and in the three ranges of T_av used by Azzoni
et al. (0-70, 70-160, 160-400 s).  Reduced process: nu sqrt(I1)/T0 at the fitted nu = 0.281 and the reference
geometry.  Full model: sd/(T_av - T_sp) at the CV-matched amplitude (results/a3_matched_long.csv, Jacobi form), with
T_sp = 16.4 model units.  Writes results/j1_linear_law.csv.
"""
import os, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "results")
S = pd.read_csv(os.path.join(HERE, "..", "data_dose_units", "results", "segment_table.csv")); S = S[S.n >= 10]
a = S.sd_detr / (S.Tav - 8.64); rows = [dict(item="data, all segments", n=len(S), alpha=a.median())]
for lo, hi in ((0, 70), (70, 160), (160, 400)):
    k = (S.Tav >= lo) & (S.Tav < hi); rows.append(dict(item=f"data, T_av {lo}-{hi} s", n=int(k.sum()), alpha=a[k].median()))
H0, Hs, nu = 0.064, 0.408, 0.281; z0, zs = 1 / (1 - H0), 1 / (1 - Hs); T0 = np.log(zs / z0); I1 = zs - z0 - T0
rows.append(dict(item="reduced process, nu = 0.281", n=np.nan, alpha=nu * np.sqrt(I1) / T0))
m = pd.read_csv(os.path.join(OUT, "a3_matched_long.csv")).set_index("form").loc["jacobi"]
rows.append(dict(item="full model, CV-matched amplitude", n=m.n, alpha=m.CV * m.Tav / (m.Tav - 16.4)))
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "j1_linear_law.csv"), index=False); print(R.round(4).to_string(index=False))
