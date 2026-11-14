"""
c3_prop41_check.py — the assumptions of Proposition 4.1 in the full deterministic model.

For pairs of stimulus levels (S1 < S2): period ratio T(S1)/T(S2) against the bound f_P(S2)/f_P(S1) of the proposition,
and the relative change of the quantities that the proposition takes as independent of the stimulus (reset and onset
values of H, recovery integral T0, offset T_sp, half-width, amplitude).

Reads results/c2_dose_axis.csv (c2_reduction_along_sweeps.py) and results/f1_curves.csv (f1_dynamic_range.py, the
curve of the reference coupling n = 1, d_H = 10^4, g_H = 30); writes results/c3_prop41_check.csv.
"""
import os, numpy as np, pandas as pd

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
PAIRS = [(0.064, 0.64), (0.05, 30.0), (0.11, 1.0), (1.0, 30.0)]   # mapped end points, whole sweep, two sub-ranges

c = pd.read_csv(os.path.join(OUT, "c2_dose_axis.csv")).set_index("S")
f = pd.read_csv(os.path.join(OUT, "f1_curves.csv")); f = f[(f.nP == 1) & (f.dH == 10000) & (f.gH == 30)].sort_values("S")
wave = lambda S, col: float(np.interp(np.log(S), np.log(f.S), f[col]))

rows = []
for S1, S2 in PAIRS:
    a, b = c.loc[S1], c.loc[S2]
    rows.append(dict(S1=S1, S2=S2, T_ratio=a.Tdet / b.Tdet, fP_ratio=b.fP / a.fP, kappa_ratio=b.kappa / a.kappa,
                     bound_holds=bool(a.Tdet / b.Tdet <= b.fP / a.fP),
                     Hoff_change=b.H_off / a.H_off - 1, Hon_change=b.H_on / a.H_on - 1, T0_change=b.T0 / a.T0 - 1,
                     Tsp_ratio=a.T_sp / b.T_sp, W50_change=wave(S2, "W50") / wave(S1, "W50") - 1,
                     amp_change=wave(S2, "amp") / wave(S1, "amp") - 1))
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "c3_prop41_check.csv"), index=False)
pd.set_option("display.width", 200); print(R.round(4).to_string(index=False))
