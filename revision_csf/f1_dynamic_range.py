"""
f1_dynamic_range.py — can a steeper stimulus-to-rate transduction recover the observed frequency range?

Rate coupling f_P = 1 + g_H P^n/(d_H^n + P^n) with Hill exponent n (manuscript: n = 1, d_H = 10^4).  For each (n, d_H):
deterministic stimulus curve (period, half-width, amplitude, onset and reset values of H), the largest period ratio
over one decade of stimulus (data: 6.95 between 3 and 30 µM), and the change of the waveform inside that decade.
Output: results/f1_curves.csv, results/f1_summary.csv
"""
import os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "round2_checks"))
import common as cm
OUT = os.path.join(HERE, "results"); BIN = os.path.join(HERE, "m1_noiseform_sde")
SV = np.geomspace(0.01, 30, 57)        # 16 points per decade
def det(a):
    nP, dH, gH, S = a
    r = cm.run(BIN, {**cm.REF, "nP": nP, "dH": dH, "gH": gH, "S": S, "sigH": 1e-9, "Tend": 60000, "Ttrans": 20000}, None).iloc[0]
    return dict(nP=nP, dH=dH, gH=gH, S=S, Tdet=r.Tdet, W50=r.W50_det, amp=r.amp_det, Hon=r.Hon_det)
if __name__ == "__main__":
    cases = [(1, 10000, 30), (1, 3000, 30), (1, 1000, 30), (2, 3000, 30), (2, 2000, 30), (2, 1500, 30), (2, 1000, 30), (3, 2000, 30), (3, 1500, 30), (2, 1500, 60), (1, 10000, 100)]
    with Pool(26) as pool: C = pd.DataFrame(pool.map(det, [(n, d, g, float(S)) for n, d, g in cases for S in SV], chunksize=4))
    C.to_csv(os.path.join(OUT, "f1_curves.csv"), index=False); rows = []
    for (n, d, g), c in C.groupby(["nP", "dH", "gH"], sort=False):
        c = c.sort_values("S").reset_index(drop=True); best = None
        for i in range(len(c) - 16):
            w = c.iloc[i:i + 17]
            if (w.Tdet <= 0).any(): continue
            r = w.Tdet.iloc[0] / w.Tdet.iloc[-1]
            if best is None or r > best[0]: best = (r, w)
        if best is None: rows.append(dict(nP=n, dH=d, gH=g, ratio_per_decade=np.nan)); continue
        r, w = best
        rows.append(dict(nP=n, dH=d, gH=g, ratio_per_decade=r, S_lo=w.S.iloc[0], S_hi=w.S.iloc[-1], T_lo=w.Tdet.iloc[0], T_hi=w.Tdet.iloc[-1], W50_ratio=w.W50.iloc[-1] / w.W50.iloc[0], amp_ratio=w.amp.iloc[-1] / w.amp.iloc[0],
                         Hon_lo=w.Hon.iloc[0], Hon_hi=w.Hon.iloc[-1], duty_hi=w.W50.iloc[-1] / w.Tdet.iloc[-1], n_osc=int((c.Tdet > 0).sum()), ratio_full_range=c[c.Tdet > 0].Tdet.max() / c[c.Tdet > 0].Tdet.min()))
    Sm = pd.DataFrame(rows); Sm.to_csv(os.path.join(OUT, "f1_summary.csv"), index=False); pd.set_option("display.width", 220); print(Sm.round(3).to_string(index=False))
