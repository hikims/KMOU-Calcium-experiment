"""
g1_noise_location_sensitivity.py — the noise-location comparison away from the single reference point, with seed
variability.

For each reference point the three noise amplitudes (H: Jacobi, calcium, IP3) are re-tuned so that the segmented,
detrended CV equals the data median, then 8 seeds give the five statistics of Table 5: median over segments, the sd of
that median across seeds, and a bootstrap interval over segments.  Reference points: the manuscript's (S = 1,
eps_H = 1000), two mapped doses (S = 0.32 ~ 15 µM, S = 0.11 ~ 5 µM), two recovery time scales (eps_H = 500, 2000) and
the cooperative rate coupling of f1 (n = 2, d_H = 2000) at S = 0.32.
Output: results/g1_grid.csv, results/g1_table.csv
"""
import os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "round2_checks"))
import common as cm
OUT = os.path.join(HERE, "results"); BIN = os.path.join(HERE, "m1_noiseform_sde")
POINTS = {"REF (S=1, epsH=1000)": {}, "S=0.32": dict(S=0.32), "S=0.11": dict(S=0.11), "epsH=500": dict(epsH=500), "epsH=2000": dict(epsH=2000), "n=2, dH=2000, S=0.32": dict(nP=2, dH=2000, S=0.32)}
KEY = {"H": "sigH", "C": "sigC", "P": "sigP"}
def tdet(extra): return float(cm.run(BIN, {**cm.REF, **extra, "sigH": 1e-9, "Tend": 40000, "Ttrans": 20000}, None).iloc[0].Tdet)
def job(a):
    pt, k, v, seed, T = a
    lst = cm.run(BIN, {**cm.REF, **POINTS[pt], KEY[k]: v, "seed": seed, "Tend": 3000 + 1450 * T}, "isi=1")
    if len(lst) < 200: return a, None, len(lst)
    return a, cm.segmented(lst, seed=seed), len(lst)
if __name__ == "__main__":
    target = cm.DMED["CV"]; pd.set_option("display.width", 250)
    with Pool(14) as pool:
        TD = dict(zip(POINTS, pool.map(tdet, list(POINTS.values())))); print("deterministic periods:", {k: round(v, 1) for k, v in TD.items()}, flush=True)
        # scale the starting grids with the recovery time: sigma_H ∝ T^-1/2 keeps nu, the others are searched widely
        grids = {pt: {"H": np.geomspace(0.012, 0.06, 8) * np.sqrt(83.76 / TD[pt]), "C": np.geomspace(0.12, 1.5, 9), "P": np.geomspace(0.006, 0.3, 10)} for pt in POINTS}
        jobs = [(pt, k, float(v), 1, TD[pt]) for pt in POINTS for k in KEY for v in grids[pt][k]]; rows = []
        for (pt, k, v, seed, T), seg, n in pool.map(job, jobs, chunksize=1):
            med = cm.seg_medians(seg) if seg is not None else {s: np.nan for s in cm.STATS}; rows.append(dict(point=pt, noise=k, sigma=v, n=n, **med))
        G = pd.DataFrame(rows); G.to_csv(os.path.join(OUT, "g1_grid.csv"), index=False)
        star, note = {}, {}
        for pt in POINTS:
            for k in KEY:
                g = G[(G.point == pt) & (G.noise == k) & G.CV.notna()].sort_values("sigma"); cv, sg = g.CV.values, g.sigma.values
                if len(cv) == 0 or np.nanmax(cv) < target: star[pt, k] = np.nan; note[pt, k] = f"segmented CV never reaches {target:.3f} (max {np.nanmax(cv) if len(cv) else float('nan'):.3f})"; continue
                i = int(np.argmax(cv >= target)); star[pt, k] = sg[0] if i == 0 else float(np.exp(np.interp(target, cv[i - 1:i + 1], np.log(sg[i - 1:i + 1])))); note[pt, k] = ""
        jobs = [(pt, k, star[pt, k], seed, TD[pt]) for (pt, k) in star if np.isfinite(star[pt, k]) for seed in range(1, 9)]; segs = {}
        for (pt, k, v, seed, T), seg, n in pool.map(job, jobs, chunksize=1):
            if seg is not None: seg["seed"] = seed; segs.setdefault((pt, k), []).append(seg)
    out = []
    for (pt, k) in star:
        rec = dict(point=pt, noise=k, sigma=star[pt, k], note=note[pt, k])
        if (pt, k) in segs:
            S = pd.concat(segs[pt, k])
            for s in cm.STATS:
                lo, hi = cm.boot_median_ci(S[s].values); rec[s] = S[s].median(); rec[s + "_seed_sd"] = S.groupby("seed")[s].median().std(ddof=1); rec[s + "_lo"] = lo; rec[s + "_hi"] = hi
        out.append(rec)
    T = pd.DataFrame(out); T.to_csv(os.path.join(OUT, "g1_table.csv"), index=False)
    show = ["point", "noise", "sigma"] + [c for s in cm.STATS for c in (s, s + "_seed_sd")]
    print(T[show].round(3).to_string(index=False)); print("data:", {k: round(v, 3) for k, v in cm.DMED.items()}); print(T[T.note != ""][["point", "noise", "note"]].to_string(index=False))
