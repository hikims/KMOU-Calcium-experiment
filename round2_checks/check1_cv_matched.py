"""Check 1: tune each noise amplitude so that the *segmented, detrended* CV equals the data median (0.196), then
compare the other statistics at that amplitude (3 seeds; bootstrap CI over segments)."""
from common import *
KEY = {"H": "sigH", "C": "sigC", "P": "sigP"}
GRID = {"H": np.geomspace(0.012, 0.06, 8), "C": np.geomspace(0.2, 1.5, 8), "P": np.geomspace(0.008, 0.25, 9)}
rows = []
for k, grid in GRID.items():
    for v in grid:
        lst = run(BIN_ISI, {KEY[k]: v}, "isi=1"); st = run(BIN_ISI, {KEY[k]: v}, None).iloc[0]
        seg = seg_medians(segmented(lst)) if len(lst) > 50 else {s: np.nan for s in STATS}
        rows.append(dict(noise=k, sigma=v, n=len(lst), CV_unseg=st.CV, **{f"{s}_seg": seg[s] for s in STATS})); print(rows[-1], flush=True)
G = pd.DataFrame(rows); G.to_csv(os.path.join(OUT, "check1_grid.csv"), index=False)
out = []
for k in GRID:
    g = G[(G.noise == k) & G.CV_seg.notna()].sort_values("sigma"); cv = g.CV_seg.values; sg = g.sigma.values
    target = DMED["CV"]
    if cv.max() < target:
        out.append(dict(noise=k, sigma_matched=np.nan, note=f"segmented CV never reaches {target:.3f}; max {cv.max():.3f} at sigma {sg[cv.argmax()]:.4g}")); continue
    i = int(np.argmax(cv >= target)); s_star = sg[0] if i == 0 else float(np.exp(np.interp(target, cv[i-1:i+1], np.log(sg[i-1:i+1]))))
    segs = []
    for seed in (1, 2, 3):
        lst = run(BIN_ISI, {KEY[k]: s_star, "seed": seed}, "isi=1"); d = segmented(lst, seed=seed); d["seed"] = seed; segs.append(d)
    S = pd.concat(segs); rec = dict(noise=k, sigma_matched=s_star, note="")
    for s in STATS:
        lo, hi = boot_median_ci(S[s].values); per_seed = S.groupby("seed")[s].median()
        rec[f"{s}"] = S[s].median(); rec[f"{s}_ci_lo"] = lo; rec[f"{s}_ci_hi"] = hi; rec[f"{s}_seed_sd"] = per_seed.std(ddof=1); rec[f"{s}_data"] = DMED[s]
    out.append(rec); print(rec, flush=True)
pd.DataFrame(out).to_csv(os.path.join(OUT, "check1_cv_matched_table.csv"), index=False); print("check1 done")
