"""Check 6: (a) bootstrap CI of the 3/30 uM train-median ratio (train-level and record-level resampling);
(b) slope of log sd vs log Tav over segments with segment-, ROI- and record-level cluster bootstrap."""
from common import *
rng = np.random.default_rng(1); NB = 5000
T = pd.read_csv(os.path.join(CODE, "data_dose_units", "results", "train_table.csv")); T = T[T.n >= 10]
t3, t30 = T[T.dose == 3], T[T.dose == 30]
def ratio_boot(level):
    out = []
    for _ in range(NB):
        parts = []
        for g in (t3, t30):
            if level == "train": s = g.isi_median.values[rng.integers(0, len(g), len(g))]
            else:
                recs = g.record.unique(); pick = recs[rng.integers(0, len(recs), len(recs))]
                s = np.concatenate([g[g.record == r].isi_median.values for r in pick])
            parts.append(np.median(s))
        out.append(parts[0] / parts[1])
    out = np.array(out); return np.median(out), np.quantile(out, .025), np.quantile(out, .975)
res = []
for lev in ("train", "record"):
    m, lo, hi = ratio_boot(lev); res.append(dict(quantity="ratio 3/30 uM (train medians)", level=lev, point=np.median(t3.isi_median) / np.median(t30.isi_median), boot_median=m, ci_lo=lo, ci_hi=hi, n3=len(t3) if lev == "train" else t3.record.nunique(), n30=len(t30) if lev == "train" else t30.record.nunique()))
S = pd.read_csv(os.path.join(CODE, "data_dose_units", "results", "segment_table.csv")); S = S[(S.n >= 10) & (S.sd_detr > 0)].copy()
S["ROI"] = S.record.astype(str) + ":" + S.ST.astype(str)
x, y = np.log(S.Tav.values), np.log(S.sd_detr.values); slope0 = np.polyfit(x, y, 1)[0]
def slope_boot(level):
    out = []
    for _ in range(NB):
        if level == "segment": idx = rng.integers(0, len(S), len(S))
        else:
            keys = S[level].unique(); pick = keys[rng.integers(0, len(keys), len(keys))]
            idx = np.concatenate([np.flatnonzero(S[level].values == k) for k in pick])
        out.append(np.polyfit(x[idx], y[idx], 1)[0])
    out = np.array(out); return np.quantile(out, .025), np.quantile(out, .975)
for lev in ("segment", "ROI", "record"):
    lo, hi = slope_boot(lev); res.append(dict(quantity="slope log sd vs log Tav (segments, detrended)", level=lev, point=slope0, boot_median=np.nan, ci_lo=lo, ci_hi=hi, n3=S[lev].nunique() if lev != "segment" else len(S), n30=np.nan))
R = pd.DataFrame(res); R.to_csv(os.path.join(OUT, "check6_bootstrap.csv"), index=False); print(R.round(3).to_string(index=False)); print("check6 done")
