import io, os, subprocess, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.abspath(__file__)); CODE = os.path.dirname(ROOT)
OUT = os.path.join(ROOT, "results"); os.makedirs(OUT, exist_ok=True)
REF = dict(mech=2, rho=0.024, gH=30, dH=10000, epsH=1000, nA=2, kA=0.4, Jin=0.02, kS=0.25, beta=4, alpha=20, qoff=0, rC=300, kR=1, S=1)
NUM = dict(Tend=120000, Ttrans=3000, dt=0.01)
BIN_ISI = os.path.join(CODE, "m1_isi_variability", "m1_isi_sde")
def run(binary, extra, mode=None, base=None):
    kv = {**(base if base is not None else {**REF, **NUM}), **extra}
    args = [binary] + ([mode] if mode else []) + [f"{k}={v}" for k, v in kv.items()]
    p = subprocess.run(args, capture_output=True, text=True)
    if p.returncode: raise SystemExit(p.stderr)
    return pd.read_csv(io.StringIO(p.stdout))
# --- data segments (same as plot_m1_isi_variability.py) ---
SEGT = pd.read_csv(os.path.join(CODE, "data_dose_units", "results", "segment_table.csv")); SEGT = SEGT[SEGT.n >= 10].copy()
SEG_N = SEGT.n.values
STATS = ["CV", "skew", "amp_cv", "corr_isi_amp", "serial"]
DMED = dict(CV=SEGT.CV_detr.median(), skew=SEGT.skew_detr.median(), amp_cv=SEGT.amp_cv.median(), corr_isi_amp=SEGT.corr_detr.median(), serial=SEGT.serial_detr.median())
def detrended(x, t):
    if len(x) >= 4: b = np.polyfit(t, x, 1); return x - np.polyval(b, t) + x.mean()
    return x
def unit_stats(x, a, t):
    xd = detrended(x, t); m = xd.mean(); r = xd - m
    return dict(CV=xd.std(ddof=1) / m, skew=float(np.mean((r / r.std()) ** 3)) if r.std() > 0 else np.nan, amp_cv=a.std(ddof=1) / a.mean(),
                corr_isi_amp=np.corrcoef(xd, a)[0, 1], serial=np.corrcoef(r[:-1], r[1:])[0, 1])
def segmented(lst, nrep=3, seed=0):
    rng = np.random.default_rng(seed); x, a, t = lst.ISI.values, lst.amp.values, lst.t.values; rows = []
    for _ in range(nrep):
        i = 0
        for n_ in rng.permutation(SEG_N):
            if i + n_ > len(x): break
            rows.append(unit_stats(x[i:i + n_], a[i:i + n_], t[i:i + n_])); i += n_
    return pd.DataFrame(rows)
def seg_medians(df): return {s: float(df[s].median()) for s in STATS}
def boot_median_ci(v, nb=2000, seed=0):
    v = np.asarray(v, float); v = v[~np.isnan(v)]; rng = np.random.default_rng(seed)
    bs = np.array([np.median(v[rng.integers(0, len(v), len(v))]) for _ in range(nb)]); return float(np.quantile(bs, .025)), float(np.quantile(bs, .975))
