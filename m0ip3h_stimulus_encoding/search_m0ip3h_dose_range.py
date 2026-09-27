"""
search_m0ip3h_dose_range.py — grid search: which (epsH, gH, dH, rho) of the IP3->H *rate* coupling
(mech=2) give a stimulus-encoded period that (a) oscillates over the whole stimulus range,
(b) decreases with stimulus, and (c) spans a ratio comparable to the data (ISI 206 s -> 31 s, 6.6x)?

Runs ./m0ip3h_stimulus_sweep for each config with an S-sweep (log-spaced), and writes
    results/m0ip3h_dose_range_search.csv      one row per config
    results/m0ip3h_dose_range_best_curves.csv  T(S) curves of the top configs
"""
import subprocess, itertools, io, os, sys
import numpy as np, pandas as pd

BIN = "./m0ip3h_stimulus_sweep"
OUT = "results"
S_LO, S_HI, S_N = 0.1, 30.0, 9          # stimulus range (nondim), log-spaced
GRID = dict(
    epsH=[70, 150, 300, 600],
    gH=[3, 10, 30, 100],
    dH=[1000, 3000, 10000, 30000],
    rho=[2.4, 0.24, 0.024],
)
FIXED = dict(mech=2, lambda_=0.001)      # lambda_ -> lambda (python keyword)

def run(cfg):
    args = [BIN, f"sweep=S", f"lo={S_LO}", f"hi={S_HI}", f"n={S_N}", "log=1"]
    for k, v in cfg.items():
        args.append(f"{k.rstrip('_')}={v}")
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    return pd.read_csv(io.StringIO(out))

rows, curves = [], []
combos = list(itertools.product(*GRID.values()))
print(f"{len(combos)} configs x {S_N} stimulus values")
for i, vals in enumerate(combos):
    cfg = dict(zip(GRID.keys(), vals)); cfg.update(FIXED)
    d = run(cfg)
    osc = d[d.osc == 1]
    n_osc = len(osc)
    if n_osc >= 2:
        T = osc.Tosc.values; S = osc.value.values
        ratio = T.max() / T.min()
        decreasing = bool(np.all(np.diff(T) < 0))
        # log-log slope of T vs S (sensitivity)
        slope = np.polyfit(np.log(S), np.log(T), 1)[0]
        duty = float(osc.duty.mean()); amp = float((osc.Cmax / osc.Cmin).mean())
    else:
        ratio, decreasing, slope, duty, amp = np.nan, False, np.nan, np.nan, np.nan
    rows.append(dict(**{k: v for k, v in cfg.items() if k != "lambda_"}, n_osc=n_osc, n_S=S_N,
                     T_min=osc.Tosc.min() if n_osc else np.nan, T_max=osc.Tosc.max() if n_osc else np.nan,
                     T_ratio=ratio, decreasing=decreasing, logslope_T_S=slope, duty_mean=duty, amp_ratio=amp))
    d["cfg"] = i; curves.append(d)
    if (i + 1) % 20 == 0: print(f"  {i+1}/{len(combos)}", file=sys.stderr)

R = pd.DataFrame(rows); R.index.name = "cfg"
os.makedirs(OUT, exist_ok=True)
R.to_csv(os.path.join(OUT, "m0ip3h_dose_range_search.csv"))

# ranking: oscillate everywhere, decreasing, largest ratio
ok = R[(R.n_osc == S_N) & (R.decreasing)].sort_values("T_ratio", ascending=False)
print("\nconfigs oscillating over the full range AND decreasing with S:", len(ok), "of", len(R))
print(ok.head(12).to_string())
best = ok.head(6).index.tolist()
C = pd.concat([c for c in curves if c.cfg.iloc[0] in best])
C.to_csv(os.path.join(OUT, "m0ip3h_dose_range_best_curves.csv"), index=False)
print("\nwrote", os.path.join(OUT, "m0ip3h_dose_range_search.csv"), "and best curves")
