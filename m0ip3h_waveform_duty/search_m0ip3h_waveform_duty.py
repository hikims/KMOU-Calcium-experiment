"""
search_m0ip3h_waveform_duty.py — random (Latin-hypercube) search for a narrow-spike regime.

Model: AIMS reaction terms + IP3->H rate coupling (from m0ip3h_stimulus_encoding, fixed) +
       three waveform ingredients, each switched off at its default:
         nA, kA  cooperative Ca activation  C^nA/(kA^nA + C^nA)   (nA=1,kA=1 original)
         Jin     basal Ca influx                                     (0 original)
         rC      Ca-dependent speed-up of H kinetics                 (0 original)
Question: is there a region where the spike is narrow (W50/T <= 0.15, data 0.05-0.15) while the period
still decreases >= 4x from S̄=0.1 to S̄=30?  (second pass: box narrowed to the relaxation regime found by hand)

Each sample is run at S̄ = S_LO and S̄ = S_HI.  Writes
    results/m0ip3h_waveform_duty_search.csv   one row per sample (both stimuli)
Usage:  python3 search_m0ip3h_waveform_duty.py [N_SAMPLES] [SEED]   (uses all CPU cores; ~3 min for 600 on 4 cores)
"""
import subprocess, io, os, sys
from multiprocessing import Pool, cpu_count
import numpy as np, pandas as pd

BIN, OUT = "./m0ip3h_waveform_sweep", "results"
N = int(sys.argv[1]) if len(sys.argv) > 1 else 600
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1
S_LO, S_HI = 0.1, 30.0
FIXED = dict(mech=2, rho=0.024, Tend=15000, Ttrans=9000, dt=0.02)
# log-uniform ranges (nA discrete)
RANGES = dict(kA=(0.2, 1.0), Jin=(0.005, 0.05), rC=(5.0, 60.0), beta=(1.0, 8.0), kS=(0.1, 0.5),
              alpha=(5.0, 40.0), epsH=(300.0, 5000.0), gH=(10.0, 300.0), dH=(3000.0, 30000.0))
NA_CHOICES = [2, 3]          # nA=1 (original, non-cooperative) never gave duty < 0.35 in the first pass

rng = np.random.default_rng(SEED)
def lhs(n, k):
    u = (np.argsort(rng.random((n, k)), axis=0) + rng.random((n, k))) / n
    return u
U = lhs(N, len(RANGES))
samples = []
for i in range(N):
    cfg = {k: float(lo * (hi / lo) ** U[i, j]) for j, (k, (lo, hi)) in enumerate(RANGES.items())}
    cfg["nA"] = int(rng.choice(NA_CHOICES)); samples.append(cfg)

def run(cfg, S):
    args = [BIN, f"S={S}"] + [f"{k}={v}" for k, v in {**FIXED, **cfg}.items()]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    return pd.read_csv(io.StringIO(out)).iloc[0]

def evaluate(cfg):
    lo, hi = run(cfg, S_LO), run(cfg, S_HI)
    r = dict(cfg)
    for tag, f in (("lo", lo), ("hi", hi)):
        for col in ("osc", "Tosc", "Cmax", "Cmin", "W50", "duty"): r[f"{col}_{tag}"] = f[col]
    both = bool(lo.osc == 1 and hi.osc == 1)
    r["both_osc"] = both
    r["T_ratio"] = lo.Tosc / hi.Tosc if both else np.nan
    r["duty_max"] = max(lo.duty, hi.duty) if both else np.nan
    r["amp_lo"] = lo.Cmax / max(lo.Cmin, 1e-9) if both else np.nan
    return r

if __name__ == "__main__":
    with Pool(cpu_count()) as pool:
        rows = pool.map(evaluate, samples, chunksize=4)

    R = pd.DataFrame(rows); R.index.name = "sample"
    os.makedirs(OUT, exist_ok=True); R.to_csv(os.path.join(OUT, "m0ip3h_waveform_duty_search.csv"))
    ok = R[R.both_osc]
    good = ok[(ok.duty_max <= 0.15) & (ok.T_ratio >= 4)]
    print(f"\n{len(ok)} of {N} oscillate at both stimuli; {len(good)} satisfy duty<=0.15 & ratio>=4")
    cols = ["nA", "kA", "Jin", "rC", "beta", "kS", "alpha", "epsH", "gH", "dH", "Tosc_lo", "Tosc_hi", "T_ratio", "duty_lo", "duty_hi", "amp_lo"]
    print("\nlowest duty (both oscillating):"); print(ok.sort_values("duty_max").head(8)[cols].round(3).to_string())
    if len(good):
        print("\nbest by T_ratio among good:"); print(good.sort_values("T_ratio", ascending=False).head(8)[cols].round(3).to_string())
    print("wrote", os.path.join(OUT, "m0ip3h_waveform_duty_search.csv"))
