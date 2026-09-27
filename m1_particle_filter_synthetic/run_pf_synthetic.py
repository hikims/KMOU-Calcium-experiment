"""
run_pf_synthetic.py — gate 2: can (eps_H, sigma_H, beta) be recovered from a synthetic fluorescence record?

1. simulate a record from the M1 SDE (stereotyped reference, S̄=1, sigma_H=0.026) with observation F = b + a C + e,
   dt_obs = 2 model units (≈ 1 s), T_obs = 2000 (≈ 18 min, ~24 spikes)          -> results/synthetic_truth.csv
2. profile log-likelihood (bootstrap PF, N particles) over each of eps_H, sigma_H, beta (others at truth),
   NSEED independent PF seeds per point                                             -> results/pf_profile_<par>.csv
3. filtered means at the truth                                                      -> results/pf_filtered_truth.csv
Usage: python3 run_pf_synthetic.py [NPART=1000] [NPOINTS=9] [NSEED=2]   (~10 s per likelihood at N=1000; all cores used)
"""
import io, os, subprocess, sys
from multiprocessing import Pool, cpu_count
import numpy as np, pandas as pd
BIN, OUT = "./m1_pf_sde", "results"; os.makedirs(OUT, exist_ok=True)
NPART = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
NPOINTS = int(sys.argv[2]) if len(sys.argv) > 2 else 9
NSEED = int(sys.argv[3]) if len(sys.argv) > 3 else 2
TRUTH = dict(epsH=1000.0, sigH=0.026, beta=4.0, a=0.11, b=0.0, sobs=0.05)
OBS = os.path.join(OUT, "synthetic_truth.csv")
def run(args): return subprocess.run([BIN] + args, capture_output=True, text=True, check=True).stdout
def loglik(cfg, seed):
    out = run([f"obs={OBS}", f"npart={NPART}", f"seed={seed}"] + [f"{k}={v}" for k, v in {**TRUTH, **cfg}.items()])
    return float(out.strip().split(",")[0])
def _eval(j): return loglik({j[0]: j[1]}, j[2])
if __name__ == "__main__":
    if not os.path.exists(OBS):
        open(OBS, "w").write(run(["simulate=1", "Tobs=2000", "dtobs=2", "seed=7"] + [f"{k}={v}" for k, v in TRUTH.items()]))
        print("simulated", OBS)
    jobs = []
    for par in ("epsH", "sigH", "beta"):
        grid = TRUTH[par] * np.logspace(-np.log10(2), np.log10(2), NPOINTS)
        for v in grid:
            for s in range(NSEED): jobs.append((par, float(v), s + 1))
    with Pool(cpu_count()) as pool:
        ll = pool.map(_eval, jobs)
    R = pd.DataFrame([dict(par=j[0], value=j[1], seed=j[2], loglik=l) for j, l in zip(jobs, ll)])
    for par, g in R.groupby("par"):
        g.to_csv(os.path.join(OUT, f"pf_profile_{par}.csv"), index=False)
        m = g.groupby("value").loglik.agg(["mean", "std"]); print(par, "\n", m.round(1).to_string())
    filt = run([f"obs={OBS}", "filt=1", f"npart={NPART}", "seed=1"] + [f"{k}={v}" for k, v in TRUTH.items()])
    open(os.path.join(OUT, "pf_filtered_truth.csv"), "w").write(filt)
    print("wrote results/pf_profile_*.csv, results/pf_filtered_truth.csv")
