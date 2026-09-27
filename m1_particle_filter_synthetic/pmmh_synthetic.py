"""
pmmh_synthetic.py — particle marginal Metropolis–Hastings on the synthetic record (optional; the calibration step
will reuse it per cell).  Random-walk proposals on log(eps_H, sigma_H, beta), flat priors on the log scale within
[truth/10, truth*10]; the PF marginal likelihood is re-estimated at every iteration (fresh seed).
Usage: python3 pmmh_synthetic.py [NITER=300] [NPART=2000] [STEP=0.08]      -> results/pmmh_chain.csv
Cost: NITER x (one PF likelihood: ~10 s at N=1000, ~20 s at N=2000).  For a usable chain aim at NITER >= 2000 and
NPART >= 4000 (log-likelihood sd across PF seeds ~4 at N=1000, ~2 at N=4000 on this record).
"""
import os, subprocess, sys
import numpy as np, pandas as pd
BIN, OUT = "./m1_pf_sde", "results"
NITER = int(sys.argv[1]) if len(sys.argv) > 1 else 300
NPART = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
STEP = float(sys.argv[3]) if len(sys.argv) > 3 else 0.08
FIXED = dict(a=0.11, b=0.0, sobs=0.05); TRUTH = dict(epsH=1000.0, sigH=0.026, beta=4.0)
OBS = os.path.join(OUT, "synthetic_truth.csv")
def loglik(theta, seed):
    out = subprocess.run([BIN, f"obs={OBS}", f"npart={NPART}", f"seed={seed}"] + [f"{k}={v}" for k, v in {**FIXED, **theta}.items()],
                         capture_output=True, text=True, check=True).stdout
    return float(out.strip().split(",")[0])
rng = np.random.default_rng(0)
theta = {k: v * np.exp(rng.normal(0, 0.3)) for k, v in TRUTH.items()}        # start away from the truth
ll = loglik(theta, 1); rows = []; acc = 0
for it in range(NITER):
    prop = {k: v * np.exp(rng.normal(0, STEP)) for k, v in theta.items()}
    if all(TRUTH[k] / 10 <= prop[k] <= TRUTH[k] * 10 for k in prop):
        llp = loglik(prop, it + 2)
        if np.log(rng.random()) < llp - ll: theta, ll, acc = prop, llp, acc + 1
    rows.append(dict(iter=it, loglik=ll, **theta))
    if it % 20 == 0: print(f"it {it}: ll {ll:.1f} acc {acc/(it+1):.2f} " + " ".join(f"{k}={v:.4g}" for k, v in theta.items()), flush=True)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "pmmh_chain.csv"), index=False); print("wrote results/pmmh_chain.csv; acceptance", acc / NITER)
