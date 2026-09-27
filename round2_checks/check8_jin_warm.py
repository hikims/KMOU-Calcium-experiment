"""Check 8: J_in = 0 population with a warm start (stimulus S_E already on during the pre-phase, cells oscillating)
versus the cold start of the paper (rest at S=0 for 3000 units, then switch)."""
from common import *
D = os.path.join(CODE, "emrich_jin_prediction"); BIN = os.path.join(D, "emrich_jin_sde")
cal = pd.read_csv(os.path.join(D, "results", "jin_calibration.csv")).iloc[0]; S_E, TWIN, NUH = cal.S_E, cal.Twin_units, cal.nuH
REFJ = [f"{k}={v}" for k, v in REF.items() if k not in ("Jin", "S")]
rows = []
for f in (0.0, 0.01, 0.075):
    for start, Spre in (("cold", 0.0), ("warm", S_E)):
        for seed in (7, 8):
            p = subprocess.run([BIN, "cells=400", f"Twin={TWIN:.2f}", "Tpre=3000", f"Spre={Spre}", "epsH_cv=0.45", f"S={S_E}", f"nuH={NUH}", "dt=0.01", "Ttrans=3000", f"Jin={0.02*f}", f"seed={seed}"] + REFJ, capture_output=True, text=True)
            d = pd.read_csv(io.StringIO(p.stdout))
            rows.append(dict(factor=f, start=start, Spre=Spre, seed=seed, median=d.n_spikes.median(), q1=d.n_spikes.quantile(.25), q3=d.n_spikes.quantile(.75), mean=d.n_spikes.mean(), frac_silent=(d.n_spikes == 0).mean())); print(rows[-1], flush=True)
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "check8_jin_warm.csv"), index=False); print(R.to_string(index=False))
tr = pd.read_csv(os.path.join(D, "results", "jin_trace_f0.00.csv")); print("stored cold-start trace at J_in=0: C max in window", tr.C.max(), "H range", tr.H.min(), tr.H.max())
