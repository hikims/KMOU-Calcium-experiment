"""
run_emrich_jin.py — driver: calibrate the stimulus to Emrich's WT spike count, then sweep the basal influx Jin.

Steps (each call of the binary is OpenMP-parallel; the driver itself is sequential)
  1. deterministic S̄ sweep at REF_ST (Jin = 0.02) -> period T(S̄)
  2. S̄_E: the stimulus whose deterministic period gives the WT median count (6 spikes / 14 min at 10 µM CCh, Emrich
     Fig. 1D; count ≈ 1 + Twin/T with the first spike shortly after the switch) with the time unit T0 (s per model
     unit; ../data_dose_units/results/dose_map_refit.csv, else 0.646)
  3. Jin sweep at S̄_E with H-noise (ν common, nuH = 0.822 = 0.026·sqrt(1000)): T, CV, amplitude, existence
  4. cell populations (N cells, ε_H lognormal with epsH_cv) under the step protocol for Jin factors
     f ∈ {1, 0.75, 0.5, 0.35, 0.25, 0.1, 0}: spike counts in Twin = 14 min -> results/jin_cells_f<f>.csv
  5. one trace per factor for the figure
Usage: python3 run_emrich_jin.py [NCELLS=400] [EPSH_CV=0.45]
"""
import os, sys, io, subprocess
import numpy as np, pandas as pd
NCELLS = int(sys.argv[1]) if len(sys.argv) > 1 else 400
EPSH_CV = float(sys.argv[2]) if len(sys.argv) > 2 else 0.45
BIN, OUT = "./emrich_jin_sde", "results"; os.makedirs(OUT, exist_ok=True)
REF = "mech=2 rho=0.024 gH=30 dH=10000 epsH=1000 nA=2 kA=0.4 Jin=0.02 kS=0.25 beta=4 alpha=20 qoff=0 rC=300 kR=1".split()
NUH = 0.822; JIN0 = 0.02
T0 = 0.646
fn = "../data_dose_units/results/dose_map_refit.csv"
if os.path.exists(fn): T0 = float(pd.read_csv(fn).T0_s.iloc[0])
WT_COUNT, TWIN_S = 6.0, 14 * 60.0
fn = "../emrich_replication/results/emrich_counts_by_genotype.csv"
if os.path.exists(fn):
    c = pd.read_csv(fn); WT_COUNT = float(c[c.genotype.str.startswith("WT")]["median"].iloc[0])
TWIN = TWIN_S / T0
def run(args):
    p = subprocess.run([BIN] + args, capture_output=True, text=True)
    if p.returncode: raise SystemExit(p.stderr)
    return pd.read_csv(io.StringIO(p.stdout))

# 1. deterministic period vs stimulus (sigH = 0 -> the stochastic part is skipped quickly: Tend small)
sw = run(["sweep=S", "lo=0.003", "hi=30", "n=41", "log=1", "Tend=4000", "Ttrans=3000", "dt=0.01"] + REF)
sw.to_csv(os.path.join(OUT, "jin_S_sweep_det.csv"), index=False)
# 2. calibrate S_E on the population: first guess from count ≈ 1 + Twin/T, then a scan of the population median
T_target = TWIN / (WT_COUNT - 1); o = sw[sw.Tdet > 0].sort_values("Tdet")
S_guess = float(np.exp(np.interp(np.log(T_target), np.log(o.Tdet.values), np.log(o.value.values))))
def cells(S, Jin, n, seed=7):
    return run([f"cells={n}", f"Twin={TWIN:.2f}", "Tpre=3000", "Spre=0", f"epsH_cv={EPSH_CV}", f"S={S}", f"nuH={NUH}", "dt=0.01", "Ttrans=3000", f"Jin={Jin}", f"seed={seed}"] + [a for a in REF if not a.startswith("Jin=")])
cal_rows = []
for S in S_guess * np.logspace(-0.6, 0.9, 9):
    d = cells(S, JIN0, max(60, NCELLS // 3)); cal_rows.append(dict(S=S, median=d.n_spikes.median(), mean=d.n_spikes.mean(), q1=d.n_spikes.quantile(.25), q3=d.n_spikes.quantile(.75)))
cal = pd.DataFrame(cal_rows); cal.to_csv(os.path.join(OUT, "jin_S_calibration_scan.csv"), index=False)
S_E = float(np.exp(np.interp(WT_COUNT, cal["mean"].values, np.log(cal.S.values))))          # mean count is monotone and smooth
pd.DataFrame([dict(T0_s=T0, Twin_units=TWIN, WT_median_count=WT_COUNT, T_target_units=T_target, T_target_s=T_target * T0, S_guess=S_guess, S_E=S_E, nuH=NUH, epsH_cv=EPSH_CV, ncells=NCELLS)]).to_csv(os.path.join(OUT, "jin_calibration.csv"), index=False)
print(f"T0 = {T0:.3f} s/unit, Twin = {TWIN:.0f} units, WT median count {WT_COUNT:.0f} -> S_guess = {S_guess:.4f} -> population-calibrated S_E = {S_E:.4f}")
print(cal.to_string(index=False))

# 3. Jin sweep at S_E with H-noise
js = run(["sweep=Jin", "lo=0", "hi=0.05", "n=26", f"S={S_E}", f"nuH={NUH}", "Tend=60000", "Ttrans=3000", "dt=0.01"] + [a for a in REF if not a.startswith("Jin=")])
js.to_csv(os.path.join(OUT, "jin_sweep_S_E.csv"), index=False); print(js[["value", "n", "Tdet", "Tav", "CV", "amp_det", "amp_cv", "W50_det"]].to_string(index=False))

# 4. cell populations under the step protocol, 5. traces
FACTORS = [1.0, 0.75, 0.5, 0.35, 0.25, 0.1, 0.0]
rows = []
for f in FACTORS:
    d = cells(S_E, JIN0 * f, NCELLS); d["factor"] = f; d.to_csv(os.path.join(OUT, f"jin_cells_f{f:.2f}.csv"), index=False)
    rows.append(dict(factor=f, Jin=JIN0 * f, median=d.n_spikes.median(), q1=d.n_spikes.quantile(.25), q3=d.n_spikes.quantile(.75), mean=d.n_spikes.mean(), frac_silent=(d.n_spikes == 0).mean(), amp_median=d[d.amp_mean > 0].amp_mean.median()))
    tr = run(["trace=1", f"Twin={TWIN:.2f}", "Tpre=3000", "Spre=0", f"S={S_E}", f"nuH={NUH}", "dt=0.01", "Ttrans=3000", f"Jin={JIN0*f}", "seed=3"] + [a for a in REF if not a.startswith("Jin=")])
    tr.to_csv(os.path.join(OUT, f"jin_trace_f{f:.2f}.csv"), index=False)
    print(f"Jin factor {f:.2f}: count median {rows[-1]['median']:.0f} [{rows[-1]['q1']:.0f}, {rows[-1]['q3']:.0f}], silent {100*rows[-1]['frac_silent']:.0f} %")
pd.DataFrame(rows).to_csv(os.path.join(OUT, "jin_counts_summary.csv"), index=False)
print("done -> results/")
