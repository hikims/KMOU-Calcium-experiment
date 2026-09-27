"""Check 3: direct population simulation at the interpolated basal-influx factor 0.075 (and 0.05, 0.1), plus the stationary ratios at that Jin."""
from common import *
D = os.path.join(CODE, "emrich_jin_prediction"); BIN = os.path.join(D, "emrich_jin_sde")
cal = pd.read_csv(os.path.join(D, "results", "jin_calibration.csv")).iloc[0]; S_E, TWIN, NUH = cal.S_E, cal.Twin_units, cal.nuH
REFJ = [f"{k}={v}" for k, v in REF.items() if k not in ("Jin", "S")]
def runj(args):
    p = subprocess.run([BIN] + args + REFJ, capture_output=True, text=True)
    if p.returncode: raise SystemExit(p.stderr)
    return pd.read_csv(io.StringIO(p.stdout))
rows = []
for f in (0.1, 0.075, 0.05):
    for seed in (7, 8):
        d = runj([f"cells=400", f"Twin={TWIN:.2f}", "Tpre=3000", "Spre=0", "epsH_cv=0.45", f"S={S_E}", f"nuH={NUH}", "dt=0.01", "Ttrans=3000", f"Jin={0.02*f}", f"seed={seed}"])
        rows.append(dict(factor=f, seed=seed, median=d.n_spikes.median(), q1=d.n_spikes.quantile(.25), q3=d.n_spikes.quantile(.75), mean=d.n_spikes.mean(), frac_silent=(d.n_spikes == 0).mean())); print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "check3_jin_population.csv"), index=False)
st = {}
for lab, f in (("ref", 1.0), ("f0.075", 0.075)):
    st[lab] = runj([f"S={S_E}", f"nuH={NUH}", "Tend=120000", "Ttrans=3000", "dt=0.01", f"Jin={0.02*f}"]).iloc[0]
rat = {c: st["f0.075"][c] / st["ref"][c] for c in ("Tdet", "Tav", "amp_det", "amp_mean", "W50_det", "CV")}
pd.DataFrame([dict(S_E=S_E, **{f"ratio_{c}": v for c, v in rat.items()}, **{f"ref_{c}": st["ref"][c] for c in rat}, **{f"f0075_{c}": st["f0.075"][c] for c in rat})]).to_csv(os.path.join(OUT, "check3_jin_ratios.csv"), index=False)
print("direct ratios TKO-like/WT at factor 0.075:", {c: round(v, 3) for c, v in rat.items()}); print("check3 done")
