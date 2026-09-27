"""Check 5: noise-induced mean shift, full model vs reduced theory, from the stored sweeps (no new runs)."""
from common import *
F = pd.read_csv(os.path.join(CODE, "m1_first_passage_theory", "results", "m1_full_vs_sigH.csv"))
Q = pd.read_csv(os.path.join(CODE, "m1_first_passage_theory", "results", "m1_fpt_vs_nu.csv")); print(Q.columns.tolist())
fP, fC, Qr = 6.22, 1.115, 0.983; RATE = fP * fC * Qr; eps_eff = 1000 / RATE; H0, HS = 0.064, 0.408; tau0 = np.log((1 - H0) / (1 - HS))
mcol = [c for c in Q.columns if c.lower().startswith("mean")][0]; nucol = "nu" if "nu" in Q.columns else Q.columns[0]
rows = []
for r in F.itertuples():
    nu = r.value * np.sqrt(eps_eff); m_th = float(np.interp(nu, Q[nucol].values, Q[mcol].values))
    rows.append(dict(sigH=r.value, nu=nu, Tdet=r.Tdet, Tav_full=r.Tav, dT_full=r.Tav - r.Tdet, dT_theory=eps_eff * (m_th - tau0), rel_shift_full=(r.Tav - r.Tdet) / r.Tdet, rel_shift_theory=eps_eff * (m_th - tau0) / r.Tdet))
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "check5_mean_shift.csv"), index=False); print(R.round(4).to_string(index=False)); print("check5 done")
