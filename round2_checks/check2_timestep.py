"""Check 2: time-step refinement (dt = 0.01, 0.005, 0.0025) at the tuned amplitudes; 3 seeds each; raw statistics."""
from common import *
CASES = {"H": {"sigH": 0.0232}, "C": {"sigC": 0.362}, "P": {"sigP": 0.0151}}
rows = []
for k, ex in CASES.items():
    for dt in (0.01, 0.005, 0.0025):
        for seed in (1, 2, 3):
            st = run(BIN_ISI, {**ex, "dt": dt, "seed": seed}, None).iloc[0]
            rows.append(dict(noise=k, dt=dt, seed=seed, n=st.n, Tav=st.Tav, sd=st.sd, CV=st.CV, skew=st["skew"], amp_cv=st.amp_cv, corr_isi_amp=st.corr_isi_amp, serial=st.corr_isi_next, Tdet=st.Tdet, amp_det=st.amp_det)); print(rows[-1], flush=True)
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "check2_timestep_raw.csv"), index=False)
A = R.groupby(["noise", "dt"]).agg(n=("n", "mean"), Tav=("Tav", "mean"), Tav_sd=("Tav", "std"), sd=("sd", "mean"), sd_sd=("sd", "std"), CV=("CV", "mean"), CV_sd=("CV", "std"),
                                    skew=("skew", "mean"), skew_sd=("skew", "std"), amp_cv=("amp_cv", "mean"), amp_cv_sd=("amp_cv", "std"), corr=("corr_isi_amp", "mean"), corr_sd=("corr_isi_amp", "std"), serial=("serial", "mean"), serial_sd=("serial", "std"), Tdet=("Tdet", "mean"), amp_det=("amp_det", "mean")).reset_index()
A.to_csv(os.path.join(OUT, "check2_timestep_summary.csv"), index=False); print(A.round(4).to_string(index=False)); print("check2 done")
