"""Check 4: the 'best joint candidate' of the M0 scan (Fig. S5d): oscillation per stimulus level and the trace at a grid level."""
from common import *
D = os.path.join(CODE, "m0_joint_scan"); BIN = os.path.join(D, "m0_joint_scan")
S = pd.read_csv(os.path.join(D, "results", "m0_joint_scan.csv")); osc = S[S.n_osc >= 3].copy()
osc["C2"] = (osc.duty_max <= 0.15) & (osc.i_lo >= 0); c = osc[osc.C2].sort_values("T_ratio", ascending=False).iloc[0]
PN = ["alpha", "beta", "gamma", "rho", "mu", "lambda", "kS", "kC", "kH", "epsH"]; kv = {k: float(c[k]) for k in PN}
print("candidate idx", int(c.idx), {k: round(v, 4) for k, v in kv.items()}, "T_ratio", c.T_ratio, "duty_max", c.duty_max, "i_lo,i_hi", c.i_lo, c.i_hi, "S_lo,S_hi", c.S_lo, c.S_hi)
def runm(mode, extra):
    p = subprocess.run([BIN, mode] + [f"{k}={v}" for k, v in {**kv, **extra}.items()], capture_output=True, text=True)
    if p.returncode or not p.stdout.strip() or "\n" not in p.stdout.strip(): return pd.DataFrame()
    return pd.read_csv(io.StringIO(p.stdout))
cur = runm("curve=1", {"ns": 9}); cur.to_csv(os.path.join(OUT, "check4_candidate_curve.csv"), index=False); print(cur.to_string(index=False))
rows = []
for Sv in cur.S.values:
    tr = runm("trace=1", {"S": Sv})
    if len(tr) == 0: rows.append(dict(S=Sv, osc_flag=int(cur[cur.S == Sv].osc.iloc[0]), amp_first_half=np.nan, amp_last_quarter=np.nan, trace_len=0)); continue
    x = tr.C.values; n = len(x)
    q = n // 4; amps = [x[j*q:(j+1)*q].max() - x[j*q:(j+1)*q].min() for j in range(4)]
    rows.append(dict(S=Sv, osc_flag=int(cur[cur.S == Sv].osc.iloc[0]), amp_q1=amps[0], amp_q2=amps[1], amp_q3=amps[2], amp_q4=amps[3], trace_len=n)); continue
    rows.append(dict(S=Sv, osc_flag=int(cur[cur.S == Sv].osc.iloc[0]), amp_first_half=x[:n//2].max() - x[:n//2].min(), amp_last_quarter=x[3*n//4:].max() - x[3*n//4:].min(), trace_len=n))
T = pd.DataFrame(rows); T.to_csv(os.path.join(OUT, "check4_candidate_traces.csv"), index=False); print(T.to_string(index=False)); print("check4 done")
