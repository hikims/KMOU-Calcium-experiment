"""Check 4b: how many of the 232 'oscillating' M0 scan points are damped transients? For each point, trace at its
oscillating mid-level and compare the C amplitude in the first and last quarter of the recorded window."""
from common import *
from multiprocessing import Pool, cpu_count
D = os.path.join(CODE, "m0_joint_scan"); BIN = os.path.join(D, "m0_joint_scan")
S = pd.read_csv(os.path.join(D, "results", "m0_joint_scan.csv")); osc = S[S.n_osc >= 3].copy()
PN = ["alpha", "beta", "gamma", "rho", "mu", "lambda", "kS", "kC", "kH", "epsH"]
LEVELS = np.geomspace(0.03, 30, 9)
def one(args):
    idx, kv, Sv = args
    p = subprocess.run([BIN, "trace=1"] + [f"{k}={v}" for k, v in {**kv, "S": Sv}.items()], capture_output=True, text=True)
    if p.returncode or "\n" not in p.stdout.strip(): return dict(idx=idx, S=Sv, ok=0)
    x = pd.read_csv(io.StringIO(p.stdout)).C.values; n = len(x); q = n // 4
    a = [x[j*q:(j+1)*q].max() - x[j*q:(j+1)*q].min() for j in range(4)]
    return dict(idx=idx, S=Sv, ok=1, amp_q1=a[0], amp_q2=a[1], amp_q3=a[2], amp_q4=a[3], ratio_q4_q1=a[3] / a[0] if a[0] > 0 else np.nan)
def main():
  jobs = []
  for r in osc.itertuples():
      kv = {k: float(osc.loc[r.Index, k]) for k in PN}
      i_mid = int((r.i_lo + r.i_hi) // 2) if r.i_lo >= 0 else 4
      jobs.append((int(r.idx), kv, float(LEVELS[i_mid])))
  with Pool(min(cpu_count(), 16)) as pool: res = pool.map(one, jobs)
  R = pd.DataFrame(res); R.to_csv(os.path.join(OUT, "check4b_transient_screen.csv"), index=False)
  ok = R[R.ok == 1]; damped = ok[ok.ratio_q4_q1 < 0.5]
  print(f"{len(R)} oscillating points traced ({int((R.ok==0).sum())} without a usable trace); {len(damped)} have last-quarter amplitude < 50% of the first quarter (damped transients); {(ok.ratio_q4_q1 < 0.1).sum()} below 10%")
  print(ok.ratio_q4_q1.describe())
if __name__ == "__main__": main()
