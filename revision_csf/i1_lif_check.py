"""
i1_lif_check.py — the bound R >= 3 against the leaky integrate-and-fire (Ornstein-Uhlenbeck) model, for which Ostojic
(J. Neurophysiol. 106, 361, 2011) notes that the region in the CV-skewness plane lies close to the gamma line SK = 2 CV.

dV = (mu - V) dt + nu dW, reset 0, threshold 1.  Supra-threshold (mu > 1): b > 0 on [0, 1] and the small-noise theory
applies with D = 1:  R = 3 T0 J2 / J1^2,  T0 = int dx/(mu-x), J1 = int dx/(mu-x)^3, J2 = int dx/(mu-x)^5.
Exact cumulants from the linear ODE system of a2 on [-L, 1] with reflection at -L (L = 4, far below reset).
Output: results/i1_lif_check.csv
"""
import os, sys, numpy as np, pandas as pd
from scipy.integrate import quad
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import a2_exact_vs_expansion as A
L = 4.0; rows = []
for mu in (1.1, 1.2, 1.5, 2.0, 5.0):
    A.FORMS["lif"] = (lambda y, mu=mu: mu - (y - L), lambda y: 1.0 + 0 * y)
    T0, J1, J2 = (quad(lambda x, p=p: (mu - x) ** (-p), 0, 1)[0] for p in (1, 3, 5)); Rlim = 3 * T0 * J2 / J1**2
    for cv_target in (0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7):
        nu = cv_target * T0 / np.sqrt(J1)                              # leading-order noise for this CV
        try: e = A.exact("lif", nu, h0=L, hs=L + 1.0)
        except Exception as ex: continue
        rows.append(dict(mu=mu, nu=nu, T0=T0, mean=e["mean"], CV=e["CV"], skew=e["skew"], skew_over_CV=e["R"], R_limit=Rlim))
R = pd.DataFrame(rows); R.to_csv(os.path.join(HERE, "results", "i1_lif_check.csv"), index=False)
pd.set_option("display.width", 200); print(R.round(4).to_string(index=False))
print("\npivot: exact skew/CV (rows mu, columns target CV)"); R["cv_t"] = R.groupby("mu").cumcount(); print(R.pivot(index="mu", columns="cv_t", values="skew_over_CV").round(3).to_string()); print(R.groupby("mu").R_limit.first().round(3).to_dict())
