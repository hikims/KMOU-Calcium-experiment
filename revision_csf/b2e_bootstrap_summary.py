"""
b2e_bootstrap_summary.py — summary of the parametric bootstrap of b2_family_comparison.py from its stored replicates.

Reads results/b2_E_parametric_bootstrap_raw.csv (total log-likelihood of every family in every replicate) and the
observed totals of results/b2_A_primary.csv; writes results/b2_E_parametric_bootstrap.csv with the rule of
isifam.boot_summary: replicates in which a fit was infeasible are excluded pairwise and counted (n_failed).
"""
import os, sys, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); OUT = os.path.join(HERE, "results")
import isifam as F
R = pd.read_csv(os.path.join(OUT, "b2_E_parametric_bootstrap_raw.csv")); names = [c for c in R.columns if c != "truth"]
obs = pd.read_csv(os.path.join(OUT, "b2_A_primary.csv")).set_index("family").logL.to_dict()
Eb = F.boot_summary(R, obs, names); Eb.to_csv(os.path.join(OUT, "b2_E_parametric_bootstrap.csv"), index=False)
pd.set_option("display.width", 220); print(Eb.round(2).to_string(index=False))
print("\nreplicates with at least one infeasible fit, by generating family:")
print(R[(R[names].abs() > 1e9).any(axis=1)].groupby("truth").size().reindex(names, fill_value=0).to_string())
