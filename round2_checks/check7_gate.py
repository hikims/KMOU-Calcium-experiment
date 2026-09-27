"""Check 7: IP3 gate saturation error 1 - P/(1+P) at the mapped doses: analytic bounds and the actual mean P from the stimulus line."""
from common import *
g, l, rho = 4.2, 0.001, 0.024
L = pd.read_csv(os.path.join(CODE, "m0ip3h_bifurcation_structure", "results", "bif_line_S_ref_st.csv")); print(L.columns.tolist()[:12])
pc = [c for c in L.columns if c.lower().startswith("pmean") or c == "P_fp" or c.lower() == "p"]; print("P columns:", pc)
rows = []
for dose, S in zip((3, 5, 10, 15, 30), (0.0635, 0.1059, 0.2118, 0.3177, 0.6353)):
    Pmin = g * S / (1 + S) / l; Pmax = Pmin + rho / l
    rec = dict(dose_uM=dose, S=S, P_lo=Pmin, P_hi=Pmax, gate_err_lo=1 / (1 + Pmax), gate_err_hi=1 / (1 + Pmin))
    if pc:
        i = int(np.argmin(np.abs(np.log(L.x.values) - np.log(S)))); rec["S_grid"] = L.x.values[i]; rec["P_line"] = L[pc[0]].values[i]; rec["gate_err_line"] = 1 / (1 + rec["P_line"])
    rows.append(rec)
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "check7_gate.csv"), index=False); print(R.to_string(index=False)); print("check7 done")
