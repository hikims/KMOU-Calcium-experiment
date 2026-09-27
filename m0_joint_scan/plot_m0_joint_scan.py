"""
plot_m0_joint_scan.py — what the original AIMS equations can and cannot do jointly (period ratio, duty, amplitude).

Reads results/m0_joint_scan.csv (from run_m0_joint_scan.sh) and the data targets from
../data_dose_units/results/{fold_change_targets,dose_table}.csv (fallback: hard-coded values of 2026-09-22).
Re-runs ./m0_joint_scan curve=1 / trace=1 for the best candidates (in parallel, one process per candidate).

Criteria (train-level data, 3 -> 30 µM):
  C1  T_ratio >= R_data (6.9); lenient: >= R_ci (4.6, ratio of the CI ends)
  C2  duty_max <= 0.15 inside the window that attains T_ratio (data W50/ISI 0.05–0.15)
  C3  amp_ratio in [0.5, 2] (indicator amplitude ratio 3/30 µM is 0.5–0.8 in the data, order one)
Main   (a) T_ratio vs duty_max of every oscillating point with the target box; (b) period curves of the best
       candidates vs the baseline; (c) pass rates of C1–C3 and their conjunctions
Supp   (a) T_ratio vs lambda (the P gate saturates unless lambda is large); (b) duty_min vs epsH;
       (c) number of oscillating stimulus points; (d) C, P, H trace of the best joint candidate
Writes results/m0_joint_scan_summary.csv, results/m0_joint_best_curves.csv, results/m0_joint_best_trace.csv
"""
import os, io, subprocess
from multiprocessing import get_context, cpu_count
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 22, "axes.titlesize": 24, "axes.labelsize": 24,
    "xtick.labelsize": 20, "ytick.labelsize": 20, "legend.fontsize": 19,
    "axes.grid": True, "grid.color": "#c8c8c8", "grid.linestyle": "--", "grid.linewidth": 0.8,
    "lines.linewidth": 2.5, "lines.markersize": 6,
    "mathtext.fontset": "stix", "font.family": "STIXGeneral",
    "savefig.dpi": 300, "legend.framealpha": 1.0,
})
NAME, OUT, FIG, BIN = "m0_joint_scan", "results", "figures", "./m0_joint_scan"
# the panels call the compiled kernel; build it if it is missing (Makefile if present, else a direct cc call)
def _build(binary):
    src = os.path.basename(binary) + ".c"
    if os.path.exists("Makefile"): subprocess.run(["make", "-s"], check=True); return
    for flags in (["-O2", "-fopenmp"], ["-O2", "-Xpreprocessor", "-fopenmp", "-lomp"], ["-O2"]):
        if subprocess.run(["cc"] + flags + ["-o", binary, src, "-lm"], capture_output=True).returncode == 0: return
    raise SystemExit(f"could not compile {src} (need a C compiler; on macOS: xcode-select --install, brew install libomp)")
if not os.path.exists(BIN): _build(BIN)

os.makedirs(FIG, exist_ok=True)
C_ALL, C_PASS, C_BASE, C_BEST, C_DATA, C_TINT = "#9ecae1", "#c0392b", "#888888", "#2ca02c", "#c0392b", "#f5d0cc"
PN = ["alpha", "beta", "gamma", "rho", "mu", "lambda", "kS", "kC", "kH", "epsH"]

# --- data targets -------------------------------------------------------------------------------------------
R_DATA, R_CI, DUTY_MAX, AMP_LO, AMP_HI = 6.9, 4.6, 0.15, 0.5, 2.0
ft = "../data_dose_units/results/fold_change_targets.csv"; dt_ = "../data_dose_units/results/dose_table.csv"
if os.path.exists(ft) and os.path.exists(dt_):
    F = pd.read_csv(ft); D = pd.read_csv(dt_)
    R_DATA = float(F[(F.kind.str.startswith("between")) & (F.lo == 3) & (F.hi == 30)].ratio.iloc[0])
    lo = D[(D.dose == 3) & (D.subset == "constant")].iloc[0]; hi = D[(D.dose == 30) & (D.subset == "constant+pre")].iloc[0]
    R_CI = float(lo.isi_train_ci_lo / hi.isi_train_ci_hi)
    DUTY_MAX = float(max(0.15, hi.W_train_median / hi.isi_train_median))
print(f"targets: T ratio {R_DATA:.2f} (CI-end ratio {R_CI:.2f}), duty <= {DUTY_MAX:.2f}, amp ratio in [{AMP_LO}, {AMP_HI}]")

# --- scan ---------------------------------------------------------------------------------------------------
S = pd.read_csv(os.path.join(OUT, "m0_joint_scan.csv"))
osc = S[S.n_osc >= 3].copy()
# damped transients identified by the trace screen of round2_checks (last-quarter amplitude < 50 % of the first quarter)
_TS = "../round2_checks/results/check4b_transient_screen.csv"
if os.path.exists(_TS):
    _ts = pd.read_csv(_TS); damped = set(_ts[(_ts.ok == 1) & (_ts.ratio_q4_q1 < 0.5)].idx)
else:
    damped = set()
osc["transient"] = osc.idx.isin(damped)
osc["C1"] = osc.T_ratio >= R_DATA; osc["C1_ci"] = osc.T_ratio >= R_CI
osc["C2"] = (osc.duty_max <= DUTY_MAX) & (osc.i_lo >= 0); osc["C3"] = osc.amp_ratio.between(AMP_LO, AMP_HI) & (osc.i_lo >= 0)
rows = [dict(criterion=k, n=int(v.sum()), frac_of_osc=float(v.mean())) for k, v in
        [("oscillating (>=3 S points)", pd.Series(True, index=osc.index)), ("C1 T_ratio >= data", osc.C1), ("C1 (CI) T_ratio >= CI ends", osc.C1_ci),
         ("C2 duty <= 0.15", osc.C2), ("C3 amp ratio in [0.5,2]", osc.C3), ("C1 & C2", osc.C1 & osc.C2), ("C1(CI) & C2", osc.C1_ci & osc.C2),
         ("C2 & C3", osc.C2 & osc.C3), ("C1 & C2 & C3", osc.C1 & osc.C2 & osc.C3), ("C1(CI) & C2 & C3", osc.C1_ci & osc.C2 & osc.C3)]]
summ = pd.DataFrame(rows); summ["n_points"] = len(S); summ.to_csv(os.path.join(OUT, "m0_joint_scan_summary.csv"), index=False)
print(f"{len(S)} points, {len(osc)} oscillate at >= 3 stimulus points; clipping events in {int((S.nclip > 0).sum())} points")
with open(os.path.join(OUT, "m0_joint_scan_summary.md"), "w") as fh:            # the README quotes this file, not hand-typed numbers
    fh.write(f"Scan: {len(S)} points, {len(osc)} oscillate at >= 3 stimulus points; max T_ratio {osc.T_ratio.max():.2f}; targets: T ratio >= {R_DATA:.2f} (CI-end ratio {R_CI:.2f}), duty <= {DUTY_MAX:.2f}\n\n")
    fh.write("| criterion | n | % of oscillating |\n|---|---:|---:|\n" + "".join(f"| {r.criterion} | {r.n} | {100*r.frac_of_osc:.2f} |\n" for r in summ.itertuples()))
print(summ.to_string(index=False))

# best candidates: (i) highest T_ratio with C2, (ii) highest T_ratio overall, (iii) lowest duty_max with T_ratio >= 2
def pick(mask, col, asc): m = osc[mask]; return None if len(m) == 0 else m.sort_values(col, ascending=asc).iloc[0]
cands = [("best joint (C2, max T ratio)", pick(osc.C2 & (osc.i_lo >= 0), "T_ratio", False)),
         ("max T ratio", pick(osc.i_lo >= 0, "T_ratio", False)),
         ("narrowest with T ratio >= 2", pick((osc.T_ratio >= 2) & (osc.i_lo >= 0), "duty_max", True))]
cands = [(lab, r) for lab, r in cands if r is not None]
_seen = set(); cands = [(lab, r) for lab, r in cands if not (int(r.idx) in _seen or _seen.add(int(r.idx)))]   # distinct points only
top = osc[osc.i_lo >= 0].sort_values("T_ratio", ascending=False).head(20); top.to_csv(os.path.join(OUT, "m0_joint_top20_by_T_ratio.csv"), index=False)
for lab, r in cands: print(f"{lab}: idx {int(r.idx)} T_ratio {r.T_ratio:.2f} duty_max {r.duty_max:.3f} amp_ratio {r.amp_ratio:.2f} " + " ".join(f"{k}={r[k]:.3g}" for k in PN))

def run_bin(args):
    lab, kv, mode = args
    out = subprocess.run([BIN, mode] + [f"{k}={v}" for k, v in kv.items()], capture_output=True, text=True).stdout
    d = pd.read_csv(io.StringIO(out)); d["cand"] = lab; return d
jobs = [(lab, {k: float(r[k]) for k in PN} | {"ns": 25}, "curve=1") for lab, r in cands] + \
       [(cands[0][0], {k: float(cands[0][1][k]) for k in PN} | {"S": float(np.sqrt(cands[0][1].S_lo * cands[0][1].S_hi)) if cands[0][1].i_lo >= 0 else 1.0}, "trace=1")] if cands else []
if jobs and os.path.exists(BIN):
    with get_context("fork").Pool(min(cpu_count(), len(jobs))) as pool: res = pool.map(run_bin, jobs)
    curves = pd.concat(res[:len(cands)]); curves.to_csv(os.path.join(OUT, "m0_joint_best_curves.csv"), index=False)
    trace = res[-1]; trace.to_csv(os.path.join(OUT, "m0_joint_best_trace.csv"), index=False)
else:
    curves = pd.read_csv(os.path.join(OUT, "m0_joint_best_curves.csv")) if os.path.exists(os.path.join(OUT, "m0_joint_best_curves.csv")) else None
    trace = pd.read_csv(os.path.join(OUT, "m0_joint_best_trace.csv")) if os.path.exists(os.path.join(OUT, "m0_joint_best_trace.csv")) else None
base_curve = pd.read_csv(os.path.join(OUT, "m0_joint_baseline_curve.csv")) if os.path.exists(os.path.join(OUT, "m0_joint_baseline_curve.csv")) else None

# --- figure assembly ---------------------------------------------------------------------------------------
def letter(a, s): a.text(-0.13, 1.06, s, transform=a.transAxes, fontsize=28, fontweight="bold", va="bottom", ha="left")
def compose(panels, ncol, name, w=8.0, h=7.5):
    n = len(panels); nrow = -(-n // ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(w * ncol, h * nrow), squeeze=False)
    for k, fn in enumerate(panels):
        a = axs[k // ncol, k % ncol]; fn(a); letter(a, f"({chr(97 + k)})")
    for k in range(n, nrow * ncol): axs[k // ncol, k % ncol].axis("off")
    plt.tight_layout(w_pad=2.5, h_pad=2.5)
    for ext in ("png", "pdf", "eps"): fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig); print(f"wrote figures/{name}.{{png,pdf,eps}}")

def p_ratio_vs_duty(a):
    w = osc[osc.i_lo >= 0]; w0 = w[~w.transient]; w1 = w[w.transient]
    a.plot(w0.duty_max, w0.T_ratio, "o", color=C_ALL, ms=6, label=f"sustained oscillations ({len(w0)})")
    if len(w1): a.plot(w1.duty_max, w1.T_ratio, "x", color="k", ms=9, mew=2, label=f"damped transients ({len(w1)})")
    a.plot(w0[w0.C1_ci & w0.C2].duty_max, w0[w0.C1_ci & w0.C2].T_ratio, "o", color=C_PASS, ms=8, label="C1(CI) & C2")
    a.axhspan(R_CI, R_DATA * 1.6, xmin=0, xmax=DUTY_MAX / 0.7, color=C_TINT, zorder=0)
    a.axhline(R_DATA, color=C_DATA, ls="-", lw=2, label=r"data $T_{3}/T_{30}$"); a.axhline(R_CI, color=C_DATA, ls=":", lw=2)
    a.axvline(DUTY_MAX, color="k", ls="--", lw=2, label="duty limit")
    if base_curve is not None:
        b = base_curve[base_curve.osc == 1]; a.plot([b.duty.max()], [b["T"].max() / b["T"].min()], "s", color=C_BASE, ms=12, label="Table 2 baseline")
    a.set_xlim(0, 0.7); a.set_yscale("log"); a.set_xlabel(r"max duty $W_{50}/T$ in window"); a.set_ylabel(r"period ratio $T_{\rm lo}/T_{\rm hi}$")
    a.set_ylim(0, 30);
    a.set_title("Joint scan of M0"); a.legend(loc="upper right", fontsize=14)

def p_curves(a):
    if base_curve is not None:
        b = base_curve[base_curve.osc == 1]; a.plot(b.S, b["T"] / b["T"].iloc[0], "s-", color=C_BASE, ms=7, label="Table 2 baseline")
    if curves is not None:
        cols = [C_BEST, "#e69f00", "#9467bd"]
        for k, (lab, g) in enumerate(curves.groupby("cand", sort=False)):
            g = g[g.osc == 1]
            if len(g): a.plot(g.S, g["T"] / g["T"].iloc[0], "o-", color=cols[k % 3], ms=7, label=lab)
    a.axhline(1 / R_DATA, color=C_DATA, ls="-", lw=2, label=r"data $T_{30}/T_{3}$")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(r"$\bar S$"); a.set_ylabel(r"$T(\bar S)\,/\,T(\bar S_{\min})$"); a.set_title("Period vs stimulus"); a.legend(loc="center right", fontsize=14)

def p_pass(a):
    keys = ["C1 T_ratio >= data", "C1 (CI) T_ratio >= CI ends", "C2 duty <= 0.15", "C3 amp ratio in [0.5,2]", "C1(CI) & C2", "C2 & C3", "C1(CI) & C2 & C3"]
    labs = [r"C1: $T$ ratio", r"C1 (CI)", r"C2: duty", r"C3: amp", "C1(CI)&C2", "C2&C3", "all"]
    v = [float(summ[summ.criterion == k].frac_of_osc.iloc[0]) * 100 for k in keys]; n = [int(summ[summ.criterion == k].n.iloc[0]) for k in keys]
    a.bar(range(len(v)), v, color=[C_ALL] * 4 + [C_PASS] * 3)
    for i, (vv, nn) in enumerate(zip(v, n)): a.text(i, vv, str(nn), ha="center", va="bottom", fontsize=16)
    a.set_xticks(range(len(v))); a.set_xticklabels(labs, rotation=35, ha="right", fontsize=15); a.set_ylabel("% of oscillating points"); a.set_title(f"Pass rates ({len(osc)} oscillating)")

def p_ratio_vs_lambda(a):
    w = osc[osc.i_lo >= 0]; a.plot(w["lambda"], w.T_ratio, "o", color=C_ALL, ms=6); a.plot(w[w.transient]["lambda"], w[w.transient].T_ratio, "x", color="k", ms=9, mew=2, label="damped transients"); a.plot(w[w.C2]["lambda"], w[w.C2].T_ratio, "o", color=C_PASS, ms=7, label="C2 flag (transient)")
    a.axhline(R_DATA, color=C_DATA, lw=2); a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(r"$\lambda$"); a.set_ylabel(r"$T_{\rm lo}/T_{\rm hi}$"); a.set_title("Period ratio vs IP$_3$ decay"); a.legend(loc="center left", fontsize=16)

def p_duty_vs_epsH(a):
    a.plot(osc.epsH, osc.duty_min, "o", color=C_ALL, ms=6); a.plot(osc[osc.C1_ci].epsH, osc[osc.C1_ci].duty_min, "o", color=C_PASS, ms=7, label="C1 (CI)")
    a.axhline(DUTY_MAX, color="k", ls="--", lw=2); a.set_xscale("log"); a.set_xlabel(r"$\varepsilon_H$"); a.set_ylabel("min duty over stimulus"); a.set_title("Spike width vs slow time scale"); a.legend(loc="lower left", fontsize=16)

def p_nosc(a):
    a.hist(S.n_osc, bins=np.arange(-0.5, S.n_osc.max() + 1.5), color=C_ALL, edgecolor="k")
    a.set_xlabel("oscillating stimulus points (of 9)"); a.set_ylabel("sample points"); a.set_title("Oscillation over the grid")

def p_trace(a):
    """the single point flagged by C2: C(t) and the peak-to-peak amplitude of C in 8 consecutive windows (log scale)
    -> the amplitude decays by four orders of magnitude over the recording window: a damped transient, not a
    sustained oscillation (round2_checks/check4_scan_candidate.py)"""
    if trace is None: a.text(0.5, 0.5, "no trace", ha="center", transform=a.transAxes); return
    t = (trace.t - trace.t.iloc[0]).values; C = trace.C.values
    a.plot(t, C, "-", color=C_BEST, lw=1.5, label=r"$C(t)$")
    a.set_xlabel(r"$t$ after the transient (model units)"); a.set_ylabel(r"$C$", color=C_BEST)
    nw = 8; edges = np.linspace(0, len(C), nw + 1).astype(int)
    amp = np.array([C[edges[k]:edges[k + 1]].max() - C[edges[k]:edges[k + 1]].min() for k in range(nw)])
    tm = np.array([0.5 * (t[edges[k]] + t[edges[k + 1] - 1]) for k in range(nw)])
    a2 = a.twinx(); a2.grid(False); a2.semilogy(tm, np.maximum(amp, 1e-7), "s--", color="k", ms=8, label="peak-to-peak amplitude per window")
    a2.set_ylabel(r"$C_{\max}-C_{\min}$ per window (log)")
    a.set_title("Point flagged by C2: damped transient")
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels(); a.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=14)

MAIN = [p_ratio_vs_duty, p_curves, p_pass]              # the original equations: where the joint criteria fail
SUPP = [p_ratio_vs_lambda, p_duty_vs_epsH, p_nosc, p_trace]
if __name__ == "__main__":
    compose(MAIN, 3, f"{NAME}_main")
    compose(SUPP, 2, f"{NAME}_supp")
