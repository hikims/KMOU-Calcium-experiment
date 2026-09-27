"""
plot_m0ip3h_amplitude_invariance.py — spike stereotypy: the stimulus changes the period, not the spike.

Facts to reconcile (processed spike tables, 3–30 µM carbachol)
  data : ISI median falls ~7x (206 -> 31 s), half-width falls ~2x (10.5 -> 4.5 s), amplitude ΔF/F0 ~ flat;
         within a cell, spike amplitude is almost uncorrelated with the preceding ISI (corr ≈ 0.09)
  model: in the previous reference regime IP3 rescales the whole H clock (recovery AND inactivation), so
         period, width and amplitude all fall together (~4x)
Structure tested here: the H gate written with separate on/off rates; IP3 accelerates recovery, and a
fraction qoff of that acceleration also acts on Ca-inactivation (qoff=1: previous model; qoff=0: the spike
itself is IP3-independent).  Ca-inactivation gets its own threshold kR (off at rest) and speed rC.

Panels (functions; MAIN / SUPP at the end)
  data_norm     data: ISI, half-width, amplitude vs dose, each normalised to its 3 µM median (median, IQR)
  model_norm    model, stereotyped reference (qoff=0): T, W50, amplitude vs S̄, normalised at S̄=0.1
                (the previous reference, qoff=1, scales all three together: see results/m0ip3h_inv_vs_stimulus_ref_prev.csv)
  ratios_qoff   ratios (value at S̄=0.1 over value at S̄=30) vs qoff, with the data ratios (3 µM / 30 µM)
  recovery      recovery function: amplitude of a spike evoked by a Ca pulse at delay τ after a free spike
                (what a Ca-triggered early spike would look like) vs τ/T_wait, both regimes
  amp_vs_isi    data: within-cell amplitude vs preceding ISI (both normalised by the cell median), corr statistics

Run after ./run_m0ip3h_invariance_sweeps.sh .  Writes figures/m0ip3h_amplitude_invariance_main.{png,pdf,eps}
(data_norm, model_norm, ratios_qoff: IP3 acts on recovery, not on the spike),
figures/m0ip3h_amplitude_invariance_supp.{png,pdf,eps} (recovery, amp_vs_isi: the Ca-pulse argument for H-triggered spikes),
results/m0ip3h_inv_ratios_vs_qoff.csv, results/m0ip3h_inv_data_by_dose.csv
"""
import os
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
C_DATA, C_PREV, C_ST, C_GREY = "#c0392b", "#e69f00", "#2ca02c", "#888888"
COL = {"T": "#1f77b4", "W": "#9467bd", "A": "#2ca02c"}       # period / width / amplitude
OUT, FIG = "results", "figures"
os.makedirs(FIG, exist_ok=True)

# --- data: canonical tables of data_dose_units (train = record:ST:dose, union of all files) -----------------------
TT, DT, EV = ("../data_dose_units/results/" + f for f in ("train_table.csv", "dose_table.csv", "events_canonical.csv"))
if not all(os.path.exists(f) for f in (TT, DT, EV)): raise SystemExit("run ../data_dose_units/data_dose_units.py first")
T_ = pd.read_csv(TT); T_ = T_[T_.n >= 10]; Dt = pd.read_csv(DT); Dt = Dt[Dt.subset == "all"].set_index("dose")
D = pd.DataFrame({"n": Dt.n_trains, "ISI": Dt.isi_train_median, "ISI_q1": Dt.isi_train_q1, "ISI_q3": Dt.isi_train_q3,
                  "W": Dt.W_train_median, "W_q1": Dt.W_train_q1, "W_q3": Dt.W_train_q3,
                  "A": Dt.amp_train_median, "A_q1": Dt.amp_train_q1, "A_q3": Dt.amp_train_q3})         # medians of train medians
D.index = D.index.astype(int); D.to_csv(os.path.join(OUT, "m0ip3h_inv_data_by_dose.csv"))
doses = np.array(sorted(D.index)); d0 = int(doses[0])
# within-train correlation of the amplitude with the preceding (detrended) ISI; normalised spike table for the scatter
cc = T_.corr_detr.dropna().values; trains = T_
CORR_MED, CORR_Q1, CORR_Q3 = np.median(cc), np.quantile(cc, .25), np.quantile(cc, .75)
E = pd.read_csv(EV); E = E[E.train.isin(T_.train) & E.isi_detr.notna()]
norm = pd.concat([pd.DataFrame({"isi": g.isi_detr / g.isi_detr.median(), "amp": g.amplitude / g.amplitude.median()}) for _, g in E.groupby("train")])
# amplitude (ΔF/F0) differs between records by 3x (indicator loading, F0), so between-dose amplitude medians are not
# comparable; the amplitude is compared WITHIN cells across the step levels (same ROI, lower/higher dose, >= 5 spikes each)
Tw = pd.read_csv(TT); Tw = Tw[Tw.n >= 5]; STEP = []
for lo_d, hi_d in [(5, 10), (10, 15), (5, 30), (5, 15)]:
    a_ = Tw[Tw.dose == lo_d].set_index(["record", "ST"]); b_ = Tw[Tw.dose == hi_d].set_index(["record", "ST"]); c_ = a_.index.intersection(b_.index)
    if len(c_): STEP.append(dict(lo=lo_d, hi=hi_d, n=len(c_), isi_ratio=float((a_.loc[c_].isi_median / b_.loc[c_].isi_median).median()),
                                 W_ratio=float((a_.loc[c_].W_median / b_.loc[c_].W_median).median()), amp_ratio=float((a_.loc[c_].amp_median / b_.loc[c_].amp_median).median())))
STEP = pd.DataFrame(STEP); STEP.to_csv(os.path.join(OUT, "m0ip3h_inv_within_cell_steps.csv"), index=False)
AMP_530 = float(STEP[(STEP.lo == 5) & (STEP.hi == 30)].amp_ratio.iloc[0]) if len(STEP) else np.nan
W_530 = float(STEP[(STEP.lo == 5) & (STEP.hi == 30)].W_ratio.iloc[0]) if len(STEP) else np.nan
print(f"data: {len(trains)} trains; corr(ISI, amp) median {CORR_MED:.2f}, IQR [{CORR_Q1:.2f}, {CORR_Q3:.2f}]")
print("within-cell step ratios (lower/higher dose):\n" + STEP.round(2).to_string(index=False))
print("data ratios 3 uM / 30 uM (train medians): ISI %.2f, W %.2f, amp %.2f" % (D.ISI[d0] / D.ISI[30], D.W[d0] / D.W[30], D.A[d0] / D.A[30]))

# --- model -----------------------------------------------------------------------------------------
def load(name): return pd.read_csv(os.path.join(OUT, name))
prev = load("m0ip3h_inv_vs_stimulus_ref_prev.csv"); prev = prev[prev.osc == 1]
st = load("m0ip3h_inv_vs_stimulus_ref_st.csv"); st = st[st.osc == 1]

def normalised(m, Sref):
    """T, W50, amp of a stimulus sweep divided by their values at S̄ = Sref (log-interpolated)"""
    out = {}
    for key, col in (("T", "Tosc"), ("W", "W50"), ("A", "amp")):
        ref = np.exp(np.interp(np.log(Sref), np.log(m.value.values), np.log(m[col].values)))
        out[key] = m[col].values / ref
    return out

# ratios vs qoff
lo, hi = load("m0ip3h_inv_vs_qoff_S0.1.csv"), load("m0ip3h_inv_vs_qoff_S30.csv")
R = pd.DataFrame({"qoff": lo.value, "T_ratio": lo.Tosc / hi.Tosc, "W_ratio": lo.W50 / hi.W50, "amp_ratio": lo.amp / hi.amp,
                  "duty_lo": lo.duty, "duty_hi": hi.duty})
R.to_csv(os.path.join(OUT, "m0ip3h_inv_ratios_vs_qoff.csv"), index=False)
rec_prev, rec_st = load("m0ip3h_inv_recovery_ref_prev.csv"), load("m0ip3h_inv_recovery_ref_st.csv")

# --- panels ----------------------------------------------------------------------------------------
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
LAB = {"T": r"ISI, $T_{\mathrm{osc}}$", "W": r"half-width, $W_{50}$", "A": "amplitude"}
MK = {"T": "o", "W": "s", "A": "^"}

def data_norm(a):
    for key, col in (("T", "ISI"), ("W", "W")):
        med = D[col].loc[doses].values / D[col][d0]
        err = np.array([med - D[col + "_q1"].loc[doses].values / D[col][d0], D[col + "_q3"].loc[doses].values / D[col][d0] - med])
        a.errorbar(doses, med, yerr=err, fmt=MK[key] + "-", color=COL[key], ms=10, capsize=5, elinewidth=1.8, label=LAB[key] + " (train medians)")
    # amplitude: within-cell step ratios chained from 5 µM (ΔF/F0 is not comparable between records)
    r = {(int(x.lo), int(x.hi)): x.amp_ratio for x in STEP.itertuples()}
    chain = {5: 1.0}; chain[10] = 1 / r.get((5, 10), np.nan); chain[15] = chain[10] / r.get((10, 15), np.nan); chain[30] = 1 / r.get((5, 30), np.nan)
    a.plot([5, 10, 15], [chain[5], chain[10], chain[15]], MK["A"] + "-", color=COL["A"], ms=11, label=r"amplitude (within-cell steps, 1 at 5 $\mu$M)")
    a.plot([5, 30], [chain[5], chain[30]], MK["A"] + "--", color=COL["A"], ms=11)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xticks(doses); a.set_xticklabels([str(int(d)) for d in doses]); a.minorticks_off()
    a.set_ylim(0.08, 3); a.set_yticks([0.1, 0.2, 0.5, 1, 2]); a.set_yticklabels(["0.1", "0.2", "0.5", "1", "2"])
    a.set_xlabel(r"carbachol ($\mu$M)"); a.set_ylabel("value / value at the lowest dose")
    a.set_title("Data: the spike hardly changes"); a.legend(loc="lower left", fontsize=17)

for m, name in ((prev, "prev"), (st, "st")):
    nm = normalised(m, 0.1)
    print(f"{name}: ratios S̄=0.1 / S̄=30 -> T {nm['T'][0]/nm['T'][-1]:.2f}, W50 {nm['W'][0]/nm['W'][-1]:.2f}, amp {nm['A'][0]/nm['A'][-1]:.2f}")

def model_norm(a):
    m = st; nm = normalised(m, 0.1)
    for key in ("T", "W", "A"):
        a.plot(m.value, nm[key], MK[key] + "-", color=COL[key], ms=8, label=LAB[key] if key != "A" else r"amplitude $C_{\max}-C_{\min}$")
    a.set_xscale("log"); a.set_yscale("log"); a.set_ylim(0.08, 3); a.set_yticks([0.1, 0.2, 0.5, 1, 2]); a.set_yticklabels(["0.1", "0.2", "0.5", "1", "2"])
    a.set_xlabel(r"stimulus $\bar S$"); a.set_ylabel(r"value / value at $\bar S=0.1$")
    a.legend(loc="lower left", fontsize=17); a.set_title("Model, stereotyped: only $T$ scales")

def ratios_qoff(a):
    for key, col in (("T", "T_ratio"), ("W", "W_ratio"), ("A", "amp_ratio")):
        a.plot(R.qoff, R[col], MK[key] + "-", color=COL[key], ms=9, label=LAB[key] + " (model)")
    a.axhline(D.ISI[d0] / D.ISI[30], color=COL["T"], ls=":", lw=2.5)
    a.axhline(W_530, color=COL["W"], ls=":", lw=2.5); a.axhline(AMP_530, color=COL["A"], ls=":", lw=2.5)
    a.text(0.02, D.ISI[d0] / D.ISI[30] * 1.08, r"data: $T$ 3/30 $\mu$M between cells; $W_{50}$, amp 5/30 $\mu$M within cells", color=C_GREY, fontsize=13)
    a.set_yscale("log"); a.set_ylim(0.4, 30); a.set_yticks([0.5, 1, 2, 5, 10, 20]); a.set_yticklabels(["0.5", "1", "2", "5", "10", "20"])
    a.set_xlabel(r"$q_{\mathrm{off}}$ (IP$_3$ acceleration on inactivation)"); a.set_ylabel(r"value($\bar S$=0.1) / value($\bar S$=30)")
    a.set_title(r"Ratios across the dose range"); a.legend(loc="upper right", fontsize=16)

def recovery(a):
    for rec, cm, lab in ((rec_prev, C_PREV, r"previous ($q_{\mathrm{off}}=1$)"), (rec_st, C_ST, r"stereotyped ($q_{\mathrm{off}}=0$)")):
        ev = rec[rec.evoked == 1]; ne = rec[rec.evoked == 0]
        a.plot(ev.tau_over_wait, ev.A_over_Afree, "o-", color=cm, ms=9, label=lab)
        a.plot(ne.tau_over_wait, ne.A_over_Afree, "x", color=cm, ms=9, mew=2)
    a.axhline(1.0, color="k", ls="--", lw=2, label=r"$H$-triggered spike ($H=H_{\mathrm{on}}$)")
    a.plot([], [], "kx", ms=9, mew=2, label="pulse fails to evoke")
    a.set_xlim(0, 1.0); a.set_ylim(0, 1.2); a.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.set_xlabel(r"delay after a spike, $\tau / T_{\mathrm{wait}}$"); a.set_ylabel(r"evoked amplitude $A/A_{\mathrm{free}}$")
    a.set_title(r"Ca-pulse-evoked spike at $\bar S=1$"); a.legend(loc="center left", fontsize=16)

def amp_vs_isi(a):
    sub = norm.sample(min(3000, len(norm)), random_state=0)
    a.plot(sub.isi, sub.amp, ".", color="#f0a8a0", ms=6, label="spikes (3000 shown)")
    edges = np.linspace(0.2, 2.2, 11); mids = 0.5 * (edges[1:] + edges[:-1])
    bm = [norm.amp[(norm.isi >= e0) & (norm.isi < e1)].median() for e0, e1 in zip(edges[:-1], edges[1:])]
    a.plot(mids, bm, "s-", color="k", ms=10, label="binned median")
    a.set_xlim(0, 2.5); a.set_ylim(0, 2.5); a.set_xlabel("preceding ISI / cell median"); a.set_ylabel("amplitude / cell median")
    a.set_title("Data: amplitude vs preceding ISI")
    a.text(0.03, 0.05, rf"corr per train: median {CORR_MED:.2f}" + "\n" + rf"IQR [{CORR_Q1:.2f}, {CORR_Q3:.2f}], $n$={len(trains)}", transform=a.transAxes, fontsize=16)
    a.legend(loc="upper right", fontsize=16)

# --- figures ----------------------------------------------------------------------------------------
MAIN = [data_norm, model_norm, ratios_qoff]     # IP3 acts on recovery, not on the spike
SUPP = [recovery, amp_vs_isi]                   # a Ca-triggered early spike would be smaller; the data show no such coupling
compose(MAIN, 3, "m0ip3h_amplitude_invariance_main")
compose(SUPP, 2, "m0ip3h_amplitude_invariance_supp")
print(R.round(3).to_string(index=False))
