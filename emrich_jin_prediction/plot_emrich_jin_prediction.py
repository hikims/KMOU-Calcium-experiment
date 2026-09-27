"""
plot_emrich_jin_prediction.py — what M1 predicts when the basal influx Jin (SOCE) is reduced, against the
spike counts of Emrich et al. 2021 (WT, STIM1/2 KO, ORAI-TKO; 10 µM CCh, 14 min).

Inputs: results/jin_sweep_S_E.csv, jin_cells_f*.csv, jin_counts_summary.csv, jin_trace_f*.csv, jin_calibration.csv
        ../emrich_replication/results/emrich_counts_by_genotype.csv (fallback: WT 6 [3, 8], ORAI-TKO 3 [2, 5])
Main  (a) period, CV and amplitude vs Jin at the calibrated stimulus (stationary statistics; ratio to Jin = 0.02)
      (b) spike count per 14 min: model populations for Jin factors (median, IQR) vs Emrich genotypes
      (c) the Jin factor that reproduces the ORAI-TKO median, and the predicted amplitude and ISI-CV there
Supp  (a) traces (C(t)) at three Jin factors after the stimulus switch; (b) count histograms WT-like vs TKO-like;
      (c) first-spike latency vs Jin factor; (d) fraction of silent cells vs Jin factor
Writes results/jin_prediction_table.csv
"""
import os, glob
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
NAME, OUT, FIG = "emrich_jin_prediction", "results", "figures"; os.makedirs(FIG, exist_ok=True)
C_T, C_CV, C_A, C_MODEL, C_DATA, C_TINT, C_GREY = "#1f77b4", "#2ca02c", "#9467bd", "#2ca02c", "#c0392b", "#f5d0cc", "#888888"
JIN0 = 0.02
cal = pd.read_csv(os.path.join(OUT, "jin_calibration.csv")).iloc[0]; T0 = float(cal.T0_s)
JS = pd.read_csv(os.path.join(OUT, "jin_sweep_S_E.csv")); CS = pd.read_csv(os.path.join(OUT, "jin_counts_summary.csv")).sort_values("factor")
CELLS = {float(os.path.basename(f)[12:-4]): pd.read_csv(f) for f in glob.glob(os.path.join(OUT, "jin_cells_f*.csv"))}
TR = {float(os.path.basename(f)[12:-4]): pd.read_csv(f) for f in glob.glob(os.path.join(OUT, "jin_trace_f*.csv"))}
EM = pd.DataFrame([dict(genotype="WT HEK293", median=6, q1=3, q3=8), dict(genotype="STIM1 KO", median=6, q1=4, q3=7), dict(genotype="STIM2 KO", median=7, q1=4, q3=9),
                   dict(genotype="STIM1/2 KO", median=6, q1=4, q3=8), dict(genotype="ORAI-TKO", median=3, q1=2, q3=5)])
fn = "../emrich_replication/results/emrich_counts_by_genotype.csv"
if os.path.exists(fn): EM = pd.read_csv(fn)
TKO = float(EM[EM.genotype.str.contains("ORAI")]["median"].iloc[0]); WT = float(EM[EM.genotype.str.startswith("WT")]["median"].iloc[0])

# the factor reproducing the TKO median (interpolate the population median vs factor)
f_grid = np.linspace(0, 1, 401); med_i = np.interp(f_grid, CS.factor, CS["median"])
idx = np.where(med_i <= TKO)[0]; F_TKO = float(f_grid[idx.max()]) if len(idx) else np.nan
ref = JS.iloc[(JS.value - JIN0).abs().idxmin()]; at = JS.iloc[(JS.value - JIN0 * F_TKO).abs().idxmin()] if np.isfinite(F_TKO) else ref
table = pd.DataFrame([dict(S_E=cal.S_E, T0_s=T0, WT_median=WT, TKO_median=TKO, Jin_factor_for_TKO=F_TKO, Jin_TKO=JIN0 * F_TKO,
                           T_ratio_TKO_over_WT=at.Tdet / ref.Tdet, amp_ratio_TKO_over_WT=at.amp_det / ref.amp_det, CV_WT=ref.CV, CV_TKO=at.CV, W50_ratio=at.W50_det / ref.W50_det,
                           Jin_osc_min=float(JS[JS.n > 0].value.min()) if (JS.n > 0).any() else np.nan)])
table.to_csv(os.path.join(OUT, "jin_prediction_table.csv"), index=False); print(table.T.to_string())

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

def p_sweep(a):
    o = JS[JS.n > 0]; f = o.value / JIN0
    a.plot(f, o.Tdet / ref.Tdet, "o-", color=C_T, label=r"$T/T_{\rm ref}$"); a.plot(f, o.amp_det / ref.amp_det, "s-", color=C_A, label=r"amp$/$amp$_{\rm ref}$")
    a.plot(f, o.CV / ref.CV, "^-", color=C_CV, label=r"CV$/$CV$_{\rm ref}$")
    if np.isfinite(F_TKO): a.axvline(F_TKO, color=C_DATA, ls=":", lw=2, label=rf"interpolated factor {F_TKO:.3f} (illustrative)")
    a.axhline(1, color="k", lw=1, ls=":"); a.set_xlabel(r"$J_{\rm in}/J_{\rm in,ref}$"); a.set_ylabel("ratio to reference"); a.set_title(rf"Influx sweep at $\bar S_E={cal.S_E:.2f}$"); a.legend(loc="upper right", fontsize=15)

# direct population runs at the small factors (round2_checks/check3_jin.py): the interpolation of the median is not valid
_DIRECT = "../round2_checks/results/check3_jin_population.csv"
DIRECT = pd.read_csv(_DIRECT).groupby("factor").agg(median=("median", "median"), q1=("q1", "min"), q3=("q3", "max")).reset_index() if os.path.exists(_DIRECT) else None
def p_counts(a):
    """population median count vs J_in factor on a symlog axis: the stored sweep (circles), the direct runs at small
    factors (squares) and the Emrich WT / ORAI-TKO medians; the median never reaches the TKO value for J_in > 0"""
    a.errorbar(CS.factor, CS["median"], yerr=[CS["median"] - CS.q1, CS.q3 - CS["median"]], fmt="o", color=C_MODEL, ms=11, capsize=6, elinewidth=2, label="model sweep (median, IQR)")
    if DIRECT is not None:
        a.errorbar(DIRECT.factor, DIRECT["median"], yerr=[DIRECT["median"] - DIRECT.q1, DIRECT.q3 - DIRECT["median"]], fmt="s", color=C_T, ms=10, capsize=5, elinewidth=2, label="direct runs, small factors (median, IQR)")
    a.set_xscale("symlog", linthresh=0.01, linscale=0.5); a.set_xlim(-0.002, 1.2); a.set_xticks([0, 0.01, 0.03, 0.1, 0.3, 1]); a.set_xticklabels(["0", "0.01", "0.03", "0.1", "0.3", "1"], fontsize=17)
    a.set_xlabel(r"$J_{\rm in}/J_{\rm in,ref}$ (symlog)")
    for g, col, ls in (("WT", C_DATA, "-"), ("ORAI", C_DATA, "--")):
        r = EM[EM.genotype.str.contains(g)].iloc[0]; a.axhspan(r.q1, r.q3, color=C_TINT, alpha=0.5 if g == "WT" else 0.3, zorder=0); a.axhline(r["median"], color=col, ls=ls, lw=2.5, label=f"Emrich {r.genotype} (median, IQR)")
    a.set_ylabel("spikes per 14 min"); a.set_title("Counts after the stimulus switch"); a.legend(loc="lower right", fontsize=14)

_RATIOS = "../round2_checks/results/check3_jin_ratios.csv"
RATIOS = pd.read_csv(_RATIOS).iloc[0] if os.path.exists(_RATIOS) else None
def p_prediction(a):
    """ratios at the illustrative factor 0.075 from the direct stationary runs (round2_checks/check3_jin.py);
    falls back to the sweep interpolation if the direct file is absent"""
    labs = [r"$T$", "amplitude", r"$W_{50}$", "ISI CV"]
    if RATIOS is not None: vals = [RATIOS.ratio_Tdet, RATIOS.ratio_amp_det, RATIOS.ratio_W50_det, RATIOS.ratio_CV]; src = "direct runs"
    else: vals = [at.Tdet / ref.Tdet, at.amp_det / ref.amp_det, at.W50_det / ref.W50_det, at.CV / ref.CV]; src = "sweep"
    a.bar(range(4), vals, color=[C_T, C_A, C_GREY, C_CV]); a.axhline(1, color="k", lw=1.5, ls=":")
    for i, v in enumerate(vals): a.text(i, v + 0.03, f"{v:.2f}", ha="center", va="bottom", fontsize=18)
    a.set_xticks(range(4)); a.set_xticklabels(labs); a.set_ylabel("ratio to reference"); a.set_ylim(0, max(vals) * 1.25)
    a.set_title(rf"Ratios at $J_{{\rm in}}$ factor 0.075 ({src}; illustrative)")

def p_traces(a):
    fs = [f for f in (1.0, 0.5, 0.25) if f in TR]; off = 0
    for f in fs:
        t = TR[f]; a.plot(t.t * T0 / 60, t.C / t.C.max() + off, "-", color=C_MODEL if f == 1 else C_T, lw=1.5); a.text(3.0, off + 0.75, rf"$J_{{\rm in}}\times{f:g}$", fontsize=16); off += 1.15
    a.set_xlabel(r"$t$ after CCh (min)"); a.set_ylabel(r"$C/C_{\max}$ (offset)"); a.set_yticks([]); a.set_title("Model traces"); a.grid(axis="y", visible=False)

def p_hist(a):
    bins = np.arange(-0.5, 20.5)
    for f, col, lab in ((1.0, C_MODEL, "WT-like (factor 1)"), (round(F_TKO, 2) if np.isfinite(F_TKO) else None, C_DATA, "TKO-like")):
        key = min(CELLS, key=lambda k: abs(k - f)) if f is not None else None
        if key is not None: a.hist(CELLS[key].n_spikes, bins=bins, color=col, alpha=0.55, label=lab + rf" (factor {key:g})")
    a.set_xlabel("spikes per 14 min"); a.set_ylabel("cells"); a.set_title("Count distributions"); a.legend(loc="upper right", fontsize=15)

def p_latency(a):
    fs = sorted(CELLS); med = [CELLS[f][CELLS[f].t_first >= 0].t_first.median() * T0 for f in fs]; q1 = [CELLS[f][CELLS[f].t_first >= 0].t_first.quantile(.25) * T0 for f in fs]; q3 = [CELLS[f][CELLS[f].t_first >= 0].t_first.quantile(.75) * T0 for f in fs]
    a.errorbar(fs, med, yerr=[np.array(med) - q1, np.array(q3) - med], fmt="o-", color=C_T, ms=10, capsize=5); a.set_xlabel(r"$J_{\rm in}/J_{\rm in,ref}$"); a.set_ylabel("first-spike latency (s)"); a.set_title("Latency")

def p_silent(a):
    a.plot(CS.factor, 100 * CS.frac_silent, "o-", color=C_DATA, ms=10); a.set_xlabel(r"$J_{\rm in}/J_{\rm in,ref}$"); a.set_ylabel("silent cells (%)"); a.set_title("Cells without spikes")

MAIN = [p_sweep, p_counts, p_prediction]
SUPP = [p_traces, p_hist, p_latency, p_silent]
if __name__ == "__main__":
    compose(MAIN, 3, f"{NAME}_main")
    compose(SUPP, 2, f"{NAME}_supp")
