"""
step_fold_change.py — fold-change encoding within cells (Proposition 2, gate 3): step records 5→10→15 µM
(03-07 E1, 03-14 E2, 03-19 E1), 5→10 (03-12 E1) and 5→30 (04-03 E1/E2), from the authors' spike tables (base files).

For each train (file, ST) with >= 5 ISIs in two consecutive stimulus levels:
  global ratio  r = T_av(level 1) / T_av(level 2)
  local ratio   r_loc = mean of the last 5 ISIs of level 1 / mean of the first 5 ISIs of level 2  (drift-minimising)
  within-level trend: relative ISI change per hour inside each level (desensitisation drift)
Proposition 2 predicts r independent of the cell's own T_av (all cells scale by the same factor); the constant-dose
medians (ISI 76 / 60 / 57 / 31 s at 5 / 10 / 15 / 30 µM) give the population prediction 1.27 (5→10), 1.05 (10→15), 2.5 (5→30).

Panels (functions; MAIN / SUPP at the end)
  example_cell   example cell: ISI vs time with the level boundaries
  per_cell       T_av at the higher level vs T_av at the lower level per cell, lines of constant ratio
  ratios         per-cell ratios (global and local) for 5→10, 10→15, 5→30 with the constant-dose prediction
Outputs: results/step_fold_change_cells.csv, results/step_fold_change_summary.csv,
figures/step_fold_change_main.{png,pdf,eps} (per_cell, ratios), figures/step_fold_change_supp.{png,pdf,eps} (example_cell)
"""
import glob, os
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
OUT, FIG = "results", "figures"; os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
C_DATA, C_LOC, C_GREY, TINT = "#c0392b", "#1f77b4", "#888888", "#f5d0cc"
STEP_COL = {"5→10": "#1f77b4", "10→15": "#2ca02c", "5→30": "#e69f00"}
ISI_MED = {3: 206.0, 5: 76.0, 10: 60.0, 15: 57.0, 30: 30.6}                 # constant-dose medians (m0 step)
NLOC = 5

ROOTS = ["../..", "../../data", "../data"]
def find(pattern):
    for r in ROOTS:
        hits = sorted(glob.glob(os.path.join(r, "*Spike Trains_data", pattern)))
        if hits: return hits
    return []
csvs = [f for f in find(os.path.join("2", "*.csv")) if "_" not in os.path.basename(f)]
if not csvs: raise SystemExit("spike tables (2/*.csv) not found under " + ", ".join(ROOTS))
P = pd.concat([pd.read_csv(f, sep=";", decimal=",", encoding="utf-8-sig").assign(file=os.path.basename(f)) for f in csvs])
step_files = [f for f, d in P.groupby("file") if d.stimulus.nunique() >= 2]
print("step records:", step_files)

rows, example = [], None
for (f, st), d in P[P.file.isin(step_files)].groupby(["file", "ST"]):
    d = d.sort_values("time"); levels = sorted(d.stimulus.unique())
    for lo, hi in zip(levels[:-1], levels[1:]):
        a, b = d[d.stimulus == lo], d[d.stimulus == hi]
        if len(a) < NLOC or len(b) < NLOC: continue
        Ta, Tb = a.ISI.mean(), b.ISI.mean()
        tr_a = np.polyfit(a.time, a.ISI, 1)[0] * 3600 / Ta if len(a) >= 3 else np.nan
        tr_b = np.polyfit(b.time, b.ISI, 1)[0] * 3600 / Tb if len(b) >= 3 else np.nan
        rows.append(dict(file=f, ST=st, step=f"{int(lo)}→{int(hi)}", lo=int(lo), hi=int(hi), n_lo=len(a), n_hi=len(b), T_lo=Ta, T_hi=Tb,
                         ratio=Ta / Tb, ratio_local=a.ISI.values[-NLOC:].mean() / b.ISI.values[:NLOC].mean(),
                         trend_lo=tr_a, trend_hi=tr_b, CV_lo=a.ISI.std(ddof=1) / Ta, CV_hi=b.ISI.std(ddof=1) / Tb))
    if example is None and len(d) >= 40 and d.stimulus.nunique() == 3: example = d
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "step_fold_change_cells.csv"), index=False)
# between-cell reference: constant-dose trains only (data_dose_units, fold_change_targets.csv); the 10 uM level exists only in step
# records, so 5->10 and 10->15 have no between-cell reference; fallback: the old pooled medians
BETWEEN = {(5, 10): ISI_MED[5] / ISI_MED[10], (10, 15): ISI_MED[10] / ISI_MED[15], (5, 30): ISI_MED[5] / ISI_MED[30]}
_ft = "../data_dose_units/results/fold_change_targets.csv"
if os.path.exists(_ft):
    _t = pd.read_csv(_ft); _t = _t[_t.kind.str.startswith("between")]
    BETWEEN = {(int(r.lo), int(r.hi)): float(r.ratio) for r in _t.itertuples()}
summ = []
for step, g in R.groupby("step"):
    pred = BETWEEN.get((int(g.lo.iloc[0]), int(g.hi.iloc[0])), np.nan)      # NaN when a level has no constant-dose records (10 uM)
    summ.append(dict(step=step, n_cells=len(g), ratio_median=g.ratio.median(), ratio_q1=g.ratio.quantile(.25), ratio_q3=g.ratio.quantile(.75),
                     ratio_cv=g.ratio.std(ddof=1) / g.ratio.mean(), Tlo_cv=g.T_lo.std(ddof=1) / g.T_lo.mean(),
                     local_median=g.ratio_local.median(), local_q1=g.ratio_local.quantile(.25), local_q3=g.ratio_local.quantile(.75),
                     pred_constant_dose=pred, trend_lo_median=g.trend_lo.median(), trend_hi_median=g.trend_hi.median(),
                     corr_ratio_Tlo=np.corrcoef(np.log(g.T_lo), np.log(g.ratio))[0, 1] if len(g) > 3 else np.nan))
S = pd.DataFrame(summ)
S["lo"] = [int(x.split("→")[0]) for x in S.step]; S["hi"] = [int(x.split("→")[1]) for x in S.step]
S = S.sort_values(["hi", "lo"]).reset_index(drop=True)      # 5→10, 10→15, 5→30
S.to_csv(os.path.join(OUT, "step_fold_change_summary.csv"), index=False); print(S.round(3).to_string(index=False))

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

def example_cell(a):
    if example is None: return
    for lv, col in zip(sorted(example.stimulus.unique()), ("#1f77b4", "#2ca02c", "#e69f00")):
        e = example[example.stimulus == lv]; a.plot(e.time / 60, e.ISI, "o-", color=col, ms=7, label=rf"{int(lv)} $\mu$M")
    a.set_xlabel(r"$t$ (min)"); a.set_ylabel("ISI (s)"); a.set_title(f"{example.file.iloc[0][:-4]}, ST {int(example.ST.iloc[0])}")
    a.legend(loc="upper right", fontsize=16)

def per_cell(a):
    for step, col in STEP_COL.items():
        g = R[R.step == step]
        if len(g): a.loglog(g.T_lo, g.T_hi, "o", color=col, ms=9, label=f"{step} (n={len(g)})")
    leg_pts = a.legend(loc="upper left", fontsize=15)          # first legend: the step points
    a.add_artist(leg_pts)                                        # keep it when the second legend is drawn
    xx = np.array([10, 600]); lines = []
    refs = [(1.0, ":", "ratio 1")] + [(BETWEEN[k], ls, f"between cells,\nconstant dose {k[0]}$\\to${k[1]}: {BETWEEN[k]:.2f}")
                                     for k, ls in (((5, 30), "-."), ((5, 10), "--"), ((10, 15), "--")) if k in BETWEEN and np.isfinite(BETWEEN[k])]
    for r_, ls, lab in refs:
        lines += a.loglog(xx, xx / r_, "k", ls=ls, lw=1.8, label=lab)
    a.legend(handles=lines, loc="lower right", fontsize=15)      # second legend: the constant-ratio lines
    a.set_xlim(10, 600); a.set_ylim(10, 600); a.set_xlabel(r"$T_{\mathrm{av}}$ at the lower dose (s)"); a.set_ylabel(r"$T_{\mathrm{av}}$ at the higher dose (s)")
    a.set_title("Fold change per cell")

def ratios(a):
    xs = np.arange(len(S)); w = 0.3
    a.bar(xs - w / 2, S.ratio_median, w, color=TINT, edgecolor=C_DATA, yerr=[S.ratio_median - S.ratio_q1, S.ratio_q3 - S.ratio_median], capsize=6, label="global ratio (median, IQR)")
    a.bar(xs + w / 2, S.local_median, w, color="white", edgecolor=C_LOC, lw=2.5, yerr=[S.local_median - S.local_q1, S.local_q3 - S.local_median], capsize=6, label="local ratio (5 ISIs each side)")
    a.plot(xs, S.pred_constant_dose, "k_", ms=40, mew=4, label="between cells (constant dose)")
    a.set_xticks(xs); a.set_xticklabels([s + r" $\mu$M" for s in S.step]); a.set_ylabel(r"$T_{\mathrm{av}}$(low) / $T_{\mathrm{av}}$(high)")
    a.set_ylim(0, max(4, S.ratio_q3.max() * 1.15)); a.set_title("Step ratios vs dose–response"); a.legend(loc="upper left", fontsize=15)

# --- figures ----------------------------------------------------------------------------------------
MAIN = [per_cell, ratios]      # the step records do not show fold-change encoding (gate 3)
SUPP = [example_cell]          # what a 3-h step record looks like (levels separated by silent intervals)
compose(MAIN, 2, "step_fold_change_main")
compose(SUPP, 1, "step_fold_change_supp", w=8.5)
