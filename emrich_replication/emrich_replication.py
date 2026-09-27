"""
emrich_replication.py — independent replication with the Emrich et al. 2021 (Cell Reports) source data
(HEK293 WT and STIM/Orai knockouts, Fura-2 ratio, 10 µM carbachol, 2 s sampling, 14–25 min records).

What the workbook contains (Source_Data_Emrich.xlsx, one sheet per figure):
  - representative single-cell Fura-2 traces: 5 cells per panel, stored as rows without a time column
    (451 points = 15 min at 2 s; 631 = 21 min (3→10→30→90 µM CCh, 5 min each); 781 = 26 min)
  - per-cell counts "total oscillations in 14 min" (Figure 1I: WT 189, STIM1 KO 230, STIM2 KO 156,
    STIM1/2 KO 152, Orai-TKO 96 cells) — population frequency, no ISI variability
  - SOCE (thapsigargin) and NFAT time courses (Figure 1B, 3A/F/K, 7): not oscillations, not used
Only the representative traces carry ISI information, so everything below is a replication on
~230 selected cells (selection bias acknowledged), not an unbiased population.

Analyses
  1. within-cell statistics per train (>= 8 spikes after the initial transient): CV, skewness, amplitude CV,
     corr(preceding ISI, amplitude), serial ISI correlation — raw and after removing a linear trend in time,
     because these 14-min records are strongly non-stationary (oscillations slow down or accelerate)
  2. dose steps 3→10→30→90 µM within single cells (Figure 1M–P, 20 cells): fold change 10→30 µM without washout
  3. population frequency by genotype from the per-cell oscillation counts (Figure 1I)
Comparison values from our HEK293 analysis (Azzoni 2025 data; m1_isi_variability): CV 0.22, skew 0.52,
amplitude CV 0.15, corr(ISI, amp) 0.10, serial corr 0.03 (detrended); dose ratio 10→30 µM: 60/31 = 1.9.

Figures
  main (1x3): (a) representative WT / STIM1/2-KO / Orai-TKO traces, (b) within-cell statistics (detrended)
              vs our HEK293 values, (c) within-cell 10→30 µM fold change vs the population prediction
  supp (1x2): (a) oscillations per cell in 14 min by genotype (Figure 1I), (b) raw vs detrended statistics
Usage: python3 emrich_replication.py [path/to/Source_Data_Emrich.xlsx]
       (default: looks for ../../*signaling_data/Source_Data_Emrich.xlsx, i.e. the Mendeley folder next to code/)
Outputs: results/emrich_traces.csv, results/emrich_train_stats.csv, results/emrich_dose_steps.csv,
         results/emrich_counts_by_genotype.csv, figures/emrich_replication_{main,supp}.{png,pdf,eps}
"""
import glob, os, re, sys
import numpy as np, pandas as pd
import openpyxl
from scipy.signal import find_peaks
from scipy.stats import skew
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
NAME, OUT, FIG = "emrich_replication", "results", "figures"
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
C_OURS, C_EMR, C_GREY, TINT = "#c0392b", "#1f77b4", "#888888", "#f5d0cc"
DT, T_CCH, MIN_SPIKES = 2.0, 60.0, 8                     # s; CCh added at 1 min; minimum spikes per train
OURS = dict(CV=0.22, skew=0.52, amp_cv=0.15, corr=0.10, serial=0.03)      # our HEK293 (detrended serial corr)
OURS_RATIO_10_30 = 60.0 / 30.6
STATS = ["CV", "skew", "amp_cv", "corr", "serial"]
LABELS = ["CV", "skew", "amp CV", "corr(ISI, amp)", "serial corr"]

# --- locate the workbook ---------------------------------------------------------------------------------
# layout: <project>/code/emrich_replication/  (this script)
#         <project>/Omnitemporal choreographies ... Ca2+ signaling_data/Source_Data_Emrich.xlsx
ROOTS = ["../..", "../../data", "../data", "."]
FN = "Source*Data*Emrich*.xlsx"          # the Mendeley file is "Source Data_Emrich.xlsx" (with a space); "Source_Data_Emrich.xlsx" also accepted
def find_xlsx():
    for r in ROOTS:
        for pat in (os.path.join(r, "*signaling_data", FN), os.path.join(r, "*Emrich*", FN), os.path.join(r, "*Omnitemporal*", FN), os.path.join(r, FN)):
            hits = sorted(glob.glob(pat))
            if hits: return hits[0]
    for r in ROOTS:                                                    # last resort: recursive search
        hits = sorted(glob.glob(os.path.join(r, "**", FN), recursive=True))
        if hits: return hits[0]
    return None
XLSX = sys.argv[1] if len(sys.argv) > 1 else find_xlsx()
if XLSX is None or not os.path.exists(XLSX):
    raise SystemExit("Source Data_Emrich.xlsx not found (searched " + ", ".join(os.path.abspath(r) for r in ROOTS) +
                     ").  Run from the emrich_replication folder, or pass the path:  python3 emrich_replication.py /path/to/Source_Data_Emrich.xlsx")
print("workbook:", XLSX)

# --- panel -> condition map (figure legends of Emrich et al. 2021) ---------------------------------------
PANELS = {
    "Figure 1D": ("WT", "10uM"), "Figure 1E": ("STIM1 KO", "10uM"), "Figure 1F": ("STIM2 KO", "10uM"), "Figure 1G": ("STIM1/2 KO", "10uM"), "Figure 1H": ("ORAI-TKO", "10uM"),
    "Figure 1M": ("WT", "steps"), "Figure 1N": ("STIM1 KO", "steps"), "Figure 1O": ("STIM2 KO", "steps"), "Figure 1P": ("STIM1/2 KO", "steps"),
    "Figure 2A": ("STIM1/2 KO", "10uM"), "Figure 2B": ("tk-STIM1 rescue", "10uM"), "Figure 2C": ("tk-STIM2 rescue", "10uM"), "Figure 2D": ("tk-S1N-S2C", "10uM"), "Figure 2E": ("tk-S2N-S1C", "10uM"),
    "Figure 3C": ("ORAI2,3 DKO", "10uM"), "Figure 3D": ("ORAI2,3-STIM1 KO", "10uM"), "Figure 3E": ("ORAI2,3-STIM2 KO", "10uM"),
    "Figure 3H": ("ORAI1,2 DKO", "10uM"), "Figure 3I": ("ORAI1,2-STIM1 KO", "10uM"), "Figure 3J": ("ORAI1,2-STIM2 KO", "10uM"),
    "Figure 3M": ("ORAI1,3 DKO", "10uM"), "Figure 3N": ("ORAI1,3-STIM1 KO", "10uM"), "Figure 3O": ("ORAI1,3-STIM2 KO", "10uM"),
    "Figure 4A": ("WT +1mM Gd", "10uM"), "Figure 4B": ("STIM1 KO +1mM Gd", "10uM"), "Figure 4C": ("STIM2 KO +1mM Gd", "10uM"), "Figure 4D": ("STIM1/2 KO +1mM Gd", "10uM"),
    "Figure 4L": ("STIM1/2 KO", "10uM"), "Figure 4M": ("STIM2.1 rescue", "10uM"), "Figure 4N": ("STIM1-F394H rescue", "10uM"), "Figure 4O": ("STIM2.1+F394H rescue", "10uM"),
    "Figure 4S": ("STIM1/2 KO", "10uM"), "Figure 4T": ("STIM2.1 EF mutant", "10uM"), "Figure 4U": ("STIM1-F394H EF mutant", "10uM"), "Figure 4V": ("combined EF mutants", "10uM"),
    "Sup Figure 1L": ("WT", "10uM"), "Sup Figure 1M": ("STIM1 KO", "10uM"), "Sup Figure 1N": ("STIM2 KO", "10uM"), "Sup Figure 1O": ("STIM1/2 KO", "10uM"),
    "Sup Figure 2A": ("WT +5uM Gd", "10uM"), "Sup Figure 2B1": ("STIM1 KO cl1 +5uM Gd", "10uM"), "Sup Figure 2B2": ("STIM1 KO cl4 +5uM Gd", "10uM"),
    "Sup Figure 2C1": ("STIM2 KO cl6 +5uM Gd", "10uM"), "Sup Figure 2C2": ("STIM2 KO cl14 +5uM Gd", "10uM"), "Sup Figure 2D1": ("STIM1/2 KO cl1 +5uM Gd", "10uM"), "Sup Figure 2D2": ("STIM1/2 KO cl20 +5uM Gd", "10uM"),
}
STEP_LEVELS = [(3, 60, 360), (10, 360, 660), (30, 660, 960), (90, 960, 1260)]      # µM, start s, end s (5 min each after 1 min)

# --- parse the workbook ----------------------------------------------------------------------------------
def isnum(v): return isinstance(v, (int, float)) and not isinstance(v, bool)
wb = openpyxl.load_workbook(XLSX, read_only=True)
traces, counts = [], {}
for ws in wb.worksheets:
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    labs = [(i, v.strip()) for i, r in enumerate(rows) for v in r if isinstance(v, str) and re.match(r"(Sup )?Figure\s*\d", v.strip())]
    for n, (i, lab) in enumerate(labs):
        end = labs[n + 1][0] if n + 1 < len(labs) else len(rows)
        if lab in PANELS:
            k = 0
            for r in rows[i + 1:end]:
                y = np.array([v for v in r if isnum(v)], float)
                if len(y) >= 200: traces.append(dict(panel=lab, genotype=PANELS[lab][0], protocol=PANELS[lab][1], cell=k, y=y)); k += 1
        if lab == "Figure 1I":                          # per-cell oscillation counts by genotype
            for r in rows[i + 1:end]:
                name = next((v for v in r if isinstance(v, str)), None); vals = [v for v in r if isnum(v)]
                if name and vals: counts[name.strip()] = np.array(vals, float)
print(f"representative traces: {len(traces)} from {len(set(t['panel'] for t in traces))} panels; count groups: {list(counts)}")
pd.DataFrame([dict(panel=t["panel"], genotype=t["genotype"], protocol=t["protocol"], cell=t["cell"], t_s=k * DT, ratio=v)
              for t in traces for k, v in enumerate(t["y"])]).to_csv(os.path.join(OUT, "emrich_traces.csv"), index=False)
CNT = pd.DataFrame([dict(genotype=g, n_cells=len(v), median=np.median(v), q1=np.quantile(v, .25), q3=np.quantile(v, .75), mean=v.mean()) for g, v in counts.items()])
CNT.to_csv(os.path.join(OUT, "emrich_counts_by_genotype.csv"), index=False)

# --- spike detection on ΔR/R0 ----------------------------------------------------------------------------------
def spikes(y, t0=T_CCH, t1=None):
    """peaks of ΔR/R0 (R0 = median of the first 50 s, before CCh) inside [t0, t1]; returns t_peak, amplitude"""
    t = np.arange(len(y)) * DT; R0 = np.median(y[:int(50 / DT)]); r = (y - R0) / R0
    d = np.diff(r); sig = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)
    pk, pr = find_peaks(r, prominence=max(0.15, 5 * sig), distance=max(1, int(10 / DT)))
    m = (t[pk] > t0) & (t[pk] < (t1 if t1 else t[-1]))
    return t[pk][m], pr["prominences"][m], t, r

def train_stats(tp, amp):
    x = np.diff(tp); a = amp[1:]
    if len(x) < MIN_SPIKES - 1: return None
    tm = tp[1:]; res = x - np.polyval(np.polyfit(tm, x, 1), tm); ares = a - np.polyval(np.polyfit(tm, a, 1), tm)
    return dict(n=len(tp), Tav=x.mean(), CV=x.std(ddof=1) / x.mean(), skew=skew(x), amp_cv=a.std(ddof=1) / a.mean(),
                corr=np.corrcoef(x, a)[0, 1], serial=np.corrcoef(x[:-1], x[1:])[0, 1], trend_per_min=np.polyfit(tm, x, 1)[0] * 60 / x.mean(),
                CV_detr=res.std(ddof=1) / x.mean(), skew_detr=skew(res), amp_cv_detr=ares.std(ddof=1) / a.mean(),
                corr_detr=np.corrcoef(res, ares)[0, 1], serial_detr=np.corrcoef(res[:-1], res[1:])[0, 1])

rows, steps = [], []
for tr in traces:
    if tr["protocol"] == "10uM":
        tp, amp, t, r = spikes(tr["y"], t0=T_CCH + 30)                    # skip the initial transient (first 30 s after CCh)
        st = train_stats(tp, amp)
        if st: rows.append(dict(panel=tr["panel"], genotype=tr["genotype"], cell=tr["cell"], **st))
    else:                                                                # dose steps: per level, skip the first 30 s
        d = dict(panel=tr["panel"], genotype=tr["genotype"], cell=tr["cell"])
        for dose, a0, a1 in STEP_LEVELS:
            tp, amp, t, r = spikes(tr["y"], t0=a0 + 30, t1=a1)
            d[f"n_{dose}"] = len(tp); d[f"Tav_{dose}"] = np.diff(tp).mean() if len(tp) >= 3 else np.nan
            d[f"amp_{dose}"] = amp.mean() if len(tp) >= 3 else np.nan
        steps.append(d)
S = pd.DataFrame(rows); S.to_csv(os.path.join(OUT, "emrich_train_stats.csv"), index=False)
ST = pd.DataFrame(steps); ST["ratio_10_30"] = ST.Tav_10 / ST.Tav_30; ST["ratio_3_10"] = ST.Tav_3 / ST.Tav_10; ST["amp_ratio_10_30"] = ST.amp_10 / ST.amp_30
ST.to_csv(os.path.join(OUT, "emrich_dose_steps.csv"), index=False)
med = {s: S[s].median() for s in STATS}; medd = {s: S[s + "_detr"].median() for s in STATS}
q = {s: (S[s + "_detr"].quantile(.25), S[s + "_detr"].quantile(.75)) for s in STATS}
print(f"trains with >= {MIN_SPIKES} spikes: {len(S)}; median trend {S.trend_per_min.median():+.2f} (rel. ISI change per min)")
print("raw      :", {s: round(med[s], 3) for s in STATS}); print("detrended:", {s: round(medd[s], 3) for s in STATS}); print("ours     :", OURS)
r1030 = ST.ratio_10_30.dropna(); print(f"dose steps: {len(r1030)} cells with 10 and 30 uM trains; ratio 10/30 median {r1030.median():.2f} [{r1030.quantile(.25):.2f}, {r1030.quantile(.75):.2f}] vs our population {OURS_RATIO_10_30:.2f}; amp ratio median {ST.amp_ratio_10_30.median():.2f}")
print(CNT.round(2).to_string(index=False))

# --- figures ----------------------------------------------------------------------------------------------
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

def p_examples(a):
    picks = [("Figure 1D", 1, "WT"), ("Figure 1G", 1, "STIM1/2 KO"), ("Figure 1H", 1, "Orai-TKO")]
    for j, (pan, c, lab) in enumerate(picks):
        tr = next(t for t in traces if t["panel"] == pan and t["cell"] == c); tp, amp, t, r = spikes(tr["y"], t0=T_CCH + 30)
        a.plot(t / 60, r + 3.5 * (2 - j), color=C_GREY, lw=1.5); a.plot(tp / 60, r[np.searchsorted(t, tp)] + 3.5 * (2 - j), "v", color=C_EMR, ms=8)
        a.text(14.8, 3.5 * (2 - j) + 2.7, lab + r", 10 $\mu$M CCh", ha="right", fontsize=17)
    a.set_yticks([]); a.set_xlim(0, 15); a.set_xlabel(r"$t$ (min)"); a.set_ylabel(r"$\Delta R/R_0$ (offset)"); a.set_title("Emrich 2021: Fura-2 traces")

def p_stats(a):
    xs = np.arange(len(STATS)); w = 0.35
    a.bar(xs - w / 2, [medd[s] for s in STATS], w, color=C_EMR, yerr=[[medd[s] - q[s][0] for s in STATS], [q[s][1] - medd[s] for s in STATS]], capsize=5, label=f"Emrich, detrended (n={len(S)})")
    a.bar(xs + w / 2, [OURS[s] for s in STATS], w, color=TINT, edgecolor=C_OURS, lw=2, label="our HEK293 (Azzoni data)")
    a.axhline(0, color="k", lw=1); a.set_xticks(xs); a.set_xticklabels(LABELS, fontsize=17); a.set_ylim(-0.6, 1.6); a.set_ylabel("median (IQR)")
    a.set_title("Within-cell statistics"); a.legend(loc="upper right", fontsize=15)

def p_steps(a):
    ok = ST.dropna(subset=["Tav_10", "Tav_30"])
    for g, col in zip(sorted(ok.genotype.unique()), ("#1f77b4", "#2ca02c", "#e69f00", "#9467bd")):
        s = ok[ok.genotype == g]; a.plot(s.Tav_10, s.Tav_30, "o", color=col, ms=10, label=f"{g} (n={len(s)})")
    xx = np.array([5, 120]); a.plot(xx, xx, "k:", lw=1.8, label="ratio 1"); a.plot(xx, xx / OURS_RATIO_10_30, "k--", lw=1.8, label=f"ratio {OURS_RATIO_10_30:.1f} (our medians)")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(5, 120); a.set_ylim(5, 120)
    a.set_xlabel(r"$T_{\mathrm{av}}$ at 10 $\mu$M (s)"); a.set_ylabel(r"$T_{\mathrm{av}}$ at 30 $\mu$M (s)"); a.set_title(r"Step 10$\to$30 $\mu$M within cells"); a.legend(loc="upper left", fontsize=14)

def p_counts(a):
    order = [g for g in ["WT HEK293", "STIM1 KO", "STIM2 KO", "STIM1/2 KO", "ORAI-TKO"] if g in counts]
    a.boxplot([counts[g] for g in order], showfliers=False, widths=0.6)
    a.set_xticks(range(1, len(order) + 1)); a.set_xticklabels(order, rotation=20, fontsize=15); a.set_ylabel("oscillations per cell in 14 min")
    a.set_title(f"Population frequency (n = {', '.join(str(len(counts[g])) for g in order)})")

def p_raw_vs_detr(a):
    xs = np.arange(len(STATS)); w = 0.35
    a.bar(xs - w / 2, [med[s] for s in STATS], w, color=C_GREY, label="raw")
    a.bar(xs + w / 2, [medd[s] for s in STATS], w, color=C_EMR, label="linear trend removed")
    a.axhline(0, color="k", lw=1); a.set_xticks(xs); a.set_xticklabels(LABELS, fontsize=17); a.set_ylabel("median over trains")
    a.set_title(f"Non-stationarity (trend {S.trend_per_min.median():+.0%}/min)"); a.legend(loc="upper right", fontsize=15)

MAIN = [p_examples, p_stats, p_steps]
SUPP = [p_counts, p_raw_vs_detr]
compose(MAIN, 3, f"{NAME}_main"); compose(SUPP, 2, f"{NAME}_supp")
