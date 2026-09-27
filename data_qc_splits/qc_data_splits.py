"""
qc_data_splits.py — data dictionary, QC rules, independent spike detection vs the authors' spike tables, and the
fixed 15 µM development / lock split (plan v3.1, sections 4–5).  Python only (no model).

Inputs (Zenodo 15705673 raw TXT, UTF-16, 'Time [s]' + ROI columns '#N (...)'; Zenodo 14779474 spike tables, base files only)
Steps
  1. record table: sampling interval, gaps > 10 s (count, max), duration, ROIs, stimulus segments (from the spike
     tables: per stimulus level [first spike − 60 s, last spike + 60 s] over all trains of the record)
  2. per ROI: ΔF/F0 with F0 = median of the pre-stimulus window (≥ 60 s) else the 10th percentile; independent spike
     detection (find_peaks, prominence ≥ max(0.15, 5 σ_noise), distance ≥ 10 s; peaks that straddle a gap are dropped);
     SNR = median prominence / σ_noise; baseline drift = |slope of the 10th-percentile envelope| in %/h, evaluated
     inside the stimulus segments after their first 120 s (onset transient excluded); saturation flag
  3. ROI ↔ ST: the authors' ST index equals the ROI number; verified by the spike-time overlap (±5 s) inside the
     authors' time window (their tables keep only regular trains and often only a segment in time)
  4. QC flags (independent of the authors' selection): n_spikes < 10, SNR < 5, drift > 10 %/h, saturation;
     table_mismatch flags ROIs whose table spikes are not found by the detector
  5. 15 µM split: development = 07-11 E1, 07-18 E1, 08-08 E1/E2; lock = 08-09 E1/E2/E3 (fixed before any fit)
Outputs: results/qc_records.csv, results/qc_rois.csv, results/qc_split_15uM.csv, results/qc_exclusions.csv,
         figures/data_qc_splits_main.{png,pdf,eps} (example trace, QC outcome by date),
         figures/data_qc_splits_supp.{png,pdf,eps} (detection vs authors' tables, SNR–drift thresholds)
"""
import glob, io, os
import numpy as np, pandas as pd
from scipy.signal import find_peaks
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
C_DATA, C_ME, C_GREY, TINT = "#c0392b", "#1f77b4", "#888888", "#f5d0cc"
GAP_S, MIN_SPIKES, MIN_SNR, MAX_DRIFT, MATCH_TOL, MIN_MATCH = 10.0, 10, 5.0, 10.0, 5.0, 0.8
DEV = ["CChCB_2024_07_11_E1", "CChCB_2024_07_18_E1", "CChCB_2024_08_08_E1", "CChCB_2024_08_08_E2"]
LOCK = ["CChCB_2024_08_09_E1", "CChCB_2024_08_09_E2", "CChCB_2024_08_09_E3"]

ROOTS = ["../..", "../../data", "../data"]
def find(pattern):
    for r in ROOTS:
        hits = sorted(glob.glob(os.path.join(r, "*Spike Trains_data", pattern)))
        if hits: return hits
    return []
raw_files = find(os.path.join("1", "*.txt")) or find(os.path.join("CChCB", "*.txt"))
csv_files = [f for f in find(os.path.join("2", "*.csv")) if "_" not in os.path.basename(f)]
if not raw_files or not csv_files: raise SystemExit("raw TXT (1/ or CChCB/) and spike tables (2/) not found under " + ", ".join(ROOTS))
def rec_of_csv(f):        # 20240307E1.csv -> CChCB_2024_03_07_E1
    b = os.path.basename(f)[:-4]; return f"CChCB_{b[:4]}_{b[4:6]}_{b[6:8]}_{b[8:]}"
TAB = {rec_of_csv(f): pd.read_csv(f, sep=";", decimal=",", encoding="utf-8-sig") for f in csv_files}

def load_raw(f):
    df = pd.read_csv(io.StringIO(open(f, "rb").read().decode("utf-16")), sep="\t").dropna(axis=1, how="all")
    t = df["Time [s]"].values.astype(float); rois = [c for c in df.columns if c.startswith("#")]
    return t, {c.split(" ")[0]: pd.to_numeric(df[c], errors="coerce").values for c in rois}

def segments(tab):
    """stimulus segments of a record from its spike table"""
    seg = []
    for dose, d in tab.groupby("stimulus"): seg.append((int(dose), d.time.min() - 60, d.time.max() + 60))
    return sorted(seg, key=lambda s: s[1])

def detect(t, F, seg, gaps):
    """ΔF/F0 and independent spike detection; returns y, peak indices, prominences, sigma_noise, F0"""
    t0 = seg[0][1] if seg else 900.0
    pre = t < t0 - 30
    F0 = np.nanmedian(F[pre]) if pre.sum() * np.median(np.diff(t)) >= 60 else np.nanpercentile(F, 10)
    y = (F - F0) / F0
    d = np.diff(y); sig = 1.4826 * np.nanmedian(np.abs(d - np.nanmedian(d))) / np.sqrt(2)      # noise from first differences
    dt = np.median(np.diff(t)); dist = max(1, int(10.0 / dt))
    pk, pr = find_peaks(np.nan_to_num(y), prominence=max(0.15, 5 * sig), distance=dist)
    prom = pr["prominences"]
    # drop peaks whose flanks straddle a gap (within one half-width window of a gap edge)
    ok = np.ones(len(pk), bool)
    for g0, g1 in gaps:
        ok &= ~((t[pk] > g0 - 15) & (t[pk] < g1 + 15))
    return y, pk[ok], prom[ok], sig, F0

def drift_pct_per_h(t, y):
    """|slope| of the running 10th-percentile baseline (windows of 300 s), in % of (1+F/F0) per hour"""
    edges = np.arange(t[0], t[-1], 300.0)
    if len(edges) < 3: return np.nan
    base = [np.nanpercentile(y[(t >= e) & (t < e + 300)], 10) for e in edges[:-1]]
    mids = edges[:-1] + 150; ok = ~np.isnan(base)
    if ok.sum() < 3: return np.nan
    return abs(np.polyfit(mids[ok], np.array(base)[ok], 1)[0]) * 3600 * 100

rec_rows, roi_rows, examples = [], [], {}
for f in raw_files:
    rec = os.path.basename(f)[:-4]; t, R = load_raw(f); tab = TAB.get(rec)
    g = np.diff(t); gaps = [(t[i], t[i + 1]) for i in np.where(g > GAP_S)[0]]
    seg = segments(tab) if tab is not None else []
    sat = max(np.nanmax(v) for v in R.values())
    rec_rows.append(dict(record=rec, date=rec[6:16], n_samples=len(t), dt_s=np.median(g), duration_s=t[-1], n_gaps=len(gaps),
                         max_gap_s=g.max(), n_rois=len(R), doses="+".join(str(s[0]) for s in seg), n_trains=tab.ST.nunique() if tab is not None else 0,
                         segments=";".join(f"{s[0]}:{s[1]:.0f}-{s[2]:.0f}" for s in seg), F_max=sat))
    # authors' trains of this record
    trains = {st: d.sort_values("time") for st, d in tab.groupby("ST")} if tab is not None else {}
    for roi, F in R.items():
        y, pk, prom, sig, F0 = detect(t, F, seg, gaps)
        tp = t[pk]
        # mapping: the authors' ST index equals the ROI number (verified below by the spike-time overlap inside the
        # authors' time window; their tables keep only regular trains, so a missing ST is not a QC failure)
        st = int(roi[1:]); d = trains.get(st); bst = st if d is not None else None
        m_mine = m_auth = np.nan; n_win = 0
        if d is not None and len(tp):
            ta = d.time.values; win = (tp >= ta.min() - MATCH_TOL) & (tp <= ta.max() + MATCH_TOL); tw = tp[win]; n_win = int(win.sum())
            m_mine = np.mean([np.min(np.abs(ta - x)) <= MATCH_TOL for x in tw]) if len(tw) else 0.0
            m_auth = np.mean([np.min(np.abs(tp - x)) <= MATCH_TOL for x in ta])
        # drift evaluated inside the stimulus segments only, skipping the first 120 s of each (onset transient)
        msk = np.zeros(len(t), bool)
        for _, s0, s1 in seg: msk |= (t >= s0 + 120) & (t <= s1)
        drift = drift_pct_per_h(t[msk], y[msk]) if msk.sum() > 10 else drift_pct_per_h(t, y)
        snr = float(np.median(prom) / sig) if len(prom) else 0.0
        reasons = []
        if len(pk) < MIN_SPIKES: reasons.append("n_spikes<10")
        if snr < MIN_SNR: reasons.append("SNR<5")
        if not np.isnan(drift) and drift > MAX_DRIFT: reasons.append("drift>10%/h")
        if np.nanmax(F) >= 65000: reasons.append("saturation")
        if bst is not None and m_auth < MIN_MATCH: reasons.append("table_mismatch")
        dose = "+".join(str(s[0]) for s in seg)
        roi_rows.append(dict(record=rec, date=rec[6:16], roi=roi, dose=dose, n_spikes_mine=len(pk), in_author_table=bst is not None,
                             n_spikes_author=len(d) if d is not None else 0, n_mine_in_window=n_win, match_mine=m_mine, match_author=m_auth,
                             SNR=snr, sigma_noise=sig, F0=F0, drift_pct_h=drift, include=len(reasons) == 0, reasons=";".join(reasons)))
        if rec == "CChCB_2024_08_08_E1" and roi == "#18": examples["trace"] = (t, y, pk, d)
REC = pd.DataFrame(rec_rows); ROI = pd.DataFrame(roi_rows)
REC.to_csv(os.path.join(OUT, "qc_records.csv"), index=False); ROI.to_csv(os.path.join(OUT, "qc_rois.csv"), index=False)
EXC = ROI[~ROI.include][["record", "roi", "reasons"]]; EXC.to_csv(os.path.join(OUT, "qc_exclusions.csv"), index=False)

# 15 uM split
r15 = ROI[ROI.record.isin(DEV + LOCK)].copy(); r15["set"] = np.where(r15.record.isin(DEV), "development", "lock")
SPL = r15.groupby(["set", "record"]).agg(n_rois=("roi", "size"), n_included=("include", "sum"), n_in_author_table=("in_author_table", "sum")).reset_index()
SPL.to_csv(os.path.join(OUT, "qc_split_15uM.csv"), index=False)
print(REC[["record", "dt_s", "n_gaps", "max_gap_s", "n_rois", "doses", "n_trains"]].to_string(index=False))
print(f"\nROIs: {len(ROI)}, included {ROI.include.sum()}; exclusion reasons:\n{ROI[~ROI.include].reasons.value_counts().to_string()}")
print("\nROIs in the authors' tables: %d; table spikes found here (median) %.2f; my peaks inside their window found in the table (median) %.2f"
      % (ROI.in_author_table.sum(), ROI.match_author.median(), ROI.match_mine.median()))
print("\n15 uM split:\n" + SPL.to_string(index=False))

# --- figures: MAIN (example trace, QC outcome) / SUPP (detection cross-check, thresholds) ---------------------
def letter(a, s): a.text(-0.13, 1.06, s, transform=a.transAxes, fontsize=28, fontweight="bold", va="bottom", ha="left")
def compose(panels, ncol, name, w=8.0, h=7.5):
    """panels: list of functions f(ax); drawn in order and lettered (a), (b), ...; saved as figures/name.{png,pdf,eps}"""
    n = len(panels); nrow = -(-n // ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(w * ncol, h * nrow), squeeze=False)
    for k, fn in enumerate(panels):
        a = axs[k // ncol, k % ncol]; fn(a); letter(a, f"({chr(97 + k)})")
    for k in range(n, nrow * ncol): axs[k // ncol, k % ncol].axis("off")
    plt.tight_layout(w_pad=2.5, h_pad=2.5)
    for ext in ("png", "pdf", "eps"): fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig); print(f"wrote figures/{name}.{{png,pdf,eps}}")

def example_trace(a):
    if "trace" not in examples: return
    t, y, pk, d = examples["trace"]; m = (t > 1500) & (t < 3500)
    a.plot(t[m], y[m], color=C_GREY, lw=1.5, label="ΔF/F0")
    a.plot(t[pk][(t[pk] > 1500) & (t[pk] < 3500)], y[pk][(t[pk] > 1500) & (t[pk] < 3500)], "v", color=C_ME, ms=11, label="detected here")
    if d is not None:
        ta = d.time.values; ta = ta[(ta > 1500) & (ta < 3500)]; a.plot(ta, np.interp(ta, t, y) + 0.25, "^", color=C_DATA, ms=11, label="authors' table (+0.25)")
    a.set_xlabel(r"$t$ (s)"); a.set_ylabel(r"$\Delta F/F_0$"); a.set_title(r"08-08 E1, ROI #18 (15 $\mu$M)"); a.legend(loc="upper right", fontsize=15)
    a.set_ylim(0, 2.85); a.set_yticks([0, 0.5, 1.0, 1.5, 2.0, 2.5])

def match_hist(a):
    a.hist(ROI.match_author.dropna(), bins=np.linspace(0, 1, 21), color=TINT, edgecolor=C_DATA, label="table spikes found here")
    a.hist(ROI.match_mine.dropna(), bins=np.linspace(0, 1, 21), histtype="step", color=C_ME, lw=3, label="my peaks (in their window) in table")
    a.axvline(MIN_MATCH, color="k", ls="--", lw=2); a.set_xlabel("match rate per ROI (±5 s)"); a.set_ylabel("ROIs"); a.set_title("Detection vs authors' tables")
    a.legend(loc="upper left", fontsize=15)

def qc_by_date(a):
    cnt = ROI.groupby("date").agg(inc=("include", "sum"), tot=("include", "size")); cnt["exc"] = cnt.tot - cnt.inc
    x = np.arange(len(cnt)); a.bar(x, cnt.inc, color="#2ca02c", label="included"); a.bar(x, cnt.exc, bottom=cnt.inc, color=TINT, edgecolor=C_DATA, label="excluded")
    a.set_xticks(x); a.set_xticklabels([d[5:] for d in cnt.index], rotation=60, fontsize=14); a.set_ylabel("ROIs"); a.set_title("QC outcome by date"); a.legend(loc="upper left", fontsize=15)

def thresholds(a):
    a.semilogy(ROI.SNR.clip(lower=0.3), ROI.drift_pct_h.clip(lower=0.05), "o", color=C_GREY, ms=6, mec="none", label="ROI")
    inc = ROI[ROI.include]; a.semilogy(inc.SNR, inc.drift_pct_h.clip(lower=0.05), "o", color="#2ca02c", ms=6, mec="none", label="included")
    a.axvline(MIN_SNR, color="k", ls="--", lw=2); a.axhline(MAX_DRIFT, color="k", ls="--", lw=2)
    a.set_xscale("log"); a.set_xlabel("SNR (median prominence / noise)"); a.set_ylabel("baseline drift (%/h)"); a.set_title("QC thresholds"); a.legend(loc="upper left", fontsize=15)

MAIN = [example_trace, qc_by_date]      # what the records look like and what survives QC (Fig 1 material)
SUPP = [match_hist, thresholds]         # the detection cross-check and the thresholds behind the QC
compose(MAIN, 2, "data_qc_splits_main")
compose(SUPP, 2, "data_qc_splits_supp")
