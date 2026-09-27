"""
data_dose_units.py — one canonical event table, train-level dose statistics, and the protocol/analysis-window table
for the HEK293 spike-train data set (Azzoni 2025, Zenodo 14779474 processed tables; 15705673 raw description).

Why this folder exists (audit items 3.1–3.2)
  * The processed folder "2/" holds 50 files.  A base file <record>.csv and its <record>_2.csv, _3.csv ... are NOT
    duplicates: the sub-files continue the same trains (spike numbering continues, e.g. base ST 2 spikes 3–13,
    _2 spikes 14–33).  Only 34 events are exact duplicates.  "Base files only" therefore truncated the long
    15 µM and 3 µM trains; the union is the right event set.
  * Pooling all events and taking a median per dose weights long trains more.  Here every train
    (record, ST, dose) is one unit: median of the train medians, bootstrap CI over trains.
  * 10 µM exists only inside 5→10→15 step records; 30 µM in 03_28 and 04_03_E1–E3 follows a 10–12 min 5 µM
    pre-exposure (Zenodo description).  Each (record, dose) is labelled constant / pre-exposed / step.
  * The analysis windows used by the QC ([first spike − 60 s, last spike + 60 s]) are compared with the nominal
    stimulus times of the Zenodo description (record_protocol.csv).

Outputs (results/)
  events_canonical.csv   record, ST, dose, protocol, spike, time, ISI (authors), isi_diff (recomputed), amplitude, W, ...
  segment_table.csv      one row per authors' segment (file; base = segment 1): the same statistics, raw — the unit
                         for within-cell statistics (m1_isi_variability, m1_sigma_tav_slope); unit_comparison.csv
                         compares segment / detrended-train / raw-train medians
  train_table.csv        one row per train: n, median/mean ISI, sd, CV, skew (raw and after removing a linear trend in
                         time: sd_detr, CV_detr, skew_detr, corr_detr, serial_detr), amp CV, corr(ISI, amp), serial,
                         trend %/min, first/last spike time, protocol label, file(s)
                         -> this table is the single source of the data statistics for m1_isi_variability,
                            m1_sigma_tav_slope, m0ip3h_amplitude_invariance and m0ip3h_stimulus_encoding
  dose_table.csv         per dose: n_records, n_trains, n_events, pooled median, median of train medians with
                         bootstrap CI, IQR of train medians, same for half-width and amplitude; constant-dose subset
  record_protocol.csv    per record and dose: nominal CCh window (Zenodo) vs analysis window (spike table)
  union_vs_base.csv      events per record: base file only vs union of all files
  dose_map_refit.csv     dose -> S̄ map of the stereotyped reference: endpoint convention (3 vs 30 µM, linear map on the
                         steepest decade of T(S̄); the middle doses are predictions) + the least-squares Hill map for reference;
                         uses ../m0ip3h_bifurcation_structure/results/bif_line_S_ref_st.csv (else the amplitude_invariance sweep)
  fold_change_targets.csv  between-cell dose ratios (constant-dose) and within-cell step ratios, train-level

Figures: figures/data_dose_units_main.{png,pdf,eps}, figures/data_dose_units_supp.{png,pdf,eps}
Usage:  python3 data_dose_units.py [NBOOT]        (bootstrap over trains runs on all cores)
"""
import glob, os, sys, re
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
NAME, OUT, FIG = "data_dose_units", "results", "figures"
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
NBOOT = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
C_POOL, C_TRAIN, C_STEP, C_PRE, C_MODEL, C_GREY = "#888888", "#c0392b", "#1f77b4", "#e69f00", "#2ca02c", "#bbbbbb"
MIN_SPIKES = 10

# --- nominal stimulus protocol (Zenodo 15705673 description) ---------------------------------------------------
def nominal(record):
    """list of (dose µM, t_start s, t_end s or None) for a record name CChCB_YYYY_MM_DD_EN"""
    d = record.replace("CChCB_", "")
    if d in ("2024_03_07_E1", "2024_03_12_E1", "2024_03_14_E2", "2024_03_19_E1"): return [(5, 600, 5400), (10, 5400, 8400), (15, 8400, None)]
    if d == "2024_03_28_E1": return [(5, 600, 1200), (30, 1200, None)]
    if d.startswith("2024_04_03"): return [(5, 300, 1020), (30, 1020, None)]
    if d == "2024_06_05_E1": return [(30, 600, None)]
    if d[:7] in ("2024_07", "2024_08") and d < "2024_08_21": return [(15, 900, None)]
    if d.startswith("2024_08_2"): return [(3, 900, None)]
    if d.startswith("2024_09"): return [(5, 600, None)]
    return []

def protocol_label(record, dose):
    nom = nominal(record); doses = [x[0] for x in nom]
    if len(doses) == 1: return "constant"
    if len(doses) == 2 and doses[0] == 5 and dose == 30: return "pre-exposed"       # brief 5 µM before 30 µM
    return "step"

# --- load: union of all files of a record, duplicates removed ----------------------------------------------------
ROOTS = ["../..", "../../data", "../data", "."]
def find_tables():
    for r in ROOTS:
        hits = sorted(glob.glob(os.path.join(r, "*Spike Trains_data", "2", "*.csv")))
        if hits: return hits
    raise SystemExit("processed spike tables (folder 2/*.csv) not found under " + ", ".join(os.path.abspath(r) for r in ROOTS))

def rec_name(fn):
    b = re.sub(r"(_\d+)?\.csv$", "", os.path.basename(fn))       # 20240711E1
    return f"CChCB_{b[:4]}_{b[4:6]}_{b[6:8]}_{b[8:]}"

files = find_tables()
raw = []
for f in files:
    d = pd.read_csv(f, sep=";", decimal=",", encoding="utf-8-sig").rename(columns={"half-width": "W", "stimulus": "dose"})
    d["file"] = os.path.basename(f); d["record"] = rec_name(f); d["is_base"] = ("_" not in os.path.basename(f))
    raw.append(d)
raw = pd.concat(raw, ignore_index=True)
raw["t_key"] = np.round(raw.time, 1)
E = raw.sort_values(["record", "ST", "dose", "time", "is_base"], ascending=[True, True, True, True, False])
E = E.drop_duplicates(["record", "ST", "dose", "t_key"], keep="first").copy()
E["protocol"] = [protocol_label(r, int(d)) for r, d in zip(E.record, E.dose)]
E["train"] = E.record + ":" + E.ST.astype(int).astype(str) + ":" + E.dose.astype(int).astype(str)
E["isi_diff"] = E.groupby("train").time.diff()
def _detrend(d):
    x, t = d.isi_diff.values, d.time.values; ok = ~np.isnan(x); out = np.full(len(x), np.nan)
    if ok.sum() >= 4: b = np.polyfit(t[ok], x[ok], 1); out[ok] = x[ok] - np.polyval(b, t[ok]) + x[ok].mean()
    elif ok.sum() > 0: out[ok] = x[ok]
    return pd.Series(out, index=d.index)
E["isi_detr"] = E.groupby("train", group_keys=False).apply(_detrend)      # slow drift removed (trains up to 2 h long)
E = E.sort_values(["record", "ST", "dose", "time"]).reset_index(drop=True)
cols = ["record", "ST", "dose", "protocol", "train", "spike", "time", "ISI", "isi_diff", "isi_detr", "amplitude", "W", "rise", "baseline", "decay50", "file"]
E[cols].to_csv(os.path.join(OUT, "events_canonical.csv"), index=False)
n_dup = len(raw) - len(E)
mism = E.isi_diff.notna() & (np.abs(E.isi_diff - E.ISI) > 2.5)
print(f"events: {len(raw)} rows in {len(files)} files -> {len(E)} unique ({n_dup} duplicates); "
      f"authors' ISI != time difference in {mism.sum()} events ({100*mism.mean():.1f} %)")
ub = raw.groupby("record").agg(n_files=("file", "nunique"), events_all=("time", "size"))
ub["events_base"] = raw[raw.is_base].groupby("record").size(); ub["events_union"] = E.groupby("record").size()
ub["gained_by_union"] = ub.events_union - ub.events_base
ub.to_csv(os.path.join(OUT, "union_vs_base.csv"))
print(f"union vs base files: {int(ub.events_base.sum())} -> {int(ub.events_union.sum())} events (+{int(ub.gained_by_union.sum())})")

# --- train table ------------------------------------------------------------------------------------------------
def train_stats(d):
    d = d.sort_values("time"); x = d.isi_diff.dropna().values if len(d) > 1 else np.array([])
    A = d.amplitude.values[1:] if len(d) > 1 else np.array([])
    r = dict(train=d.train.iloc[0], record=d.record.iloc[0], ST=int(d.ST.iloc[0]), dose=int(d.dose.iloc[0]), protocol=d.protocol.iloc[0],
             n=len(d), t_first=d.time.min(), t_last=d.time.max(), files=";".join(sorted(d.file.unique())),
             isi_median=np.nan, Tav=np.nan, sd=np.nan, CV=np.nan, skew=np.nan, sd_detr=np.nan, CV_detr=np.nan, skew_detr=np.nan,
             amp_median=np.nanmedian(d.amplitude), amp_cv=np.nan, W_median=np.nanmedian(d.W), corr_isi_amp=np.nan, corr_detr=np.nan,
             serial=np.nan, serial_detr=np.nan, trend_pct_per_min=np.nan)
    if len(x) >= 3:
        m = x.mean(); s = x.std(ddof=1); r.update(isi_median=np.median(x), Tav=m, sd=s, CV=s / m)
        r["skew"] = float(np.mean(((x - m) / x.std()) ** 3)) if x.std() > 0 else np.nan
        r["amp_cv"] = A.std(ddof=1) / A.mean() if A.mean() > 0 else np.nan
        if len(x) >= 4:
            r["corr_isi_amp"] = np.corrcoef(x, A)[0, 1]
            r["serial"] = np.corrcoef(x[:-1], x[1:])[0, 1]
            t = d.time.values[1:]; b = np.polyfit(t, x, 1); res = x - np.polyval(b, t)
            r["serial_detr"] = np.corrcoef(res[:-1], res[1:])[0, 1]
            r["trend_pct_per_min"] = 100 * 60 * b[0] / m
            sdd = res.std(ddof=1); r.update(sd_detr=sdd, CV_detr=sdd / m, skew_detr=float(np.mean((res / res.std()) ** 3)) if res.std() > 0 else np.nan,
                                            corr_detr=np.corrcoef(res, A)[0, 1])
    return r
T = pd.DataFrame([train_stats(d) for _, d in E.groupby("train", sort=False)])
T.to_csv(os.path.join(OUT, "train_table.csv"), index=False)
Tn = T[T.n >= MIN_SPIKES].copy()
# --- segment table: the authors' files are stationary segments of a train (base file = segment 1, _2, _3, ... = later
#     segments of the same ROI); statistics inside a segment need no detrending.  Unit for within-cell statistics.
E["segment"] = E.file.str.replace(r"\.csv$", "", regex=True)
def seg_stats(d):
    d = d.sort_values("time"); r = train_stats(d.assign(isi_diff=d.time.diff())); r["segment"] = d.segment.iloc[0]; r["seg_index"] = int(d.file.iloc[0].split("_")[1][:-4]) if "_" in d.file.iloc[0] else 1
    return r
SEG = pd.DataFrame([seg_stats(d) for _, d in E.groupby(["train", "segment"], sort=False)])
SEG = SEG[["segment", "seg_index"] + [c for c in SEG.columns if c not in ("segment", "seg_index")]]
SEG.to_csv(os.path.join(OUT, "segment_table.csv"), index=False)
SEGn = SEG[SEG.n >= MIN_SPIKES]
cmp_ = pd.DataFrame({"unit": ["segment (authors', raw)", "segment (authors', detrended)", "train (union, detrended)", "train (union, raw)"],
                     "n_units": [len(SEGn), len(SEGn), len(Tn), len(Tn)],
                     "CV": [SEGn.CV.median(), SEGn.CV_detr.median(), Tn.CV_detr.median(), Tn.CV.median()],
                     "skew": [SEGn["skew"].median(), SEGn.skew_detr.median(), Tn.skew_detr.median(), Tn["skew"].median()],
                     "amp_cv": [SEGn.amp_cv.median()] * 2 + [Tn.amp_cv.median()] * 2,
                     "corr_isi_amp": [SEGn.corr_isi_amp.median(), SEGn.corr_detr.median(), Tn.corr_detr.median(), Tn.corr_isi_amp.median()],
                     "serial": [SEGn.serial.median(), SEGn.serial_detr.median(), Tn.serial_detr.median(), Tn.serial.median()],
                     "median_n_spikes": [SEGn.n.median()] * 2 + [Tn.n.median()] * 2})
cmp_.to_csv(os.path.join(OUT, "unit_comparison.csv"), index=False)
print(f"segments: {len(SEG)} total, {len(SEGn)} with >= {MIN_SPIKES} spikes"); print(cmp_.round(3).to_string(index=False))
print(f"trains: {len(T)} total, {len(Tn)} with >= {MIN_SPIKES} spikes")

# --- dose table: pooled vs train-level, bootstrap over trains (parallel) -------------------------------------------
def boot_median(args):
    vals, nb, seed = args; rng = np.random.default_rng(seed); n = len(vals)
    return np.array([np.median(vals[rng.integers(0, n, n)]) for _ in range(nb)])
def boot_ci(vals, nb=NBOOT):
    vals = np.asarray(vals, float); vals = vals[~np.isnan(vals)]
    if len(vals) < 3: return (np.nan, np.nan)
    nw = max(1, cpu_count()); per = -(-nb // nw)
    with get_context("fork").Pool(nw) as pool: bs = np.concatenate(pool.map(boot_median, [(vals, per, 1000 + k) for k in range(nw)]))
    return (float(np.quantile(bs, .025)), float(np.quantile(bs, .975)))

rows = []
for dose, g in Tn.groupby("dose"):
    ev = E[(E.dose == dose) & E.train.isin(g.train)]
    for subset, gg in (("all", g), ("constant", g[g.protocol == "constant"]), ("constant+pre", g[g.protocol.isin(["constant", "pre-exposed"])])):
        if len(gg) == 0: continue
        evs = ev[ev.train.isin(gg.train)]
        lo, hi = boot_ci(gg.isi_median.values)
        rows.append(dict(dose=dose, subset=subset, n_records=gg.record.nunique(), n_trains=len(gg), n_events=len(evs),
                         protocols="+".join(sorted(gg.protocol.unique())),
                         isi_pooled_median=np.nanmedian(evs.isi_diff), isi_pooled_q1=np.nanquantile(evs.isi_diff, .25), isi_pooled_q3=np.nanquantile(evs.isi_diff, .75),
                         isi_train_median=np.median(gg.isi_median), isi_train_ci_lo=lo, isi_train_ci_hi=hi,
                         isi_train_q1=np.quantile(gg.isi_median, .25), isi_train_q3=np.quantile(gg.isi_median, .75),
                         isi_record_median=float(np.median(gg.groupby("record").isi_median.median())),
                         W_train_median=np.median(gg.W_median), W_train_q1=np.quantile(gg.W_median, .25), W_train_q3=np.quantile(gg.W_median, .75),
                         amp_train_median=np.median(gg.amp_median), amp_train_q1=np.quantile(gg.amp_median, .25), amp_train_q3=np.quantile(gg.amp_median, .75),
                         CV_train_median=np.nanmedian(gg.CV), serial_detr_median=np.nanmedian(gg.serial_detr), trend_median=np.nanmedian(gg.trend_pct_per_min)))
Dt = pd.DataFrame(rows); Dt.to_csv(os.path.join(OUT, "dose_table.csv"), index=False)
print(Dt[Dt.subset == "all"][["dose", "n_records", "n_trains", "n_events", "isi_pooled_median", "isi_train_median", "isi_train_ci_lo", "isi_train_ci_hi", "isi_record_median", "protocols"]].to_string(index=False))

# --- record protocol vs analysis windows ------------------------------------------------------------------------
prow = []
for rec, g in E.groupby("record"):
    for dose, t0, t1 in nominal(rec):
        gd = g[g.dose == dose]
        prow.append(dict(record=rec, dose=dose, protocol=protocol_label(rec, dose), nominal_start=t0, nominal_end=t1 if t1 else np.nan,
                         n_events=len(gd), n_trains=gd.train.nunique(),
                         first_spike=gd.time.min() if len(gd) else np.nan, last_spike=gd.time.max() if len(gd) else np.nan,
                         window_start=gd.time.min() - 60 if len(gd) else np.nan, window_end=gd.time.max() + 60 if len(gd) else np.nan,
                         first_spike_minus_nominal=(gd.time.min() - t0) if len(gd) else np.nan))
PR = pd.DataFrame(prow); PR.to_csv(os.path.join(OUT, "record_protocol.csv"), index=False)
outside = PR[(PR.first_spike < PR.nominal_start - 1)]
print(f"record_protocol: {len(PR)} (record, dose) rows; first spike before the nominal start in {len(outside)} rows")

# --- fold-change targets: between-cell (constant dose) and within-cell (steps) ----------------------------------
def med_const(d):
    r = Dt[(Dt.dose == d) & (Dt.subset == ("constant+pre" if d == 30 else "constant"))]
    return float(r.isi_train_median.iloc[0]) if len(r) else np.nan
ft = []
for lo_d, hi_d in [(3, 5), (3, 15), (3, 30), (5, 15), (5, 30), (15, 30)]:
    ft.append(dict(kind="between-cell (constant dose)", lo=lo_d, hi=hi_d, ratio=med_const(lo_d) / med_const(hi_d), n_lo=int(Tn[(Tn.dose == lo_d) & (Tn.protocol != "step")].shape[0]), n_hi=int(Tn[(Tn.dose == hi_d) & (Tn.protocol != "step")].shape[0])))
# within-cell: same record & ST, both levels with >= 5 spikes each; ratio of train medians
Tw = T[T.n >= 5]
for lo_d, hi_d in [(5, 10), (10, 15), (5, 30), (5, 15)]:
    a = Tw[(Tw.dose == lo_d)].set_index(["record", "ST"]); b = Tw[(Tw.dose == hi_d)].set_index(["record", "ST"])
    common = a.index.intersection(b.index)
    if len(common) == 0: continue
    rr = (a.loc[common].isi_median / b.loc[common].isi_median).values
    ft.append(dict(kind="within-cell (step)", lo=lo_d, hi=hi_d, ratio=float(np.median(rr)), n_lo=len(common), n_hi=len(common)))
FT = pd.DataFrame(ft); FT.to_csv(os.path.join(OUT, "fold_change_targets.csv"), index=False); print(FT.to_string(index=False))

# --- dose map of the stereotyped reference: (a) endpoint convention and (b) least squares ------------------------------
# (a) Only the 3 and 30 µM levels are used quantitatively (the middle levels sit on a plateau dominated by
#     between-record variation).  The dose enters through a monotone map S̄ = c U/(K+U); across one decade of dose the
#     model period can fall at most by max_S T(S)/T(10 S) (attained for the linear map, K -> inf).  The endpoint map
#     takes the decade of steepest descent, c = S̄_lo/3, and T0 from the geometric mean of the two endpoint ratios, so the
#     3 and 30 µM errors are equal and opposite in log; the middle doses are then predictions, not fits.
# (b) The least-squares Hill map on the four constant-dose levels is kept for reference.
REFIT = None
cand = [p_ for p_ in ["../m0ip3h_bifurcation_structure/results/bif_line_S_ref_st.csv", "bif_line_S_ref_st.csv",
                      "../m0ip3h_amplitude_invariance/results/m0ip3h_inv_vs_stimulus_ref_st.csv", "m0ip3h_inv_vs_stimulus_ref_st.csv"] if os.path.exists(p_)]
if cand:
    sw = pd.read_csv(cand[0]); sw = sw[sw.osc == 1]; Sgrid = sw["x"].values if "x" in sw else sw.value.values; Tgrid = sw["T"].values if "T" in sw else sw.Tosc.values
    lT = lambda s: np.interp(np.log(s), np.log(Sgrid), np.log(Tgrid))                      # log T as a function of S̄
    isi3, isi30 = med_const(3), med_const(30); R_data = isi3 / isi30
    # (a) endpoints: steepest decade
    ls_ = np.linspace(np.log(Sgrid.min()), np.log(Sgrid.max()) - np.log(10), 600)
    dec = np.array([lT(np.exp(l)) - lT(np.exp(l) * 10) for l in ls_]); k = int(np.argmax(dec)); S_lo = float(np.exp(ls_[k])); R_model = float(np.exp(dec[k]))
    c_end = S_lo / 3.0; T0_end = float(np.sqrt(isi3 * isi30 / np.exp(lT(S_lo) + lT(10 * S_lo))))
    S_end = {d: c_end * d for d in [3, 5, 10, 15, 30]}
    pred = {d: T0_end * np.exp(lT(S_end[d])) for d in S_end}
    # (b) least squares on the four constant-dose levels
    dd = np.array([3, 5, 15, 30]); isi = np.array([med_const(d) for d in dd]); best = None
    for lc in np.linspace(np.log(0.01), np.log(300), 80):
        for lK in np.linspace(np.log(0.3), np.log(3000), 80):
            m = np.exp(np.array([lT(np.exp(lc) * d / (np.exp(lK) + d)) for d in dd]))
            T0 = np.exp(np.mean(np.log(isi) - np.log(m))); err = np.sum((np.log(isi) - np.log(T0 * m)) ** 2)
            if best is None or err < best[0]: best = (err, np.exp(lc), np.exp(lK), T0)
    err, c, K, T0 = best
    REFIT = dict(method="endpoints", c=c_end, K=np.inf, T0=T0_end, S_dose=S_end, pred=pred, R_data=R_data, R_model=R_model, S_lo=S_lo,
                 ls=dict(c=c, K=K, T0=T0, log_rms=np.sqrt(err / len(dd))), Sgrid=Sgrid, Tgrid=Tgrid, dd=dd, isi=isi, source=cand[0])
    pd.DataFrame([dict(method="endpoints (3 vs 30 uM, linear map on the steepest decade)", c=c_end, K=np.inf, T0_s=T0_end,
                       **{f"S_{d}uM": v for d, v in S_end.items()}, **{f"pred_ISI_{d}uM": v for d, v in pred.items()},
                       data_ratio_3_30=R_data, model_max_ratio_per_decade=R_model, S_lo_of_steepest_decade=S_lo,
                       ls_c=c, ls_K=K, ls_T0_s=T0, ls_log_rms=np.sqrt(err / len(dd)), period_curve=cand[0])]).to_csv(os.path.join(OUT, "dose_map_refit.csv"), index=False)
    print(f"dose map (endpoints): data 3/30 ratio {R_data:.2f}; model max ratio per decade {R_model:.2f} at S_lo = {S_lo:.3g}; "
          f"T0 = {T0_end:.3f} s/unit; S(3..30 uM) = " + ", ".join(f"{d}:{v:.3g}" for d, v in S_end.items())
          + "; predicted ISI 5/10/15 uM = " + ", ".join(f"{pred[d]:.0f}" for d in (5, 10, 15)) + f" s (observed {med_const(5):.0f}/{float(Dt[(Dt.dose == 10) & (Dt.subset == 'all')].isi_train_median.iloc[0]):.0f}/{med_const(15):.0f})")
    print(f"dose map (least squares, reference): c={c:.2f}, K={K:.0f} uM, T0={T0:.3f} s, log-RMS {np.sqrt(err / len(dd)):.3f}")
else:
    print("dose map skipped (no REF_ST period curve found)")

# --- figure assembly -------------------------------------------------------------------------------------------
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

DOSES = [3, 5, 10, 15, 30]
def p_dose_units(a):
    """pooled-event median vs median of train medians (CI over trains); step-only levels hollow"""
    A = Dt[Dt.subset == "all"].set_index("dose")
    a.plot(A.index, A.isi_pooled_median, "s--", color=C_POOL, ms=11, label="pooled events (previous)")
    for d in A.index:
        r = A.loc[d]; hollow = (r.protocols == "step")
        a.errorbar([d], [r.isi_train_median], yerr=[[r.isi_train_median - r.isi_train_ci_lo], [r.isi_train_ci_hi - r.isi_train_median]],
                   fmt="o", color=C_TRAIN, mfc="white" if hollow else C_TRAIN, ms=12, capsize=6, elinewidth=2, mew=2.5)
        a.text(d, r.isi_train_ci_hi * 1.0, f"{int(r.n_trains)}", ha="center", va="bottom", fontsize=17, color=C_TRAIN)
    a.plot([], [], "o", color=C_TRAIN, ms=12, label="train medians (95 % CI)"); a.plot([], [], "o", color=C_TRAIN, mfc="white", ms=12, mew=2.5, label="step level only")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xticks(DOSES); a.set_xticklabels([str(d) for d in DOSES]); a.minorticks_off()
    a.set_xlabel(r"carbachol ($\mu$M)"); a.set_ylabel("median ISI (s)"); a.set_title("Aggregation unit"); a.legend(loc="upper right", fontsize=16)

def p_train_strip(a):
    """every train's median ISI by dose, coloured by protocol"""
    rng = np.random.default_rng(3)
    for lab, col, mk in (("constant", C_TRAIN, "o"), ("pre-exposed", C_PRE, "D"), ("step", C_STEP, "^")):
        g = Tn[Tn.protocol == lab]
        a.plot(g.dose * np.exp(rng.normal(0, 0.05, len(g))), g.isi_median, mk, color=col, ms=7, alpha=0.75, label=lab)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xticks(DOSES); a.set_xticklabels([str(d) for d in DOSES]); a.minorticks_off()
    a.set_xlabel(r"carbachol ($\mu$M)"); a.set_ylabel("train median ISI (s)"); a.set_title(rf"Trains ($n \geq {MIN_SPIKES}$ spikes)"); a.legend(loc="upper right", fontsize=16)

def p_dose_map(a):
    if REFIT is None: a.text(0.5, 0.5, "REF_ST period curve\nnot found", ha="center", va="center", transform=a.transAxes); return
    U = np.logspace(np.log10(2), np.log10(40), 200); lT = lambda s: np.interp(np.log(s), np.log(REFIT["Sgrid"]), np.log(REFIT["Tgrid"]))
    T_end = REFIT["T0"] * np.exp(np.array([lT(REFIT["c"] * u) for u in U]))
    ls = REFIT["ls"]; T_ls = ls["T0"] * np.exp(np.array([lT(ls["c"] * u / (ls["K"] + u)) for u in U]))
    A = Dt[(Dt.subset == "all")].set_index("dose")
    for d in REFIT["dd"]:
        r = Dt[(Dt.dose == d) & (Dt.subset == ("constant+pre" if d == 30 else "constant"))].iloc[0]
        a.errorbar([d], [r.isi_train_median], yerr=[[r.isi_train_median - r.isi_train_ci_lo], [r.isi_train_ci_hi - r.isi_train_median]], fmt="o", color=C_TRAIN, ms=12, capsize=6, elinewidth=2)
    a.plot([], [], "o", color=C_TRAIN, ms=12, label="constant dose, train medians")
    r10 = A.loc[10]; a.plot([10], [r10.isi_train_median], "o", color=C_TRAIN, mfc="white", ms=12, mew=2.5, label=r"10 $\mu$M (step only)")
    a.plot(U, T_end, "-", color=C_MODEL, label=rf"endpoints: $t_{{\rm unit}}={REFIT['T0']:.2f}$ s, $\bar S \propto U$")
    a.plot(U, T_ls, "--", color=C_GREY, label=rf"least squares ($K={ls['K']:.0f}\,\mu$M)")
    a.text(0.03, 0.05, rf"3/30 $\mu$M ratio: data {REFIT['R_data']:.1f}, model $\leq$ {REFIT['R_model']:.1f} per decade", transform=a.transAxes, fontsize=15)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xticks(DOSES); a.set_xticklabels([str(d) for d in DOSES]); a.minorticks_off()
    a.set_xlabel(r"carbachol ($\mu$M)"); a.set_ylabel("ISI (s)"); a.set_title("Dose map (endpoint convention)"); a.legend(loc="upper right", fontsize=14)

def p_protocol(a):
    """nominal CCh windows (bars) vs first/last spike (markers) per record"""
    recs = sorted(PR.record.unique()); y = {r: i for i, r in enumerate(recs)}
    cmap = {3: "#a6cee3", 5: "#b2df8a", 10: "#fdbf6f", 15: "#fb9a99", 30: "#cab2d6"}
    for _, r in PR.iterrows():
        t1 = r.nominal_end if not np.isnan(r.nominal_end) else max(PR.last_spike.max(), r.nominal_start + 600)
        a.barh(y[r.record], t1 - r.nominal_start, left=r.nominal_start, height=0.7, color=cmap[int(r.dose)], edgecolor="none")
        if not np.isnan(r.first_spike): a.plot([r.first_spike, r.last_spike], [y[r.record]] * 2, "k|-", lw=1.5, ms=9)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    hd = [Patch(color=c, label=rf"{d} $\mu$M") for d, c in cmap.items()] + [Line2D([], [], color="k", marker="|", ms=9, lw=1.5, label="first–last spike")]
    a.set_yticks(range(len(recs))); a.set_yticklabels([r.replace("CChCB_2024_", "") for r in recs], fontsize=13); a.invert_yaxis()
    a.set_xlabel(r"$t$ (s)"); a.set_title("Nominal stimulus vs spikes"); a.legend(handles=hd, loc="lower right", fontsize=13, ncol=2); a.grid(axis="y", visible=False)

def p_union(a):
    u = ub.sort_index()
    x = np.arange(len(u)); a.bar(x - 0.2, u.events_base, 0.4, color=C_GREY, label="base file only")
    a.bar(x + 0.2, u.events_union, 0.4, color=C_TRAIN, label="union of files")
    a.set_xticks(x); a.set_xticklabels([r.replace("CChCB_2024_", "") for r in u.index], rotation=90, fontsize=12)
    a.set_ylabel("events"); a.set_title("Events recovered per record"); a.legend(loc="upper right", fontsize=16)

def p_trend(a):
    """train-level ISI drift (% of mean per min) by dose"""
    data = [Tn[Tn.dose == d].trend_pct_per_min.dropna().values for d in DOSES]
    bp = a.boxplot(data, positions=range(len(DOSES)), widths=0.6, patch_artist=True, showfliers=True)
    for b in bp["boxes"]: b.set(facecolor="#f5d0cc", edgecolor=C_TRAIN, linewidth=2)
    for k in ("whiskers", "caps", "medians"): [l.set(color=C_TRAIN, linewidth=2) for l in bp[k]]
    a.axhline(0, color="k", lw=1.5, ls=":")
    a.set_xticks(range(len(DOSES))); a.set_xticklabels([str(d) for d in DOSES]); a.set_xlabel(r"carbachol ($\mu$M)")
    a.set_ylabel("ISI trend (% of mean / min)"); a.set_title("Within-train drift")

def p_stats_units(a):
    """CV and serial correlation per train: median with IQR, by dose"""
    for col, lab, c, off in (("CV", "CV", C_TRAIN, -0.12), ("serial_detr", "serial corr (detrended)", C_STEP, 0.12)):
        med = [np.nanmedian(Tn[Tn.dose == d][col]) for d in DOSES]; q1 = [np.nanquantile(Tn[Tn.dose == d][col], .25) for d in DOSES]; q3 = [np.nanquantile(Tn[Tn.dose == d][col], .75) for d in DOSES]
        a.errorbar(np.arange(len(DOSES)) + off, med, yerr=[np.array(med) - q1, np.array(q3) - med], fmt="o", color=c, ms=10, capsize=5, elinewidth=2, label=lab)
    a.axhline(0, color="k", lw=1.5, ls=":"); a.set_xticks(range(len(DOSES))); a.set_xticklabels([str(d) for d in DOSES])
    a.set_xlabel(r"carbachol ($\mu$M)"); a.set_ylabel("per-train statistic"); a.set_title("ISI statistics by dose"); a.legend(loc="upper right", fontsize=16)

MAIN = [p_dose_units, p_train_strip, p_dose_map]        # the aggregation unit changes the dose curve; the refit
SUPP = [p_protocol, p_union, p_trend]                   # protocol windows, union vs base, drift (p_stats_units is Fig 7d)
if __name__ == "__main__":
    compose(MAIN, 3, f"{NAME}_main")
    compose(SUPP, 3, f"{NAME}_supp")
