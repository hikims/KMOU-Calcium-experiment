"""Check 11: (a) the 24 cross-segment intervals by dose and their effect on the dose medians; (b) acquisition gaps
(>10 s in the raw time axis) that fall inside an author source segment; (c) the 141 / 146 ROI counts."""
from common import *
import glob
DD = os.path.join(CODE, "data_dose_units", "results"); QC = os.path.join(CODE, "data_qc_splits", "results")
E = pd.read_csv(os.path.join(DD, "events_canonical.csv")).sort_values(["train", "time"])
E["prev_file"] = E.groupby("train").file.shift(); E["cross"] = E.isi_diff.notna() & (E.file != E.prev_file)
# (a) by dose, and train medians without the cross intervals
cr = E[E.cross]; print("cross intervals by dose:", cr.groupby("dose").size().to_dict(), "; trains affected:", cr.train.nunique())
T = pd.read_csv(os.path.join(DD, "train_table.csv")); T = T[T.n >= 10]
med_all = T.groupby("dose").isi_median.median()
Ex = E[~E.cross & E.isi_diff.notna()]; med_tr = Ex.groupby("train").isi_diff.median(); tr = T.set_index("train")
med_x = pd.Series({d: np.median([med_tr[t] for t in tr[tr.dose == d].index if t in med_tr]) for d in (3, 5, 10, 15, 30)})
rows = [dict(dose=d, median_all=med_all[d], median_no_cross=med_x[d]) for d in (3, 5, 10, 15, 30)]
A = pd.DataFrame(rows); A["ratio_3_30_all"] = med_all[3] / med_all[30]; A["ratio_3_30_no_cross"] = med_x[3] / med_x[30]; print(A.round(2).to_string(index=False))
A.to_csv(os.path.join(OUT, "check11_cross_by_dose.csv"), index=False)
# (b) gaps inside source segments: raw time axes from the Zenodo txt files
RAW = glob.glob(os.path.join(CODE, "..", "Fluorescence Intensity Data from HEK293 Ca²⁺ Spike Trains_data", "CChCB", "*.txt"))
S = pd.read_csv(os.path.join(DD, "segment_table.csv")); S = S[S.n >= 10]
gap_rows = []
for f in RAW:
    rec = os.path.basename(f).replace(".txt", "")
    try: t = pd.read_csv(io.StringIO(open(f, "rb").read().decode("utf-16")), sep="\t")["Time [s]"].values.astype(float)
    except Exception as e: print("skip", rec, e); continue
    g = np.diff(t); gaps = [(t[i], t[i + 1]) for i in np.where(g > 10)[0]]
    for _, s in S[S.record == rec].iterrows():
        inside = [(a, b) for a, b in gaps if a >= s.t_first and b <= s.t_last]
        if inside: gap_rows.append(dict(record=rec, segment=s.segment, ST=s.ST, dose=s.dose, n_gaps_inside=len(inside), max_gap=max(b - a for a, b in inside)))
G = pd.DataFrame(gap_rows); G.to_csv(os.path.join(OUT, "check11_gaps_in_segments.csv"), index=False)
print(f"raw records read: {len(RAW)}; segments with an acquisition gap inside: {len(G)} of {len(S)}"); print(G.head(20).to_string(index=False) if len(G) else "")
# (c) ROI counts
q = pd.read_csv(os.path.join(QC, "qc_rois.csv")); n141 = int(q.in_author_table.sum())
S["roi"] = S.record + ":" + S.ST.astype(str); n146 = S.roi.nunique(); nT = (T.record + ":" + T.ST.astype(str)).nunique()
print(f"ROIs in author tables matched to raw files: {n141}; ROI units among the {len(S)} analysed segments: {n146}; among trains: {nT}")
allE = pd.read_csv(os.path.join(DD, "events_canonical.csv")); print("distinct (record, ST) in all events:", (allE.record + ":" + allE.ST.astype(str)).nunique())
pd.DataFrame([dict(rois_matched_to_raw=n141, roi_units_segments=n146, roi_units_trains=nT, roi_units_all_events=(allE.record + ":" + allE.ST.astype(str)).nunique())]).to_csv(os.path.join(OUT, "check11_roi_counts.csv"), index=False)
