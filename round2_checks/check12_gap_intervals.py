"""Check 12: intervals of the analysed segments that span an acquisition gap (>10 s) of the raw time axis; segment
statistics with those intervals removed; and the 141 vs 146 ROI counts."""
from common import *
import glob
DD = os.path.join(CODE, "data_dose_units", "results"); QC = os.path.join(CODE, "data_qc_splits", "results")
E = pd.read_csv(os.path.join(DD, "events_canonical.csv")); E["segment"] = E.file.str.replace(r"\.csv$", "", regex=True)
S = pd.read_csv(os.path.join(DD, "segment_table.csv")); S = S[S.n >= 10]
RAW = {os.path.basename(f).replace(".txt", ""): f for f in glob.glob(os.path.join(CODE, "..", "Fluorescence Intensity Data from HEK293 Ca²⁺ Spike Trains_data", "CChCB", "*.txt"))}
gaps = {}
for rec, f in RAW.items():
    t = pd.read_csv(io.StringIO(open(f, "rb").read().decode("utf-16")), sep="\t")["Time [s]"].values.astype(float)
    g = np.diff(t); gaps[rec] = [(t[i], t[i + 1]) for i in np.where(g > 10)[0]]
def stats(x, t, keep_pair=None):
    """CV, skew on the kept intervals; serial correlation only over pairs of originally consecutive intervals"""
    xd = detrended(x, t); r = xd - xd.mean()
    if keep_pair is None: keep_pair = np.ones(len(r) - 1, bool)
    return dict(CV=xd.std(ddof=1) / xd.mean(), skew=float(np.mean((r / r.std()) ** 3)), serial=np.corrcoef(r[:-1][keep_pair], r[1:][keep_pair])[0, 1])
rows = []; n_span = 0; n_int = 0
for _, s in S.iterrows():
    e = E[(E.train == s.train) & (E.segment == s.segment)].sort_values("time"); tt = e.time.values
    if len(tt) < 10: continue
    x = np.diff(tt); tm = tt[1:]; n_int += len(x)
    span = np.array([any(tt[i] < a and tt[i + 1] > b for a, b in gaps.get(s.record, [])) for i in range(len(x))]); n_span += int(span.sum())
    a = stats(x, tm); keep = ~span
    if keep.sum() >= 4:
        idx = np.flatnonzero(keep); pair_ok = np.diff(idx) == 1          # pairs that were consecutive before the removal
        b = stats(x[keep], tm[keep], keep_pair=pair_ok)
    else: b = a
    rows.append(dict(segment=s.segment, train=s.train, dose=s.dose, n_int=len(x), n_span=int(span.sum()), **{k + "_all": a[k] for k in a}, **{k + "_nogap": b[k] for k in b}))
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "check12_gap_intervals.csv"), index=False)
print(f"segments {len(R)}, intervals {n_int}, spanning a gap {n_span} ({100*n_span/n_int:.2f} %), segments affected {(R.n_span>0).sum()}")
for k in ("CV", "skew", "serial"): print(f"median {k}: all {R[k+'_all'].median():.3f} -> without gap-spanning intervals {R[k+'_nogap'].median():.3f}")
# ROI counts
q = pd.read_csv(os.path.join(QC, "qc_rois.csv")); ev = set(zip(E.record, E.ST.astype(int)))
qq = q[q.in_author_table]; qs = set(zip(qq.record, qq.roi.str.replace("#", "").astype(int)))
print("events ROI units", len(ev), "; QC ROIs in author table", len(qs), "; in events but not QC:", sorted(ev - qs)[:10], "; in QC not events:", sorted(qs - ev)[:10])
