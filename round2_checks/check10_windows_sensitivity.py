"""Check 10: (a) intervals spanning a boundary between two source segments in the union trains;
(b) sensitivity of the dose summaries, the segment statistics and the slope to excluding the three record-dose windows
whose first tabulated spike precedes the nominal start (2024-03-07-E1 at 10 and 15 uM, 2024-08-21-E1 at 3 uM)."""
from common import *
DD = os.path.join(CODE, "data_dose_units", "results")
E = pd.read_csv(os.path.join(DD, "events_canonical.csv")); T = pd.read_csv(os.path.join(DD, "train_table.csv")); S = pd.read_csv(os.path.join(DD, "segment_table.csv"))
# (a) cross-boundary intervals
E = E.sort_values(["train", "time"]); E["prev_file"] = E.groupby("train").file.shift(); E["cross"] = E.isi_diff.notna() & (E.file != E.prev_file)
n_isi = int(E.isi_diff.notna().sum()); n_cross = int(E.cross.sum())
print(f"union trains: {n_isi} intervals, {n_cross} span a boundary between two source files ({100*n_cross/n_isi:.1f} %)")
# (b) exclusion of the flagged windows
flag = {("CChCB_2024_03_07_E1", 10), ("CChCB_2024_03_07_E1", 15), ("CChCB_2024_08_21_E1", 3)}
T = T[T.n >= 10].copy(); T["flag"] = [(r, d) in flag for r, d in zip(T.record, T.dose)]
S = S[(S.n >= 10) & (S.sd_detr > 0)].copy(); S["flag"] = [(r, d) in flag for r, d in zip(S.record, S.dose)]
print(f"flagged trains {int(T.flag.sum())} of {len(T)}, flagged segments {int(S.flag.sum())} of {len(S)}")
rows = []
for lab, Tm, Sm in (("all", T, S), ("flagged windows excluded", T[~T.flag], S[~S.flag])):
    med = Tm.groupby("dose").isi_median.median(); n = Tm.groupby("dose").size()
    x, y = np.log(Sm.Tav.values), np.log(Sm.sd_detr.values); slope = np.polyfit(x, y, 1)[0]
    rows.append(dict(subset=lab, **{f"ISI_{d}uM": med.get(d, np.nan) for d in (3, 5, 10, 15, 30)}, **{f"n_{d}": int(n.get(d, 0)) for d in (3, 5, 10, 15, 30)},
                     ratio_3_30=med[3] / med[30], seg_n=len(Sm), seg_CV=Sm.CV_detr.median(), seg_skew=Sm.skew_detr.median(), seg_serial=Sm.serial_detr.median(), slope=slope))
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "check10_windows_sensitivity.csv"), index=False); print(R.round(3).T.to_string())
pd.DataFrame([dict(n_isi=n_isi, n_cross=n_cross)]).to_csv(os.path.join(OUT, "check10_cross_boundary.csv"), index=False)
