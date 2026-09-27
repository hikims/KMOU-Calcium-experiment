# data_qc_splits — data dictionary, QC rules, detection cross-check, 15 µM development/lock split (plan §4–5)

Python only.  Run:  `python3 qc_data_splits.py`  (~10 s; needs `../../*Spike Trains_data/1/*.txt` (or `CChCB/`) and `2/*.csv`).

## What it fixes before any fitting
- Real time axis `Time [s]` (never `ND.T`).  Every record has gaps > 10 s (1–6 per record, max 36–125 s); spikes within
  15 s of a gap edge are dropped and intervals across gaps are not used.
- Stimulus segments per record are taken from the authors' spike tables (per level: first spike − 60 s … last spike + 60 s)
  and written to `results/qc_records.csv` (column `segments`).  Base csv files only; `*_2.csv` … are segments of the same trains.
- ΔF/F0 with F0 = median of the pre-stimulus window (≥ 60 s) else the 10th percentile of the trace.
- Independent detection: `find_peaks`, prominence ≥ max(0.15, 5 σ_noise), distance ≥ 10 s.
- ROI ↔ ST: the authors' ST index is the ROI number.  Verified: of the 141 ROIs present in the tables, the table spikes
  are found by the detector with median match 1.00 (±5 s), and my peaks inside their time window are in the table with
  median match 1.00.  The authors keep only regular trains and often only a time segment, so 247 ROIs are absent from
  the tables; they are QC'd here on their own.
- QC (independent of the authors' selection): n_spikes ≥ 10, SNR ≥ 5, baseline drift ≤ 10 %/h (10th-percentile envelope,
  inside the stimulus segments after their first 120 s), no saturation, table_mismatch (< 0.8 of table spikes found).
  Result: 239 / 388 ROIs included; reasons in `results/qc_exclusions.csv`.
  By dose: 15 µM 128/170, 3 µM 20/39, 5 µM 25/30, 30 µM 27/52, step 5→10(→15) 33/43, step 5→30 6/54.
  The 5→30 records (04-03, 0.1 s sampling, 2000 s) fail mostly on drift (photobleaching at the high frame rate);
  a per-segment detrending rule for those two records is a decision to take before calibration, not a default here.
- 15 µM split fixed by file metadata: development = 07-11 E1, 07-18 E1, 08-08 E1, 08-08 E2 (106 ROIs, 84 after QC);
  lock = 08-09 E1/E2/E3 (64 ROIs, 44 after QC).  `results/qc_split_15uM.csv`.

## Outputs
`results/qc_records.csv` (23 records), `results/qc_rois.csv` (388 ROIs: counts, mapping, SNR, drift, include, reasons),
`results/qc_exclusions.csv`, `results/qc_split_15uM.csv`, `figures/data_qc_splits_{main,supp}.{png,pdf,eps}` (2×2:
example trace with both detections; match-rate histograms; QC outcome by date; SNR–drift plane with thresholds).

## Figures
- `figures/data_qc_splits_main.{png,pdf,eps}` (paper): example trace; QC outcome by date
- `figures/data_qc_splits_supp.{png,pdf,eps}` (supplement): detection vs authors' tables; SNR–drift thresholds
