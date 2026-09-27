# data_dose_units — aggregation unit, event union, protocol windows (audit items 3.1–3.2)

## What changed and why
1. **Event set.** The processed folder `2/` has 50 files; `<record>_2.csv, _3.csv, ...` are *continuations* of the
   trains in `<record>.csv` (spike numbering continues), not duplicates.  Union of all files: 6 094 rows → 6 060
   unique events (34 exact duplicates).  "Base files only" (used by `m1_isi_variability`, `m1_sigma_tav_slope`,
   `data_qc_splits`) had 5 529 events, i.e. it truncated the long 15 µM and 3 µM trains.  The canonical table
   is `results/events_canonical.csv`; every train is `record:ST:dose`.
2. **Unit of aggregation.** Per dose the previous scripts pooled events; here one train = one unit
   (median of train medians, bootstrap CI over trains, plus the record-level median).
3. **Protocol.** `results/record_protocol.csv` lists the nominal CCh windows of the Zenodo description
   (15705673) next to the first/last spike of the spike tables.  10 µM exists only as a step level;
   30 µM in 03_28 and 04_03_E1–E3 follows a 10–12 min 5 µM pre-exposure (label `pre-exposed`); the
   remaining (record, dose) pairs are `constant`.  Three (record, dose) rows have their first spike before
   the nominal start (03_07_E1 at 10 and 15 µM by ~420 s, 08_21_E1 at 3 µM by 645 s) — either a clock
   offset in the description or a labelling of spikes by the segment they belong to.

## Results (train-level, n ≥ 10 spikes; 172 trains)
| dose (µM) | records | trains | pooled median (previous) | median of train medians [95 % CI] | protocol |
|---|---|---|---|---|---|
| 3 | 5 | 16 | 207 | 226 [171, 319] | constant |
| 5 | 8 | 39 | 76 | 92 [72, 116] | constant + step |
| 10 | 4 | 14 | 56 | 128 [76, 172] | step only |
| 15 | 10 | 71 | 56 | 94 [80, 110] | constant + step |
| 30 | 5 | 32 | 30 | 33 [28, 37] | constant + pre-exposed |

* 5 and 15 µM are indistinguishable at the train level (ratio 1.02 between cells, 1.00 within cells);
  the dose response is 3 → (5 ≈ 15) → 30, a 6.9-fold range.  Between-record spread at one dose
  (15 µM: 40–174 s) is as large as the dose effect.
* Within-cell step ratios (train medians, ≥ 5 spikes each level): 5→10 0.85 (n = 10), 10→15 1.08 (12),
  5→30 1.79 (4), 5→15 1.00 (9).
* Detrended statistics (linear trend in time removed per train; the union makes trains up to 2 h long and the
  drift inflates raw values): CV 0.213 [0.14, 0.34] (raw 0.253), skew 0.55 (raw 0.77), amp CV 0.145,
  corr(ISI, amp) 0.05, serial corr −0.02 (raw 0.23).  These are the targets of the stationary model.
* **Dose map, endpoint convention (a).**  Only 3 and 30 µM are used quantitatively (ratio 6.95 [CI-end ratio 4.6]);
  the middle levels sit on a plateau dominated by between-record variation.  With any monotone map
  S̄ = cU/(K+U) the stereotyped reference can lower its period by at most 2.8× per decade of dose (linear map on
  the steepest decade, S̄_lo = 0.064), so it reproduces 2.8 of the 6.95-fold change.  T0 = 0.527 s per model
  unit (geometric mean of the two endpoint ratios), S̄(3, 5, 10, 15, 30 µM) = 0.064, 0.11, 0.21, 0.32, 0.64;
  predicted 5/10/15 µM ISI 116/84/69 s vs observed 94/128/92 s.  The least-squares Hill map (K → 3000 µM,
  log-RMS 0.37) is kept in the same file for reference.  `results/dose_map_refit.csv`.

## Files
- `data_dose_units.py [NBOOT]` — everything (bootstrap over trains on all cores, fork pool); ~6 s
- `results/` — `events_canonical.csv`, `train_table.csv`, `dose_table.csv`, `record_protocol.csv`,
  `union_vs_base.csv`, `fold_change_targets.csv`, `dose_map_refit.csv`
- `figures/data_dose_units_main.*` — (a) pooled vs train-level dose curve, (b) all train medians by protocol,
  (c) dose-map refit; `figures/data_dose_units_supp.*` — (a) nominal windows vs spikes, (b) union vs base
  events per record, (c) within-train drift, (d) per-train CV and detrended serial correlation by dose

## Downstream
`m0ip3h_stimulus_encoding` (dose medians, CI), `m0ip3h_amplitude_invariance` (dose table, within-cell steps,
train correlations), `m1_isi_variability` and `m1_sigma_tav_slope` (train statistics, detrended) read
`train_table.csv` / `dose_table.csv` / `events_canonical.csv`; `emrich_jin_prediction`, `m0_joint_scan` and
`m0ip3h_bifurcation_structure` read `dose_table.csv` / `dose_map_refit.csv`.  Run this folder first.

## Update 2026-09-22 (segments)
`results/segment_table.csv`: the authors' files are stationary segments of a train (base file = segment 1; `_2`, `_3`
= later segments of the same ROI).  Within-cell statistics use segments (detrended inside the segment);
dose medians use trains; between-record comparisons use records.  `results/unit_comparison.csv` shows that the
segment and detrended-train medians agree (CV 0.196 vs 0.213, skew 0.49 vs 0.55, serial −0.03 vs −0.02).
