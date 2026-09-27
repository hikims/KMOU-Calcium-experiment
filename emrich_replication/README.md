# emrich_replication — independent replication with the Emrich et al. 2021 (Cell Reports) source data

Python only.  Run: `python3 emrich_replication.py [path/to/Source_Data_Emrich.xlsx]`
(default: `../../*signaling_data/Source_Data_Emrich.xlsx`, i.e. the Mendeley folder `Omnitemporal chore…Ca2+ signaling_data` next to `code/`; run from this folder or pass the path).

## What the source-data workbook actually contains
- Representative single-cell Fura-2 traces only: **5 cells per panel**, 46 panels = 230 traces (rows without a time
  column; 2 s sampling; 451 points = 15 min, 631 = 21 min dose steps 3→10→30→90 µM CCh (5 min each), 781 = 26 min).
  Genotypes: WT, STIM1/2/1+2 KO, Orai-TKO, Orai double KO × STIM KO, STIM rescues/mutants, ± Gd³⁺.  10 µM CCh, 2 mM Ca²⁺.
- Per-cell counts "oscillations in 14 min" (Figure 1I): WT 189, STIM1 KO 230, STIM2 KO 156, STIM1/2 KO 152, Orai-TKO 96 cells.
- The large time-stamped blocks (Figure 1B, 3A/F/K, 7, Sup 5A) are thapsigargin-SOCE and NFAT time courses, not oscillations.
So the ISI information comes from 230 *selected* cells (selection bias) in 14-min, strongly non-stationary records
(median ISI trend +15 % per minute).  σ_H, hold-out predictions and the σ–T_av law over a wide T_av range are not possible here.

## Results (this run)
- Within-cell statistics (128 trains ≥ 8 spikes, linear trend removed): CV 0.42 (ours 0.22), skew 0.17 (0.52),
  amplitude CV 0.25 (0.15), corr(ISI, amp) 0.32 (0.10), serial corr 0.09 (0.03).  Raw values are dominated by the
  trend (CV 0.87, serial 0.64).  Same ordering as ours (serial ≈ 0 after detrending, positive weak ISI–amplitude
  coupling) but everything is noisier: 8–15 spikes per train and a nonlinear slowing that linear detrending only partly removes.
- Dose steps within single cells (14 cells with trains at both 10 and 30 µM): T_av(10)/T_av(30) median 1.63
  [IQR 1.54–2.10] — our population medians give 1.96 — with amplitude ratio 1.01.  A clean within-cell fold change with an
  invariant spike amplitude: the stereotypy result (m0ip3h_amplitude_invariance) and Proposition 2 in an independent lab.
- Population frequency (Figure 1I): WT 6, STIM1 KO 6, STIM2 KO 7, STIM1/2 KO 6, Orai-TKO 3 oscillations per 14 min
  (medians): losing STIM does not slow the oscillations; losing all Orai halves the count.  Usable as a J_in perturbation
  test of the model (to be added: J_in sweep of the stereotyped regime → predicted period and CV).

## Outputs
`results/emrich_traces.csv` (all representative traces, long format), `results/emrich_train_stats.csv`,
`results/emrich_dose_steps.csv`, `results/emrich_counts_by_genotype.csv`,
`figures/emrich_replication_main.{png,pdf,eps}` (traces; within-cell statistics vs ours; 10→30 µM step),
`figures/emrich_replication_supp.{png,pdf,eps}` (counts by genotype; raw vs detrended statistics)
