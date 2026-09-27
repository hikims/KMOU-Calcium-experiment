# emrich_jin_prediction — what M1 predicts when the basal influx (SOCE) is reduced; Emrich 2021 counts

## Question
Emrich et al. (2021) count CCh-evoked spikes over 14 min (10 µM): WT median 6 [3, 8], STIM1 KO 6, STIM2 KO 7,
STIM1/2 KO 6, ORAI-TKO 3 [2, 5] (`../emrich_replication`).  In the stereotyped reference the only influx term
is the basal J_in.  What does M1 (H-noise, ν common) predict for the count, period, amplitude, width and ISI-CV
when J_in is reduced, and which reduction reproduces the ORAI-TKO count?

## Method (all binary calls OpenMP-parallel: sweep values / cells on separate threads with their own RNG)
1. deterministic period vs S̄ at J_in = 0.02
2. stimulus calibration on a population: N cells with ε_H lognormal (cv 0.45) and ν = 0.822 common, rest at
   S̄ = 0 for 3000 units, then the stimulus is switched on and spikes are counted for 14 min = 1300 model units
   (T0 = 0.527 s/unit from `../data_dose_units`); S̄_E is the stimulus whose population mean count equals the WT
   mean (6.0) → S̄_E = 0.13 (the Azzoni endpoint map gives S̄(10 µM) = 0.21; same order)
3. J_in sweep at S̄_E: stationary period, CV, amplitude, width (60 000 units per value)
4. populations for J_in factors 1, 0.75, 0.5, 0.35, 0.25, 0.1, 0 (300 cells each) and one trace per factor

## Result
- The count is insensitive to a partial loss of influx: factor 1 → 0.1 moves the median from 5 [4, 7] to 4 [3, 5];
  at J_in = 0 every cell is silent within 14 min (deterministic period ≈ 13 min plus latency).
- The ORAI-TKO median (3) is reached only at a J_in factor ≈ 0.07 (`results/jin_prediction_table.csv`), i.e. the
  model needs an almost complete loss of basal influx, with a graded but steep transition between 0.1 and 0.
- Predicted TKO/WT ratios at that factor: period 1.39, amplitude 1.22, W50 1.10, ISI-CV 1.04 — larger and
  slightly wider spikes at a lower rate, same relative variability.  STIM-single KOs (count unchanged) correspond
  to factors ≥ 0.35.  The population ISI-CV at S̄_E is 0.31–0.4 (Emrich WT trains, detrended: 0.42).
- Not in the model: store depletion (a lower ER content in TKO cells would reduce the amplitude, opposite to the
  prediction) — the amplitude ratio is therefore the discriminating observable.

## Files
- `emrich_jin_sde.c` — M1 kernel with `sweep=`, `cells=N` (step protocol, ε_H spread), `trace=1`
- `run_emrich_jin.py [NCELLS] [EPSH_CV]` → `results/jin_S_sweep_det.csv`, `jin_S_calibration_scan.csv`,
  `jin_calibration.csv`, `jin_sweep_S_E.csv`, `jin_cells_f*.csv`, `jin_trace_f*.csv`, `jin_counts_summary.csv`
- `plot_emrich_jin_prediction.py` → `results/jin_prediction_table.csv`; `figures/emrich_jin_prediction_main.*`
  (a: sweep ratios, b: counts vs Emrich, c: predicted ratios), `figures/emrich_jin_prediction_supp.*`
  (a: traces, b: count histograms, c: latency, d: silent fraction)

## Run
    make && python3 run_emrich_jin.py 400 0.45 && python3 plot_emrich_jin_prediction.py     (~5 min on 2 threads)
