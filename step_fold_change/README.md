# step_fold_change — fold-change encoding inside single cells (Proposition 2, gate 3)

Python only.  Run: `python3 step_fold_change.py` (needs the spike tables `../../*Spike Trains_data/2/*.csv`).

## What it computes
Step records (authors' tables, base files): 5→10→15 µM (03-07 E1, 03-14 E2, 03-19 E1), 5→10 (03-12 E1), 5→30 (04-03 E1, E2).
Per train with ≥ 5 ISIs at two consecutive levels: global ratio T_av(low)/T_av(high), local ratio (last 5 ISIs / first 5
ISIs across the step), within-level ISI trend (relative change per hour), and the comparison with the constant-dose
medians (prediction 1.27 for 5→10, 1.05 for 10→15, 2.5 for 5→30).

## Result (gate 3 outcome)
| step  | cells | global ratio median [IQR] | local ratio | constant-dose prediction |
|-------|------:|---------------------------|-------------|--------------------------|
| 5→10  | 10    | 0.89 [0.75, 1.13]         | 0.91        | 1.27 |
| 10→15 | 12    | 1.11 [1.02, 1.30]         | 1.04        | 1.05 |
| 5→30  | 4     | 1.78 [1.73, 1.89]         | 1.20        | 2.48 |

- In the 3-h records the levels are separated by 30–40 min without spikes (e.g. 03-07 E1: 5 µM spikes end at 65 min,
  10 µM spikes start at 105 min), so these are not clean steps; the 5→10 step does not shorten the period at all, and
  the ratio correlates with the cell's own T_av (corr(log T_low, log ratio) 0.55–0.9) — the opposite of Proposition 2.
- The 5→30 records (0.1 s sampling, 2000 s) show strong within-level trends (−5 to −7 relative change per hour ≈ fast
  drift/bleaching) and only 4 usable cells; the local ratio (1.20) is far below the global one (1.78).
- Consequence for the plan: gate 3 cannot be passed with these step records.  As foreseen in §12, Proposition 2 is to be
  tested on the constant-dose populations (dose→T_av curve, m0ip3h_stimulus_encoding) and on the digitised Thurley 2014
  fold-change data; the HEK293 step records enter the paper as the desensitisation/limitation result (M1d).

## Outputs
`results/step_fold_change_cells.csv`, `results/step_fold_change_summary.csv`, `figures/step_fold_change_{main,supp}.{png,pdf,eps}` (1×3)

## Figures
- `figures/step_fold_change_main.{png,pdf,eps}` (paper): fold change per cell; step ratios
- `figures/step_fold_change_supp.{png,pdf,eps}` (supplement): example cell
