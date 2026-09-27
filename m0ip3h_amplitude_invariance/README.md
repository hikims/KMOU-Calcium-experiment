# m0ip3h_amplitude_invariance — spike stereotypy (the stimulus changes the period, not the spike)

## Question
In the data the ISI falls ~7x from 3 to 30 µM carbachol, but the spike itself hardly changes:
half-width falls only ~2x and the ΔF/F0 amplitude is flat; within a cell the amplitude is almost
uncorrelated with the preceding ISI (corr ≈ 0.09).  In the reference regime of the previous steps IP3
accelerates the whole H clock, so period, width and amplitude all fall together (~4x).
Which structure of the H kinetics makes the spike stereotyped?

## Model change (switchable; qoff = 1 reproduces the previous model exactly)
The H equation is written as a two-rate gate (identical to (Q − H)/ε_H):

    dH/dt = (1/ε_H) [ fP · fC · Q(C) (1 − H)  −  fPq · fC · (1 − Q(C)) H ]

    fP  = 1 + g_H P/(d_H + P)          IP3 accelerates recovery
    fPq = 1 + qoff · g_H P/(d_H + P)   fraction qoff of it also accelerates Ca-inactivation
    fC  = 1 + r_C C²/(k_R² + C²)       Ca accelerates the H kinetics; k_R is its own threshold
                                       (k_R = k_H reproduces the previous steps; k_R = 1 switches it off at rest)

Stereotyped reference `REF_ST`: previous reference with `qoff=0 kR=1 rC=300`
(rC is raised because with k_R = 1 the speed-up no longer acts at rest, and the inactivation alone has to
make the spike narrow).  `kP` (half-saturation of the IP3 gain of release) is also available; it does not
change the amplitude in this regime and is not used in the figure.

## Files
- `m0ip3h_invariance_sweep.c` — RK4 kernel; modes: single run / `sweep=` / `trace=1` / `recovery=1`
  (`recovery`: after a free spike, a Ca pulse `pulse` is applied at delay τ and the evoked amplitude is measured
  in a window that ends before the free spike would arrive)
- `run_m0ip3h_invariance_sweeps.sh` — stimulus sweeps for both regimes, qoff sweep at S̄ = 0.1 and 30,
  rC sweep at qoff = 0, recovery functions at S̄ = 1 → `results/`
- `plot_m0ip3h_amplitude_invariance.py` — two figures → `figures/m0ip3h_amplitude_invariance_{main,supp}.{png,pdf,eps}` (1x3, main)
  and `figures/m0ip3h_amplitude_invariance_recovery.{png,pdf,eps}` (1x2, to be merged into the m1_isi_variability figure);
  also `results/m0ip3h_inv_ratios_vs_qoff.csv`, `results/m0ip3h_inv_data_by_dose.csv`

## Run
    make && ./run_m0ip3h_invariance_sweeps.sh && python3 plot_m0ip3h_amplitude_invariance.py
(~15 s + a few seconds for the figure; needs the processed spike tables in `../../*Spike Trains_data/2/`)

## Panels
Main (one message: IP3 acts on recovery, not on the spike)
  (a) data: ISI, half-width, amplitude vs dose normalised to 3 µM
  (b) stereotyped regime: only T scales (the previous regime scales all three together, ratios printed by the script)
  (c) ratios S̄=0.1 / S̄=30 vs qoff with the data ratios (dotted)
Recovery (kept separate; belongs with the noise-location argument of m1_isi_variability)
  (a) recovery function A(τ)/A_free for both regimes (Ca-pulse-evoked spike) and the H-triggered level
  (b) data: within-cell amplitude vs preceding ISI

## Findings
- qoff = 0 (IP3 acts on recovery only) gives a dose-invariant spike: T ratio 4.0, W50 ratio 0.9, amplitude ratio 0.8
  (previous regime: 4.4 / 4.2 / 4.0).  The duty then rises with dose (0.024 → 0.10; data 0.053 → 0.14).
- Intermediate qoff interpolates: qoff ≈ 0.3–0.4 reproduces the data half-width ratio (2.3) but then the
  amplitude falls with it (ratio ≈ 2): width and amplitude cannot be separated by qoff alone.
  Release gain (alpha, kP) does not change the amplitude in the relaxation regime.
- Recovery function: a spike evoked by a Ca pulse before the free time has a smaller amplitude in BOTH regimes
  (graded in H: A/A_free ≈ 0.5 → 1 over τ/T_wait 0.7 → 1).  Hence a within-cell corr(ISI, amp) ≈ 0 (data)
  is expected only if early spikes are triggered through H (H reaches its threshold H_on at a random time,
  so every spike starts from H ≈ H_on), not through Ca fluctuations.  This is the argument for putting the
  M1 noise in H (next computation: m1_isi_variability).

## Figures
- `figures/m0ip3h_amplitude_invariance_main.{png,pdf,eps}` (paper): data normalised; stereotyped model; ratios vs q_off
- `figures/m0ip3h_amplitude_invariance_supp.{png,pdf,eps}` (supplement): recovery function; amplitude vs preceding ISI

## Update 2026-09-22 (data source)
Dose medians are medians of train medians (`../data_dose_units/results/dose_table.csv`); train correlations use
the detrended ISI.  ΔF/F0 amplitude differs 3× between records (indicator loading, F0), so the amplitude is now
compared **within cells across the step levels** (`results/m0ip3h_inv_within_cell_steps.csv`): 5→10 0.94 (n = 10),
10→15 1.08 (12), 5→30 1.26 (4), 5→15 1.06 (9), with W50 ratios 1.02–1.10 and ISI ratios 0.85–1.79.  Between doses:
ISI 3/30 µM 6.95, W50 2.30.  Panel (a) shows ISI and W50 between doses and the amplitude chain within cells;
panel (c) uses the within-cell 5/30 amplitude ratio.  corr(ISI, amp) within trains 0.05 [−0.09, 0.20].
