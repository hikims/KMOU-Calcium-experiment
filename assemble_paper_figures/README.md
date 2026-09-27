# assemble_paper_figures — paper figures 1–8 and supplement S1–S24 from the folders' panel functions

`assemble_paper_figures.py [main|supp|all] [--serial]` executes each folder's plot script inside its folder
(`compose` replaced by a no-op, `__name__ != "__main__"`), takes the panel functions from the namespace and
composes the figures below; one process per figure (fork pool, all cores), a folder's script is loaded once per
process.  Output: `figures/fig<k>_main_<topic>.{png,pdf,eps}` and `figures/figS<k>_supp_<folder>_<main|supp>.*`.
The panel code lives only in the computation folders; this script contains no plotting of its own.

| Fig | panels (folder : function) | message |
|---|---|---|
| 1 | data_dose_units:p_dose_units, m0_period_waveform:p_period_vs_stimulus, m0_period_waveform:p_one_cycle, m0_joint_scan:p_ratio_vs_duty | data facts; the original equations fail jointly |
| 2 | m0ip3h_stimulus_encoding:p_direction, m0ip3h_waveform_duty:p_duty_vs_beta, p_ref_vs_stimulus, m0ip3h_bifurcation_structure:p_map_S_epsH | rate coupling, pump, relaxation regime |
| 3 | m0ip3h_amplitude_invariance:data_norm, model_norm, ratios_qoff, m0ip3h_bifurcation_structure:p_map_beta_qoff | stereotypy: q_off |
| 4 | m1_first_passage_theory:moments, reduction, m1_small_nu_asymptotics:p_mean_sd, p_skew | reduction theorem, closed moments, small-ν corollary |
| 5 | m1_isi_variability:trace_panel, isi_dist, noise_loc | H-noise and the discriminating statistics |
| 6 | m1_sigma_tav_slope:law, slopes | variability law |
| 7 | data_qc_splits:example_trace, data_dose_units:p_train_strip, step_fold_change:per_cell, data_dose_units:p_stats_units | Azzoni data |
| 8 | emrich_replication:p_examples, p_stats, p_steps, emrich_jin_prediction:p_counts | independent data and the J_in prediction |

Supplement: every folder's SUPP list, plus the MAIN lists of m0_joint_scan, m0ip3h_amplitude_observation,
m1_particle_filter_synthetic, data_qc_splits, step_fold_change, emrich_jin_prediction, m0ip3h_bifurcation_structure.
To change a figure edit MAIN_FIGS / SUPP_ORDER at the top of the script.

## Run
    python3 assemble_paper_figures.py all          (~1 min on 2 processes; the folders' results/ must exist)
Requires the sibling folders under `../` with their `results/` (and the compiled binaries where a panel re-runs one).
