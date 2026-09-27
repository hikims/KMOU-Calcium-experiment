#!/bin/bash
# Everything from scratch, in dependency order (~2 h on 8 cores).  Run from code/.
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)}
step(){ echo; echo "=== $1 ==="; }
build(){ ( cd "$1" && if [ -f Makefile ]; then make -s; else cc -O2 -fopenmp -o "$2" "$2.c" -lm 2>/dev/null || cc -O2 -Xpreprocessor -fopenmp -o "$2" "$2.c" -lm -lomp 2>/dev/null || cc -O2 -o "$2" "$2.c" -lm; fi ); }
[ -d data_dose_units ] || { echo "run this from the code/ directory"; exit 1; }
step data_dose_units;               (cd data_dose_units && python3 data_dose_units.py)
step data_qc_splits;                (cd data_qc_splits && python3 qc_data_splits.py)
step m0_period_waveform;            build m0_period_waveform m0_period_waveform_sweep; (cd m0_period_waveform && ./run_m0_sweeps.sh && python3 plot_m0_period_waveform.py)
step m0_joint_scan;                 build m0_joint_scan m0_joint_scan; (cd m0_joint_scan && ./run_m0_joint_scan.sh 4000 && python3 plot_m0_joint_scan.py)
step m0ip3h_stimulus_encoding;      build m0ip3h_stimulus_encoding m0ip3h_stimulus_sweep; (cd m0ip3h_stimulus_encoding && ./run_m0ip3h_sweeps.sh && python3 search_m0ip3h_dose_range.py && python3 plot_m0ip3h_stimulus_encoding.py)
step m0ip3h_waveform_duty;          build m0ip3h_waveform_duty m0ip3h_waveform_sweep; (cd m0ip3h_waveform_duty && ./run_m0ip3h_waveform_sweeps.sh && python3 search_m0ip3h_waveform_duty.py && python3 plot_m0ip3h_waveform_duty.py)
step m0ip3h_amplitude_observation;  build m0ip3h_amplitude_observation m0ip3h_amplitude_sweep; (cd m0ip3h_amplitude_observation && ./run_m0ip3h_amplitude_sweeps.sh && python3 plot_m0ip3h_amplitude_observation.py)
step m0ip3h_amplitude_invariance;   build m0ip3h_amplitude_invariance m0ip3h_invariance_sweep; (cd m0ip3h_amplitude_invariance && ./run_m0ip3h_invariance_sweeps.sh && python3 plot_m0ip3h_amplitude_invariance.py)
step m0ip3h_bifurcation_structure;  build m0ip3h_bifurcation_structure m0ip3h_bifurcation_grid; (cd m0ip3h_bifurcation_structure && ./run_m0ip3h_bifurcation.sh && python3 plot_m0ip3h_bifurcation_structure.py)
step "data_dose_units again (dose map needs the bifurcation line)"; (cd data_dose_units && python3 data_dose_units.py)
step m1_first_passage_theory;       build m1_first_passage_theory m1_fpt_jacobi; build m1_first_passage_theory m1_full_sde; (cd m1_first_passage_theory && ./run_m1_fpt.sh && python3 plot_m1_first_passage_theory.py)
step m1_small_nu_asymptotics;       build m1_small_nu_asymptotics m1_fpt_asymptotics; (cd m1_small_nu_asymptotics && ./run_m1_small_nu.sh && python3 plot_m1_small_nu_asymptotics.py)
step m1_isi_variability;            build m1_isi_variability m1_isi_sde; (cd m1_isi_variability && ./run_m1_isi_sweeps.sh && python3 plot_m1_isi_variability.py)
step m1_sigma_tav_slope;            build m1_sigma_tav_slope m1_slope_sde; (cd m1_sigma_tav_slope && ./run_m1_slope_sweeps.sh && python3 plot_m1_sigma_tav_slope.py)
step m1_particle_filter_synthetic;  build m1_particle_filter_synthetic m1_pf_sde; (cd m1_particle_filter_synthetic && python3 run_pf_synthetic.py && python3 plot_m1_particle_filter_synthetic.py)
step step_fold_change;              (cd step_fold_change && python3 step_fold_change.py)
step emrich_replication;            (cd emrich_replication && python3 emrich_replication.py)
step emrich_jin_prediction;         build emrich_jin_prediction emrich_jin_sde; (cd emrich_jin_prediction && python3 run_emrich_jin.py 400 0.45 && python3 plot_emrich_jin_prediction.py)
step assemble_paper_figures;        (cd assemble_paper_figures && rm -f figures/*.png figures/*.pdf figures/*.eps && python3 assemble_paper_figures.py all)
echo; echo "done -> assemble_paper_figures/figures/"
