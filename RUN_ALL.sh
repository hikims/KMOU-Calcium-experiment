#!/bin/bash
# One pass over everything that changed on 2026-09-22.  Run from code/ (the directory that holds the folders).
# Order matters: data_dose_units feeds the others; assemble_paper_figures comes last.
# Needs: a C compiler (macOS: xcode-select --install; OpenMP via brew install libomp); python3 with numpy, pandas,
# scipy, matplotlib, openpyxl.
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)}
step(){ echo; echo "=== $1 ==="; }
# build DIR TARGET : compile DIR/TARGET.c -> DIR/TARGET if the binary is missing (Makefile if present, else cc directly)
build(){ ( cd "$1" && if [ -x "$2" ]; then exit 0; fi
           if [ -f Makefile ]; then make -s; exit 0; fi
           cc -O2 -fopenmp -o "$2" "$2.c" -lm 2>/dev/null || cc -O2 -Xpreprocessor -fopenmp -o "$2" "$2.c" -lm -lomp 2>/dev/null || { echo "  (no OpenMP: serial build of $2)"; cc -O2 -o "$2" "$2.c" -lm; } ); }
[ -d data_dose_units ] || { echo "run this from the code/ directory (data_dose_units not found here)"; exit 1; }
step "data_dose_units (tables, segments, dose map)";        (cd data_dose_units && python3 data_dose_units.py)
step "m0_joint_scan (plot only; scan results kept)";       build m0_joint_scan m0_joint_scan; (cd m0_joint_scan && python3 plot_m0_joint_scan.py)
step "m0ip3h_stimulus_encoding (plot)";                     build m0ip3h_stimulus_encoding m0ip3h_stimulus_sweep; (cd m0ip3h_stimulus_encoding && python3 plot_m0ip3h_stimulus_encoding.py)
step "kernels of waveform_duty / amplitude_observation";    build m0ip3h_waveform_duty m0ip3h_waveform_sweep; build m0ip3h_amplitude_observation m0ip3h_amplitude_sweep
step "m0ip3h_amplitude_invariance (plot)";                  build m0ip3h_amplitude_invariance m0ip3h_invariance_sweep; (cd m0ip3h_amplitude_invariance && python3 plot_m0ip3h_amplitude_invariance.py)
step "m0ip3h_bifurcation_structure (plot)";                 (cd m0ip3h_bifurcation_structure && python3 plot_m0ip3h_bifurcation_structure.py)
step "m1_first_passage_theory (plot + error table)";        build m1_first_passage_theory m1_fpt_jacobi; (cd m1_first_passage_theory && python3 plot_m1_first_passage_theory.py)
step "m1_isi_variability (re-tunes sigma, ~5 min)";         build m1_isi_variability m1_isi_sde; (cd m1_isi_variability && python3 plot_m1_isi_variability.py)
step "m1_sigma_tav_slope (dose-axis sweeps, ~3 min)";       build m1_sigma_tav_slope m1_slope_sde; (cd m1_sigma_tav_slope && ./run_m1_slope_sweeps.sh && python3 plot_m1_sigma_tav_slope.py)
step "emrich_jin_prediction (new T0, ~3 min)";              build emrich_jin_prediction emrich_jin_sde; (cd emrich_jin_prediction && python3 run_emrich_jin.py 400 0.45 && python3 plot_emrich_jin_prediction.py)
step "assemble_paper_figures (Fig 1-7 + supplement)";       (cd assemble_paper_figures && rm -f figures/*.png figures/*.pdf figures/*.eps && python3 assemble_paper_figures.py all)
echo; echo "done.  Figures: assemble_paper_figures/figures/ ; failures (if any): assemble_paper_figures/figures/assemble_log.txt"
