#!/bin/bash
# Amplitude of the relaxation regime vs stimulus (C_peak, C_rest) and gain sweeps that change amplitude.
REF="mech=2 rho=0.024 gH=30 dH=10000 epsH=1000 nA=2 kA=0.4 Jin=0.02 kS=0.25 rC=30 beta=4 alpha=20"   # kR is not a key of this binary: it uses kR = kH = 0.15 (REF_PREV)
NUM="Tend=15000 Ttrans=9000 dt=0.02"
B=./m0ip3h_amplitude_sweep
set -e; mkdir -p results
$B sweep=S     lo=0.1 hi=30 n=13 log=1 $REF $NUM           > results/m0ip3h_amp_period_vs_stimulus_ref.csv
$B sweep=alpha lo=8   hi=40 n=13 log=1 $REF $NUM S=1       > results/m0ip3h_amp_vs_alpha_S1.csv
$B sweep=beta  lo=2   hi=8  n=13 log=1 $REF $NUM S=1       > results/m0ip3h_amp_vs_beta_S1.csv
$B sweep=kS    lo=0.1 hi=1  n=13 log=1 $REF $NUM S=1       > results/m0ip3h_amp_vs_kS_S1.csv
echo "done -> results/"
