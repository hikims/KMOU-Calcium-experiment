#!/bin/bash
# M0 (deterministic AIMS reaction system, Table 1 baseline): parameter sweeps.
# Question: can M0 produce (i) a stimulus-dependent period and (ii) a spike-like waveform?
# Results -> results/ ; then:  python3 plot_m0_period_waveform.py
set -e
mkdir -p results
B=./m0_period_waveform_sweep
$B sweep=S      lo=0     hi=10  n=21              > results/m0_period_vs_stimulus_lambda0.001.csv
$B sweep=S      lo=0     hi=10  n=21 lambda=1     > results/m0_period_vs_stimulus_lambda1.csv
$B sweep=lambda lo=0.001 hi=10  n=13 log=1 S=1    > results/m0_period_vs_lambda_S1.csv
$B sweep=epsH   lo=10    hi=400 n=25 log=1        > results/m0_period_vs_epsH.csv
$B sweep=beta   lo=0.05  hi=0.35 n=31             > results/m0_period_vs_beta.csv
$B trace=1 Tend=4000 Ttrans=3000                  > results/m0_waveform_trace_baseline_epsH70.csv
$B trace=1 Tend=4000 Ttrans=3000 epsH=20          > results/m0_waveform_trace_epsH20.csv
echo "done -> results/"
