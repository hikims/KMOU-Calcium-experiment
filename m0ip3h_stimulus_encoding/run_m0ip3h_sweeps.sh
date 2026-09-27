#!/bin/bash
# M0 + IP3->H coupling: does the stimulus encode the period?
#   mech=1 threshold coupling  kH(P) = kH (dH + gH P)/(dH + P)   (neutral gH=1)
#   mech=2 rate coupling       dH/dt = (1 + gH P/(dH+P)) (Q-H)/epsH (neutral gH=0)
# Results -> results/ ; then: python3 search_m0ip3h_dose_range.py ; python3 plot_m0ip3h_stimulus_encoding.py
set -e
mkdir -p results
B=./m0ip3h_stimulus_sweep
SW="sweep=S lo=0.1 hi=30 n=13 log=1"
# reference: original model (no coupling)
$B $SW                                        > results/m0ip3h_period_vs_stimulus_mech0.csv
# threshold coupling, Table-1 lambda (P ~ 10^3) and lambda=1 (P ~ 1)
$B $SW mech=1 gH=5  dH=3000                   > results/m0ip3h_period_vs_stimulus_mech1_gH5_dH3000.csv
$B $SW mech=1 gH=5  dH=3    lambda=1          > results/m0ip3h_period_vs_stimulus_mech1_gH5_dH3_lambda1.csv
# rate coupling: Table-1 rho, then reduced rho (wider IP3 range), then larger epsH (stay oscillatory)
$B $SW mech=2 gH=5  dH=3000                   > results/m0ip3h_period_vs_stimulus_mech2_gH5_dH3000.csv
$B $SW mech=2 gH=10 dH=3000 rho=0.024         > results/m0ip3h_period_vs_stimulus_mech2_gH10_dH3000_rho0.024.csv
$B $SW mech=2 gH=10 dH=3000 rho=0.024 epsH=300 > results/m0ip3h_period_vs_stimulus_mech2_gH10_dH3000_rho0.024_epsH300.csv
# gH sweep at fixed low/high stimulus (rate coupling)
$B sweep=gH lo=0 hi=30 n=16 mech=2 dH=3000 rho=0.024 epsH=300 S=0.3  > results/m0ip3h_period_vs_gH_S0.3.csv
$B sweep=gH lo=0 hi=30 n=16 mech=2 dH=3000 rho=0.024 epsH=300 S=30   > results/m0ip3h_period_vs_gH_S30.csv
echo "done -> results/"
