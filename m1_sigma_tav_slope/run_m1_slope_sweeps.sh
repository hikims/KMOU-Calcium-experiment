#!/bin/bash
# σ–T_av law: cells differ in the recovery time scale ε_H; which noise convention gives sd ∝ T_av (log-log slope 1)?
#   M1 (ν common):   σ_H = ν/√ε_H  with ν = σ_H*·√1000 (σ_H* = 0.026 tuned in m1_isi_variability at ε_H = 1000)
#   M1 (σ_H common): σ_H = 0.026 for every cell
#   M2:              σ_P = 0.0167 (tuned to the same CV at ε_H = 1000)
#   Ca-noise:        σ_C = 0.39
REF="mech=2 rho=0.024 gH=30 dH=10000 nA=2 kA=0.4 Jin=0.02 kS=0.25 beta=4 alpha=20 qoff=0 rC=300 kR=1 S=1"
NUM="Tend=200000 Ttrans=3000 dt=0.01"      # >= 500 ISIs at the longest period; ~20 s per point
SW="sweep=epsH lo=200 hi=8000 n=9 log=1"
B=./m1_slope_sde
set -e; mkdir -p results
$B $SW $REF $NUM nuH=0.822    > results/m1_slope_M1_nu_common.csv
$B $SW $REF $NUM sigH=0.026   > results/m1_slope_M1_sig_common.csv
$B $SW $REF $NUM sigP=0.0167  > results/m1_slope_M2.csv
$B $SW $REF $NUM sigC=0.39    > results/m1_slope_Ca.csv
echo "done -> results/"
# the other axis of the data: the same cell at different stimuli (dose), ε_H = 1000 fixed; the first-passage CV
# depends on H* (S̄), so the slope along this axis is not 1 even with ν common
REFS="mech=2 rho=0.024 gH=30 dH=10000 nA=2 kA=0.4 Jin=0.02 kS=0.25 beta=4 alpha=20 qoff=0 rC=300 kR=1 epsH=1000"
$B sweep=S lo=0.05 hi=30 n=9 log=1 $REFS $NUM nuH=0.822  > results/m1_slope_M1_nu_common_dose.csv
$B sweep=S lo=0.05 hi=30 n=9 log=1 $REFS $NUM sigP=0.0167 > results/m1_slope_M2_dose.csv
echo "done (dose axis) -> results/"
