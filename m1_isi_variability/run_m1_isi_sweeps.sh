#!/bin/bash
# M1: ISI statistics vs noise amplitude for three noise locations (H: M1, P: M2, Ca: comparison), at S̄ = 1,
# on the stereotyped reference regime of m0ip3h_amplitude_invariance.
REF="mech=2 rho=0.024 gH=30 dH=10000 epsH=1000 nA=2 kA=0.4 Jin=0.02 kS=0.25 beta=4 alpha=20 qoff=0 rC=300 kR=1"
NUM="Tend=120000 Ttrans=3000 dt=0.01"     # ~1400 ISIs per point; ~10 s per point
B=./m1_isi_sde
set -e; mkdir -p results
$B sweep=sigH lo=0.003 hi=0.3   n=13 log=1 $REF $NUM S=1 > results/m1_isi_vs_sigH.csv
$B sweep=sigC lo=0.03  hi=1.0   n=13 log=1 $REF $NUM S=1 > results/m1_isi_vs_sigC.csv
$B sweep=sigP lo=0.002 hi=0.05  n=13 log=1 $REF $NUM S=1 > results/m1_isi_vs_sigP.csv
echo "done -> results/  (the plot script picks the noise amplitude matching the data CV and writes the ISI lists / traces)"
