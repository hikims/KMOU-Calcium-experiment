#!/bin/bash
# First-passage moments of the reduced recovery process (Dynkin quadrature), Monte-Carlo checks of the 1-D SDE,
# and the full 3-variable SDE for the reduction check.  H0 = H_min = 0.064, H* = H_on = 0.408 (stereotyped regime, S̄ = 1).
H0=0.064; HS=0.408
F=./m1_fpt_jacobi; B=./m1_full_sde
REF="mech=2 rho=0.024 gH=30 dH=10000 nA=2 kA=0.4 Jin=0.02 kS=0.25 beta=4 alpha=20 qoff=0 rC=300 kR=1 S=1"
set -e; mkdir -p results
$F sweep=nu lo=0.02 hi=3 n=25 log=1 H0=$H0 Hs=$HS                 > results/m1_fpt_vs_nu.csv
for nu in 0.1 0.3 1.0; do $F sweep=Hs lo=0.15 hi=0.85 n=29 H0=$H0 nu=$nu > results/m1_fpt_vs_Hs_nu$nu.csv; done
echo "kind,N,mean,sd,CV,skew,nu" > results/m1_fpt_mc.csv
for nu in 0.1 0.3 1.0 2.0; do $F mc=20000 nu=$nu H0=$H0 Hs=$HS dt=1e-4 | tail -1 | sed "s/\$/,$nu/" >> results/m1_fpt_mc.csv; done
# full model: sigma_H sweep at eps_H = 1000 and eps_H sweep with nu common (as in m1_isi_variability / m1_sigma_tav_slope)
$B sweep=sigH lo=0.003 hi=0.3 n=13 log=1 $REF epsH=1000 Tend=120000 Ttrans=3000 dt=0.01 > results/m1_full_vs_sigH.csv
$B sweep=epsH lo=200 hi=8000 n=9 log=1 $REF nuH=0.822 Tend=200000 Ttrans=3000 dt=0.01  > results/m1_full_vs_epsH.csv
echo "done -> results/"
