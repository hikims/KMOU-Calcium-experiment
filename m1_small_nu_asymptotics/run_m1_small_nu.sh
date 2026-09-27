#!/bin/bash
# Small-ν asymptotics of the first-passage moments: high-resolution quadrature, grid convergence, Monte-Carlo,
# FPT density, and the limit ratio skew/CV.  H0 = H_min = 0.064, H* = H_on = 0.408 (stereotyped regime, S̄ = 1).
H0=0.064; HS=0.408; B=./m1_fpt_asymptotics
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)}
set -e; mkdir -p results; echo "threads: $OMP_NUM_THREADS"
python3 derive_small_nu.py
$B sweep=nu lo=0.01 hi=3 n=31 log=1 H0=$H0 Hs=$HS                       > results/m1_asym_vs_nu.csv
echo "ny,nw,nu,H0,Hs,mean,sd,CV,skew,mean_det,mean_as,sd_as,skew_as,ratio_as" > results/m1_asym_convergence.csv
for g in "2001 1000" "4001 2000" "8001 4000" "16001 8000"; do set -- $g
  $B sweep=nu lo=0.02 hi=0.2 n=3 log=1 H0=$H0 Hs=$HS ny=$1 nw=$2 | tail -n +2 | sed "s/^/$1,$2,/" >> results/m1_asym_convergence.csv; done
echo "kind,N,mean,sd,CV,skew,nu" > results/m1_asym_mc.csv
for nu in 0.05 0.1 0.2 0.3 0.5 1.0 2.0; do $B mc=200000 nu=$nu H0=$H0 Hs=$HS dt=1e-4 | sed -n 2p | sed "s/\$/,$nu/" >> results/m1_asym_mc.csv; done
$B mc=400000 nu=0.3 H0=$H0 Hs=$HS dt=1e-4 hist=1 nbins=80 | tail -n +3 > results/m1_asym_hist_nu0.3.csv
$B mc=400000 nu=1.0 H0=$H0 Hs=$HS dt=1e-4 hist=1 nbins=80 | tail -n +3 > results/m1_asym_hist_nu1.0.csv
for h0 in 0.02 0.064 0.15; do $B ratio=1 H0=$h0 Hslo=0.15 Hshi=0.9 n=41 > results/m1_asym_ratio_vs_Hs_H0$h0.csv; done
echo "done -> results/"
