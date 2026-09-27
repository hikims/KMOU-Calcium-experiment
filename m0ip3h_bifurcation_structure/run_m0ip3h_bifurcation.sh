#!/bin/bash
# Bifurcation structure of the stereotyped reference REF_ST (and the previous reference for comparison):
# oscillation region, fixed-point stability and cycle features on (S̄, ε_H) and (β, qoff) grids; 1-D lines.
COMMON="mech=2 rho=0.024 gH=30 dH=10000 nA=2 kA=0.4 Jin=0.02 kS=0.25 alpha=20"
REF_ST="$COMMON epsH=1000 beta=4 qoff=0 rC=300 kR=1"
REF_PREV="$COMMON epsH=1000 beta=4 qoff=1 rC=30"
NUM="Tend=15000 Ttrans=9000 dt=0.02"
B=./m0ip3h_bifurcation_grid
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)}
set -e; mkdir -p results; echo "threads: $OMP_NUM_THREADS"
$B grid=1 px=S xlo=0.01 xhi=100 nx=41 logx=1 py=epsH ylo=10 yhi=10000 ny=41 logy=1 $REF_ST $NUM   > results/bif_grid_S_epsH_ref_st.csv
$B grid=1 px=beta xlo=0.5 xhi=40 nx=41 logx=1 py=qoff ylo=0 yhi=1 ny=21 logy=0 $REF_ST $NUM S=1   > results/bif_grid_beta_qoff_S1.csv
$B line=1 px=S xlo=0.001 xhi=1000 nx=121 logx=1 $REF_ST $NUM                                   > results/bif_line_S_ref_st.csv
$B line=1 px=S xlo=0.001 xhi=1000 nx=121 logx=1 $REF_PREV $NUM                                 > results/bif_line_S_ref_prev.csv
$B line=1 px=epsH xlo=10 xhi=10000 nx=61 logx=1 $REF_ST $NUM S=1                               > results/bif_line_epsH_ref_st.csv
$B line=1 px=beta xlo=0.5 xhi=40 nx=61 logx=1 $REF_ST $NUM S=1                                 > results/bif_line_beta_ref_st.csv
echo "done -> results/"
