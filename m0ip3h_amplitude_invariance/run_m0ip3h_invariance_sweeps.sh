#!/bin/bash
# Spike stereotypy: does the stimulus change only the period (data) or the whole spike (model)?
#   REF_PREV : reference regime of the previous steps (IP3 rescales the whole H clock, qoff=1)
#   REF_ST   : stereotyped-spike regime: IP3 accelerates recovery only (qoff=0); Ca-dependent
#              inactivation with its own threshold kR=1 (off at rest) and speed rC=300
COMMON="mech=2 rho=0.024 gH=30 dH=10000 epsH=1000 nA=2 kA=0.4 Jin=0.02 kS=0.25 beta=4 alpha=20"
REF_PREV="$COMMON qoff=1 rC=30 kR=0.15"
REF_ST="$COMMON qoff=0 rC=300 kR=1"
BASE_ST="$COMMON rC=300 kR=1"          # qoff varied
NUM="Tend=15000 Ttrans=9000 dt=0.02"
B=./m0ip3h_invariance_sweep
set -e; mkdir -p results
# period / width / amplitude vs stimulus
$B sweep=S lo=0.1 hi=30 n=13 log=1 $REF_PREV $NUM            > results/m0ip3h_inv_vs_stimulus_ref_prev.csv
$B sweep=S lo=0.1 hi=30 n=13 log=1 $REF_ST   $NUM            > results/m0ip3h_inv_vs_stimulus_ref_st.csv
# the split parameter qoff at the two ends of the stimulus range
$B sweep=qoff lo=0 hi=1 n=11 $BASE_ST $NUM S=0.1             > results/m0ip3h_inv_vs_qoff_S0.1.csv
$B sweep=qoff lo=0 hi=1 n=11 $BASE_ST $NUM S=30              > results/m0ip3h_inv_vs_qoff_S30.csv
# inactivation speed rC (sets the spike width, hence the duty) at qoff=0
$B sweep=rC lo=30 hi=3000 n=11 log=1 $COMMON qoff=0 kR=1 $NUM S=0.1 > results/m0ip3h_inv_vs_rC_S0.1.csv
$B sweep=rC lo=30 hi=3000 n=11 log=1 $COMMON qoff=0 kR=1 $NUM S=30  > results/m0ip3h_inv_vs_rC_S30.csv
# recovery function: amplitude of a spike evoked by a Ca pulse at delay tau after a free spike
$B recovery=1 ntau=25 pulse=1 $REF_PREV $NUM S=1             > results/m0ip3h_inv_recovery_ref_prev.csv
$B recovery=1 ntau=25 pulse=1 $REF_ST   $NUM S=1             > results/m0ip3h_inv_recovery_ref_st.csv
echo "done -> results/"
