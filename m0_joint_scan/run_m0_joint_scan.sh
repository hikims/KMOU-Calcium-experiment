#!/bin/bash
# Joint Latin-hypercube scan of the original AIMS M0 equations (all cores via OpenMP).
#   NPTS  sample points (default 4000; ~0.5 s per point per thread)
#   NS    stimulus grid points in [0.03, 30] (log-spaced)
NPTS=${1:-4000}; NS=${2:-9}
set -e; mkdir -p results
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)}
echo "threads: $OMP_NUM_THREADS, points: $NPTS"
./m0_joint_scan npts=$NPTS ns=$NS seed=1 > results/m0_joint_scan.csv
# reference curve of the Table 1 baseline on the same stimulus grid, and its trace
./m0_joint_scan curve=1 ns=25 > results/m0_joint_baseline_curve.csv
echo "done -> results/m0_joint_scan.csv"
