#!/bin/bash
cd "$(dirname "$0")"
for c in check5_mean_shift check6_bootstrap check7_gate check4_scan_candidate check3_jin check2_timestep check1_cv_matched; do
  echo "===== $c $(date +%H:%M:%S)"; python3 $c.py > results/$c.log 2>&1 && tail -3 results/$c.log || { echo "FAILED $c"; tail -20 results/$c.log; }
done
echo "ALL DONE $(date +%H:%M:%S)"
