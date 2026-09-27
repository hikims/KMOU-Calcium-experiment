#!/bin/bash
# Waveform (duty cycle) of the AIMS reaction system + IP3->H rate coupling.
# Three ingredients, each off at its default (so the binary reproduces earlier results exactly):
#   nA,kA  cooperative Ca activation C^nA/(kA^nA+C^nA)   (nA=1,kA=1 original)
#   Jin    basal Ca influx                                 (0 original)
#   rC     Ca-dependent speed-up of H kinetics             (0 original)
# Reference relaxation regime found by hand (2026-09-21):
REF="mech=2 rho=0.024 gH=30 dH=10000 epsH=1000 nA=2 kA=0.4 Jin=0.02 kS=0.25 rC=30 beta=4 alpha=20"   # kR is not a key of this binary: it uses kR = kH = 0.15 (REF_PREV)
NUM="Tend=15000 Ttrans=9000 dt=0.02"
B=./m0ip3h_waveform_sweep
set -e; mkdir -p results
# which ingredient narrows the spike? (S=1, one ingredient at a time from the encoding config)
ENC="mech=2 rho=0.024 gH=30 dH=10000 epsH=1000 S=1 $NUM"
$B sweep=rC   lo=0   hi=60 n=13       $ENC                       > results/m0ip3h_duty_vs_rC_original_activation.csv
$B sweep=rC   lo=0   hi=60 n=13       $ENC nA=2 kA=0.4 Jin=0.02  > results/m0ip3h_duty_vs_rC_cooperative.csv
$B sweep=beta lo=0.3 hi=8  n=16 log=1 $ENC nA=2 kA=0.4 Jin=0.02 rC=30 alpha=20 > results/m0ip3h_duty_vs_beta.csv
$B sweep=alpha lo=3 hi=40  n=13 log=1 $ENC nA=2 kA=0.4 Jin=0.02 rC=30 beta=4   > results/m0ip3h_duty_vs_alpha.csv
$B sweep=Jin  lo=0.003 hi=0.1 n=13 log=1 $ENC nA=2 kA=0.4 rC=30 beta=4 alpha=20 > results/m0ip3h_duty_vs_Jin.csv
# reference regime: stimulus sweep (encoding preserved?) and gain sweep
$B sweep=S  lo=0.1 hi=30 n=13 log=1 $REF $NUM                    > results/m0ip3h_ref_period_vs_stimulus.csv
$B sweep=gH lo=3   hi=300 n=13 log=1 $REF $NUM S=0.1             > results/m0ip3h_ref_period_vs_gH_S0.1.csv
$B sweep=gH lo=3   hi=300 n=13 log=1 $REF $NUM S=30              > results/m0ip3h_ref_period_vs_gH_S30.csv
# reference regime traces at low / high stimulus
$B trace=1 $REF $NUM S=0.1                                       > results/m0ip3h_ref_trace_S0.1.csv
$B trace=1 $REF $NUM S=30                                        > results/m0ip3h_ref_trace_S30.csv
echo "done -> results/"
