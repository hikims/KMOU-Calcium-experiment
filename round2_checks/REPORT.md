# Round-2 review: computational checks (2026-09-23)

Scripts `check1`–`check7` in this folder; outputs in `results/`. Binaries are the existing kernels
(`m1_isi_variability/m1_isi_sde`, `emrich_jin_prediction/emrich_jin_sde`, `m0_joint_scan/m0_joint_scan`);
no stored result of the paper was overwritten.

## 1. Noise comparison with the CV matched *after* segmentation and detrending (`check1`)
Amplitude tuned so that the segmented, detrended CV equals the data median 0.196 (grid + log-interpolation),
then 3 seeds at the matched amplitude; median over segments with a bootstrap 95 % CI over segments.

| noise | σ matched | CV | skew | amp CV | ρ(ISI,amp) | ρ₁ |
|---|---|---|---|---|---|---|
| H | 0.0244 | 0.191 | 0.52 [0.48, 0.59] | 0.089 | −0.10 [−0.12, −0.07] | −0.07 |
| Ca | 0.423 | 0.209 | 1.11 [1.03, 1.20] | 0.264 | +0.36 [0.34, 0.38] | −0.01 |
| P | 0.0335 | 0.192 | 0.45 [0.40, 0.51] | 0.037 | −0.71 [−0.73, −0.67] | +0.65 [0.60, 0.68] |
| data | — | 0.196 | 0.49 | 0.135 | 0.05 | −0.03 |

Conclusion: with a CV-matched comparison the picture of Table 4 is unchanged. Ca-noise fails on skewness, amplitude CV
and ρ(ISI,amp); P-noise now matches the skewness but fails on amplitude CV, ρ(ISI,amp) and ρ₁. H-noise matches CV,
skewness and ρ₁, underestimates the amplitude CV (0.09 vs 0.135) and gives a small negative ρ(ISI,amp) where the data
give a small positive one. The paper can keep a *moderated* discrimination statement and should cite these numbers
(or replace Table 4 by them). Note σ_H matched after segmentation is 0.0244 (ν = 0.295), not 0.0232.

## 2. Time-step refinement (`check2`)
dt = 0.01, 0.005, 0.0025 at the Table-4 amplitudes; 3 seeds each (≈1400 ISIs per run).

| noise | dt | T_av | sd | CV | skew | amp CV | ρ(ISI,amp) | ρ₁ |
|---|---|---|---|---|---|---|---|---|
| H | 0.01 / 0.005 / 0.0025 | 85.5 / 85.4 / 85.5 (±0.4) | 16.5 / 16.4 / 16.5 | 0.193 / 0.192 / 0.193 | 0.78 / 0.79 / 0.77 | 0.086 / 0.086 / 0.086 | −0.06 / −0.03 / −0.06 (±0.03) | 0.03 / −0.02 / 0.01 (±0.04) |
| Ca | 0.01 / 0.005 / 0.0025 | 72.7 / 73.9 / 74.6 (±0.5) | 13.1 / 13.5 / 13.4 | 0.181 / 0.182 / 0.180 | 1.76 / 1.73 / 1.87 | 0.249 / 0.249 / 0.249 | 0.43 / 0.44 / 0.45 | 0.12 / 0.13 / 0.11 |
| P | 0.01 / 0.005 / 0.0025 | 84.9 / 84.9 / 85.7 (±2) | 17.2 / 16.2 / 16.3 | 0.203 / 0.190 / 0.190 | 0.42 / 0.41 / 0.34 | 0.033 / 0.030 / 0.030 | −0.94 / −0.94 / −0.94 | 0.93 / 0.92 / 0.92 |

Conclusion: the H-noise model (the paper's model) is converged at dt = 0.01 within the seed spread. Ca-noise shows a
≈1.3 % drift of the mean interval per halving of dt and P-noise a ≈6 % change of the sd between 0.01 and 0.005; the
statistics used to reject the two controls (amplitude CV, ρ(ISI,amp), ρ₁) are insensitive to dt.

## 3. Basal-influx scenario: direct population runs (`check3`)
400 cells, 2 seeds, same protocol as S9.3.

| J_in factor | median count [IQR] | mean | silent |
|---|---|---|---|
| 0.10 | 4 [3, 6] | 4.5–4.7 | 0.3–0.5 % |
| 0.075 | 4 [3, 6] | 4.5–4.6 | 0.3–0.5 % |
| 0.05 | 4 [3, 5–6] | 4.4–4.5 | 0–0.3 % |
| 0.03 | 4 [3, 5] | 4.3 | 0–0.5 % |
| 0.02 | 4 [3, 5] | 4.2–4.3 | 0.3–0.8 % |
| 0.01 | 4 [3, 5] | 4.1–4.2 | 0.3–0.5 % |
| 0 (stored) | 0 | 0 | 100 % |

Conclusion: the interpolated "TKO count of 3 at factor ≈ 0.075" is **not** reproduced. The median count stays at 4
for every partial reduction down to 1 % of the reference influx and drops to 0 only at J_in = 0; the mean falls
smoothly from 4.5 to 4.1. A reduction of the basal influx alone therefore does not reproduce the Orai-TKO median in
this model. Stationary ratios at factor 0.075 (direct run): period 1.41, amplitude 1.24 (1.27 mean), W50 1.11, CV 1.10.
S9.3 and Fig. S25 must be rewritten as a negative/partial result.

## 4. M0 joint scan: the C2 point is a damped transient (`check4`, `check4b`)
The single point flagged C2 (idx 601: β = 0.20, ε_H = 76, T ≈ 850–900) has a C amplitude of 1.96 → 0.055 → 0.0014 →
0.00004 over the four quarters of the recording window at every "oscillating" level. It is a slowly damped transient
that the amplitude criterion (window max − min ≥ 0.05·max) does not catch. Consequently no sampled point satisfies C2;
Table S2 rows "C2" and "C2 & C3" become 0 and Fig. S5(d) must be replaced or dropped. `check4b` screens all 232
oscillating points the same way (see `results/check4b_transient_screen.log`).

## 5. Noise-induced mean shift: sign differs between reduced process and full model (`check5`, `check5b`)
Reduced theory (fixed reset/threshold): mean recovery time decreases, −1.7 % at ν = 0.36. Full model: mean ISI
*increases*, +3.3 % at ν = 0.36 (+1.3 % at ν = 0.25). The increase is quantitatively explained by the upward drift of
the H value at spike onset (H_on 0.4085 → 0.4204 at σ_H = 0.03): ε_eff Δln[(1−H₀)/(1−H_on)] = 2.98 vs observed 2.79
model units at ν = 0.36 (agreement within 10 % for ν ≤ 0.4; at larger ν premature spikes counteract). The sentence
"the mean interval is shifted below the deterministic value by Var/2" applies to the reduced process only; the
Discussion sentence on "noise-induced early crossings" must be restricted accordingly.

## 6. Bootstrap (`check6`)
- Ratio of train-median ISIs 3/30 µM: 6.95; 95 % CI [4.75, 10.4] (train resampling, n = 16/32),
  [4.72, 9.72] (record resampling, n = 5/5). Replace "CI-end ratio 4.6" by "lower CI bound ≈ 4.7 (ratio bootstrap)".
- Slope log σ vs log T_av (196 detrended segments): 1.33; 95 % CI [1.23, 1.42] (segments), [1.22, 1.43] (146 ROIs),
  [1.18, 1.44] (23 records). Clustering barely widens the interval.

## 7. IP₃ gate saturation (`check7`)
1 − P̄/(1+P̄) at the mapped doses: 3 µM 0.36–0.40 %, 5 µM 0.23–0.25 %, 10 µM 0.13 %, 15 µM 0.10 %, 30 µM 0.06 %.
"Within 10⁻³ of one over the whole dose range" holds only for S̄ ≳ 0.3; replace by the range above.

## 4b. Transient screen of all oscillating scan points (`check4b`)
231 of the 232 oscillating points traced at the middle of their monotone window (one without a usable trace).
8 points have a last-quarter amplitude below 50 % of the first quarter (4 below 10 %): damped transients. None of them
passes C1 (max period ratio 1.07). Among the 223 sustained points: C2 = 0, C1(CI) = 1, C3 = 170.

## Manuscript revision (2026-09-23)
Applied to `manuscript/ecrc/main.tex`, `supplement.tex`, `references.bib` (git commit "Round-2 revision"). New
supplementary Tables S7 (protocol windows), S8 (CV-matched noise comparison), S9 (time-step refinement),
S10 (direct influx runs) in a final section S12; S2.2, S5.3, S7, S8, S9.2, S9.3 rewritten as described above.

## Round 3 (2026-09-23): additional checks
- `check8`: J_in = 0 silence is a cold-start artefact (C decays to 1e-6 in the rest phase; C = 0 absorbing for n_A = 2).
  Warm start (S_E on during the pre-phase): J_in = 0 -> median 2 [0,3], 32 % silent; factors 0.01 / 0.075 -> 5 / 6.
- `check9`: dt refinement at the re-tuned amplitudes (Table S8): H stable (Tav 85.7/85.6/85.7, CV 0.202/0.203/0.202);
  Ca mean +3.3 % over two halvings; P CV -6 % at the first halving; correlations unchanged.
- `check10`: 24 of 5888 union-train intervals (0.4 %) span a source-segment boundary. Excluding the three
  protocol-discrepant windows (16 trains/segments): 3 uM median 226 -> 206 s, ratio 6.95 -> 6.30, segment CV 0.196 -> 0.207,
  skew 0.49 -> 0.51, slope 1.33 -> 1.36.
- Manuscript: main abstract/intro/2.3/4.1/4.2/5.3/6/7/Discussion wording as in the round-3 report; supplement S1, S2.2,
  S5.3, S6, S7, S9 rewritten; new Tables S11 (warm start), S12 (dt at calibration B), S13 (window sensitivity);
  Section S12 renamed "Additional validation and sensitivity checks".

## Round 4 (2026-09-23)
- `check11`: the 24 cross-segment intervals all lie in 15 uM trains (20 trains); removing them changes only the 15 uM
  train median (93.5 -> 92.0 s). 40 of 196 source segments contain a raw acquisition gap > 10 s (max 54 s).
- `check12`: 37 of 5864 segment intervals (0.63 %, 36 segments) span a gap; removing them: CV 0.196 -> 0.193,
  skew 0.493 -> 0.487, serial -0.029 -> -0.024. ROI counts: 146 ROI units in the tables, 141 matched to a raw column
  (5 table ROIs of 3 records have no raw column).
- Manuscript: warm-start "success interval" claim removed; C_c = mu k_A^2/alpha = 8e-6 sufficient condition added to
  S9.3; abstract R limit restricted to nu -> 0; 35/32 windows in main; wording per the round-4 list.

## Round 5 (2026-09-23)
- REF_PREV in the code = COMMON + qoff=1, rC=30 with kR at its kernel default kH = 0.15 (REF_ST uses kR = 1); Table 2 and
  Fig. 3 captions now state this. Duty at S = 1: REF_PREV 0.13 (Fig. 3c), REF_ST 0.069 (0.024 at 0.1, 0.10 at 30).
  The observation-model folder (Figs S9-S10) also uses REF_PREV with its own LS map T0 = 0.72 s; stated in S3.3.
- check12: serial correlation recomputed without joining the neighbours of a removed gap-spanning interval into a new
  lag-one pair; value unchanged at three decimals (-0.024).
- Rounding unified to one decimal for the dose medians (226.5 / 92.5 / 128.0 / 93.5 / 32.6 s) in main text and Table S13.

## Round 6 (2026-09-23)
- Exact S = 1 duties with the author's kernel (dt 0.02, transient 9000): REF_ST T = 83.53, W50 = 5.64, duty 0.0675;
  REF_PREV (kR = 0.15) T = 60.11, W50 = 7.92, duty 0.1318. Both agree with the reviewer's independent ODE values.
- M0 scan window rule written into the main text (adjacent grid levels, all oscillating, non-increasing T with 1e-4
  tolerance; 51 points without a qualifying pair screened at the central level S ~ 0.95).
- Symbols: time unit -> t_unit, moment functions -> M_n (M_0 = 1), T_0 kept for the deterministic recovery time.
- "Hopf curve" -> "candidate Hopf locus Re lambda = 0 of the linearised rest state" throughout.

## Round 7 (2026-09-23, pasted review)
- Window rule unified: both documents now say "qualifying contiguous run with the largest end-point period ratio"
  (what m0_joint_scan.c does); "of these" sentence split so that the C1 count refers to the 181 windowed sets.
- Table S8 rows regenerated at three decimals from check1 results; caption restored (medians, 95 % segment bootstrap, 2000 resamples).
- Fig. 6(a) caption: H_on shown as deterministic reference, not a firing rule. Legends of Figs S10(b), S26(a,c): T_0 -> t_unit
  (plot scripts of data_dose_units, m0ip3h_stimulus_encoding, m0ip3h_amplitude_observation edited; figures regenerated).
- Discussion \paragraph headings without trailing periods; Proposition 7.1 kept on one page (needspace).
