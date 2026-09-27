# First-passage analysis of stochastic recovery in a calcium relaxation oscillator

Code and stored results for

> M. Kim, *First-passage analysis of stochastic recovery in a calcium relaxation oscillator calibrated on public HEK293 spike trains*, submitted to Communications in Nonlinear Science and Numerical Simulation (2026).

Every figure of the main text (Figs. 1–7) and of the supplement (Figs. S1–S26, Tables S1–S13) is produced by the scripts in this repository. Figures are not stored here; the `results/` folders hold the CSV output that the figures and the supplementary tables are built from, so that the reported numbers can be checked without re-running the long sweeps. The C kernels integrate the deterministic and stochastic models (OpenMP), the Python scripts run the sweeps, analyse the public data and draw the figures.

## Layout

Each computation lives in one folder with its own `README.md`, a C kernel (`*.c` + `Makefile`) where needed, run scripts (`run_*.sh`, `*.py`), a `results/` folder with the stored CSV output and a `plot_*.py` that draws the panels. `assemble_paper_figures/` composes the paper figures from those panels.

| folder | content | outputs |
|---|---|---|
| `data_dose_units` | union of the authors' spike tables, train and segment statistics, detrending, dose tables, end-point dose map (the reference tables read by the other folders) | Figs. 1b, 1d, 2a, S1, S26 |
| `data_qc_splits` | raw-trace peak detection, cross-check against the tables, ROI quality control | Figs. 1a, S2, S3 |
| `step_fold_change` | within-cell step ratios | Figs. 1c, S20, S21 |
| `m0_period_waveform` | period and waveform of the base model M0 | Figs. 2b, 2c, S4 |
| `m0_joint_scan` | Latin-hypercube scan of the ten M0 parameters | Figs. 2d, S5, S6 |
| `m0ip3h_stimulus_encoding` | IP3–receptor coupling and stimulus encoding | Figs. 3a, S7 |
| `m0ip3h_waveform_duty` | spike shape and duty cycle | Figs. 3b, 3c, S8 |
| `m0ip3h_amplitude_observation` | fluorescence observation model | Figs. S9, S10 |
| `m0ip3h_amplitude_invariance` | stereotypy of the spike, q_off | Figs. 4a–c, S11 |
| `m0ip3h_bifurcation_structure` | rest-state stability, Re λ = 0 locus, oscillation region | Figs. 3d, 4d, S12, S13 |
| `m1_first_passage_theory` | reduced recovery process, quadrature of the moments, full-model check | Figs. 5a–c, S14 |
| `m1_small_nu_asymptotics` | symbolic small-noise expansion, convergence, densities | Figs. 5d, 5e, S15 |
| `m1_isi_variability` | noise location (H, IP3, Ca) with the same-pipeline comparison | Figs. 6, S16 |
| `m1_sigma_tav_slope` | variability law and its slopes | Figs. 7, S17 |
| `m1_particle_filter_synthetic` | conditional likelihood sensitivity on a synthetic record | Figs. S18, S19 |
| `emrich_replication` | re-analysis of the Emrich et al. source data | Figs. S22, S23 |
| `emrich_jin_prediction` | basal-influx scenario | Figs. S24, S25 |
| `round2_checks` | review-stage checks: CV-matched noise comparison, time-step refinement, direct influx runs, transient screen, bootstrap, window and gap sensitivity (`REPORT.md` summarises them) | Tables S7–S13 |
| `assemble_paper_figures` | figure assembly (writes `figures/`, not tracked) | Figs. 1–7, S1–S26 |

Model definitions, parameter sets (`REF_ST`, `REF_PREV`) and the analysis conventions are those of the paper; the header of each C file states the equations it integrates.

## Requirements

- A C compiler. OpenMP is used when available (macOS: `xcode-select --install` and `brew install libomp`; Linux: gcc). Without OpenMP the kernels build serially.
- Python 3.10 or later with `numpy`, `pandas`, `scipy`, `matplotlib`, `sympy`, `openpyxl` (`pip install -r requirements.txt`). The stored results were produced with numpy 2.3, pandas 2.2, scipy 1.17, matplotlib 3.10, sympy 1.14.

## Data

The public data are not included. Download the two Zenodo records and the Emrich source data and place them next to the repository (or in a sibling folder named `data/`):

```
<project>/
  calcium-fpt-code/                                                  (this repository)
  Fluorescence Intensity Data from HEK293 Ca²⁺ Spike Trains_data/   one folder holding both Zenodo records:
      CChCB/*.txt   raw fluorescence traces      (Zenodo 10.5281/zenodo.15705673)
      2/*.csv       processed spike tables        (Zenodo 10.5281/zenodo.14779474)
  Omnitemporal choreographies ... Ca2+ signaling_data/               Emrich et al. 2021, Source_Data_Emrich.xlsx
```

The scripts locate the folders by glob (`../../*Spike Trains_data/`, `../../*signaling_data/`), so the folder names may keep the Zenodo titles; `round2_checks/check11_gaps_rois.py` and `check12_gap_intervals.py` use the folder name shown above. A sibling `data/` folder holding the same subfolders also works. See `data_dose_units/README.md` for the exact patterns.

## Running

From the repository root:

```bash
./RUN_ALL.sh      # regenerate all figures from the stored results/ (~20 min on 8 cores; compiles the kernels first)
./RUN_FULL.sh     # recompute everything from the raw data (~2 h on 8 cores)
```

`RUN_ALL.sh` re-runs only the light steps and reuses the stored scans (joint scan, bifurcation grids, asymptotics, particle filter). Each folder can also be run on its own (`make && ./run_*.sh && python3 plot_*.py`; see its `README.md`). The assembled figures appear in `assemble_paper_figures/figures/`; failures, if any, are listed in `assemble_paper_figures/figures/assemble_log.txt`.

The review-stage checks are run with `cd round2_checks && ./run_all.sh`; they use the compiled kernels of the other folders and write to `round2_checks/results/`.

## Citation

If you use this code, please cite the paper above and this repository (see `CITATION.cff`).

## License

MIT License; see `LICENSE`. The data sets have their own licenses (the Zenodo records are CC-BY-4.0).
