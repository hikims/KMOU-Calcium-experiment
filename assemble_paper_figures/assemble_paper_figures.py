"""
assemble_paper_figures.py — build the paper's main figures (Fig 1–7) and the supplementary figures (S1–…) from the
panel functions of the computation folders, without duplicating any plotting code.

How it works
  * Every plot script defines panel functions `def p_xxx(a)` (or `xxx(a)`) and lists MAIN / SUPP.  This script
    executes a folder's plot script inside that folder (so its relative paths work) with `compose` replaced by a
    no-op and __name__ != "__main__", then takes the panel functions from the resulting namespace.
  * A figure is a list of (folder, panel-function name); panels are lettered (a), (b), ... and saved as
    figures/<name>.{png,pdf,eps} with "main" or "supp" in the name.
  * Figures are rendered in parallel, one process per figure (fork pool; a folder's script is loaded at most once
    per process).  Scripts that use multiprocessing themselves run serially inside a worker.
Usage:  python3 assemble_paper_figures.py [main|supp|all] [--serial]
Expects the sibling folders under ../ (the code/ directory of the project).
"""
import os, re, sys, io, contextlib, traceback
import multiprocessing as mp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures"); os.makedirs(FIG, exist_ok=True)
SCRIPTS = {                                            # folder -> plot script
    "m0_period_waveform": "plot_m0_period_waveform.py",
    "m0ip3h_stimulus_encoding": "plot_m0ip3h_stimulus_encoding.py",
    "m0ip3h_waveform_duty": "plot_m0ip3h_waveform_duty.py",
    "m0ip3h_amplitude_observation": "plot_m0ip3h_amplitude_observation.py",
    "m0ip3h_amplitude_invariance": "plot_m0ip3h_amplitude_invariance.py",
    "m0ip3h_bifurcation_structure": "plot_m0ip3h_bifurcation_structure.py",
    "m0_joint_scan": "plot_m0_joint_scan.py",
    "m1_isi_variability": "plot_m1_isi_variability.py",
    "m1_sigma_tav_slope": "plot_m1_sigma_tav_slope.py",
    "m1_first_passage_theory": "plot_m1_first_passage_theory.py",
    "m1_small_nu_asymptotics": "plot_m1_small_nu_asymptotics.py",
    "m1_particle_filter_synthetic": "plot_m1_particle_filter_synthetic.py",
    "data_qc_splits": "qc_data_splits.py",
    "data_dose_units": "data_dose_units.py",
    "step_fold_change": "step_fold_change.py",
    "emrich_replication": "emrich_replication.py",
    "emrich_jin_prediction": "plot_emrich_jin_prediction.py",
}
# ---- paper figures: (name, ncol, [(folder, panel), ...]) --------------------------------------------------------
MAIN_FIGS = [
    ("fig1_main_data_and_M0", 2, [("data_dose_units", "p_dose_units"), ("m0_period_waveform", "p_period_vs_stimulus"), ("m0_period_waveform", "p_one_cycle"), ("m0_joint_scan", "p_ratio_vs_duty")]),
    ("fig2_main_rate_coupling_relaxation", 2, [("m0ip3h_stimulus_encoding", "p_direction"), ("m0ip3h_waveform_duty", "p_duty_vs_beta"), ("m0ip3h_waveform_duty", "p_ref_vs_stimulus"), ("m0ip3h_bifurcation_structure", "p_map_S_epsH")]),
    ("fig3_main_stereotypy_qoff", 2, [("m0ip3h_amplitude_invariance", "data_norm"), ("m0ip3h_amplitude_invariance", "model_norm"), ("m0ip3h_amplitude_invariance", "ratios_qoff"), ("m0ip3h_bifurcation_structure", "p_map_beta_qoff")]),
    ("fig4_main_first_passage_theory", 3, [("m1_first_passage_theory", "moments"), ("m1_first_passage_theory", "reduction"), ("m1_first_passage_theory", "reduction_error"), ("m1_small_nu_asymptotics", "p_mean_sd"), ("m1_small_nu_asymptotics", "p_skew")]),
    ("fig5_main_noise_model_comparison", 3, [("m1_isi_variability", "trace_panel"), ("m1_isi_variability", "isi_dist"), ("m1_isi_variability", "noise_loc")]),
    ("fig6_main_variability_law", 2, [("m1_sigma_tav_slope", "law"), ("m1_sigma_tav_slope", "slopes")]),
    ("fig7_main_azzoni_data", 2, [("data_qc_splits", "example_trace"), ("data_dose_units", "p_train_strip"), ("step_fold_change", "per_cell"), ("data_dose_units", "p_stats_units")]),
]
# Emrich (independent-source consistency check, WT only) and the J_in scenario are supplementary material
# supplementary: every folder's SUPP list (read from the script) plus MAIN lists not used above
SUPP_ORDER = ["data_dose_units", "data_qc_splits", "m0_period_waveform", "m0_joint_scan", "m0ip3h_stimulus_encoding", "m0ip3h_waveform_duty", "m0ip3h_amplitude_observation",
              "m0ip3h_amplitude_invariance", "m0ip3h_bifurcation_structure", "m1_first_passage_theory", "m1_small_nu_asymptotics", "m1_isi_variability", "m1_sigma_tav_slope",
              "m1_particle_filter_synthetic", "step_fold_change", "emrich_replication", "emrich_jin_prediction"]
EXTRA_MAIN_AS_SUPP = ["m0_joint_scan", "m0ip3h_amplitude_observation", "m1_particle_filter_synthetic", "data_qc_splits", "step_fold_change", "emrich_replication", "emrich_jin_prediction", "m0ip3h_bifurcation_structure"]

# ---- loading a folder's panel functions -------------------------------------------------------------------------
_CACHE = {}
class _SerialPool:
    def __init__(self, *a, **k): pass
    def __enter__(self): return self
    def __exit__(self, *a): return False
    def map(self, f, it): return list(map(f, it))
    def close(self): pass
    def join(self): pass
class _SerialCtx:
    def Pool(self, *a, **k): return _SerialPool()
def _patch_multiprocessing():
    mp.get_context = lambda *a, **k: _SerialCtx(); mp.Pool = lambda *a, **k: _SerialPool()

def load_folder(folder):
    if folder in _CACHE: return _CACHE[folder]
    path = os.path.join(CODE, folder); script = os.path.join(path, SCRIPTS[folder])
    if not os.path.exists(script): raise FileNotFoundError(script)
    src = open(script, encoding="utf-8").read()
    src = re.sub(r"^def compose\(.*?(?=^\S)", "def compose(*a, **k): pass\n", src, flags=re.S | re.M)   # no figure files from the script itself
    ns = {"__name__": "__assemble__", "__file__": script}
    cwd, argv = os.getcwd(), sys.argv
    os.chdir(path); sys.argv = [script]; sys.path.insert(0, path)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(src, script, "exec"), ns)
    finally:
        os.chdir(cwd); sys.argv = argv; sys.path.pop(0); plt.close("all")
    _CACHE[folder] = (path, ns); return _CACHE[folder]

def panel_list(folder, which):
    """names in the script's MAIN = [...] / SUPP = [...] line (read from the source, the script is not executed here)"""
    src = open(os.path.join(CODE, folder, SCRIPTS[folder]), encoding="utf-8").read()
    m = re.search(rf"^{which}\s*=\s*\[([^\]]*)\]", src, flags=re.M)
    return [(folder, n.strip()) for n in m.group(1).split(",") if n.strip()] if m else []

# ---- rendering ---------------------------------------------------------------------------------------------------
def letter(a, s): a.text(-0.13, 1.06, s, transform=a.transAxes, fontsize=28, fontweight="bold", va="bottom", ha="left")
def render(spec, w=8.0, h=7.5):
    name, ncol, panels = spec
    try:
        _patch_multiprocessing()
        for folder, _ in panels: load_folder(folder)           # scripts set rcParams: load them before the axes are created
        n = len(panels); nrow = -(-n // ncol)
        fig, axs = plt.subplots(nrow, ncol, figsize=(w * ncol, h * nrow), squeeze=False)
        for k, (folder, fname) in enumerate(panels):
            path, ns = load_folder(folder); a = axs[k // ncol, k % ncol]
            cwd = os.getcwd(); os.chdir(path)
            try:
                with contextlib.redirect_stdout(io.StringIO()): ns[fname](a)
            finally: os.chdir(cwd)
            letter(a, f"({chr(97 + k)})")
        for k in range(n, nrow * ncol): axs[k // ncol, k % ncol].axis("off")
        plt.tight_layout(w_pad=2.5, h_pad=2.5)
        for ext in ("png", "pdf", "eps"): fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
        plt.close(fig); return f"wrote figures/{name}.{{png,pdf,eps}}  [{', '.join(f'{f}:{p}' for f, p in panels)}]"
    except Exception:
        return f"FAILED {name}:\n" + traceback.format_exc()

def supp_specs():
    specs, k = [], 1
    for folder in SUPP_ORDER:
        try:
            for which in (["SUPP", "MAIN"] if folder in EXTRA_MAIN_AS_SUPP else ["SUPP"]):
                pl = panel_list(folder, which)
                if not pl: continue
                ncol = 3 if len(pl) in (3, 5, 6, 9) else (2 if len(pl) in (2, 4, 7, 8) else 1)
                specs.append((f"figS{k}_supp_{folder}_{which.lower()}", ncol, pl)); k += 1
        except Exception as e:
            print(f"skip {folder}: {e}")
    # MAIN panels of the remaining folders that no main figure uses (e.g. the dose maps, the M0ip3h one-cycle match)
    used = {p for _, _, pl in MAIN_FIGS + specs for p in pl}
    left = [p for folder in SUPP_ORDER if folder not in EXTRA_MAIN_AS_SUPP for p in panel_list(folder, "MAIN") if p not in used]
    if left: specs.append((f"figS{k}_supp_main_panels_not_in_fig1to7", 3 if len(left) in (3, 5, 6) else 2, left))
    return specs

if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "all"
    serial = "--serial" in sys.argv
    specs = []
    if what in ("main", "all"): specs += MAIN_FIGS
    if what in ("supp", "all"): specs += supp_specs()
    print(f"{len(specs)} figures on {1 if serial else mp.cpu_count()} processes")
    if serial: out = [render(s) for s in specs]
    else:
        with mp.get_context("fork").Pool(mp.cpu_count()) as pool: out = pool.map(render, specs, chunksize=1)
    print("\n".join(out))
    with open(os.path.join(FIG, "assemble_log.txt"), "w") as fh: fh.write("\n".join(out) + "\n")     # failures with tracebacks land here
    nfail = sum(o.startswith("FAILED") for o in out); print(f"{len(out) - nfail} figures written, {nfail} failed -> figures/assemble_log.txt")
