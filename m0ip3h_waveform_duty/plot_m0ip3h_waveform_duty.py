"""
plot_m0ip3h_waveform_duty.py — what makes the Ca spike narrow (duty W50/T ~ 0.1) while keeping stimulus
encoding?  Model = AIMS reaction terms + IP3->H rate coupling + {cooperative activation (nA,kA), basal influx Jin,
Ca-dependent H speed-up rC}.

Main figure  (figures/m0ip3h_waveform_duty_main.*, 1x3)
  (a) duty vs SERCA rate β in the cooperative/relaxation regime (the decisive knob) + data band
  (b) reference regime: period and duty vs stimulus (encoding preserved with a narrow spike)
  (c) one normalized cycle of the reference regime vs a measured spike
Supplement   (figures/m0ip3h_waveform_duty_supp.*, 1x3)
  (a) ablation: duty of the reference regime when each ingredient is switched back to the original
  (b) random search: period ratio vs duty (both stimuli oscillating); data targets; feasible set
  (c) reference regime: C(t) at low and high stimulus (a few cycles)

Run after ./run_m0ip3h_waveform_sweeps.sh and python3 search_m0ip3h_waveform_duty.py
"""
import glob, io, os, subprocess
import numpy as np, pandas as pd
from scipy.signal import find_peaks
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 22, "axes.titlesize": 24, "axes.labelsize": 24,
    "xtick.labelsize": 20, "ytick.labelsize": 20, "legend.fontsize": 19,
    "axes.grid": True, "grid.color": "#c8c8c8", "grid.linestyle": "--", "grid.linewidth": 0.8,
    "lines.linewidth": 2.5, "lines.markersize": 6,
    "mathtext.fontset": "stix", "font.family": "STIXGeneral",
    "savefig.dpi": 300, "legend.framealpha": 1.0,
})
C_DATA, C_LO, C_HI, C_REF, C_GREY = "#c0392b", "#1f77b4", "#e69f00", "#2ca02c", "#888888"
TINT_DATA = "#f5d0cc"
NAME, OUT, FIG, BIN = "m0ip3h_waveform_duty", "results", "figures", "./m0ip3h_waveform_sweep"
# the panels call the compiled kernel; build it if it is missing (Makefile if present, else a direct cc call)
def _build(binary):
    src = os.path.basename(binary) + ".c"
    if os.path.exists("Makefile"): subprocess.run(["make", "-s"], check=True); return
    for flags in (["-O2", "-fopenmp"], ["-O2", "-Xpreprocessor", "-fopenmp", "-lomp"], ["-O2"]):
        if subprocess.run(["cc"] + flags + ["-o", binary, src, "-lm"], capture_output=True).returncode == 0: return
    raise SystemExit(f"could not compile {src} (need a C compiler; on macOS: xcode-select --install, brew install libomp)")
if not os.path.exists(BIN): _build(BIN)

os.makedirs(FIG, exist_ok=True)

REF = dict(mech=2, rho=0.024, gH=30, dH=10000, epsH=1000, nA=2, kA=0.4, Jin=0.02, kS=0.25, rC=30, beta=4, alpha=20)
NUM = dict(Tend=15000, Ttrans=9000, dt=0.02)

# --- figure assembly (MAIN / SUPP lists at the end of the file) ---------------------------------------------------------
def letter(a, s): a.text(-0.13, 1.06, s, transform=a.transAxes, fontsize=28, fontweight="bold", va="bottom", ha="left")
def compose(panels, ncol, name, w=8.0, h=7.5):
    """panels: list of functions f(ax); drawn in order and lettered (a), (b), ...; saved as figures/name.{png,pdf,eps}"""
    n = len(panels); nrow = -(-n // ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(w * ncol, h * nrow), squeeze=False)
    for k, fn in enumerate(panels):
        a = axs[k // ncol, k % ncol]; fn(a); letter(a, f"({chr(97 + k)})")
    for k in range(n, nrow * ncol): axs[k // ncol, k % ncol].axis("off")
    plt.tight_layout(w_pad=2.5, h_pad=2.5)
    for ext in ("png", "pdf", "eps"): fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig); print(f"wrote figures/{name}.{{png,pdf,eps}}")

# --- data ----------------------------------------------------------------------------------------
ROOTS = ["../..", "../../data", "../data"]
def find(pattern):
    for r in ROOTS:
        hits = sorted(glob.glob(os.path.join(r, "*Spike Trains_data", pattern)))
        if hits: return hits
    return []
ISI_MED = {3: 206.0, 5: 76.0, 10: 60.0, 15: 57.0, 30: 30.6}; DUTY_LO, DUTY_HI = 0.05, 0.15
csvs = find(os.path.join("2", "*.csv"))
if csvs:
    P = pd.concat([pd.read_csv(f, sep=";", decimal=",", encoding="utf-8-sig") for f in csvs]); g = P.groupby("stimulus")
    ISI_MED = {int(k): float(v) for k, v in g.ISI.median().items()}
    dd = (g["half-width"].median() / g.ISI.median()); DUTY_LO, DUTY_HI = float(dd.min()), float(dd.max())
DATA_RATIO = max(ISI_MED.values()) / min(ISI_MED.values())

def data_spike(dose_file="CChCB_2024_08_08_E1.txt", stim_t=900.0, target_isi=57.0, tol=0.25):
    hits = find(os.path.join("CChCB", dose_file)) or find(os.path.join("1", dose_file))
    if not hits: return None
    df = pd.read_csv(io.StringIO(open(hits[0], "rb").read().decode("utf-16")), sep="\t").dropna(axis=1, how="all")
    t = df["Time [s]"].values; dt = np.median(np.diff(t)); best = None
    for col in [c for c in df.columns if c.startswith("#")]:
        F = pd.to_numeric(df[col], errors="coerce").values; F0 = np.nanmedian(F[t < stim_t]); y = (F - F0) / F0
        ys = pd.Series(y).rolling(3, center=True, min_periods=1).mean().values
        pk, _ = find_peaks(ys, prominence=0.3, distance=max(1, int(15 / dt))); pk = pk[t[pk] > stim_t]
        if len(pk) < 4: continue
        isis = np.diff(t[pk])
        for k in range(1, len(pk) - 1):
            i0, i1 = isis[k - 1], isis[k]
            if abs(i0 - target_isi) > tol * target_isi or abs(i1 - target_isi) > tol * target_isi: continue
            seg = slice(pk[k - 1], pk[k + 1] + 1); base = np.nanmin(ys[seg]); amp = ys[pk[k]] - base
            if amp <= 0: continue
            ph = (t[seg] - t[pk[k]]) / i0; yy = (ys[seg] - base) / amp
            edge = yy[(np.abs(ph) > 0.35) & (np.abs(ph) <= 0.5)]
            if edge.size == 0 or edge.max() > 0.15: continue
            sc = abs(i0 - target_isi) + abs(i1 - target_isi)
            if best is None or sc < best[0]: best = (sc, ph, yy, i0)
    return None if best is None else best[1:]

# --- helpers ---------------------------------------------------------------------------------------
def load(name): return pd.read_csv(os.path.join(OUT, name))
def osc(d): return d[d.osc == 1]
def run_one(cfg):
    args = [BIN] + [f"{k}={v}" for k, v in {**NUM, **cfg}.items()]
    return pd.read_csv(io.StringIO(subprocess.run(args, capture_output=True, text=True, check=True).stdout)).iloc[0]
def one_cycle_from_trace(tr):
    c, t = tr.C.values, tr.t.values
    pk, _ = find_peaks(c, prominence=0.5 * (c.max() - c.min()))
    if len(pk) < 3: return None
    j = len(pk) // 2; ip = pk[j]; T = t[pk[j + 1]] - t[ip]
    seg = (t >= t[ip] - 0.5 * T) & (t <= t[ip] + 0.5 * T); base = c[seg].min()
    return (t[seg] - t[ip]) / T, (c[seg] - base) / (c[ip] - base), T

# ablation runs (used by the supplement panel and the printed summary)
variants = [("reference", {}), (r"$n_A=1$ (orig. activation)", dict(nA=1, kA=1)), (r"$r_C=0$ (orig. $H$ kinetics)", dict(rC=0)),
            (r"$\beta=0.26$ (Table 1)", dict(beta=0.26)), (r"$\alpha=3.5$ (Table 1)", dict(alpha=3.5)), (r"$J_{\mathrm{in}}=0$", dict(Jin=0))]
abl_labels, abl_vals, abl_cols = [], [], []
for lab, mod in variants:
    f = run_one({**REF, **mod, "S": 1.0}); abl_labels.append(lab)
    abl_vals.append(f.duty if f.osc == 1 else 0.0); abl_cols.append(C_REF if lab == "reference" else (C_GREY if f.osc == 1 else "white"))
# random search
RS = pd.read_csv(os.path.join(OUT, "m0ip3h_waveform_duty_search.csv")); ok = RS[RS.both_osc]
good = ok[(ok.duty_max <= 0.15) & (ok.T_ratio >= 4)]

# --- panels ---------------------------------------------------------------------------------
def p_duty_vs_beta(a):
    d = load("m0ip3h_duty_vs_beta.csv"); o = osc(d)
    a.axhspan(DUTY_LO, DUTY_HI, color=TINT_DATA, lw=0, zorder=0, label=f"data {DUTY_LO:.2f}–{DUTY_HI:.2f}")
    a.semilogx(o.value, o.duty, "o-", color=C_REF, label=r"cooperative regime ($n_A=2$, $r_C=30$)")
    lost = d[d.osc == 0]
    if len(lost): a.semilogx(lost.value, np.zeros(len(lost)), "x", color=C_REF, ms=10, mew=2)
    a.axvline(0.26, color=C_GREY, ls="--", lw=1.5); a.text(0.28, 0.62, "Table 1\n" + r"$\beta=0.26$", fontsize=16, color=C_GREY, va="top")
    a.set_ylim(0, 0.8); a.set_xlabel(r"$\beta$ (SERCA)"); a.set_ylabel(r"$W_{50}/T_{\mathrm{osc}}$")
    a.set_title("Pump strength sets the duty cycle"); a.legend(loc="upper right", fontsize=17)

def p_ablation(a):
    a.barh(range(len(abl_labels)), abl_vals, color=abl_cols, edgecolor="k")
    for i, v in enumerate(abl_vals):
        if v == 0: a.text(0.01, i, "no oscillation", va="center", fontsize=16)
    a.axvspan(DUTY_LO, DUTY_HI, color=TINT_DATA, lw=0, zorder=0)
    a.set_yticks(range(len(abl_labels))); a.set_yticklabels(abl_labels, fontsize=17); a.invert_yaxis()
    a.set_xlim(0, 0.7); a.set_xlabel(r"$W_{50}/T_{\mathrm{osc}}$ at $\bar S=1$"); a.set_title("Ablation of the reference regime")

def p_search(a):
    a.axhspan(DUTY_LO, DUTY_HI, color=TINT_DATA, lw=0, zorder=0, label="data duty")
    a.axvline(DATA_RATIO, color=C_DATA, ls="--", lw=2, label=rf"data ratio {DATA_RATIO:.1f}$\times$")
    a.plot(ok.T_ratio, ok.duty_max, "o", color=C_GREY, ms=7, label=f"both oscillate (n={len(ok)})")
    a.plot(good.T_ratio, good.duty_max, "o", color=C_REF, ms=9, label=f"duty≤0.15 & ratio≥4 (n={len(good)})")
    a.set_ylim(0, 0.7); a.set_xlim(1, None); a.set_xlabel(r"period ratio $T(\bar S=0.1)/T(\bar S=30)$"); a.set_ylabel(r"max $W_{50}/T_{\mathrm{osc}}$")
    a.set_title(f"Random search ({len(RS)} samples)"); a.legend(loc="upper right", fontsize=16)

def p_ref_vs_stimulus(a):
    d = osc(load("m0ip3h_ref_period_vs_stimulus.csv"))
    a.semilogx(d.value, d.Tosc, "o-", color=C_REF, label=r"$T_{\mathrm{osc}}$")
    a.set_ylim(0, None); a.set_xlabel(r"stimulus $\bar S$"); a.set_ylabel(r"$T_{\mathrm{osc}}$", color=C_REF)
    a2 = a.twinx(); a2.semilogx(d.value, d.duty, "s--", color=C_DATA, label=r"$W_{50}/T_{\mathrm{osc}}$"); a2.set_ylim(0, 0.5)
    a2.set_ylabel(r"$W_{50}/T_{\mathrm{osc}}$", color=C_DATA); a2.grid(False)
    a.set_title(rf"Reference regime:" + "\n" + rf"ratio {d.Tosc.max()/d.Tosc.min():.1f}$\times$, duty {d.duty.mean():.2f}")
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels(); a.legend(h1 + h2, l1 + l2, loc="center right")

def p_traces(a):
    for fn, lab, col in [("m0ip3h_ref_trace_S0.1.csv", r"$\bar S=0.1$", C_LO), ("m0ip3h_ref_trace_S30.csv", r"$\bar S=30$", C_HI)]:
        tr = load(fn); t = tr.t - tr.t.iloc[0]; m = t <= 600
        a.plot(t[m], tr.C[m], color=col, label=lab)
    a.set_xlabel(r"$t$ (model units)"); a.set_ylabel(r"$C$"); a.set_title("Reference regime:\n spike trains"); a.legend(loc="upper right")
    a.set_ylim(0, 40); a.set_yticks([0, 5, 10, 15, 20, 25, 30])

def p_one_cycle(a):
    for fn, lab, col in [("m0ip3h_ref_trace_S0.1.csv", r"model, $\bar S=0.1$", C_LO), ("m0ip3h_ref_trace_S30.csv", r"model, $\bar S=30$", C_HI)]:
        r = one_cycle_from_trace(load(fn))
        if r: a.plot(r[0], r[1], color=col, label=lab)
    sp = data_spike(target_isi=ISI_MED.get(15, 57.0))
    if sp is not None:
        ph, y, isi15 = sp; m = (ph >= -0.5) & (ph <= 0.5); a.plot(ph[m], y[m], color=C_DATA, lw=3, label=r"data, 15 $\mu$M")
    a.set_xlim(-0.5, 0.5); a.set_ylim(0, 1.2); a.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.set_xlabel(r"phase $(t-t_{\mathrm{peak}})/T$"); a.set_ylabel("normalized amplitude"); a.set_title("One cycle:\n reference regime vs measured spike")
    a.legend(loc="upper left", ncol=1, bbox_to_anchor=(0.6, 0.995), frameon=False, columnspacing=0.8, handlelength=1.2, fontsize=17)

# --- figures --------------------------------------------------------------------------------
MAIN = [p_duty_vs_beta, p_ref_vs_stimulus, p_one_cycle]  # the pump sets the duty; encoding kept; spike shape matches
SUPP = [p_ablation, p_search, p_traces]                  # ablation, random search, raw traces
compose(MAIN, 3, f"{NAME}_main")
compose(SUPP, 3, f"{NAME}_supp")
print("ablation duty:", dict(zip(abl_labels, np.round(abl_vals, 3))))
print(f"search: {len(ok)} both-oscillating, {len(good)} feasible; best ratio among feasible = {good.T_ratio.max():.2f}" if len(good) else "search: no feasible sample")
