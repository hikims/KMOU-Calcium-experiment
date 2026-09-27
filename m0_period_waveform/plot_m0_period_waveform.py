"""
plot_m0_period_waveform.py — HEK293 data vs the deterministic 0-D AIMS model (M0).

Main figure  (figures/m0_period_waveform_main.*, 1x3)
  (a) data: median ISI by carbachol dose, IQR bars, n cells            [processed CSVs]
  (b) M0: period vs stimulus S̄ (two λ) — the original model cannot encode the dose
  (c) one normalized cycle: M0 (ε_H = 70, 20) vs a measured spike — the spike is far too broad
Supplement   (figures/m0_period_waveform_supp.*, 1x3)
  (a) M0: period vs ε_H, reference T ≈ 2.8 ε_H, Hopf boundary
  (b) M0: period vs β, β*_lin marker
  (c) duty cycle W50/T: model sweeps vs data band

Run after ./run_m0_sweeps.sh .
Data locations (relative to this folder; first match wins):
  processed CSVs : ../../*Spike Trains_data/2/*.csv
  raw TXT        : ../../*Spike Trains_data/CChCB/*.txt   or  .../1/*.txt
If the data are not found, hard-coded medians (2026-09-21) are used and the data spike is omitted.
"""
import glob, io, os
import numpy as np, pandas as pd
from scipy.signal import find_peaks
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- global style (shared by all figure scripts) ------------------------------------------
plt.rcParams.update({
    "font.size": 22, "axes.titlesize": 24, "axes.labelsize": 24,
    "xtick.labelsize": 20, "ytick.labelsize": 20, "legend.fontsize": 19,
    "axes.grid": True, "grid.color": "#c8c8c8", "grid.linestyle": "--", "grid.linewidth": 0.8,
    "lines.linewidth": 2.5, "lines.markersize": 6,
    "mathtext.fontset": "stix", "font.family": "STIXGeneral",
    "savefig.dpi": 300, "legend.framealpha": 1.0,
})
C_DATA, C_M0, C_EPS, C_BETA, C_LAM1 = "#c0392b", "#1f77b4", "#e69f00", "#2ca02c", "#ff7f0e"
TINT_DATA, TINT_LAM1 = "#f5d0cc", "#ffe4cc"      # opaque tints (EPS has no transparency)

NAME, OUT, FIG = "m0_period_waveform", "results", "figures"
os.makedirs(FIG, exist_ok=True)

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

# --- data -----------------------------------------------------------------------------------
ROOTS = ["../..", "../../data", "../data"]           # project root candidates
def find(pattern):
    for r in ROOTS:
        hits = sorted(glob.glob(os.path.join(r, "*Spike Trains_data", pattern)))
        if hits: return hits
    return []

ISI_MED = {3: 206.0, 5: 76.0, 10: 60.0, 15: 57.0, 30: 30.6}   # fallback medians (s)
ISI_Q1, ISI_Q3, N_CELLS, DUTY_BY_DOSE = {}, {}, {}, {}
csvs = find(os.path.join("2", "*.csv"))
if csvs:
    P = pd.concat([pd.read_csv(f, sep=";", decimal=",", encoding="utf-8-sig").assign(file=os.path.basename(f)) for f in csvs])
    g = P.groupby("stimulus")
    ISI_MED = {int(k): float(v) for k, v in g.ISI.median().items()}
    ISI_Q1 = {int(k): float(v) for k, v in g.ISI.quantile(0.25).items()}
    ISI_Q3 = {int(k): float(v) for k, v in g.ISI.quantile(0.75).items()}
    N_CELLS = {int(k): int(v) for k, v in g.apply(lambda d: (d.file + "_" + d.ST.astype(str)).nunique()).items()}
    DUTY_BY_DOSE = {int(k): float(v) for k, v in (g["half-width"].median() / g.ISI.median()).items()}
    print("data: medians", {k: round(v, 1) for k, v in ISI_MED.items()}, "| duty by dose", {k: round(v, 3) for k, v in DUTY_BY_DOSE.items()})
else:
    print("processed CSVs not found: using hard-coded medians")
duty_lo = min(DUTY_BY_DOSE.values()) if DUTY_BY_DOSE else 0.09
duty_hi = max(DUTY_BY_DOSE.values()) if DUTY_BY_DOSE else 0.16

def data_spike(dose_file="CChCB_2024_08_08_E1.txt", stim_t=900.0, target_isi=57.0, tol=0.25):
    """One measured spike (15 uM, 1 s sampling) from any ROI of the file: both neighbouring ISIs within
    +-tol of target_isi and a clean baseline (window minimum < 10% of the peak).  Returned normalized:
    phase = (t - t_peak)/ISI_prev, amplitude 0..1.  Picks the candidate whose ISI is closest to target."""
    hits = find(os.path.join("CChCB", dose_file)) or find(os.path.join("1", dose_file))
    if not hits: return None
    txt = open(hits[0], "rb").read().decode("utf-16")
    df = pd.read_csv(io.StringIO(txt), sep="\t").dropna(axis=1, how="all")
    t = df["Time [s]"].values; dt = np.median(np.diff(t))
    best = None
    for col in [c for c in df.columns if c.startswith("#")]:
        F = pd.to_numeric(df[col], errors="coerce").values
        F0 = np.nanmedian(F[t < stim_t]); y = (F - F0) / F0
        ys = pd.Series(y).rolling(3, center=True, min_periods=1).mean().values
        pk, _ = find_peaks(ys, prominence=0.3, distance=max(1, int(15 / dt)))
        pk = pk[t[pk] > stim_t]
        if len(pk) < 4: continue
        isis = np.diff(t[pk])
        for k in range(1, len(pk) - 1):
            i0, i1 = isis[k - 1], isis[k]
            if abs(i0 - target_isi) > tol * target_isi or abs(i1 - target_isi) > tol * target_isi: continue
            a, b, ip = pk[k - 1], pk[k + 1], pk[k]
            seg = slice(a, b + 1); base = np.nanmin(ys[seg]); amp = ys[ip] - base
            if amp <= 0 or (ys[ip] - amp * 0.0) <= 0: continue
            ph = (t[seg] - t[ip]) / i0; yy = (ys[seg] - base) / amp
            edge = yy[(np.abs(ph) > 0.35) & (np.abs(ph) <= 0.5)]
            if edge.size == 0 or edge.max() > 0.15: continue
            score = abs(i0 - target_isi) + abs(i1 - target_isi)
            if best is None or score < best[0]: best = (score, ph, yy, i0, col)
    if best is None: return None
    _, ph, yy, isi, col = best
    print(f"data spike: ROI {col}, ISI {isi:.0f} s")
    return ph, yy, isi

# --- model results -------------------------------------------------------------------------
def load(name): return pd.read_csv(os.path.join(OUT, name))
def one_cycle(name):
    tr = load(name); c = tr.C.values; t = tr.t.values
    pk, _ = find_peaks(c, prominence=0.5 * (c.max() - c.min()))
    j = len(pk) // 2; ip = pk[j]; T = t[pk[j + 1]] - t[ip]
    seg = (t >= t[ip] - 0.5 * T) & (t <= t[ip] + 0.5 * T)
    base = c[seg].min()
    return (t[seg] - t[ip]) / T, (c[seg] - base) / (c[ip] - base), T

# --- panels ---------------------------------------------------------------------------------
def p_data_isi(a):
    doses = sorted(ISI_MED); med = np.array([ISI_MED[d] for d in doses])
    if ISI_Q1:
        err = np.array([[ISI_MED[d] - ISI_Q1[d] for d in doses], [ISI_Q3[d] - ISI_MED[d] for d in doses]])
        a.errorbar(doses, med, yerr=err, fmt="s-", color=C_DATA, ms=10, capsize=6, elinewidth=2)
        ticklab = [f"{d}\n(n={N_CELLS[d]})" for d in doses]
    else:
        a.plot(doses, med, "s-", color=C_DATA, ms=10); ticklab = [str(d) for d in doses]
    a.set_xscale("log"); a.set_xticks(doses); a.set_xticklabels(ticklab, fontsize=17); a.minorticks_off()
    a.set_xlabel(r"carbachol ($\mu$M)"); a.set_ylabel("ISI (s)"); a.set_ylim(0, None)
    a.set_title(rf"Data: median ISI (IQR), ratio {med.max()/med.min():.1f}$\times$")

def p_period_vs_stimulus(a):
    for fn, lab, col in [("m0_period_vs_stimulus_lambda0.001.csv", r"$\lambda=0.001$ (Table 2)", C_M0), ("m0_period_vs_stimulus_lambda1.csv", r"$\lambda=1$", C_LAM1)]:
        d = load(fn); o = d[d.osc == 1]; a.plot(o.value, o.Tosc, "o-", color=col, label=lab)
        if (d.osc == 0).any(): a.axvspan(d.value.min(), d[d.osc == 0].value.max(), color=TINT_LAM1, lw=0, zorder=0)
    a.set_ylim(0, 250); a.set_xlabel(r"stimulus $\bar{S}$"); a.set_ylabel(r"$T_{\mathrm{osc}}$")
    a.set_title("M0: period vs stimulus"); a.legend(loc="lower right")
    a.text(0.03, 0.06, "shaded: no oscillation", transform=a.transAxes, fontsize=16, color=C_LAM1)

def p_period_vs_epsH(a):
    d = load("m0_period_vs_epsH.csv"); o = d[d.osc == 1]
    k = np.median(o.Tosc / o.value)
    a.loglog(o.value, o.Tosc, "o-", color=C_EPS, label="M0")
    a.loglog(o.value, k * o.value, "k--", lw=1.5, label=rf"$T_{{\mathrm{{osc}}}}={k:.1f}\,\varepsilon_H$")
    lowfail = d[(d.osc == 0) & (d.value < o.value.min())]
    eh_hopf = lowfail.value.max() if len(lowfail) else None
    if eh_hopf: a.axvline(eh_hopf, color="k", ls=":", lw=1.5); a.text(eh_hopf * 1.05, o.Tosc.max() * 0.7, "Hopf", fontsize=18)
    a.set_xlabel(r"$\varepsilon_H$"); a.set_ylabel(r"$T_{\mathrm{osc}}$"); a.set_title(r"M0: period vs $\varepsilon_H$"); a.legend(loc="lower right")

def p_period_vs_beta(a):
    d = load("m0_period_vs_beta.csv"); o = d[d.osc == 1]
    a.plot(o.value, o.Tosc, "o-", color=C_BETA); a.axvline(0.270, color="k", ls="--", lw=1.5)
    a.text(0.265, o.Tosc.max() * 0.97, r"$\beta^*_{\mathrm{lin}}=0.270$", fontsize=18, va="center", ha="right")
    a.set_xlabel(r"$\beta$ (SERCA)"); a.set_ylabel(r"$T_{\mathrm{osc}}$"); a.set_title(r"M0: period vs $\beta$")

def p_duty(a):
    for fn, lab, col in [("m0_period_vs_stimulus_lambda1.csv", r"$\bar{S}$ sweep", C_LAM1), ("m0_period_vs_epsH.csv", r"$\varepsilon_H$ sweep", C_EPS), ("m0_period_vs_beta.csv", r"$\beta$ sweep", C_BETA)]:
        d = load(fn); o = d[d.osc == 1]; a.plot(o.Tosc, o.duty, "o", color=col, ms=8, label=lab)
    a.axhspan(duty_lo, duty_hi, color=TINT_DATA, lw=0, zorder=0, label=f"data {duty_lo:.2f}–{duty_hi:.2f}")
    a.set_ylim(0, 0.7); a.set_xlabel(r"$T_{\mathrm{osc}}$"); a.set_ylabel(r"$W_{50}/T_{\mathrm{osc}}$")
    a.set_title("Duty cycle: M0 vs data"); a.legend(loc="center right")

def p_one_cycle(a):
    for fn, lab, col in [("m0_waveform_trace_baseline_epsH70.csv", r"M0, $\varepsilon_H=70$", C_M0), ("m0_waveform_trace_epsH20.csv", r"M0, $\varepsilon_H=20$", C_EPS)]:
        ph, y, T = one_cycle(fn); a.plot(ph, y, color=col, label=lab)
    sp = data_spike(target_isi=ISI_MED.get(15))
    if sp is not None:
        ph, y, isi = sp; m = (ph >= -0.5) & (ph <= 0.5)
        a.plot(ph[m], y[m], color=C_DATA, lw=3, label=rf"data, 15 $\mu$M (ISI {isi:.0f} s)")
    a.set_xlim(-0.5, 0.5); a.set_ylim(0, 1.3); a.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.set_xlabel(r"phase $(t-t_{\mathrm{peak}})/T$"); a.set_ylabel("normalized amplitude")
    a.set_title("One cycle: M0 vs measured spike")
    a.legend(loc="upper left", ncol=1, bbox_to_anchor=(0.5, 0.995), frameon=False, columnspacing=0.8, handlelength=1.2, fontsize=17)

# --- figures --------------------------------------------------------------------------------
MAIN = [p_data_isi, p_period_vs_stimulus, p_one_cycle]   # the data fact and why M0 fails: no dose encoding, spike too broad
SUPP = [p_period_vs_epsH, p_period_vs_beta, p_duty]      # parameter dependences of M0
compose(MAIN, 3, f"{NAME}_main")
compose(SUPP, 3, f"{NAME}_supp")
