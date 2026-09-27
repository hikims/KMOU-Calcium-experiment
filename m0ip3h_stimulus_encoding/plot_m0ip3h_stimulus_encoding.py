"""
plot_m0ip3h_stimulus_encoding.py — can an IP3 -> H coupling make the stimulus encode the period?

Main figure  (figures/m0ip3h_stimulus_encoding_main.*, 1x2)
  (a) direction: T(S̄) for original (mech 0), threshold coupling (mech 1), rate coupling (mech 2)
  (b) best config mapped to data: S̄ = c·U/(K+U), t = T0·s fitted to the 5 dose medians
Supplement   (figures/m0ip3h_stimulus_encoding_supp.*, 2x2)
  (a) rate coupling: Table-1 ρ / reduced ρ / reduced ρ + larger ε_H (oscillation kept over the range)
  (b) grid search: achievable period ratio vs mean duty for all full-range, decreasing configs
  (c) gain sweep: T vs g_H at low and high stimulus
  (d) one cycle of the best config at low/high S̄ vs a measured spike (waveform still broad)

Run after ./run_m0ip3h_sweeps.sh and python3 search_m0ip3h_dose_range.py
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
C_DATA, C_M0, C_M1, C_M2, C_LO, C_HI = "#c0392b", "#1f77b4", "#2ca02c", "#e69f00", "#1f77b4", "#e69f00"
TINT_DATA = "#f5d0cc"
NAME, OUT, FIG, BIN = "m0ip3h_stimulus_encoding", "results", "figures", "./m0ip3h_stimulus_sweep"
# the panels call the compiled kernel; build it if it is missing (Makefile if present, else a direct cc call)
def _build(binary):
    src = os.path.basename(binary) + ".c"
    if os.path.exists("Makefile"): subprocess.run(["make", "-s"], check=True); return
    for flags in (["-O2", "-fopenmp"], ["-O2", "-Xpreprocessor", "-fopenmp", "-lomp"], ["-O2"]):
        if subprocess.run(["cc"] + flags + ["-o", binary, src, "-lm"], capture_output=True).returncode == 0: return
    raise SystemExit(f"could not compile {src} (need a C compiler; on macOS: xcode-select --install, brew install libomp)")
if not os.path.exists(BIN): _build(BIN)

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

# --- data (same loaders as m0_period_waveform) ------------------------------------------------
ROOTS = ["../..", "../../data", "../data"]
def find(pattern):
    for r in ROOTS:
        hits = sorted(glob.glob(os.path.join(r, "*Spike Trains_data", pattern)))
        if hits: return hits
    return []
ISI_MED = {3: 226.5, 5: 92.5, 10: 128.0, 15: 93.5, 30: 32.6}; ISI_Q1 = ISI_Q3 = None       # train-level medians (2026-09-22)
DUTY_LO, DUTY_HI = 0.05, 0.15
FIT_DOSES = [3, 5, 15, 30]                     # constant-dose levels (10 µM exists only inside step records)
dt_ = "../data_dose_units/results/dose_table.csv"
if os.path.exists(dt_):
    Dt = pd.read_csv(dt_)
    A_ = Dt[Dt.subset == "all"].set_index("dose"); C_ = Dt[Dt.subset == "constant+pre"].set_index("dose")
    ISI_MED = {int(d): float(C_.isi_train_median[d]) if d in C_.index else float(A_.isi_train_median[d]) for d in A_.index}
    ISI_Q1 = {int(d): float(C_.isi_train_ci_lo[d]) if d in C_.index else float(A_.isi_train_ci_lo[d]) for d in A_.index}
    ISI_Q3 = {int(d): float(C_.isi_train_ci_hi[d]) if d in C_.index else float(A_.isi_train_ci_hi[d]) for d in A_.index}
    dd = A_.W_train_median / A_.isi_train_median; DUTY_LO, DUTY_HI = float(dd.min()), float(dd.max())
else:
    print("data_dose_units results not found: using the hard-coded train-level medians")

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

# --- helpers ----------------------------------------------------------------------------------
def load(name): return pd.read_csv(os.path.join(OUT, name))
def osc(d): return d[d.osc == 1]
def trace(cfg, S):
    """run the binary for a trace of the last cycles at stimulus S; return one normalized cycle"""
    args = [BIN, "trace=1", "Tend=6000", "Ttrans=4000", f"S={S}"] + [f"{k}={v}" for k, v in cfg.items()]
    tr = pd.read_csv(io.StringIO(subprocess.run(args, capture_output=True, text=True, check=True).stdout))
    c, t = tr.C.values, tr.t.values
    pk, _ = find_peaks(c, prominence=0.5 * (c.max() - c.min()))
    if len(pk) < 3: return None
    j = len(pk) // 2; ip = pk[j]; T = t[pk[j + 1]] - t[ip]
    seg = (t >= t[ip] - 0.5 * T) & (t <= t[ip] + 0.5 * T); base = c[seg].min()
    return (t[seg] - t[ip]) / T, (c[seg] - base) / (c[ip] - base), T

R = pd.read_csv(os.path.join(OUT, "m0ip3h_dose_range_search.csv"), index_col="cfg")
ok = R[(R.n_osc == R.n_S) & (R.decreasing)].sort_values("T_ratio", ascending=False)
best_cfg = ok.iloc[0]
BEST = dict(mech=2, epsH=best_cfg.epsH, gH=best_cfg.gH, dH=best_cfg.dH, rho=best_cfg.rho)
curves = pd.read_csv(os.path.join(OUT, "m0ip3h_dose_range_best_curves.csv"))
bc = osc(curves[curves.cfg == ok.index[0]])

# dose mapping of the best config: S = c U/(K+U), t = T0 s  (log least squares to the 5 medians)
doses = np.array(sorted(ISI_MED)); isi_all = np.array([ISI_MED[d] for d in doses])
fit_doses = np.array([d for d in doses if d in FIT_DOSES]); isi = np.array([ISI_MED[d] for d in fit_doses])
Sgrid, Tgrid = bc.value.values, bc.Tosc.values
def model_isi(c, K, T0):
    Sd = c * fit_doses / (K + fit_doses); return T0 * np.interp(np.log(Sd), np.log(Sgrid), Tgrid, left=np.nan, right=np.nan)
best_fit = None
for lc in np.linspace(np.log(0.1), np.log(300), 60):
    for lK in np.linspace(np.log(0.3), np.log(3000), 60):
        m = model_isi(np.exp(lc), np.exp(lK), 1.0)
        if np.any(np.isnan(m)): continue
        T0 = np.exp(np.mean(np.log(isi) - np.log(m)))
        err = np.sum((np.log(isi) - np.log(T0 * m)) ** 2)
        if best_fit is None or err < best_fit[0]: best_fit = (err, np.exp(lc), np.exp(lK), T0)
err, cfit, Kfit, T0fit = best_fit

# --- panels ---------------------------------------------------------------------------------
def p_direction(a):
    for fn, lab, col in [("m0ip3h_period_vs_stimulus_mech0.csv", "original (no coupling)", C_M0),
                         ("m0ip3h_period_vs_stimulus_mech1_gH5_dH3000.csv", r"threshold: $\kappa_H(P)$", C_M1),
                         ("m0ip3h_period_vs_stimulus_mech2_gH5_dH3000.csv", r"rate: $\varepsilon_H^{-1}(1+g_H\frac{P}{d_H+P})$", C_M2)]:
        o = osc(load(fn)); a.semilogx(o.value, o.Tosc, "o-", color=col, label=lab)
    a.set_ylim(0, None); a.set_xlabel(r"stimulus $\bar{S}$"); a.set_ylabel(r"$T_{\mathrm{osc}}$")
    a.set_title(r"Which coupling shortens $T_{\mathrm{osc}}$ with $\bar S$?"); a.legend(loc="center left", fontsize=17)

def p_rate_range(a):
    for fn, lab, col in [("m0ip3h_period_vs_stimulus_mech2_gH5_dH3000.csv", r"$\rho=2.4$ (Table 1), $g_H=5$", "#999999"),
                         ("m0ip3h_period_vs_stimulus_mech2_gH10_dH3000_rho0.024.csv", r"$\rho=0.024$, $g_H=10$", C_M1),
                         ("m0ip3h_period_vs_stimulus_mech2_gH10_dH3000_rho0.024_epsH300.csv", r"$\rho=0.024$, $g_H=10$, $\varepsilon_H=300$", C_M2)]:
        d = load(fn); o = osc(d); a.semilogx(o.value, o.Tosc, "o-", color=col, label=lab)
        lost = d[d.osc == 0]
        if len(lost): a.semilogx(lost.value, np.zeros(len(lost)), "x", color=col, ms=10, mew=2)
    a.set_ylim(0, None); a.set_xlabel(r"stimulus $\bar{S}$"); a.set_ylabel(r"$T_{\mathrm{osc}}$")
    a.set_title(r"Rate coupling: range of $\bar S$ with oscillation ($\times$: none)"); a.legend(loc="upper right", fontsize=16)

def p_grid_search(a):
    a.axhspan(DUTY_LO, DUTY_HI, color=TINT_DATA, lw=0, zorder=0, label=f"data duty {DUTY_LO:.2f}–{DUTY_HI:.2f}")
    a.axvline(max(ISI_MED.values()) / min(ISI_MED.values()), color=C_DATA, ls="--", lw=2, label=rf"data ratio {max(ISI_MED.values())/min(ISI_MED.values()):.1f}$\times$")
    for eh, m in zip(sorted(ok.epsH.unique()), ["o", "s", "^", "D"]):
        sub = ok[ok.epsH == eh]; a.plot(sub.T_ratio, sub.duty_mean, m, ms=9, alpha=1, label=rf"$\varepsilon_H={int(eh)}$")
    a.plot(best_cfg.T_ratio, best_cfg.duty_mean, "*", color="k", ms=22, label="best")
    a.set_ylim(0, 0.7); a.set_xlabel(r"period ratio $T_{\max}/T_{\min}$ over $\bar S\in[0.1,30]$"); a.set_ylabel(r"mean $W_{50}/T_{\mathrm{osc}}$")
    a.set_title(f"Grid search ({len(ok)} full-range, decreasing configs)"); a.legend(loc="center", bbox_to_anchor=(0.45, 0.42), ncol=2, fontsize=15)

def p_dose_map(a):
    Ufine = np.logspace(np.log10(2), np.log10(40), 100); Sfine = cfit * Ufine / (Kfit + Ufine)
    Tfine = T0fit * np.interp(np.log(Sfine), np.log(Sgrid), Tgrid, left=np.nan, right=np.nan)
    if ISI_Q1:
        errb = np.array([[ISI_MED[d] - ISI_Q1[d] for d in fit_doses], [ISI_Q3[d] - ISI_MED[d] for d in fit_doses]])
        a.errorbar(fit_doses, isi, yerr=errb, fmt="s", color=C_DATA, ms=10, capsize=6, elinewidth=2, label="data, constant dose (train medians, 95 % CI)")
        oth = [d for d in doses if d not in FIT_DOSES]
        if oth: a.errorbar(oth, [ISI_MED[d] for d in oth], yerr=[[ISI_MED[d] - ISI_Q1[d] for d in oth], [ISI_Q3[d] - ISI_MED[d] for d in oth]], fmt="s", color=C_DATA, mfc="white", mew=2.5, ms=10, capsize=6, elinewidth=2, label="step level only (not fitted)")
    else:
        a.plot(fit_doses, isi, "s", color=C_DATA, ms=10, label="data (median)")
    a.plot(Ufine, Tfine, "-", color=C_M2, label=rf"model, $t_{{\rm unit}}={T0fit:.2f}$ s, $K={Kfit:.1f}\,\mu$M")
    a.set_xscale("log"); a.set_xticks(doses); a.set_xticklabels([str(int(d)) for d in doses]); a.minorticks_off()
    a.set_ylim(0, None); a.set_xlabel(r"carbachol ($\mu$M)"); a.set_ylabel("ISI (s)")
    a.set_title(rf"Best config mapped to data: $\bar S = c\,U/(K+U)$"); a.legend(loc="upper right", fontsize=17)

def p_gain(a):
    for fn, lab, col in [("m0ip3h_period_vs_gH_S0.3.csv", r"$\bar S=0.3$", C_LO), ("m0ip3h_period_vs_gH_S30.csv", r"$\bar S=30$", C_HI)]:
        o = osc(load(fn)); a.plot(o.value, o.Tosc, "o-", color=col, label=lab)
    a.set_ylim(0, None); a.set_xlabel(r"coupling gain $g_H$"); a.set_ylabel(r"$T_{\mathrm{osc}}$")
    a.set_title(r"Gain sets the ratio ($\rho=0.024$, $d_H=3000$, $\varepsilon_H=300$)"); a.legend(loc="upper right")

def p_waveform(a):
    for S, lab, col in [(bc.value.min(), rf"model, $\bar S={bc.value.min():.1f}$", C_LO), (bc.value.max(), rf"model, $\bar S={bc.value.max():.0f}$", C_HI)]:
        r = trace(BEST, S)
        if r: a.plot(r[0], r[1], color=col, label=lab)
    sp = data_spike(target_isi=ISI_MED.get(15, 57.0))
    if sp is not None:
        ph, y, isi15 = sp; m = (ph >= -0.5) & (ph <= 0.5); a.plot(ph[m], y[m], color=C_DATA, lw=3, label=rf"data, 15 $\mu$M")
    a.set_xlim(-0.5, 0.5); a.set_ylim(0, 1.2); a.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.set_xlabel(r"phase $(t-t_{\mathrm{peak}})/T$"); a.set_ylabel("normalized amplitude")
    a.set_title("Best config: waveform still broad")
    a.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.995), frameon=False, columnspacing=0.8, handlelength=1.2, fontsize=17)

# --- figures --------------------------------------------------------------------------------
MAIN = [p_direction, p_dose_map]                        # rate coupling gives the right direction; best config on the dose data
SUPP = [p_rate_range, p_grid_search, p_gain, p_waveform]  # range, search, gain, waveform still broad
compose(MAIN, 2, f"{NAME}_main")
compose(SUPP, 2, f"{NAME}_supp")
print("best config:", BEST, f"| T ratio {best_cfg.T_ratio:.2f}, duty {best_cfg.duty_mean:.2f}")
print(f"dose mapping fit: c={cfit:.2f}, K={Kfit:.2f} uM, T0={T0fit:.3f} s, log-RMS err={np.sqrt(err/len(isi)):.3f}")
