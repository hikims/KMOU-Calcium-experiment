"""
plot_m0ip3h_amplitude_observation.py — amplitude and the fluorescence observation model.

Facts to reconcile
  data : ΔF/F0 spike amplitude ~0.9–1.5, roughly independent of dose; within-cell corr(ISI, amplitude) ≈ 0
  model: C_peak/C_rest ~ 400–1600 in the relaxation regime, C_peak falls ~4x from S̄=0.1 to S̄=30
Observation model (single-wavelength indicator, F_min normalized to 1):
  F(C) = 1 + R · C/(K_d + C),   ΔF/F0 = (F(C) − F(C_rest)) / F(C_rest)

Main figure  (figures/m0ip3h_amplitude_observation_main.*, 1x2)
  (a) amplitude vs dose: data vs model under the chosen (K_d, R)
  (b) ΔF/F0(t): model (chosen observation model, seconds via T0) vs a measured trace
Supplement   (figures/m0ip3h_amplitude_observation_supp.*, 2x2)
  (a) data: ΔF/F0 amplitude by dose (median, IQR) and within-cell corr(ISI, amplitude)
  (b) model: C_peak, C_rest vs dose (S̄ mapped by S̄ = cU/(K+U), T0 fitted to the ISI medians)
  (c) ΔF/F0 amplitude at 15 µM vs K_d for several R; data band
  (d) F-space duty W50/T at 15 µM vs K_d for several R; data band  -> feasible window

Run after ./run_m0ip3h_amplitude_sweeps.sh .  Also writes results/m0ip3h_observation_grid.csv.
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
C_DATA, C_REF, C_GREY = "#c0392b", "#2ca02c", "#888888"
TINT_DATA = "#f5d0cc"
RCOL = {1: "#1f77b4", 3: "#2ca02c", 10: "#e69f00", 30: "#9467bd"}
NAME, OUT, FIG, BIN = "m0ip3h_amplitude_observation", "results", "figures", "./m0ip3h_amplitude_sweep"
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

# --- data ------------------------------------------------------------------------------------------
ROOTS = ["../..", "../../data", "../data"]
def find(pattern):
    for r in ROOTS:
        hits = sorted(glob.glob(os.path.join(r, "*Spike Trains_data", pattern)))
        if hits: return hits
    return []
ISI_MED = {3: 206.0, 5: 76.0, 10: 60.0, 15: 57.0, 30: 30.6}
AMP = {3: (0.87, 0.67, 1.16), 5: (1.23, 1.01, 1.56), 10: (1.50, 1.16, 2.05), 15: (0.90, 0.42, 1.25), 30: (1.41, 1.16, 1.79)}
DUTY_LO, DUTY_HI, CORR_MED = 0.05, 0.15, 0.09
csvs = find(os.path.join("2", "*.csv"))
if csvs:
    P = pd.concat([pd.read_csv(f, sep=";", decimal=",", encoding="utf-8-sig").assign(file=os.path.basename(f)) for f in csvs]); g = P.groupby("stimulus")
    ISI_MED = {int(k): float(v) for k, v in g.ISI.median().items()}
    AMP = {int(k): (float(d.amplitude.median()), float(d.amplitude.quantile(.25)), float(d.amplitude.quantile(.75))) for k, d in g}
    dd = (g["half-width"].median() / g.ISI.median()); DUTY_LO, DUTY_HI = float(dd.min()), float(dd.max())
    cc = P.groupby(["file", "ST"]).apply(lambda d: np.corrcoef(d.ISI, d.amplitude)[0, 1] if len(d) >= 10 else np.nan).dropna()
    CORR_MED, CORR_Q1, CORR_Q3 = float(cc.median()), float(cc.quantile(.25)), float(cc.quantile(.75))
else:
    CORR_Q1, CORR_Q3 = -0.23, 0.36
doses = np.array(sorted(ISI_MED)); isi = np.array([ISI_MED[d] for d in doses])
med = [AMP[d][0] for d in doses]; err = np.array([[AMP[d][0] - AMP[d][1] for d in doses], [AMP[d][2] - AMP[d][0] for d in doses]])

def data_trace(dose_file="CChCB_2024_08_08_E1.txt", roi_name="#18 (New)", t0=2400, span=400, stim_t=900.0):
    hits = find(os.path.join("CChCB", dose_file)) or find(os.path.join("1", dose_file))
    if not hits: return None
    df = pd.read_csv(io.StringIO(open(hits[0], "rb").read().decode("utf-16")), sep="\t").dropna(axis=1, how="all")
    t = df["Time [s]"].values; col = roi_name if roi_name in df.columns else [c for c in df.columns if c.startswith("#")][0]
    F = pd.to_numeric(df[col], errors="coerce").values; F0 = np.nanmedian(F[t < stim_t]); y = (F - F0) / F0
    m = (t >= t0) & (t <= t0 + span); return t[m] - t0, y[m]

# --- model: dose mapping ------------------------------------------------------------------------------
def load(name): return pd.read_csv(os.path.join(OUT, name))
ref = load("m0ip3h_amp_period_vs_stimulus_ref.csv"); ref = ref[ref.osc == 1]
Sg, Tg = ref.value.values, ref.Tosc.values
def model_T(c, K): return np.interp(np.log(c * doses / (K + doses)), np.log(Sg), Tg, left=np.nan, right=np.nan)
best = None
for lc in np.linspace(np.log(0.1), np.log(300), 60):
    for lK in np.linspace(np.log(0.3), np.log(300), 60):
        m = model_T(np.exp(lc), np.exp(lK))
        if np.any(np.isnan(m)): continue
        T0 = np.exp(np.mean(np.log(isi) - np.log(m))); e = np.sum((np.log(isi) - np.log(T0 * m)) ** 2)
        if best is None or e < best[0]: best = (e, np.exp(lc), np.exp(lK), T0)
_, cfit, Kfit, T0 = best
S_of_dose = {int(d): cfit * d / (Kfit + d) for d in doses}
print(f"dose mapping: c={cfit:.2f}, K={Kfit:.1f} uM, T0={T0:.3f} s ; S̄(15 uM)={S_of_dose[15]:.2f}")

def trace(S):
    args = [BIN, "trace=1", f"S={S}"] + [f"{k}={v}" for k, v in {**NUM, **REF}.items()]
    tr = pd.read_csv(io.StringIO(subprocess.run(args, capture_output=True, text=True, check=True).stdout))
    return tr.t.values, tr.C.values
TR = {int(d): trace(S_of_dose[int(d)]) for d in doses}

def features(t, y):
    """period, peak, rest, W50/T of a periodic signal y(t) (last cycles)"""
    pk, _ = find_peaks(y, prominence=0.5 * (y.max() - y.min()))
    if len(pk) < 3: return np.nan, np.nan, np.nan, np.nan
    T = np.mean(np.diff(t[pk])); ip = pk[len(pk) // 2]
    lo = np.searchsorted(t, t[ip] - 0.5 * T); hi = np.searchsorted(t, t[ip] + 0.5 * T)
    seg = y[lo:hi]; base = seg.min(); half = base + 0.5 * (y[ip] - base)
    a = ip
    while a > lo and y[a] > half: a -= 1
    b = ip
    while b < hi - 1 and y[b] > half: b += 1
    return T, y[ip], base, (t[b] - t[a]) / T

def dFF(C, Kd, R):
    g = C / (Kd + C); g0 = g.min(); return R * (g - g0) / (1 + R * g0)

# --- observation grid at 15 uM ----------------------------------------------------------------------------
Kd_grid = np.logspace(-2, 2, 41); R_list = [1, 3, 10, 30]
t15, C15 = TR[15]; rows = []
for R in R_list:
    for Kd in Kd_grid:
        T, pk, base, duty = features(t15, dFF(C15, Kd, R)); rows.append(dict(R=R, Kd=Kd, amp=pk - base, duty_F=duty))
G = pd.DataFrame(rows); G.to_csv(os.path.join(OUT, "m0ip3h_observation_grid.csv"), index=False)
amp_lo, amp_hi = AMP[15][1], AMP[15][2]
feas = G[(G.amp >= amp_lo) & (G.amp <= amp_hi) & (G.duty_F >= DUTY_LO) & (G.duty_F <= DUTY_HI)]
if len(feas):
    # choose the feasible point whose amplitude is closest to the data median at 15 uM
    ch = feas.iloc[(feas.amp - AMP[15][0]).abs().argmin()]; Kd_c, R_c = float(ch.Kd), int(ch.R)
else:
    ch = G.iloc[(G.amp - AMP[15][0]).abs().argmin()]; Kd_c, R_c = float(ch.Kd), int(ch.R)
print(f"feasible (K_d, R) points: {len(feas)} ; chosen K_d={Kd_c:.3g}, R={R_c}")
Cduty = features(t15, C15)[3]
mamp = [features(TR[int(d)][0], dFF(TR[int(d)][1], Kd_c, R_c))[1] - features(TR[int(d)][0], dFF(TR[int(d)][1], Kd_c, R_c))[2] for d in doses]

# --- panels ---------------------------------------------------------------------------------
def p_data_amp(a):
    a.errorbar(doses, med, yerr=err, fmt="s-", color=C_DATA, ms=10, capsize=6, elinewidth=2)
    a.set_xscale("log"); a.set_xticks(doses); a.set_xticklabels([str(int(d)) for d in doses]); a.minorticks_off()
    a.set_ylim(0, 2.5); a.set_xlabel(r"carbachol ($\mu$M)"); a.set_ylabel(r"spike amplitude $\Delta F/F_0$")
    a.set_title("Data: amplitude vs dose (median, IQR)")
    a.text(0.03, 0.05, rf"within-cell corr(ISI, amp): median {CORR_MED:.2f}" + "\n" + rf"IQR [{CORR_Q1:.2f}, {CORR_Q3:.2f}]", transform=a.transAxes, fontsize=16)

def p_model_C(a):
    Cpk = [features(*TR[int(d)])[1] for d in doses]; Crs = [features(*TR[int(d)])[2] for d in doses]
    a.semilogy(doses, Cpk, "o-", color=C_REF, label=r"$C_{\mathrm{peak}}$"); a.semilogy(doses, Crs, "s--", color=C_GREY, label=r"$C_{\mathrm{rest}}$")
    a.set_xscale("linear"); a.set_xticks(doses); a.set_xticklabels([str(int(d)) for d in doses]); a.set_xlim(0, 32)
    a.set_xlabel(r"carbachol ($\mu$M)"); a.set_ylabel(r"$C$ (model units)")
    a.set_title(rf"Model: $C_{{\mathrm{{peak}}}}/C_{{\mathrm{{rest}}}}$ = {Cpk[0]/Crs[0]:.0f}–{Cpk[-1]/Crs[-1]:.0f}"); a.legend(loc="center right")

def p_amp_vs_Kd(a):
    a.axhspan(amp_lo, amp_hi, color=TINT_DATA, lw=0, zorder=0, label=r"data IQR (15 $\mu$M)")
    for R in R_list:
        s = G[G.R == R]; a.semilogx(s.Kd, s.amp, "-", color=RCOL[R], label=rf"$R={R}$")
    a.plot(Kd_c, float(ch.amp), "*", color="k", ms=22, label="chosen")
    a.set_ylim(0, 25); a.set_xlabel(r"indicator $K_d$ (model $C$ units)"); a.set_ylabel(r"model $\Delta F/F_0$ amplitude")
    a.set_title(r"Observation model: amplitude at 15 $\mu$M"); a.legend(loc="upper left", fontsize=16)

def p_duty_vs_Kd(a):
    a.axhspan(DUTY_LO, DUTY_HI, color=TINT_DATA, lw=0, zorder=0, label="data duty")
    for R in R_list:
        s = G[G.R == R]; a.semilogx(s.Kd, s.duty_F, "-", color=RCOL[R], label=rf"$R={R}$")
    a.axhline(Cduty, color=C_REF, ls=":", lw=2, label=rf"$C$-space duty {Cduty:.2f}")
    a.plot(Kd_c, float(ch.duty_F), "*", color="k", ms=22)
    a.set_ylim(0, 0.4); a.set_xlabel(r"indicator $K_d$ (model $C$ units)"); a.set_ylabel(r"$W_{50}/T$ in $\Delta F/F_0$")
    a.set_title("Saturation broadens F-space spike"); a.legend(loc="upper right", fontsize=16)
    a.set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])

def p_amp_vs_dose(a):
    a.errorbar(doses, med, yerr=err, fmt="s", color=C_DATA, ms=10, capsize=6, elinewidth=2, label="data")
    a.plot(doses, mamp, "o-", color=C_REF, label=rf"model, $K_d={Kd_c:.2g}$, $R={R_c}$")
    a.set_xscale("log"); a.set_xticks(doses); a.set_xticklabels([str(int(d)) for d in doses]); a.minorticks_off()
    a.set_ylim(0, 3.0); a.set_yticks([0, 0.5, 1.0, 1.5, 2.0, 2.5]); a.set_xlabel(r"carbachol ($\mu$M)"); a.set_ylabel(r"$\Delta F/F_0$ amplitude")
    a.set_title("Amplitude vs dose: data vs model"); a.legend(loc="upper right")

def p_traces(a):
    tt, CC = TR[15]; y = dFF(CC, Kd_c, R_c); ts = (tt - tt[0]) * T0; m = ts <= 400
    a.plot(ts[m], y[m], color=C_REF, label=rf"model, 15 $\mu$M ($t_{{\rm unit}}={T0:.2f}$ s)")
    dtr = data_trace()
    if dtr is not None: a.plot(dtr[0], dtr[1] + 1.6, color=C_DATA, label=r"data, 15 $\mu$M, ROI #18 (offset)")
    a.set_ylim(0, 4.2); a.set_yticks([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]); a.set_xlabel(r"$t$ (s)"); a.set_ylabel(r"$\Delta F/F_0$")
    a.set_title(r"$\Delta F/F_0(t)$: model vs data"); a.legend(loc="upper right", fontsize=16)

# --- figures --------------------------------------------------------------------------------
MAIN = [p_amp_vs_dose, p_traces]                          # amplitude vs dose and ΔF/F0(t) under the chosen observation model
SUPP = [p_data_amp, p_model_C, p_amp_vs_Kd, p_duty_vs_Kd]  # how (K_d, R) was chosen
compose(MAIN, 2, f"{NAME}_main")
compose(SUPP, 2, f"{NAME}_supp")
print("model amplitude by dose under chosen observation model:", dict(zip(doses.astype(int), np.round(mamp, 2))))
