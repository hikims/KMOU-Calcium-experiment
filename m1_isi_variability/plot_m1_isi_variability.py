"""
plot_m1_isi_variability.py — M1: bounded H-noise reproduces the within-cell ISI statistics; where the noise sits.

Data statistics from ../data_dose_units/results/train_table.csv (union of all processed files; one train = one
(record, ST, dose) with >= 10 spikes):
  CV of ISI, skewness, amplitude CV, corr(preceding ISI, amplitude), serial corr(ISI_k, ISI_{k+1})
  (all statistics after removing a linear trend in time per train — the union of the processed files makes trains up to
2 h long and slow drift inflates CV, skew and serial correlation; raw values are kept in the table for reference)
Model: the stereotyped reference regime + noise in H (M1), in P (M2) or in Ca (comparison), each tuned so that
the ISI CV equals the data median.  The four other statistics then discriminate the noise location.

Panels (functions; MAIN / SUPP at the end)
  trace         trace at the tuned H-noise: C(t) and H(t), H_on of the deterministic cycle
  cv_vs_sigma   CV and skewness vs sigma_H, data bands -> tuned sigma_H*
  isi_dist      ISI distribution (ISI / mean) at sigma_H*: model vs data (per-train normalised, pooled)
  noise_loc     the five statistics for H-, Ca-, P-noise vs data (median, IQR)

Run after ./run_m1_isi_sweeps.sh (calls ./m1_isi_sde for the ISI lists and the trace).  Writes
figures/m1_isi_variability_main.{png,pdf,eps} (trace, isi_dist, noise_loc), figures/m1_isi_variability_supp.{png,pdf,eps}
(cv_vs_sigma), results/m1_data_train_stats.csv, results/m1_isi_list_{H,C,P}.csv,
results/m1_trace_H.csv, results/m1_noise_location_table.csv
"""
import io, os, subprocess
import numpy as np, pandas as pd
from scipy.stats import skew, gamma as gamma_dist
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
C_DATA, C_H, C_C, C_P, C_GREY = "#c0392b", "#2ca02c", "#e69f00", "#9467bd", "#888888"
TINT_DATA = "#f5d0cc"
OUT, FIG, BIN = "results", "figures", "./m1_isi_sde"
# the panels call the compiled kernel; build it if it is missing (Makefile if present, else a direct cc call)
def _build(binary):
    src = os.path.basename(binary) + ".c"
    if os.path.exists("Makefile"): subprocess.run(["make", "-s"], check=True); return
    for flags in (["-O2", "-fopenmp"], ["-O2", "-Xpreprocessor", "-fopenmp", "-lomp"], ["-O2"]):
        if subprocess.run(["cc"] + flags + ["-o", binary, src, "-lm"], capture_output=True).returncode == 0: return
    raise SystemExit(f"could not compile {src} (need a C compiler; on macOS: xcode-select --install, brew install libomp)")
if not os.path.exists(BIN): _build(BIN)

os.makedirs(FIG, exist_ok=True)
REF = dict(mech=2, rho=0.024, gH=30, dH=10000, epsH=1000, nA=2, kA=0.4, Jin=0.02, kS=0.25, beta=4, alpha=20, qoff=0, rC=300, kR=1, S=1)
NUM = dict(Tend=120000, Ttrans=3000, dt=0.01)

# --- data: the authors' segments (data_dose_units/results/segment_table.csv), detrended -------------------------------
# Unit = one stationary segment as delimited by the data authors (file; base file = segment 1); a linear trend in time
# is removed inside each segment before CV, skewness, corr(ISI, amp) and the serial correlation are computed
# (the same formulas are applied to the synthetic records below).  Trains (union, detrended) are the robustness check.
ST_ = "../data_dose_units/results/segment_table.csv"; EV = "../data_dose_units/results/events_canonical.csv"; UC = "../data_dose_units/results/unit_comparison.csv"
if not (os.path.exists(ST_) and os.path.exists(EV)): raise SystemExit("run ../data_dose_units/data_dose_units.py first (segment_table.csv, events_canonical.csv)")
D = pd.read_csv(ST_); D = D[D.n >= 10].copy()
D["CV_raw"], D["skew_raw"], D["corr_raw"], D["serial_raw"] = D.CV, D["skew"], D.corr_isi_amp, D.serial
D["CV"], D["skew"], D["corr_isi_amp"], D["serial"] = D.CV_detr, D.skew_detr, D.corr_detr, D.serial_detr
D.to_csv(os.path.join(OUT, "m1_data_train_stats.csv"), index=False)
E = pd.read_csv(EV); E["segment"] = E.file.str.replace(r"\.csv$", "", regex=True); E = E[E.isi_diff.notna()]
def detrended(x, t):
    """ISI with a linear trend in time removed, mean kept (identical to data_dose_units)"""
    if len(x) >= 4: b = np.polyfit(t, x, 1); return x - np.polyval(b, t) + x.mean()
    return x
pooled = np.concatenate([(lambda x: x / x.mean())(detrended(g.isi_diff.values, g.time.values)) for _, g in E.groupby(["train", "segment"]) if len(g) >= 10])
STATS = ["CV", "skew", "amp_cv", "corr_isi_amp", "serial"]
dmed = {s: D[s].median() for s in STATS + ["serial_raw", "CV_raw", "skew_raw", "corr_raw"]}
dq = {s: (D[s].quantile(.25), D[s].quantile(.75)) for s in STATS + ["serial_raw", "CV_raw", "skew_raw", "corr_raw"]}
SEG_N = D.n.values                                   # segment lengths (spikes): the synthetic records are cut the same way
print(f"data (segments, detrended): {len(D)} segments, median {int(np.median(SEG_N))} spikes; " + ", ".join(f"{s} {dmed[s]:.3f} [{dq[s][0]:.2f},{dq[s][1]:.2f}]" for s in STATS)
      + f"; raw: CV {dmed['CV_raw']:.3f}, skew {dmed['skew_raw']:.3f}, serial {dmed['serial_raw']:.3f}")
if os.path.exists(UC): print(pd.read_csv(UC).round(3).to_string(index=False))

# --- model sweeps -> tuned noise amplitudes -----------------------------------------------------------
def load(name): return pd.read_csv(os.path.join(OUT, name))
SW = {"H": load("m1_isi_vs_sigH.csv"), "C": load("m1_isi_vs_sigC.csv"), "P": load("m1_isi_vs_sigP.csv")}
KEY = {"H": "sigH", "C": "sigC", "P": "sigP"}
def tune(sw, target):
    """noise amplitude at which the model CV equals the target (log-interpolation on the increasing branch)"""
    sw = sw[sw.n > 0]; cv = sw.CV.values; v = sw.value.values
    i = np.argmax(cv >= target) if np.any(cv >= target) else len(cv) - 1
    if i == 0: return v[0]
    return float(np.exp(np.interp(target, cv[i - 1:i + 1], np.log(v[i - 1:i + 1]))))
SIG = {k: tune(SW[k], dmed["CV"]) for k in SW}
print("tuned noise amplitudes (CV = data median):", {k: round(v, 4) for k, v in SIG.items()})

def run(extra, mode):
    args = [BIN] + ([mode] if mode else []) + [f"{k}={v}" for k, v in {**REF, **NUM, **extra}.items()]
    return pd.read_csv(io.StringIO(subprocess.run(args, capture_output=True, text=True, check=True).stdout))
LISTS, STAT = {}, {}
for k in SW:
    LISTS[k] = run({KEY[k]: SIG[k]}, "isi=1"); LISTS[k].to_csv(os.path.join(OUT, f"m1_isi_list_{k}.csv"), index=False)
    STAT[k] = run({KEY[k]: SIG[k]}, None).iloc[0]
TR = run({"sigH": SIG["H"], "tspan": 700}, "trace=1"); TR.to_csv(os.path.join(OUT, "m1_trace_H.csv"), index=False)
det = run({}, None).iloc[0]                       # deterministic reference (all sigmas 0)
def unit_stats(x, a, t):
    """CV, skew, amp CV, corr(ISI, amp), serial corr of one unit after detrending (the data formulas)"""
    xd = detrended(x, t); m = xd.mean(); r = xd - m
    return dict(CV=xd.std(ddof=1) / m, skew=float(np.mean((r / r.std()) ** 3)) if r.std() > 0 else np.nan, amp_cv=a.std(ddof=1) / a.mean(),
                corr_isi_amp=np.corrcoef(xd, a)[0, 1], serial=np.corrcoef(r[:-1], r[1:])[0, 1])
def segmented(lst, nrep=3, seed=0):
    """cut the long synthetic ISI list into consecutive segments whose lengths follow the data segments; per-unit stats"""
    rng = np.random.default_rng(seed); x, a, t = lst.ISI.values, lst.amp.values, lst.t.values; rows = []
    for _ in range(nrep):
        i = 0
        for n_ in rng.permutation(SEG_N):
            if i + n_ > len(x): break
            rows.append(unit_stats(x[i:i + n_], a[i:i + n_], t[i:i + n_])); i += n_
    return pd.DataFrame(rows)
SEGM = {k: segmented(LISTS[k]) for k in SW}
tab = pd.DataFrame({k: {"sigma": SIG[k], "CV": STAT[k].CV, "skew": STAT[k]["skew"], "amp_cv": STAT[k].amp_cv,
                        "corr_isi_amp": STAT[k].corr_isi_amp, "serial": STAT[k].corr_isi_next} for k in SW}).T
for k in SW:                                                                                       # same unit length + detrending as the data
    tab.loc[f"{k} (segmented)"] = {"sigma": SIG[k], **{s: SEGM[k][s].median() for s in STATS}}
tab.loc["data (median)"] = {"sigma": np.nan, **{s: dmed[s] for s in STATS}}                       # detrended segment statistics
tab.loc["data (raw)"] = {"sigma": np.nan, "CV": dmed["CV_raw"], "skew": dmed["skew_raw"], "amp_cv": dmed["amp_cv"], "corr_isi_amp": dmed["corr_raw"], "serial": dmed["serial_raw"]}
tab.to_csv(os.path.join(OUT, "m1_noise_location_table.csv")); print(tab.round(3).to_string())

# --- panels ----------------------------------------------------------------------------------------
def letter(a, s): a.text(-0.13, 1.06, s, transform=a.transAxes, fontsize=28, fontweight="bold", va="bottom", ha="left")
def compose(panels, ncol, name, w=8.0, h=7.5):
    n = len(panels); nrow = -(-n // ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(w * ncol, h * nrow), squeeze=False)
    for k, fn in enumerate(panels):
        a = axs[k // ncol, k % ncol]; fn(a); letter(a, f"({chr(97 + k)})")
    for k in range(n, nrow * ncol): axs[k // ncol, k % ncol].axis("off")
    plt.tight_layout(w_pad=2.5, h_pad=2.5)
    for ext in ("png", "pdf", "eps"): fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig); print(f"wrote figures/{name}.{{png,pdf,eps}}")

def trace_panel(a):
    a2 = a.twinx(); a2.grid(False)
    a.plot(TR.t, TR.C, color=C_H, label=r"$C$")
    a2.plot(TR.t, TR.H, color="k", lw=2, label=r"$H$"); a2.axhline(det.Hon_det, color="k", ls="--", lw=2, label=r"$H_{\mathrm{on}}$ (deterministic)")
    a.set_xlim(0, 700); a.set_ylim(0, 2.2 * det.amp_det); a2.set_ylim(0, 1.0)
    a.set_xlabel(r"$t$ (model units)"); a.set_ylabel(r"$C$", color=C_H); a2.set_ylabel(r"$H$")
    a.set_title(rf"H-noise, $\sigma_H={SIG['H']:.3f}$, $\bar S=1$")
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels(); a.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=16)

def cv_vs_sigma(a):
    sw = SW["H"][SW["H"].n > 0]
    a.axhspan(dq["CV"][0], dq["CV"][1], color=TINT_DATA, lw=0, zorder=0, label="data CV (IQR)")
    a.semilogx(sw.value, sw.CV, "o-", color=C_H, label="CV (model)")
    a.semilogx(sw.value, sw["skew"], "s--", color=C_GREY, label="skewness (model)")
    a.axhline(dmed["skew"], color=C_DATA, ls=":", lw=2.5, label="data skewness (median)")
    a.axvline(SIG["H"], color="k", ls="-.", lw=1.8); a.text(SIG["H"] * 1.1, 1.75, rf"$\sigma_H^*={SIG['H']:.3f}$", fontsize=17)
    a.set_ylim(0, 2.2); a.set_xlabel(r"$\sigma_H$"); a.set_ylabel("CV, skewness of ISI"); a.set_title("ISI variability vs noise amplitude")
    a.legend(loc="upper left", fontsize=15)

def isi_dist(a):
    bins = np.linspace(0, 3, 37)
    a.hist(pooled, bins=bins, density=True, color=TINT_DATA, edgecolor=C_DATA, label=f"data ({len(D)} segments, detrended)")
    xm = LISTS["H"].ISI.values / LISTS["H"].ISI.mean()
    a.hist(xm, bins=bins, density=True, histtype="step", color=C_H, lw=3, label=rf"model, H-noise $\sigma_H^*$")
    xx = np.linspace(0.01, 3, 300); k_g = 1 / dmed["CV"] ** 2
    a.plot(xx, gamma_dist.pdf(xx, k_g, scale=1 / k_g), "k--", lw=2, label=rf"gamma, CV={dmed['CV']:.2f}")
    a.set_xlim(0, 3); a.set_xlabel("ISI / mean ISI"); a.set_ylabel("density"); a.set_title("ISI distribution shape")
    a.legend(loc="upper right", fontsize=15)

def noise_loc(a):
    xs = np.arange(len(STATS)); w = 0.2
    for j, (k, col, lab) in enumerate((("H", C_H, r"H-noise (M1)"), ("C", C_C, "Ca-noise"), ("P", C_P, "P-noise (M2)"))):
        med = [SEGM[k][s].median() for s in STATS]; err = np.array([[SEGM[k][s].median() - SEGM[k][s].quantile(.25) for s in STATS], [SEGM[k][s].quantile(.75) - SEGM[k][s].median() for s in STATS]])
        a.bar(xs + (j - 1) * w, med, w, color=col, label=lab); a.errorbar(xs + (j - 1) * w, med, yerr=err, fmt="none", ecolor="k", elinewidth=1.2, capsize=3)
    med = [dmed[s] for s in STATS]; err = np.array([[dmed[s] - dq[s][0] for s in STATS], [dq[s][1] - dmed[s] for s in STATS]])
    a.errorbar(xs + 2 * w, med, yerr=err, fmt="D", color=C_DATA, ms=11, capsize=5, elinewidth=2, label="data (median, IQR)")
    a.axhline(0, color="k", lw=1)
    a.set_xticks(xs + 0.5 * w); a.set_xticklabels(["CV", "skew", "amp CV", r"$\rho$(ISI, amp)", r"$\rho_1$ (serial)"], fontsize=16, rotation=12)
    a.set_ylim(-1.05, 2.2); a.set_ylabel("value (per segment, detrended)"); a.set_title("Specified noise models vs data")
    a.legend(loc="upper right", fontsize=15, ncol=1)

# --- figures ----------------------------------------------------------------------------------------
MAIN = [trace_panel, isi_dist, noise_loc]     # spikes start at H_on; the ISI shape; the noise sits in H
SUPP = [cv_vs_sigma]                          # how sigma_H* was tuned
compose(MAIN, 3, "m1_isi_variability_main")
compose(SUPP, 1, "m1_isi_variability_supp", w=8.5)
