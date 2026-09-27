"""
plot_m1_sigma_tav_slope.py — the variability law σ_ISI ∝ T_av across cells and what it says about the noise.

Data: per train (record, ST, dose; ../data_dose_units/results/train_table.csv, union of all files, >= 10 spikes)
T_av and the sd of the detrended ISI (a linear trend in time removed per train); log σ – log T_av slope with a
train-bootstrap CI, and the linear fit sd = a (T_av − T_ref).
Model: cells differ only in the recovery time scale ε_H (sweep 200–8000).  Four noise conventions:
  M1 with common rescaled intensity ν (σ_H = ν/√ε_H)  -> theory: sd ∝ T_av (slope 1, constant CV)
  M1 with common σ_H                                  -> sd ∝ ε_H^{3/2} (slope 1.5)
  M2 (multiplicative noise in the slow P)             -> slow modulation, slope ≈ 1 but serial correlation ≈ 1
  Ca-noise
Model time is converted with T0 from ../data_dose_units/results/dose_map_refit.csv (0.646 s; fallback 0.533 s).

Panels (functions; MAIN / SUPP at the end)
  law        sd vs T_av (log-log): data trains by dose + model lines
  cv_cells   CV vs T_av: data + model
  slopes     log-log slopes: data (bootstrap CI) and the four model conventions
Run after ./run_m1_slope_sweeps.sh .  Writes figures/m1_sigma_tav_slope_main.{png,pdf,eps} (law, slopes),
figures/m1_sigma_tav_slope_supp.{png,pdf,eps} (cv_cells),
results/m1_slope_data_trains.csv, results/m1_slope_summary.csv
"""
import os
import numpy as np, pandas as pd
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
OUT, FIG = "results", "figures"; os.makedirs(FIG, exist_ok=True)
T0 = 0.533
_fn = "../data_dose_units/results/dose_map_refit.csv"
if os.path.exists(_fn): T0 = float(pd.read_csv(_fn).T0_s.iloc[0])
DOSE_COL = {3: "#1f77b4", 5: "#17becf", 10: "#2ca02c", 15: "#e69f00", 30: "#c0392b"}
MODELS = [("m1_slope_M1_nu_common.csv", r"M1, $\nu$ common (cells: $\varepsilon_H$)", "#2ca02c", "-"),
          ("m1_slope_M1_sig_common.csv", r"M1, $\sigma_H$ common (cells)", "#2ca02c", "--"),
          ("m1_slope_M2.csv", "M2, P-noise (cells)", "#9467bd", "-"),
          ("m1_slope_Ca.csv", "Ca-noise (cells)", "#e69f00", "-"),
          ("m1_slope_M1_nu_common_dose.csv", r"M1, fixed $\sigma_H$ (dose: $\bar S$)", "#2ca02c", ":"),
          ("m1_slope_M2_dose.csv", "M2, P-noise (dose)", "#9467bd", ":")]
MODELS = [m for m in MODELS if os.path.exists(os.path.join("results", m[0]))]

# --- data: the authors' segments (data_dose_units/results/segment_table.csv; sd after detrending inside the segment);
#     the union trains (detrended) are the robustness check --------------------------------------------------------
ST_ = "../data_dose_units/results/segment_table.csv"; TT = "../data_dose_units/results/train_table.csv"
if not os.path.exists(ST_): raise SystemExit("run ../data_dose_units/data_dose_units.py first (segment_table.csv)")
S_ = pd.read_csv(ST_); S_ = S_[S_.n >= 10]
D = pd.DataFrame(dict(record=S_.record, ST=S_.ST, segment=S_.segment, dose=S_.dose.astype(int), n=S_.n, Tav=S_.Tav, sd=S_.sd_detr, sd_raw=S_.sd)).reset_index(drop=True)
D.to_csv(os.path.join(OUT, "m1_slope_data_trains.csv"), index=False)
slope_raw = np.polyfit(np.log(D.Tav.values), np.log(D.sd_raw.values), 1)[0]
T_ = pd.read_csv(TT); T_ = T_[T_.n >= 10]; slope_train = np.polyfit(np.log(T_.Tav.values), np.log(T_.sd_detr.values), 1)[0]
lx, ly = np.log(D.Tav.values), np.log(D.sd.values)
slope, icpt = np.polyfit(lx, ly, 1)
rng = np.random.default_rng(0); bs = np.array([np.polyfit(lx[i], ly[i], 1)[0] for i in (rng.integers(0, len(D), len(D)) for _ in range(2000))])
lin_a, lin_b = np.polyfit(D.Tav, D.sd, 1); Tref = -lin_b / lin_a
Xfe = np.column_stack([lx] + [(D.dose == d).astype(float) for d in sorted(D.dose.unique())]); slope_within = float(np.linalg.lstsq(Xfe, ly, rcond=None)[0][0])
print(f"data: {len(D)} segments; log-log slope {slope:.3f} (bootstrap 95% CI {np.quantile(bs, .025):.2f}–{np.quantile(bs, .975):.2f}); "
      f"linear sd = {lin_a:.3f} (T_av − {Tref:.1f} s); raw (undetrended) slope {slope_raw:.3f}; within-dose slope {slope_within:.3f}")

# --- model -----------------------------------------------------------------------------------------
summ = [dict(model="data", slope=slope, lo=np.quantile(bs, .025), hi=np.quantile(bs, .975), Tref_s=Tref, a=lin_a),
        dict(model="data (raw, with drift)", slope=slope_raw, lo=np.nan, hi=np.nan, Tref_s=np.nan, a=np.nan),
        dict(model="data (union trains, detrended)", slope=slope_train, lo=np.nan, hi=np.nan, Tref_s=np.nan, a=np.nan),
        dict(model="data (within dose, fixed effects)", slope=slope_within, lo=np.nan, hi=np.nan, Tref_s=np.nan, a=np.nan)]
M = {}
for fn, lab, col, ls in MODELS:
    m = pd.read_csv(os.path.join(OUT, fn)); m = m[m.n > 0]; M[fn] = m
    s = np.polyfit(np.log(m.Tav), np.log(m.sd), 1)[0]; a, b = np.polyfit(m.Tav * T0, m.sd * T0, 1)
    summ.append(dict(model=lab, slope=s, lo=np.nan, hi=np.nan, Tref_s=-b / a, a=a))
S = pd.DataFrame(summ); S.to_csv(os.path.join(OUT, "m1_slope_summary.csv"), index=False); print(S.round(3).to_string(index=False))

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

def law(a):
    for dose, col in DOSE_COL.items():
        s = D[D.dose == dose]; a.loglog(s.Tav, s.sd, "o", color=col, ms=7, mec="none", label=rf"{dose} $\mu$M")
    xx = np.array([10, 800]); a.loglog(xx, np.exp(icpt) * xx ** slope, "k-", lw=2, label=rf"data fit, slope {slope:.2f}")
    for fn, lab, col, ls in MODELS:
        m = M[fn]; a.loglog(m.Tav * T0, m.sd * T0, ls, color=col, lw=3, label=lab)
    a.set_xlim(10, 800); a.set_ylim(1, 300); a.set_xlabel(r"$T_{\mathrm{av}}$ (s)"); a.set_ylabel(r"$\sigma_{\mathrm{ISI}}$ (s)")
    a.set_title(r"Variability law $\sigma \propto T_{\mathrm{av}}$"); a.legend(loc="lower right", fontsize=10, ncol=2)

def cv_cells(a):
    for dose, col in DOSE_COL.items():
        s = D[D.dose == dose]; a.semilogx(s.Tav, s.sd / s.Tav, "o", color=col, ms=7, mec="none")
    for fn, lab, col, ls in MODELS:
        m = M[fn]; a.semilogx(m.Tav * T0, m.CV, ls, color=col, lw=3, label=lab)
    a.set_xlim(10, 800); a.set_ylim(0, 0.8); a.set_xlabel(r"$T_{\mathrm{av}}$ (s)"); a.set_ylabel("CV of ISI")
    a.set_title("CV across cells"); a.legend(loc="upper left", fontsize=15)

def slopes(a):
    labels = ["data"] + [lab for _, lab, _, _ in MODELS]; cols = ["#c0392b"] + [col for _, _, col, _ in MODELS]
    vals = S[~S.model.str.startswith("data (")].slope.values
    a.bar(range(len(vals)), vals, color=cols, edgecolor="k")
    a.errorbar([0], [vals[0]], yerr=[[vals[0] - S.lo[0]], [S.hi[0] - vals[0]]], fmt="none", ecolor="k", capsize=8, elinewidth=2.5)
    a.axhline(1.0, color="k", ls="--", lw=2); a.axhline(0.5, color="k", ls=":", lw=2)
    a.set_xticks(range(len(vals))); a.set_xticklabels(labels, fontsize=14, rotation=35, ha="right")
    a.set_ylim(0, 1.8); a.set_ylabel(r"log $\sigma$ – log $T_{\mathrm{av}}$ slope"); a.set_title("Slope of the law")

# --- figures ----------------------------------------------------------------------------------------
MAIN = [law, slopes]      # the law in the data and which noise convention reproduces its slope
SUPP = [cv_cells]         # the same information as CV vs T_av
compose(MAIN, 2, "m1_sigma_tav_slope_main")
compose(SUPP, 1, "m1_sigma_tav_slope_supp", w=8.5)
