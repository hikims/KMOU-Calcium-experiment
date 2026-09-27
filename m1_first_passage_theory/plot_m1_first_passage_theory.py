"""
plot_m1_first_passage_theory.py — the first-passage reduction (Proposition 1): closed moments vs the full SDE.

Reduced process between spikes (rescaled time s' = s/ε_eff):  dH = (1−H) ds' + ν √(H(1−H)) dW,  H0 = H_min -> H* = H_on
  ISI = T_sp + ε_eff τ_FP,   ε_eff = ε_H / [fP fC Q]_rest,   ν = σ_H √ε_eff
  T_sp = T_det − ε_eff ln((1−H0)/(1−H*))  (spike + early-recovery time, from the deterministic period of each cell)
Constants of the stereotyped regime at S̄ = 1 (m0ip3h_amplitude_invariance, S=1): C_rest = 0.0196, P̄ = 2108,
H_min = 0.064, H_on = 0.408, T_det = 83.5.

Panels (functions; MAIN / SUPP at the end)
  moments     mean and CV of τ_FP vs ν (quadrature), Monte-Carlo of the 1-D SDE (points), full model (ISI CV vs its ν)
  cv_vs_Hs    CV vs H* for three ν, H0 fixed: which (ν, H*) give the data CV (band)
  reduction   reduction check: full-model T_av and σ_ISI (σ_H sweep and ε_H sweep) vs the theory prediction, y = x
Run after ./run_m1_fpt.sh .  Writes figures/m1_first_passage_theory_main.{png,pdf,eps} (moments, reduction),
figures/m1_first_passage_theory_supp.{png,pdf,eps} (cv_vs_Hs), results/m1_fpt_reduction_check.csv
"""
import io, os, subprocess
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
OUT, FIG, FPT = "results", "figures", "./m1_fpt_jacobi"; os.makedirs(FIG, exist_ok=True)
# the panels call the compiled kernel; build it if it is missing (Makefile if present, else a direct cc call)
def _build(binary):
    src = os.path.basename(binary) + ".c"
    if os.path.exists("Makefile"): subprocess.run(["make", "-s"], check=True); return
    for flags in (["-O2", "-fopenmp"], ["-O2", "-Xpreprocessor", "-fopenmp", "-lomp"], ["-O2"]):
        if subprocess.run(["cc"] + flags + ["-o", binary, src, "-lm"], capture_output=True).returncode == 0: return
    raise SystemExit(f"could not compile {src} (need a C compiler; on macOS: xcode-select --install, brew install libomp)")
if not os.path.exists(FPT): _build(FPT)

C_TH, C_MC, C_FULL, C_DATA, TINT = "k", "#1f77b4", "#2ca02c", "#c0392b", "#f5d0cc"
# regime constants (S̄ = 1)
C_REST, PBAR, H0, HS, TDET = 0.0196, 2108.0, 0.064, 0.408, 83.5
gH, dH, rC, kR, kH = 30.0, 10000.0, 300.0, 1.0, 0.15
fP = 1 + gH * PBAR / (dH + PBAR); fC = 1 + rC * C_REST ** 2 / (kR ** 2 + C_REST ** 2); Q = kH ** 2 / (kH ** 2 + C_REST ** 2)
RATE = fP * fC * Q                                   # ε_eff = ε_H / RATE
mean_det = np.log((1 - H0) / (1 - HS)); EPS_EFF_1000 = 1000.0 / RATE
T_SP = TDET - EPS_EFF_1000 * mean_det
print(f"fP={fP:.2f}, fC={fC:.3f}, Q={Q:.3f} -> ε_eff = ε_H/{RATE:.2f}; at ε_H=1000: ε_eff={EPS_EFF_1000:.1f}, T_sp={T_SP:.1f}, ε_eff·ln = {EPS_EFF_1000*mean_det:.1f}")
DATA_CV = (0.14, 0.31, 0.196)                        # IQR and median of the per-segment CV (data_dose_units, detrended)
_uc = "../data_dose_units/results/segment_table.csv"
if os.path.exists(_uc):
    _d = pd.read_csv(_uc); _d = _d[_d.n >= 10].CV_detr; DATA_CV = (float(_d.quantile(.25)), float(_d.quantile(.75)), float(_d.median()))

def load(n): return pd.read_csv(os.path.join(OUT, n))
th = load("m1_fpt_vs_nu.csv"); mc = load("m1_fpt_mc.csv")
full_s = load("m1_full_vs_sigH.csv"); full_s = full_s[full_s.n > 0]
full_e = load("m1_full_vs_epsH.csv"); full_e = full_e[full_e.n > 0]

def theory(nu):
    out = subprocess.run([FPT, f"nu={nu}", f"H0={H0}", f"Hs={HS}"], capture_output=True, text=True, check=True).stdout
    return pd.read_csv(io.StringIO(out)).iloc[0]

# reduction check: predicted T_av, sd for every full-model point
rows = []
for name, df, eps_col in (("sigH sweep", full_s, None), ("epsH sweep", full_e, "value")):
    for _, r in df.iterrows():
        epsH = r.value if eps_col else 1000.0; sigH = r.value if not eps_col else 0.822 / np.sqrt(epsH)
        eps_eff = epsH / RATE; nu = sigH * np.sqrt(eps_eff); t = theory(nu)
        Tsp = r.Tdet - eps_eff * mean_det                  # spike + early-recovery time of that cell (deterministic period minus the recovery time)
        rows.append(dict(sweep=name, epsH=epsH, sigH=sigH, nu=nu, Tav_full=r.Tav, sd_full=r.sd, CV_full=r.CV,
                         Tav_th=Tsp + eps_eff * t["mean"], sd_th=eps_eff * t["sd"], CV_th=eps_eff * t["sd"] / (Tsp + eps_eff * t["mean"])))
R = pd.DataFrame(rows)
R["Tav_err_pct"] = 100 * (R.Tav_th - R.Tav_full).abs() / R.Tav_full; R["sd_err_pct"] = 100 * (R.sd_th - R.sd_full).abs() / R.sd_full; R["CV_err_pct"] = 100 * (R.CV_th - R.CV_full).abs() / R.CV_full
R.to_csv(os.path.join(OUT, "m1_fpt_reduction_check.csv"), index=False)
print(R.round(3).to_string(index=False))
# error summary by ν range (the approximation is a small-ν reduction; report where it holds and where it fails)
ERR = pd.DataFrame([dict(nu_range=lab, n=int(m.sum()), Tav_err_median=R[m].Tav_err_pct.median(), Tav_err_max=R[m].Tav_err_pct.max(), sd_err_median=R[m].sd_err_pct.median(), sd_err_max=R[m].sd_err_pct.max(), CV_err_max=R[m].CV_err_pct.max())
                    for lab, m in (("nu <= 0.5", R.nu <= 0.5), ("0.5 < nu <= 1", (R.nu > 0.5) & (R.nu <= 1)), ("nu > 1", R.nu > 1))])
ERR.to_csv(os.path.join(OUT, "m1_fpt_reduction_error_summary.csv"), index=False); print(ERR.round(2).to_string(index=False))
# effective vs nominal noise intensity at the reference (nu0 = sigma_H sqrt(eps_H) is what the run scripts set; nu = sigma_H sqrt(eps_eff) enters the reduced process)
print(f"nu0 = 0.822 (sigma_H sqrt(eps_H)) corresponds to nu = {0.822/np.sqrt(RATE):.3f} (sigma_H sqrt(eps_eff)) in the reduced process")

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

def moments(a):
    a.semilogx(th.nu, th["mean"] / th.mean_det, "-", color=C_TH, label=r"mean $\tau_{FP}$ / deterministic")
    a.semilogx(th.nu, th.CV, "--", color=C_TH, label=r"CV of $\tau_{FP}$")
    a.semilogx(mc.nu, mc["mean"] / mean_det, "o", color=C_MC, ms=11, mfc="white", mew=2.5, label="Monte-Carlo (1-D SDE)")
    a.semilogx(mc.nu, mc.CV, "s", color=C_MC, ms=11, mfc="white", mew=2.5)
    Rs = R[R.sweep == "sigH sweep"]
    a.semilogx(Rs.nu, Rs.CV_full, "^", color=C_FULL, ms=10, label="full model, CV of ISI")
    a.semilogx(Rs.nu, Rs.CV_th, ":", color=C_FULL, lw=3, label=r"$\varepsilon s/(T_{sp}+\varepsilon m)$")
    a.axhspan(DATA_CV[0], DATA_CV[1], color=TINT, lw=0, zorder=0, label="data CV (IQR)")
    a.set_xlim(0.02, 3); a.set_ylim(0, 1.2); a.set_xlabel(r"$\nu = \sigma_H\sqrt{\varepsilon_{\mathrm{eff}}}$"); a.set_ylabel("value")
    a.set_title(r"First-passage moments, $H_0=0.06$, $H^*=0.41$"); a.legend(loc="center left", fontsize=14)

def cv_vs_Hs(a):
    for nu, col in ((0.1, "#1f77b4"), (0.3, "#2ca02c"), (1.0, "#e69f00")):
        d = load(f"m1_fpt_vs_Hs_nu{nu}.csv"); a.plot(d.Hs, d.CV, "-", color=col, label=rf"$\nu={nu}$")
    a.axhspan(DATA_CV[0], DATA_CV[1], color=TINT, lw=0, zorder=0, label="data CV (IQR)")
    a.axvline(HS, color="k", ls="-.", lw=1.8); a.text(HS + 0.01, 1.05, r"$H^*=H_{\mathrm{on}}$", fontsize=17)
    a.set_xlim(0.15, 0.85); a.set_ylim(0, 1.2); a.set_xlabel(r"threshold $H^*$"); a.set_ylabel(r"CV of $\tau_{FP}$")
    a.set_title(r"CV depends on $(\nu, H^*)$, not on $\varepsilon_H$"); a.legend(loc="upper right", fontsize=15)

def reduction(a):
    for name, mk in (("sigH sweep", "o"), ("epsH sweep", "s")):
        d = R[R.sweep == name]
        a.loglog(d.Tav_th, d.Tav_full, mk, color=C_FULL, ms=10, mfc="white", mew=2.5, label=rf"$T_{{\mathrm{{av}}}}$, {name}")
        a.loglog(d.sd_th, d.sd_full, mk, color=C_MC, ms=10, label=rf"$\sigma_{{\mathrm{{ISI}}}}$, {name}")
    xx = np.array([1, 2000]); a.loglog(xx, xx, "k--", lw=2)
    a.set_xlim(1, 2000); a.set_ylim(1, 2000); a.set_xlabel("theory (model units)"); a.set_ylabel("full 3-variable SDE")
    a.set_title("Reduction check"); a.legend(loc="upper left", fontsize=14)

def reduction_error(a):
    """relative error of the reduced-process prediction vs the full SDE, as a function of nu (both sweeps)"""
    for name, mk in (("sigH sweep", "o"), ("epsH sweep", "s")):
        d = R[R.sweep == name].sort_values("nu")
        a.semilogx(d.nu, d.Tav_err_pct, mk + "-", color=C_FULL, ms=9, mfc="white", mew=2, label=rf"$T_{{\mathrm{{av}}}}$, {name}")
        a.semilogx(d.nu, d.sd_err_pct, mk + "-", color=C_MC, ms=9, label=rf"$\sigma_{{\mathrm{{ISI}}}}$, {name}")
    a.axhline(5, color="k", ls=":", lw=2); a.axvspan(0.02, 0.5, color=TINT, lw=0, zorder=0); a.text(0.025, 32, "data range", fontsize=16)
    a.set_xlim(0.02, 3); a.set_ylim(0, 35); a.set_xlabel(r"$\nu = \sigma_H\sqrt{\varepsilon_{\mathrm{eff}}}$"); a.set_ylabel("relative error (%)")
    a.set_title("Where the reduction holds"); a.legend(loc="upper left", fontsize=14)

# --- figures ----------------------------------------------------------------------------------------
MAIN = [moments, reduction, reduction_error]     # closed moments vs Monte-Carlo and vs the full SDE; error vs nu
SUPP = [cv_vs_Hs]               # the (nu, H*) map
compose(MAIN, 3, "m1_first_passage_theory_main")
compose(SUPP, 1, "m1_first_passage_theory_supp", w=8.5)
