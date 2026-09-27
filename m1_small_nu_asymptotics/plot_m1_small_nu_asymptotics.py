"""
plot_m1_small_nu_asymptotics.py — small-ν expansion of the first-passage moments vs quadrature and Monte-Carlo.

Closed forms (derive_small_nu.py; a = ν²/2, z0 = 1/(1−H0), zs = 1/(1−H*), T0 = ln(zs/z0), I1 = zs − z0 − T0,
K3 = zs² − z0² − 4(zs − z0) + 2T0):
   E[τ]   = T0 − a I1 + O(a²)          Var τ = 2a I1 + O(a²)          κ3 = 6a² K3 + O(a³)
   => sd ≈ ν √I1,  E[τ] ≈ T0 − Var/2,  skew ≈ (3/2) ν K3 / I1^{3/2},  skew/CV → (3/2) K3 T0 / I1²  (ν → 0)
Main  (a) mean shift T0 − E[τ] and sd vs ν: quadrature (points), MC (crosses), leading orders (lines)
      (b) skewness vs ν: quadrature, MC, leading order; gamma (2 CV) and inverse-Gaussian (3 CV) references
Supp  (a) grid convergence of skew/ν at small ν; (b) FPT density at ν = 0.3 (MC) vs gamma / inverse Gaussian
      with the same mean and variance; (c) the limit ratio skew/CV vs H* for three H0, with the data value
Writes results/m1_asym_summary.csv (leading-order coefficients and their numerical checks)
"""
import os, glob
import numpy as np, pandas as pd
from scipy.stats import gamma as gamma_dist, invgauss
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
NAME, OUT, FIG = "m1_small_nu_asymptotics", "results", "figures"; os.makedirs(FIG, exist_ok=True)
C_Q, C_MC, C_AS, C_MEAN, C_SD, C_GAM, C_IG, C_DATA = "#1f77b4", "k", "#c0392b", "#1f77b4", "#2ca02c", "#e69f00", "#9467bd", "#c0392b"
H0, HS = 0.064, 0.408
DATA_SKEW_OVER_CV = 0.52 / 0.22                # Azzoni trains: skew 0.52, CV 0.22 (m1_isi_variability)
ts = "../m1_isi_variability/results/m1_noise_location_table.csv"
if os.path.exists(ts):
    t = pd.read_csv(ts, index_col=0); r = t.loc["data (median)"]; DATA_SKEW_OVER_CV = float(r["skew"]) / float(r["CV"])

Q = pd.read_csv(os.path.join(OUT, "m1_asym_vs_nu.csv")); MC = pd.read_csv(os.path.join(OUT, "m1_asym_mc.csv"))
CV = pd.read_csv(os.path.join(OUT, "m1_asym_convergence.csv"))
HIST = pd.read_csv(os.path.join(OUT, "m1_asym_hist_nu0.3.csv"))
RAT = {os.path.basename(f).split("_H0")[1][:-4]: pd.read_csv(f) for f in sorted(glob.glob(os.path.join(OUT, "m1_asym_ratio_vs_Hs_H0*.csv")))}

z0, zs = 1 / (1 - H0), 1 / (1 - HS); T0 = np.log(zs / z0); I1 = zs - z0 - T0; K3 = zs**2 - z0**2 - 4 * (zs - z0) + 2 * T0
SD_C, SK_C, RATIO = np.sqrt(I1), 1.5 * K3 / I1**1.5, 1.5 * K3 * T0 / I1**2
small = Q[Q.nu <= 0.1]
summ = pd.DataFrame([dict(H0=H0, Hs=HS, T0=T0, I1=I1, K3=K3, sd_over_nu=SD_C, mean_shift_over_nu2=I1 / 2, skew_over_nu=SK_C, skew_over_CV_limit=RATIO,
                          quad_sd_over_nu_small=float((small.sd / small.nu).mean()), quad_shift_over_nu2_small=float(((small.mean_det - small["mean"]) / small.nu**2).mean()),
                          quad_skew_over_nu_small=float((small["skew"] / small.nu).mean()), data_skew_over_CV=DATA_SKEW_OVER_CV)])
summ.to_csv(os.path.join(OUT, "m1_asym_summary.csv"), index=False); print(summ.T.to_string())

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

def p_mean_sd(a):
    nu = np.logspace(np.log10(Q.nu.min()), np.log10(Q.nu.max()), 100)
    a.loglog(Q.nu, Q.mean_det - Q["mean"], "o", color=C_MEAN, ms=8, label=r"$T_0 - E[\tau]$ (quadrature)")
    a.loglog(Q.nu, Q.sd, "s", color=C_SD, ms=8, label=r"sd $\tau$ (quadrature)")
    a.loglog(MC.nu, T0 - MC["mean"], "x", color=C_MC, ms=12, mew=2.5, label="Monte-Carlo")
    a.loglog(MC.nu, MC.sd, "x", color=C_MC, ms=12, mew=2.5)
    a.loglog(nu, 0.5 * I1 * nu**2, "-", color=C_MEAN, lw=2, label=r"$\nu^2 I_1/2$"); a.loglog(nu, SD_C * nu, "-", color=C_SD, lw=2, label=r"$\nu\sqrt{I_1}$")
    a.set_xlabel(r"$\nu$"); a.set_ylabel("rescaled time"); a.set_title("Mean shift and sd"); a.legend(loc="lower right", fontsize=15)

def p_skew(a):
    nu = np.logspace(np.log10(Q.nu.min()), np.log10(Q.nu.max()), 100)
    q = Q[Q.nu >= 0.03]                                    # below 0.03 the third cumulant needs n_y > 16000 (see p_convergence)
    a.loglog(q.nu, q["skew"], "o", color=C_Q, ms=8, label=r"quadrature ($\nu\geq0.03$)")
    a.loglog(MC.nu, MC["skew"], "x", color=C_MC, ms=12, mew=2.5, label="Monte-Carlo")
    a.loglog(nu, SK_C * nu, "-", color=C_AS, lw=2, label=r"$\frac{3}{2}\nu K_3/I_1^{3/2}$")
    a.loglog(Q.nu, 2 * Q.CV, "--", color=C_GAM, lw=2, label=r"gamma: $2\,$CV"); a.loglog(Q.nu, 3 * Q.CV, ":", color=C_IG, lw=2.5, label=r"inv. Gaussian: $3\,$CV")
    a.set_xlabel(r"$\nu$"); a.set_ylabel("skewness"); a.set_title("Skewness"); a.legend(loc="lower right", fontsize=15)

def p_convergence(a):
    for nu, g in CV.groupby("nu"):
        a.semilogx(g.ny, g["skew"] / g.nu, "o-", ms=9, label=rf"$\nu={nu:.2g}$")
    a.axhline(SK_C, color=C_AS, lw=2, ls="--", label="leading order")
    a.set_xlabel(r"$y$-grid points $n_y$ ($n_w = n_y/2$)"); a.set_ylabel(r"skew$/\nu$"); a.set_title("Quadrature convergence"); a.legend(loc="upper right", fontsize=15)

def p_density(a):
    r = MC[np.isclose(MC.nu, 0.3)].iloc[0]; m, v = r["mean"], r.sd**2
    t = np.linspace(HIST.t.min(), HIST.t.max(), 400)
    a.bar(HIST.t, HIST.density, width=HIST.t.iloc[1] - HIST.t.iloc[0], color="#cfe3f3", edgecolor="none", label=r"MC, $\nu=0.3$")
    a.plot(t, gamma_dist.pdf(t, m * m / v, scale=v / m), "--", color=C_GAM, lw=2.5, label="gamma (same mean, var)")
    lam = m**3 / v; a.plot(t, invgauss.pdf(t, m / lam, scale=lam), ":", color=C_IG, lw=3, label="inverse Gaussian")
    a.set_xlabel(r"$\tau$ (rescaled time)"); a.set_ylabel("density"); a.set_title(rf"FPT density (CV $= {r.CV:.2f}$)"); a.legend(loc="upper right", fontsize=15)

def p_ratio(a):
    cols = {"0.02": "#9ecae1", "0.064": C_Q, "0.15": "#08306b"}
    for h0, d in RAT.items(): a.plot(d.Hs, d.skew_over_CV_limit, "-", color=cols.get(h0, "k"), label=rf"$H_0={h0}$")
    a.axhline(DATA_SKEW_OVER_CV, color=C_DATA, lw=2, ls="--", label="data (train medians)")
    a.axhline(2, color=C_GAM, lw=1.5, ls=":"); a.axhline(3, color=C_IG, lw=1.5, ls=":")
    a.plot([HS], [RATIO], "o", color=C_AS, ms=12, label=rf"$H^*={HS}$: {RATIO:.2f}")
    a.set_ylim(0, 8); a.set_xlabel(r"$H^*$"); a.set_ylabel(r"$\lim_{\nu\to0}$ skew$/$CV"); a.set_title("Shape ratio"); a.legend(loc="upper right", fontsize=14)

MAIN = [p_mean_sd, p_skew]                 # the expansion is exact to leading order; skew is O(nu) with a computable coefficient
SUPP = [p_convergence, p_density, p_ratio]
if __name__ == "__main__":
    compose(MAIN, 2, f"{NAME}_main")
    compose(SUPP, 3, f"{NAME}_supp")
