"""
plot_m1_particle_filter_synthetic.py — gate 2: synthetic record, filtered hidden state, profile likelihoods.

Main figure  (figures/m1_particle_filter_synthetic_main.*, 1x2)
  (a) hidden H(t) (truth) vs filtered H — is the recovery variable recovered from F alone?
  (b) profile log-likelihood of eps_H and beta (relative to the maximum; χ²(1) 95 % drop = 1.92 marked)
Supplement   (figures/m1_particle_filter_synthetic_supp.*, 1x2)
  (a) synthetic F(t) (first 600 model units) with the PF-filtered C (scaled to F) — does the filter track the spikes?
  (b) profile log-likelihood of sigma_H (weakly identified from one 24-spike record)
Run after run_pf_synthetic.py .  Also writes results/pf_profile_summary.csv
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
NAME, OUT, FIG = "m1_particle_filter_synthetic", "results", "figures"; os.makedirs(FIG, exist_ok=True)
C_OBS, C_PF, C_TRUE, C_EPS, C_BETA, C_SIG = "#888888", "#1f77b4", "k", "#2ca02c", "#e69f00", "#9467bd"
TRUTH = dict(epsH=1000.0, sigH=0.026, beta=4.0)
LAB = {"epsH": r"$\varepsilon_H$", "sigH": r"$\sigma_H$", "beta": r"$\beta$"}
truth = pd.read_csv(os.path.join(OUT, "synthetic_truth.csv")); filt = pd.read_csv(os.path.join(OUT, "pf_filtered_truth.csv"))
prof = {p: pd.read_csv(os.path.join(OUT, f"pf_profile_{p}.csv")) for p in TRUTH}
m = truth.t <= 600

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

# --- profile summary (computed once, independent of the figures) ---------------------------------------------------------
rows = []
PROF = {}
for p in TRUTH:
    g = prof[p].groupby("value").loglik.agg(["mean", "std"]).reset_index(); rel = g["mean"] - g["mean"].max(); PROF[p] = (g, rel)
    inside = g.value[rel >= -1.92]; rows.append(dict(par=p, truth=TRUTH[p], mle=float(g.value[rel.idxmax()]), ci_lo=float(inside.min()), ci_hi=float(inside.max())))
S = pd.DataFrame(rows); S.to_csv(os.path.join(OUT, "pf_profile_summary.csv"), index=False); print(S.to_string(index=False))

def profile(a, p, col, label):
    g, rel = PROF[p]
    a.errorbar(g.value / TRUTH[p], rel, yerr=g["std"].fillna(0), fmt="o-", color=col, ms=9, capsize=5, label=label)

# --- panels ---------------------------------------------------------------------------------
def p_record(a):
    a.plot(truth.t[m], truth.F[m], ".", color=C_OBS, ms=5, label=r"observed $F$")
    a.plot(filt.t[m], 0.0 + 0.11 * filt.Cf[m], "-", color=C_PF, lw=2, label=r"PF filtered $b + a\,C$")
    a.set_xlabel(r"$t$ (model units)"); a.set_ylabel(r"$F$"); a.set_title("Synthetic record and filter"); a.legend(loc="upper center", fontsize=14)

def p_hidden_H(a):
    a.plot(truth.t[m], truth.H[m], "-", color=C_TRUE, lw=2, label=r"true $H$")
    a.plot(filt.t[m], filt.Hf[m], "-", color=C_PF, lw=2, label=r"filtered $H$")
    a.set_ylim(0, 0.7); a.set_xlabel(r"$t$ (model units)"); a.set_ylabel(r"$H$"); a.set_title("Hidden recovery variable"); a.legend(loc="upper right", fontsize=16)

def p_profiles_eps_beta(a):
    profile(a, "epsH", C_EPS, LAB["epsH"]); profile(a, "beta", C_BETA, LAB["beta"])
    a.axhline(-1.92, color="k", ls="--", lw=2, label="$\\chi^2_1$ 95 % drop\n(reference scale only)"); a.axvline(1.0, color="k", ls=":", lw=1.5)
    a.set_xscale("log"); a.set_yscale("symlog", linthresh=10); a.set_ylim(-2e4, 10); a.set_xticks([0.5, 0.7, 1, 1.4, 2]); a.set_xticklabels(["0.5", "0.7", "1", "1.4", "2"]); a.minorticks_off()
    a.set_xlabel("parameter / truth"); a.set_ylabel(r"$\log L - \max \log L$ (symlog)"); a.set_title(r"Conditional likelihood sections: $\varepsilon_H$, $\beta$"); a.legend(loc="upper right", fontsize=14)

def p_profile_sigma(a):
    profile(a, "sigH", C_SIG, LAB["sigH"])
    a.axhline(-1.92, color="k", ls="--", lw=2, label="$\\chi^2_1$ 95 % drop\n(reference scale only)"); a.axvline(1.0, color="k", ls=":", lw=1.5)
    a.set_xscale("log"); a.set_ylim(-60, 5); a.set_xticks([0.5, 0.7, 1, 1.4, 2]); a.set_xticklabels(["0.5", "0.7", "1", "1.4", "2"]); a.minorticks_off()
    a.set_xlabel("parameter / truth"); a.set_ylabel(r"$\log L - \max \log L$"); a.set_title(r"Conditional section: $\sigma_H$ (weak)"); a.legend(loc="lower right", fontsize=16)

# --- figures --------------------------------------------------------------------------------
MAIN = [p_hidden_H, p_profiles_eps_beta]     # the hidden H is recovered; eps_H and beta are sharply identified (gate 2)
SUPP = [p_record, p_profile_sigma]           # the record itself and the weak sigma_H profile
compose(MAIN, 2, f"{NAME}_main")
compose(SUPP, 2, f"{NAME}_supp")
