"""
plot_revision_figures.py — supplementary figures added in revision (S27-S30), from the result files of this folder.
  figR1_theory_checks        skew/CV vs noise for two forms and the LIF model; bounds of Theorem 5.2; boundary check
  figR2_interval_densities   likelihood differences, offset profiles, fitted noise by dose
  figR3_law_components       components of the sigma-T_av law, CV by dose, calibration of the estimators
  figR4_dynamic_range        period curves for n = 1, 2, 3; dose response; spike shape
Output: figures/figR*.{png,pdf,eps}
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 22, "axes.titlesize": 24, "axes.labelsize": 24, "xtick.labelsize": 20, "ytick.labelsize": 20, "legend.fontsize": 16,
                     "axes.grid": True, "grid.color": "#c8c8c8", "grid.linestyle": "--", "grid.linewidth": 0.8, "lines.linewidth": 2.5, "lines.markersize": 8,
                     "mathtext.fontset": "stix", "font.family": "STIXGeneral", "savefig.dpi": 200, "legend.framealpha": 1.0})
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "results"); FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
L = lambda n: pd.read_csv(os.path.join(OUT, n))
BLUE, GREEN, RED, ORANGE, PURPLE, GREY = "#1f77b4", "#2ca02c", "#c0392b", "#e69f00", "#9467bd", "#666666"
TINTS = {BLUE: "#c6dbef", GREEN: "#c7e9c0", RED: "#f5d0cc", PURPLE: "#dcd0ea", GREY: "#d9d9d9"}     # opaque tints (EPS has no transparency)
def letter(a, s): a.text(-0.13, 1.06, s, transform=a.transAxes, fontsize=28, fontweight="bold", va="bottom", ha="left")
def compose(panels, ncol, name, w=8.0, h=7.5):
    n = len(panels); nrow = -(-n // ncol); fig, axs = plt.subplots(nrow, ncol, figsize=(w * ncol, h * nrow), squeeze=False)
    for k, fn in enumerate(panels): a = axs[k // ncol, k % ncol]; fn(a); letter(a, f"({chr(97 + k)})")
    for k in range(n, nrow * ncol): axs[k // ncol, k % ncol].axis("off")
    plt.tight_layout(w_pad=2.5, h_pad=2.5)
    for ext in ("png", "pdf", "eps"): fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig); print("wrote", name)

# ------------------------------------------------------------------ R1 theory checks
X = L("a2_exact_vs_expansion.csv"); LIF = L("i1_lif_check.csv"); HB = L("h1_remainder_bounds.csv"); E1 = L("e1_reflecting_boundary.csv"); E1b = L("e1b_boundary_exact_step.csv")
def p_ratio(a):
    for form, col, lab in (("jacobi", BLUE, r"Jacobi $H(1-H)$"), ("cle", GREEN, r"chemical Langevin")):
        d = X[(X.form == form) & (X.CV < 0.75)]; a.plot(d.CV, d.R, "o-", color=col, label=lab + ", exact"); a.plot(d.CV, d.R_next.where(d.R_next > 2.4), ":", color=col, lw=2)
        a.axhline(d.R_limit.iloc[0], color=col, ls="--", lw=1.5)
    for mu, mk in ((1.5, "s"), (2.0, "^"), (5.0, "v")):
        d = LIF[LIF.mu == mu]; a.plot(d.CV, d.skew_over_CV, mk + "-", color=GREY, mfc="white", mew=2, lw=1.5, label=rf"LIF, $\mu={mu:g}$")
    a.axhline(3, color=RED, lw=2); a.text(0.02, 2.9, "inverse Gaussian", color=RED, fontsize=16); a.axhline(2, color=ORANGE, lw=2); a.text(0.02, 2.04, "gamma", color=ORANGE, fontsize=16)
    a.set_xlim(0, 0.72); a.set_ylim(1.9, 4.3); a.set_xlabel(r"CV of $\tau$"); a.set_ylabel("skewness / CV"); a.set_title("Ratio against the noise level"); a.legend(loc="upper right", fontsize=13)
def p_bounds(a):
    d = HB[np.isclose(HB.H0, 0.064) & (HB.nu <= 0.85)]
    for q, col, lab in (("mean", GREY, "mean"), ("var", BLUE, "variance"), ("k3", GREEN, r"$\kappa_3$")):
        lead = {"mean": d["mean_lo"], "var": d["var_hi"], "k3": d["k3_hi"]}[q]
        a.fill_between(d.nu, d[q + "_lo"] / lead, d[q + "_hi"] / lead, color=TINTS[col], lw=0, zorder=1); a.plot(d.nu, d[q] / lead, "o-", color=col, label=lab, zorder=3)
    a.axvline(0.281, color="k", ls="-.", lw=1.5); a.text(0.29, 1.12, r"fitted $\nu$", fontsize=16)
    a.set_xscale("log"); a.set_xlim(0.02, 0.85); a.set_ylim(0, 1.2); a.set_xlabel(r"$\nu$"); a.set_ylabel("exact / leading term"); a.set_title("Bounds of Theorem 5.2(b)"); a.legend(loc="lower left")
def p_boundary(a):
    a.plot(E1.nu, 100 * E1.err_mean, "o-", color=RED, label="Euler–Maruyama, reflection"); a.plot(E1b.nu, 100 * E1b.err_mean, "s-", color=BLUE, label="exact square-root step")
    a.axhline(0, color="k", lw=1); a.axvline(np.sqrt(2), color="k", ls="-.", lw=1.5); a.text(1.45, -20, r"$\nu^2=2$", fontsize=16)
    a.set_xlabel(r"$\nu$"); a.set_ylabel(r"error of $\mathrm{E}[\tau]$ against (12)  (%)"); a.set_ylim(-25, 5); a.set_title("Boundary 0: two schemes"); a.legend(loc="lower left")
compose([p_ratio, p_bounds, p_boundary], 3, "figR1_theory_checks")

# ------------------------------------------------------------------ R2 interval densities
A1 = L("b2_A_primary.csv").set_index("family"); C1 = L("b3_contamination.csv"); C1 = C1[C1.eps == 0.01].set_index("family"); B = L("b2_B_offset_profile.csv"); P = L("b3_proportional_offset.csv"); SEG = L("b2_A_segments.csv"); DA = L("c2_dose_axis.csv")
ORDER = ["lognorm+off", "invgauss+off", "lognorm", "invgauss", "fpt_cle+off", "gamma+off", "gamma"]; NAME = {"lognorm+off": "log-normal+off", "invgauss+off": "inv. Gaussian+off", "lognorm": "log-normal", "invgauss": "inv. Gaussian", "fpt_cle+off": "FPT chem. Langevin+off", "gamma+off": "gamma+off", "gamma": "gamma"}
def p_dll(a):
    y = np.arange(len(ORDER))[::-1]
    for T, dy, col, lab, mfc in ((A1, 0.15, BLUE, "plain likelihood", BLUE), (C1, -0.15, RED, "with 1 % outlier component", "white")):
        v = T.loc[ORDER]; a.errorbar(v.dlogL_vs_ref, y + dy, xerr=[v.dlogL_vs_ref - v.seg_lo, v.seg_hi - v.dlogL_vs_ref], fmt="o", color=col, mfc=mfc, mew=2, capsize=5, lw=2, label=lab)
    a.axvline(0, color="k", lw=1.5); a.set_yticks(y); a.set_yticklabels([NAME[o] for o in ORDER], fontsize=17); a.set_xlabel(r"$\log L-\log L$(FPT Jacobi+off)"); a.set_title("Likelihood of the 196 segments"); a.legend(loc="lower left", fontsize=14); a.set_xlim(-300, 200)
FAM = [("gamma", ORANGE, "gamma"), ("invgauss", RED, "inverse Gaussian"), ("lognorm", PURPLE, "log-normal"), ("fpt_jacobi", BLUE, "FPT Jacobi"), ("fpt_cle", GREEN, "FPT chem. Langevin")]
def p_off(a):
    ref = B.set_index("delta").loc[0.0]
    for f, col, lab in FAM: a.plot(B.delta, B[f] - ref[f], "o-", color=col, label=lab)
    a.axvline(8.64, color="k", ls="-.", lw=1.5); a.text(8.45, -19, r"$T_{\mathrm{sp}}$ of the model", fontsize=15, ha="right")
    a.set_xlabel("common offset (s)"); a.set_ylabel(r"$\log L$ gain over no offset"); a.set_title("Constant offset"); a.legend(loc="upper left", fontsize=13)
def p_prop(a):
    Pp = P[P.phi <= 0.2]; ref = Pp.set_index("phi").loc[0.0]
    for f, col, lab in FAM: a.plot(Pp.phi, Pp[f] - ref[f], "o-", color=col, label=lab)
    a.axvspan(0.17, 0.22, color="#f5d0cc", lw=0, zorder=0); a.text(0.211, 30, r"model $T_{\mathrm{sp}}/T_{\mathrm{det}}$", fontsize=15, rotation=90, va="bottom")
    a.set_xlim(-0.005, 0.225); a.set_xlabel("offset / segment mean"); a.set_ylabel(r"$\log L$ gain over no offset"); a.set_title("Proportional offset")
def p_nu(a):
    doses = [3, 5, 10, 15, 30]; g = SEG.groupby("dose").nu_jacobi; med, lo, hi = g.median(), g.quantile(.25), g.quantile(.75)
    a.errorbar(doses, med.loc[doses], yerr=[med.loc[doses] - lo.loc[doses], hi.loc[doses] - med.loc[doses]], fmt="o", color=BLUE, capsize=6, ms=11, label="segment fits (median, IQR)")
    Sm = np.array([0.064, 0.11, 0.21, 0.32, 0.64]); nu_m = np.interp(np.log(Sm), np.log(DA.S), DA.nu); a.plot(doses, nu_m * med.loc[15] / nu_m[3], "s--", color=RED, mfc="white", mew=2, label=r"model, fixed $\sigma_H$ (scaled at 15 µM)")
    a.axhline(med.loc[15], color=GREEN, ls=":", lw=2.5, label="fixed number of units $N$")
    a.set_xscale("log"); a.set_xticks(doses); a.set_xticklabels(doses); a.set_xlabel("carbachol (µM)"); a.set_ylabel(r"fitted $\nu$ (Jacobi form)"); a.set_ylim(0, 0.55); a.set_title("Noise intensity by dose"); a.legend(loc="lower left", fontsize=14)
compose([p_dll, p_off, p_prop, p_nu], 2, "figR2_interval_densities", w=9.0)

# ------------------------------------------------------------------ R3 components of the law
SL = L("c1_slopes.csv"); SL = SL[SL.regressor == "log T_av"].set_index("component"); CAL = L("c1_calibration_raw.csv"); F2 = L("f2_by_dose.csv"); A4 = L("a4_slope_sweeps.csv")
def p_comp(a):
    comps = ["pooled", "xr", "xc", "xs", "xw"]; labs = ["pooled", "between\nrecords", "between cells\nof a record", "same cell,\nother dose", "same cell\nand dose"]; x = np.arange(len(comps)); v = SL.loc[comps]
    a.errorbar(x, v.slope, yerr=[v.slope - v.rec_lo, v.rec_hi - v.slope], fmt="o", color="k", capsize=7, ms=12, lw=2.5, label="data (record bootstrap)")
    for xi, val, col, mk, lab in ((1, 1.61, RED, "s", r"dose axis, fixed $\sigma_H$"), (3, 1.61, RED, "s", None), (1, 1.09, GREEN, "D", r"dose axis, fixed $N$"), (3, 1.09, GREEN, "D", None), (2, 1.02, BLUE, "^", r"cell axis, common $\nu$"), (2, 1.49, ORANGE, "v", r"cell axis, common $\sigma_H$")):
        a.plot(xi + 0.22, val, mk, color=col, ms=13, mfc="white", mew=2.5, label=lab)
    a.set_xticks(x); a.set_xticklabels(labs, fontsize=15); a.set_ylim(0.8, 2.25); a.set_ylabel(r"slope of $\log\sigma$ on the component of $\log T_{\mathrm{av}}$"); a.set_title(r"Components of the $\sigma$–$T_{\mathrm{av}}$ law"); a.legend(loc="upper center", fontsize=13, ncol=2)
def p_cvdose(a):
    doses = [3, 5, 10, 15, 30]; d1 = F2[F2.variant.str.startswith("n=1")]; d2 = F2[F2.variant.str.startswith("n=2")]
    a.plot(doses, d1.data_CV, "o-", color="k", ms=12, label="data (segment median)"); a.plot(doses, d1.CV, "s--", color=BLUE, mfc="white", mew=2, label=r"model $n=1$, fixed $\sigma_H$"); a.plot(doses, d2.CV, "^--", color=RED, mfc="white", mew=2, label=r"model $n=2$, fixed $\sigma_H$")
    c = A4[(A4.case.str.startswith("chemical")) & (A4.axis == "dose")]; cv = np.interp(np.log([0.064, 0.11, 0.21, 0.32, 0.64]), np.log(c.value), c.CV); a.plot(doses, cv, "D:", color=GREEN, mfc="white", mew=2, label=r"chemical Langevin, fixed $N$ (unsegmented)")
    a.set_xscale("log"); a.set_xticks(doses); a.set_xticklabels(doses); a.set_ylim(0, 0.5); a.set_xlabel("carbachol (µM)"); a.set_ylabel("CV of the ISI"); a.set_title("CV by dose (tuned at 15 µM only)"); a.legend(loc="upper right", fontsize=13)
def p_cal(a):
    comps = [("pooled", "k", "o"), ("xr", RED, "s"), ("xc", BLUE, "^"), ("xs", GREEN, "D"), ("xw", ORANGE, "v")]; names = {"pooled": "pooled", "xr": "between records", "xc": "between cells", "xs": "same cell, other dose", "xw": "same cell and dose"}
    for k, (c, col, mk) in enumerate(comps):
        g = CAL.groupby("s_true")[c]; a.errorbar(g.mean().index + (k - 2) * 0.012, g.mean().values, yerr=g.std().values, fmt=mk, color=col, capsize=5, ms=10, mfc="white", mew=2, label=names[c])
    a.plot([0.95, 1.55], [0.95, 1.55], "k--", lw=1.5); a.set_xlabel("true exponent"); a.set_ylabel("estimated slope (mean ± sd, 300 data sets)"); a.set_title("Calibration on synthetic data"); a.legend(loc="upper left", fontsize=13)
compose([p_comp, p_cvdose, p_cal], 3, "figR3_law_components")

# ------------------------------------------------------------------ R4 dynamic range
CU = L("f1_curves.csv"); SM = L("f1_summary.csv"); DM = L("f2_dose_map.csv"); DT = pd.read_csv(os.path.join(HERE, "..", "data_dose_units", "results", "dose_table.csv"))
VAR = [((1, 10000, 30), BLUE, r"$n=1$, $d_H=10^4$ (manuscript)"), ((2, 2000, 30), RED, r"$n=2$, $d_H=2000$"), ((3, 2000, 30), PURPLE, r"$n=3$, $d_H=2000$")]
sel = lambda T, k: T[(T.nP == k[0]) & (T.dH == k[1]) & (T.gH == k[2])]
def p_period(a):
    for k, col, lab in VAR:
        c = sel(CU, k); c = c[c.Tdet > 0]; s = sel(SM, k).iloc[0]; a.loglog(c.S, c.Tdet, "-", color=col, label=lab + rf": {s.ratio_per_decade:.1f} per decade"); w = c[(c.S >= s.S_lo * 0.999) & (c.S <= s.S_hi * 1.001)]; a.loglog(w.S, w.Tdet, "-", color=TINTS[col], lw=11, zorder=1, solid_capstyle="butt")
    a.set_xlabel(r"stimulus $\bar S$"); a.set_ylabel("period (model units)"); a.set_title("Period curves; thick: steepest decade"); a.legend(loc="upper right", fontsize=13)
def p_dose(a):
    doses = np.array([3, 5, 10, 15, 30.0]); med, lo, hi = [], [], []
    for d_ in doses: r = DT[DT.dose == d_].iloc[0]; med.append(r.isi_train_median); lo.append(r.isi_train_ci_lo); hi.append(r.isi_train_ci_hi)
    med, lo, hi = map(np.array, (med, lo, hi)); a.errorbar(doses, med, yerr=[med - lo, hi - med], fmt="o", color="k", capsize=7, ms=12, label="data (train medians)")
    for v, col, lab in (("n=1, dH=1e4 (manuscript)", BLUE, r"$n=1$"), ("n=2, dH=2000", RED, r"$n=2$")):
        for mp_, ls, mk in (("end points, linear", "-", "s"), ("least squares, linear", ":", "^")):
            r = DM[(DM.variant == v) & (DM["map"] == mp_)].iloc[0]; a.plot(doses, [r[f"ISI_{int(d_)}"] for d_ in doses], mk + ls, color=col, mfc="white", mew=2, label=f"{lab}, {mp_.split(',')[0]}")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xticks(doses); a.set_xticklabels([int(d_) for d_ in doses]); a.set_xlabel("carbachol (µM)"); a.set_ylabel("ISI (s)"); a.set_title("Dose response"); a.legend(loc="lower left", fontsize=13)
def p_shape(a):
    for k, col, lab in VAR[:2]:
        c = sel(CU, k); c = c[c.Tdet > 0]; s = sel(SM, k).iloc[0]; w = c[(c.S >= s.S_lo * 0.999) & (c.S <= s.S_hi * 1.001)]
        a.semilogx(w.S, w.W50 / w.W50.iloc[0], "-", color=col, label=lab.split(",")[0] + ", half-width"); a.semilogx(w.S, w.amp / w.amp.iloc[0], "--", color=col, label=lab.split(",")[0] + ", amplitude")
    a.axhline(1.10, color="k", ls=":", lw=2); a.axhline(1.26, color="k", ls="-.", lw=2); a.text(0.07, 1.105, "data: width, 5→30 µM within cells", fontsize=14); a.text(0.07, 1.265, "data: amplitude", fontsize=14)
    a.set_ylim(0.95, 1.4); a.set_xlabel(r"stimulus $\bar S$ (steepest decade)"); a.set_ylabel("relative to the low end"); a.set_title("Spike shape across the decade"); a.legend(loc="upper left", fontsize=13)
compose([p_period, p_dose, p_shape], 3, "figR4_dynamic_range")
