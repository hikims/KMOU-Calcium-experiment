"""
plot_m0ip3h_bifurcation_structure.py — where the stereotyped reference oscillates, how the cycle is born, and
how far it sits from the Hopf boundary (relaxation regime).

Inputs (run_m0ip3h_bifurcation.sh): results/bif_grid_S_epsH_ref_st.csv, bif_grid_beta_qoff_S1.csv,
bif_line_S_ref_{st,prev}.csv, bif_line_epsH_ref_st.csv, bif_line_beta_ref_st.csv
Columns: x,y,osc,T,amp,W50,duty,Cmin,Cmax,Hmin,Hmax,nfp,C_fp,ReL,ImL,n_unstable,class,hopf_ratio,nclip
Main  (a) (S̄, ε_H): period of the cycle (colour) with the Hopf boundary Re λ = 0 of the fixed point and the
          oscillation boundary; the reference and the dose range S̄(3–30 µM) marked
      (b) (β, q_off): duty W50/T of the cycle, oscillation boundary; reference marked
      (c) 1-D in S̄: period and amplitude for REF_ST and the previous reference; fixed-point stability shaded
Supp  (a) Re λ map; (b) T·Im λ/2π (Hopf ratio) map; (c) amplitude map (β, q_off); (d) ε_H line; (e) β line;
      (f) class map (steady / Hopf-born cycle / coexistence / multiple fixed points)
Writes results/bif_summary.csv (Hopf points on the lines, reference values)
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm

plt.rcParams.update({
    "font.size": 22, "axes.titlesize": 24, "axes.labelsize": 24,
    "xtick.labelsize": 20, "ytick.labelsize": 20, "legend.fontsize": 19,
    "axes.grid": True, "grid.color": "#c8c8c8", "grid.linestyle": "--", "grid.linewidth": 0.8,
    "lines.linewidth": 2.5, "lines.markersize": 6,
    "mathtext.fontset": "stix", "font.family": "STIXGeneral",
    "savefig.dpi": 300, "legend.framealpha": 1.0,
})
NAME, OUT, FIG = "m0ip3h_bifurcation_structure", "results", "figures"; os.makedirs(FIG, exist_ok=True)
C_ST, C_PREV, C_HOPF, C_OSC, C_REF, C_TINT = "#2ca02c", "#e69f00", "#c0392b", "k", "#c0392b", "#f5d0cc"
REF = dict(S=1.0, epsH=1000.0, beta=4.0, qoff=0.0)
S_DOSE = None
fn = "../data_dose_units/results/dose_map_refit.csv"
if os.path.exists(fn):
    r = pd.read_csv(fn).iloc[0]; S_DOSE = (float(r.S_3uM), float(r.S_30uM))

def load(n): return pd.read_csv(os.path.join(OUT, n))
G1 = load("bif_grid_S_epsH_ref_st.csv"); G2 = load("bif_grid_beta_qoff_S1.csv")
LS = load("bif_line_S_ref_st.csv"); LP = load("bif_line_S_ref_prev.csv"); LE = load("bif_line_epsH_ref_st.csv"); LB = load("bif_line_beta_ref_st.csv")

def pivot(G, col):
    P = G.pivot(index="y", columns="x", values=col); return P.columns.values, P.index.values, P.values
def hopf_points(L):
    """x where Re λ changes sign along a line (linear interpolation)"""
    x, r = L.x.values, L.ReL.values; pts = []
    for i in range(len(x) - 1):
        if r[i] * r[i + 1] < 0: pts.append(float(np.exp(np.interp(0, [r[i], r[i + 1]], [np.log(x[i]), np.log(x[i + 1])]) if r[i] < r[i + 1] else np.interp(0, [r[i + 1], r[i]], [np.log(x[i + 1]), np.log(x[i])]))))
    return pts
rows = []
for lab, L in (("S_ref_st", LS), ("S_ref_prev", LP), ("epsH_ref_st", LE), ("beta_ref_st", LB)):
    o = L[L.osc == 1]
    rows.append(dict(line=lab, hopf_points=";".join(f"{p:.4g}" for p in hopf_points(L)), osc_min=o.x.min() if len(o) else np.nan, osc_max=o.x.max() if len(o) else np.nan,
                     n_osc=int(len(o)), n_coexist=int((L["class"] % 10 == 2).sum()), n_multi_fp=int((L["class"] >= 10).sum()), hopf_ratio_max=float(o.hopf_ratio.max()) if len(o) else np.nan))
ref_row = G1.iloc[((np.log(G1.x / REF["S"]))**2 + (np.log(G1.y / REF["epsH"]))**2).idxmin()]
rows.append(dict(line="REF_ST grid point", hopf_points="", osc_min=ref_row.x, osc_max=ref_row.y, n_osc=int(ref_row.osc), n_coexist=int(ref_row["class"] % 10 == 2), n_multi_fp=int(ref_row["class"] >= 10), hopf_ratio_max=ref_row.hopf_ratio))
S = pd.DataFrame(rows); S.to_csv(os.path.join(OUT, "bif_summary.csv"), index=False); print(S.to_string(index=False))
print(f"grid 1: {int((G1.osc == 1).sum())}/{len(G1)} oscillate; classes: {G1['class'].value_counts().to_dict()}; clipped cells: {int((G1.nclip > 0).sum())}")
print(f"grid 2: {int((G2.osc == 1).sum())}/{len(G2)} oscillate; classes: {G2['class'].value_counts().to_dict()}")

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

def boundaries(a, G, logx=True, logy=True):
    x, y, osc = pivot(G, "osc"); _, _, rel = pivot(G, "ReL")
    X, Y = np.meshgrid(x, y)
    a.contour(X, Y, osc, levels=[0.5], colors=C_OSC, linewidths=2.5)
    a.contour(X, Y, rel, levels=[0.0], colors=C_HOPF, linewidths=2.5, linestyles="--")
    a.plot([], [], "-", color=C_OSC, label="oscillation\nboundary"); a.plot([], [], "--", color=C_HOPF, label="Hopf:\nRe$\\,\\lambda=0$")

def p_map_S_epsH(a):
    x, y, T = pivot(G1, "T"); T = np.where(pivot(G1, "osc")[2] == 1, T, np.nan)
    pc = a.pcolormesh(x, y, T, norm=LogNorm(vmin=np.nanmin(T), vmax=np.nanmax(T)), cmap="viridis", shading="nearest")
    plt.colorbar(pc, ax=a, pad=0.02).set_label(r"period $T$"); boundaries(a, G1)
    if S_DOSE: a.axvspan(S_DOSE[0], S_DOSE[1], color=C_TINT, alpha=0.35, zorder=0); a.plot([], [], "s", color=C_TINT, ms=12, label=r"$\bar S$(3–30 $\mu$M)")
    a.plot([REF["S"]], [REF["epsH"]], "*", color="w", mec="k", ms=22, label=r"REF$_{ST}$")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(r"$\bar S$"); a.set_ylabel(r"$\varepsilon_H$"); a.set_title("Oscillation region and period"); a.legend(loc="lower right", fontsize=14); a.grid(False)

def p_map_beta_qoff(a):
    x, y, D = pivot(G2, "duty"); D = np.where(pivot(G2, "osc")[2] == 1, D, np.nan)
    pc = a.pcolormesh(x, y, D, vmin=0, vmax=min(0.6, np.nanmax(D)), cmap="magma", shading="nearest")
    plt.colorbar(pc, ax=a, pad=0.02).set_label(r"duty $W_{50}/T$"); boundaries(a, G2, logy=False)
    a.plot([REF["beta"]], [REF["qoff"]], "*", color="w", mec="k", ms=22, label=r"REF$_{ST}$")
    a.set_xscale("log"); a.set_xlabel(r"$\beta$"); a.set_ylabel(r"$q_{\rm off}$"); a.set_title(r"Spike width at $\bar S=1$"); a.legend(loc="upper right", fontsize=14); a.grid(False)

def p_line_S(a):
    for L, col, lab in ((LS, C_ST, r"REF$_{ST}$"), (LP, C_PREV, "previous ref.")):
        o = L[L.osc == 1]; a.loglog(o.x, o["T"], "-", color=col, label=lab + r": $T$"); a.loglog(o.x, o.amp * 10, "--", color=col, label=lab + r": $10\times$amp")
    st, x, i = LS.ReL.values < 0, LS.x.values, 0
    while i < len(x):                                   # contiguous stretches with a stable fixed point (REF_ST)
        if st[i]:
            j = i
            while j + 1 < len(x) and st[j + 1]: j += 1
            a.axvspan(x[i], x[j], color="#dddddd", alpha=0.6, zorder=0); i = j + 1
        else: i += 1
    if S_DOSE: a.axvspan(S_DOSE[0], S_DOSE[1], color=C_TINT, alpha=0.35, zorder=0)
    if st.any(): a.plot([], [], "s", color="#dddddd", ms=12, label="fixed point stable")
    a.set_xlabel(r"$\bar S$"); a.set_ylabel(r"$T$, $10\times$amplitude"); a.set_title("Along the stimulus axis"); a.legend(loc="upper right", fontsize=13)

def p_map_ReL(a):
    x, y, R = pivot(G1, "ReL"); m = np.nanmax(np.abs(R)) if np.isfinite(R).any() else 1
    pc = a.pcolormesh(x, y, R, cmap="RdBu_r", vmin=-0.05, vmax=0.05, shading="nearest"); plt.colorbar(pc, ax=a, pad=0.02).set_label(r"max Re$\,\lambda$ (clipped)")
    boundaries(a, G1); a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(r"$\bar S$"); a.set_ylabel(r"$\varepsilon_H$"); a.set_title("Fixed-point stability"); a.grid(False)

def p_map_hopf_ratio(a):
    x, y, H = pivot(G1, "hopf_ratio"); H = np.where(H > 0, H, np.nan)
    pc = a.pcolormesh(x, y, H, norm=LogNorm(vmin=0.5, vmax=max(2, np.nanmax(H))), cmap="cividis", shading="nearest"); plt.colorbar(pc, ax=a, pad=0.02).set_label(r"$T\,{\rm Im}\lambda/2\pi$")
    boundaries(a, G1); a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(r"$\bar S$"); a.set_ylabel(r"$\varepsilon_H$"); a.set_title("Relaxation index"); a.grid(False)

def p_map_amp(a):
    x, y, A = pivot(G2, "amp"); A = np.where(pivot(G2, "osc")[2] == 1, A, np.nan)
    pc = a.pcolormesh(x, y, A, cmap="viridis", shading="nearest"); plt.colorbar(pc, ax=a, pad=0.02).set_label("amplitude")
    boundaries(a, G2, logy=False); a.set_xscale("log"); a.set_xlabel(r"$\beta$"); a.set_ylabel(r"$q_{\rm off}$"); a.set_title(r"Amplitude at $\bar S=1$"); a.grid(False)

def line_panel(a, L, xlabel, title):
    o = L[L.osc == 1]; a.loglog(o.x, o["T"], "-", color=C_ST, label=r"$T$"); a.loglog(o.x, o.amp * 10, "--", color=C_ST, label=r"$10\times$amp"); a.loglog(o.x, o.W50, ":", color=C_ST, lw=3, label=r"$W_{50}$")
    for p in hopf_points(L): a.axvline(p, color=C_HOPF, ls="--", lw=2)
    a.plot([], [], "--", color=C_HOPF, label="Hopf"); a.set_xlabel(xlabel); a.set_title(title); a.legend(loc="upper center", fontsize=14)
def p_line_epsH(a): line_panel(a, LE, r"$\varepsilon_H$", r"Along $\varepsilon_H$ ($\bar S=1$)")
def p_line_beta(a): line_panel(a, LB, r"$\beta$", r"Along $\beta$ ($\bar S=1$)")

def p_class(a):
    x, y, C = pivot(G1, "class"); codes = {0: 0, 1: 1, 2: 2, 3: 3}; M = np.vectorize(lambda c: (4 if c >= 10 else codes.get(int(c), 3)))(C)
    cmap = ListedColormap(["#dddddd", "#2ca02c", "#e69f00", "#9467bd", "#c0392b"]); norm = BoundaryNorm(np.arange(-0.5, 5.5), cmap.N)
    a.pcolormesh(x, y, M, cmap=cmap, norm=norm, shading="nearest")
    from matplotlib.patches import Patch
    a.legend(handles=[Patch(color=c, label=l) for c, l in zip(cmap.colors, ["steady", "cycle, FP unstable", "cycle + stable FP", "no cycle, FP unstable", "several FPs"])], loc="lower right", fontsize=13)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(r"$\bar S$"); a.set_ylabel(r"$\varepsilon_H$"); a.set_title("Classification"); a.grid(False)

MAIN = [p_map_S_epsH, p_map_beta_qoff, p_line_S]
SUPP = [p_map_ReL, p_map_hopf_ratio, p_map_amp, p_line_epsH, p_line_beta, p_class]
if __name__ == "__main__":
    compose(MAIN, 3, f"{NAME}_main")
    compose(SUPP, 3, f"{NAME}_supp")
