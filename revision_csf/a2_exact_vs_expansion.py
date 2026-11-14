"""
a2_exact_vs_expansion.py — exact first-passage cumulants for several noise forms vs the general expansion of a1.

Exact: u = -k1', v = -k2', w = -k3' solve the linear first-order system
    a D u' + b u = 1,   a D v' + b v = 2 a D u^2,   a D w' + b w = 6 a D u v
on (0, H*), with u = v = w = 0 at a reflecting regular boundary (D(0) > 0) and the bounded (outer) solution at an
entrance boundary (Jacobi, nu^2 <= 2).  k_n(H0) = int_{H0}^{H*}.  Radau integration, rtol 1e-10.
Monte Carlo (Euler-Maruyama with reflection at 0, two time steps) is an independent check at three noise levels.
Output: results/a2_exact_vs_expansion.csv, results/a2_mc_check.csv
"""
import os, sys, numpy as np, pandas as pd
from scipy.integrate import solve_ivp
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
HINF, H0, HS = 0.997, 0.064, 0.408
FORMS = {"jacobi": (lambda y: 1 - y, lambda y: y * (1 - y)),
         "one_minus_H": (lambda y: 1 - y, lambda y: 1 - y),
         "cle": (lambda y: HINF - y, lambda y: HINF * (1 - y) + (1 - HINF) * y)}

def exact(form, nu, h0=H0, hs=HS):
    b, D = FORMS[form]; a = nu**2 / 2
    def rhs(yv, s):
        u, v, w = s[:3]; d = a * D(yv); bb = b(yv)
        du, dv, dw = (1 - bb * u) / d, (2 * d * u * u - bb * v) / d, (6 * d * u * v - bb * w) / d
        on = 1.0 if yv >= h0 else 0.0
        return [du, dv, dw, on * u, on * v, on * w]
    if D(0.0) > 0: y0, s0 = 0.0, [0, 0, 0, 0, 0, 0]
    else:
        y0 = 1e-7; s0 = [1 / b(y0), 0, 0, 0, 0, 0]          # entrance boundary: bounded solution, u(0) = 1/b(0)
    # integrate in two legs so that the indicator switch at H0 is exact
    s1 = solve_ivp(rhs, (y0, h0), s0, method="Radau", rtol=1e-10, atol=1e-14).y[:, -1]; s1[3:] = 0
    s2 = solve_ivp(rhs, (h0, hs), s1, method="Radau", rtol=1e-10, atol=1e-14).y[:, -1]
    k1, k2, k3 = s2[3:]
    return dict(mean=k1, var=k2, k3=k3, CV=np.sqrt(k2) / k1, skew=k3 / k2**1.5, R=k3 * k1 / k2**2)

if __name__ == "__main__":
    C = pd.read_csv(os.path.join(OUT, "a1_noise_forms_reference.csv")).set_index("form")
    rows = []
    for form in FORMS:
        c = C.loc[form]
        for nu in [0.02, 0.05, 0.1, 0.15, 0.2, 0.281, 0.4, 0.5, 0.7, 1.0]:
            a = nu**2 / 2; e = exact(form, nu)
            mean2 = c.T0 + a * c.mean_a1 + a * a * c.mean_a2; var2 = a * c.var_a1 + a * a * c.var_a2
            k3_lead, k3_next = a * a * c.k3_a2, a * a * c.k3_a2 + a**3 * c.k3_a3
            rows.append(dict(form=form, nu=nu, mean=e["mean"], mean_exp=mean2, var=e["var"], var_exp=var2, k3=e["k3"], k3_lead=k3_lead, k3_next=k3_next,
                             CV=e["CV"], skew=e["skew"], R=e["R"], R_limit=c.R, R_next=c.R * (1 + c.cR_per_a * a),
                             err_mean=mean2 / e["mean"] - 1, err_sd=np.sqrt(max(var2, 0)) / np.sqrt(e["var"]) - 1, err_k3_lead=k3_lead / e["k3"] - 1, err_k3_next=k3_next / e["k3"] - 1))
    X = pd.DataFrame(rows); X.to_csv(os.path.join(OUT, "a2_exact_vs_expansion.csv"), index=False)
    pd.set_option("display.width", 250)
    print(X[["form", "nu", "mean", "CV", "skew", "R", "R_limit", "R_next", "err_mean", "err_sd", "err_k3_lead", "err_k3_next"]].round(4).to_string(index=False))

    # cross-check Jacobi against the manuscript's quadrature of the Siegert recursion, if present
    p = os.path.join(os.path.dirname(OUT), "..", "m1_small_nu_asymptotics", "results")
    if os.path.isdir(p): print("\nexisting quadrature files:", sorted(os.listdir(p))[:12])

    if "--mc" in sys.argv:
        def mc(form, nu, dt, npath=200000, seed=1):
            b, D = FORMS[form]; rng = np.random.default_rng(seed); h = np.full(npath, H0); t = np.zeros(npath); alive = np.arange(npath); tau = np.empty(npath); sq = np.sqrt(dt); n = 0
            while len(alive):
                x = h[alive]; x = x + b(x) * dt + nu * np.sqrt(np.clip(D(x), 0, None)) * sq * rng.standard_normal(len(alive)); x = np.abs(x); n += 1
                hit = x >= HS; tau[alive[hit]] = n * dt; h[alive] = x; alive = alive[~hit]
            return tau
        rows = []
        for form in FORMS:
            for nu in [0.15, 0.3, 0.6]:
                e = exact(form, nu); r = dict(form=form, nu=nu, mean=e["mean"], CV=e["CV"], skew=e["skew"])
                for dt in (4e-4, 1e-4):
                    tau = mc(form, nu, dt); m = tau.mean(); s = tau.std(); r[f"mc_mean_dt{dt:g}"] = m; r[f"mc_CV_dt{dt:g}"] = s / m; r[f"mc_skew_dt{dt:g}"] = float(np.mean(((tau - m) / s) ** 3))
                rows.append(r); print(r, flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(OUT, "a2_mc_check.csv"), index=False)
