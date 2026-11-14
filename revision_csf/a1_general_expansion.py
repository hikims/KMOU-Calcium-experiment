"""
a1_general_expansion.py — small-noise expansion of the first-passage cumulants for a general diffusion coefficient.

Reduced recovery process in rescaled time:  dH = b(H) ds + nu sqrt(D(H)) dW,  b > 0 on [H0, H*],  a = nu^2/2.
Cumulant functions k_n(x) of the passage time to H* obey (from log E exp(-lambda tau)):
    L k1 = -1,   L k2 = -2 a D (k1')^2,   L k3 = -6 a D k1' k2',      L = b d/dx + a D d^2/dx^2 .
With u = -k1', v = -k2', w = -k3' the outer (regular) expansion is
    u = 1/b - a (D/b) u',   v = 2 a D u^2 / b - a (D/b) v',   w = 6 a D u v / b - a (D/b) w'
and  k_n(H0) = int_{H0}^{H*} (u, v, w) dy.  Leading orders:
    T0 = int 1/b,   Var = nu^2 J1,  J1 = int D/b^3,   k3 = 3 nu^4 J2,  J2 = int D^2/b^5,
    R = lim skew/CV = 3 T0 J2 / J1^2 >= 3   (Cauchy-Schwarz with dmu = dy/b, g = D/b^2; equality iff D ∝ b^2).
Part 1 checks the recursion symbolically against Theorem 5.2 (Jacobi).  Part 2 tabulates coefficients for the
noise forms; Part 3 the dependence of R on the geometry.  Output: results/a1_*.csv
"""
import os, numpy as np, pandas as pd, sympy as sp
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
y, a = sp.symbols("y a", positive=True)

def series_uvw(b, D, order=3):
    """coefficients of a^j in u, v, w (lists indexed by power of a)"""
    op = lambda f: sp.simplify(-(D / b) * sp.diff(f, y))
    u = [1 / b]
    for j in range(1, order + 1): u.append(op(u[j - 1]))
    v = [sp.Integer(0)]
    for j in range(1, order + 1):
        src = sum(u[i] * u[j - 1 - i] for i in range(j))                # coefficient of a^{j-1} in u^2
        v.append(sp.simplify(2 * D * src / b + op(v[j - 1])))
    w = [sp.Integer(0), sp.Integer(0)]
    for j in range(2, order + 1):
        src = sum(u[i] * v[j - 1 - i] for i in range(j))                # coefficient of a^{j-1} in u v
        w.append(sp.simplify(6 * D * src / b + op(w[j - 1])))
    return u, v, w

# ---- Part 1: Jacobi, symbolic check against Theorem 5.2 --------------------------------------------------
H0, Hs = sp.symbols("H0 Hs", positive=True)
u, v, w = series_uvw(1 - y, y * (1 - y), order=3)
I = lambda f: sp.integrate(sp.simplify(f), (y, H0, Hs))
z0, zs = 1 / (1 - H0), 1 / (1 - Hs); T0 = sp.log(zs / z0)
paper = dict(mean1=-(zs - z0 - T0), mean2=zs**2 - z0**2 - 3 * (zs - z0) + T0,
             var1=2 * (zs - z0 - T0), var2=-4 * (zs**2 - z0**2) + 14 * (zs - z0) - 6 * T0,
             k3_2=6 * (zs**2 - z0**2 - 4 * (zs - z0) + 2 * T0))
mine = dict(mean1=I(u[1]), mean2=I(u[2]), var1=I(v[1]), var2=I(v[2]), k3_2=I(w[2]))
num = {H0: sp.Rational(64, 1000), Hs: sp.Rational(408, 1000)}
print("Part 1: Jacobi recursion vs Theorem 5.2 (difference at the reference geometry)")
for k in paper:
    d = float((mine[k] - paper[k]).subs(num)); print(f"  {k:6s} recursion {float(mine[k].subs(num)):+.8f}  paper {float(paper[k].subs(num)):+.8f}  diff {d:+.1e}")
    assert abs(d) < 1e-10
print("  third cumulant, next order (a^3 coefficient, new):", sp.simplify(I(w[3]).subs(num)).evalf(8))

# ---- Part 2: noise forms at the reference geometry --------------------------------------------------------
from scipy.integrate import quad
HINF = 0.997
FORMS = {   # name: (b, D, description)
    "jacobi":      (1 - y, y * (1 - y),                 "H(1-H)  (manuscript)"),
    "one_minus_H": (1 - y, 1 - y,                       "1-H  (chemical Langevin, inactivation rate -> 0)"),
    "cle":         (HINF - y, HINF * (1 - y) + (1 - HINF) * y, "alpha(1-H)+beta H, H_inf=0.997 (chemical Langevin)"),
    "additive":    (1 - y, sp.Integer(1),               "constant"),
    "H":           (1 - y, y,                           "H"),
}
def coeffs(b, D, h0, hs):
    u, v, w = series_uvw(b, D, order=3)
    q = lambda f: quad(sp.lambdify(y, f, "numpy"), h0, hs, epsabs=1e-13, epsrel=1e-12)[0]
    T0 = q(u[0]); m1, m2 = q(u[1]), q(u[2]); v1, v2 = q(v[1]), q(v[2]); k2, k3n = q(w[2]), q(w[3])
    J1, J2 = v1 / 2, k2 / 12
    R = 3 * T0 * J2 / J1**2
    # next-order correction of skew/CV:  skew/CV = R (1 + cR a + ...)
    #   sd^2 = a v1 + a^2 v2, k3 = a^2 k2 + a^3 k3n, mean = T0 + a m1
    cR = k3n / k2 - 2 * v2 / v1 + m1 / T0          # skew/CV = k3*mean/var^2
    return dict(T0=T0, mean_a1=m1, mean_a2=m2, var_a1=v1, var_a2=v2, k3_a2=k2, k3_a3=k3n, J1=J1, J2=J2, R=R, cR_per_a=cR)
rows = []
for nm, (b, D, desc) in FORMS.items():
    c = coeffs(b, D, 0.064, 0.408); c.update(form=nm, D=desc); rows.append(c)
T = pd.DataFrame(rows).set_index("form"); T.to_csv(os.path.join(OUT, "a1_noise_forms_reference.csv"))
pd.set_option("display.width", 250)
print("\nPart 2: coefficients at (H0, H*) = (0.064, 0.408)   [mean = T0 + a m1 + a^2 m2, var = a v1 + a^2 v2, k3 = a^2 k3_2 + a^3 k3_3]")
print(T[["T0", "mean_a1", "mean_a2", "var_a1", "var_a2", "k3_a2", "k3_a3", "R", "cR_per_a"]].round(5).to_string())

# ---- Part 3: R over the geometry ---------------------------------------------------------------------------
print("\nPart 3: R(H0, H*) for the forms (rows H0, columns H*)")
grid0, grids = [0.02, 0.064, 0.15, 0.3], [0.2, 0.408, 0.6, 0.8, 0.95]
rows = []
for nm in ("jacobi", "one_minus_H", "cle"):
    b, D, _ = FORMS[nm]
    fT = sp.lambdify(y, 1 / b); f1 = sp.lambdify(y, D / b**3); f2 = sp.lambdify(y, D**2 / b**5)
    for h0 in grid0:
        for hs in grids:
            if hs <= h0: continue
            T0, J1, J2 = (quad(f, h0, hs)[0] for f in (fT, f1, f2)); rows.append(dict(form=nm, H0=h0, Hs=hs, R=3 * T0 * J2 / J1**2))
G = pd.DataFrame(rows); G.to_csv(os.path.join(OUT, "a1_R_geometry.csv"), index=False)
for nm in ("jacobi", "one_minus_H", "cle"):
    print(nm); print(G[G.form == nm].pivot(index="H0", columns="Hs", values="R").round(3).to_string())
