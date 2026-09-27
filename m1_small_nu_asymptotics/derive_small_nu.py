"""
derive_small_nu.py — symbolic small-ν expansion of the Dynkin moments of the reduced recovery process.

T_n(x) = n ∫_x^{H*} ∫_0^1 T_{n-1}(y w^{1/k}) / (1 − y w^{1/k}) dw dy,  k = 2/ν²,  a := ν²/2 = 1/k.
With w^{a} = exp(a ln w) and ∫_0^1 (ln w)^m dw = (−1)^m m!, the inner integral is the formal series
   ∫_0^1 F(y w^a) dw = Σ_m (−a)^m (y ∂_y)^m F(y),      F = T_{n-1}(y)/(1−y)
(an asymptotic series in a).  In z = 1/(1−y):  y ∂_y = (z² − z) ∂_z,  dy = dz/z², 1/(1−y) = z, so every term is a
polynomial in z times a power of ln z and the integrals are elementary.  Output: results/asymptotic_coefficients.csv
and results/asymptotic_coefficients.txt (z0 = 1/(1−H0), zs = 1/(1−H*), T0 = ln(zs/z0)).
"""
import os, sympy as sp
z, zs, zx, a = sp.symbols('z z_s z_x a', positive=True)
ORD = 2
L = lambda f: sp.expand((z**2 - z) * sp.diff(f, z))
def opseries(F):
    out, g = 0, F
    for m in range(ORD + 1): out += (-a)**m * g; g = L(g)
    return sp.expand(out)
def integ_y(expr):
    tot = 0
    for t in sp.Add.make_args(sp.expand(expr / z**2)):
        c, rest = t.as_independent(z); F = sp.integrate(rest, z); tot += c * (F.subs(z, zs) - F.subs(z, zx))
    return sp.expand(tot)
def Tn(n, Tprev_z):
    integrand = sp.expand(sp.series(opseries(sp.expand(Tprev_z * z)), a, 0, ORD + 1).removeO())
    return sp.expand(n * integ_y(integrand))
T1 = Tn(1, sp.Integer(1)); T2 = Tn(2, T1.subs(zx, z)); T3 = Tn(3, T2.subs(zx, z))
tr = lambda e: sp.expand(sp.series(sp.expand(e), a, 0, ORD + 1).removeO())
var, k3 = tr(T2 - T1**2), tr(T3 - 3*T1*T2 + 2*T1**3)
os.makedirs("results", exist_ok=True)
rows, lines = [], ["a = nu^2/2, z_x = 1/(1-H0), z_s = 1/(1-H*)"]
for nm, e in (("mean", T1), ("var", var), ("kappa3", k3)):
    for n in range(ORD + 1):
        c = sp.simplify(e.coeff(a, n)); rows.append(dict(moment=nm, order=n, coefficient=str(c))); lines.append(f"{nm}: a^{n} : {c}")
import pandas as pd; pd.DataFrame(rows).to_csv("results/asymptotic_coefficients.csv", index=False)
open("results/asymptotic_coefficients.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
