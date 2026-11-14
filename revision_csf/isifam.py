"""
isifam.py — interval-density families with a fixed offset, maximum-likelihood fits per segment (library for b2).

Every family describes the shifted interval y = x - delta (delta >= 0 fixed) with two free parameters, written as the
mean m and the coefficient of variation cv of y:
    gamma     shape 1/cv^2, scale m cv^2            invgauss   mean m, shape m/cv^2
    lognorm   sigma^2 = ln(1+cv^2)                  fpt_<form> y = scale * tau, tau = first-passage time of the reduced
                                                    recovery process (b1 tables); cv <-> nu, scale = m / E[tau]
Observed intervals are differences of spike times recorded on a frame grid of step q, so the likelihood of an interval
is the density averaged with the triangular kernel of half-width q (9 nodes); q = 0 gives the density itself.
"""
import os, pickle, numpy as np
from scipy.special import gammaln, logsumexp
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "results")
NODES = np.linspace(-0.8, 0.8, 9); LOGW = np.log((1 - np.abs(NODES)) / np.sum(1 - np.abs(NODES)))
FLOOR = -700.0

# ---- standardised first-passage tables -------------------------------------------------------------------
ZG = np.arange(-12.0, 16.0 + 1e-9, 1 / 160)
class FPT:
    def __init__(self, form):
        res = pickle.load(open(os.path.join(OUT, f"b1_fpt_density_{form}.pkl"), "rb")); res = sorted(res, key=lambda r: r["mean"] and r["sd"] / r["mean"])
        self.cv = np.array([r["sd"] / r["mean"] for r in res]); self.nu = np.array([r["nu"] for r in res]); self.mean = np.array([r["mean"] for r in res])
        assert np.all(np.diff(self.cv) > 0)
        self.logg = np.empty((len(res), len(ZG))); self.slope = np.empty(len(res)); self.cdf = []
        for k, r in enumerate(res):
            s = r["ds"] * np.arange(len(r["f"])); z = (s - r["mean"]) / r["sd"]; g = r["sd"] * r["f"]
            with np.errstate(divide="ignore"): lg = np.where(g > 0, np.log(np.maximum(g, 1e-300)), FLOOR)
            self.logg[k] = np.maximum(np.interp(ZG, z, lg, left=FLOOR, right=FLOOR), FLOOR)
            self.slope[k] = (self.logg[k, -1] - self.logg[k, -161]) / 1.0                        # d log g / dz in the far tail
            c = np.cumsum(g) * (z[1] - z[0]); self.cdf.append((z, c / c[-1]))
    def _one(self, k, z):
        lg = np.interp(z, ZG, self.logg[k], left=FLOOR)
        hi = z > ZG[-1]
        if np.any(hi): lg = np.where(hi, self.logg[k, -1] + self.slope[k] * (z - ZG[-1]), lg)
        return lg
    def logpdf(self, y, m, cv):
        """log density of y = scale*tau with mean m and coefficient of variation cv (cv inside the table range)"""
        k = int(np.clip(np.searchsorted(self.cv, cv) - 1, 0, len(self.cv) - 2)); w = (np.log(cv) - np.log(self.cv[k])) / (np.log(self.cv[k + 1]) - np.log(self.cv[k]))
        sd = m * cv; z = (y - m) / sd
        lg = (1 - w) * self._one(k, z) + w * self._one(k + 1, z) - np.log(sd)
        return np.where(y > 0, np.maximum(lg, FLOOR), -np.inf)
    def nu_of_cv(self, cv): return float(np.exp(np.interp(np.log(cv), np.log(self.cv), np.log(self.nu))))
    def mean_of_cv(self, cv): return float(np.interp(np.log(cv), np.log(self.cv), self.mean))
    def sample(self, rng, n, m, cv):
        k = int(np.clip(np.searchsorted(self.cv, cv) - 1, 0, len(self.cv) - 2)); w = (np.log(cv) - np.log(self.cv[k])) / (np.log(self.cv[k + 1]) - np.log(self.cv[k]))
        kk = np.where(rng.random(n) < w, k + 1, k); u = rng.random(n); z = np.empty(n)
        for j in (k, k + 1):
            i = kk == j; z[i] = np.interp(u[i], self.cdf[j][1], self.cdf[j][0])
        return m + m * cv * z
_FPT = {}
def fpt(form):
    if form not in _FPT: _FPT[form] = FPT(form)
    return _FPT[form]

# ---- log densities in (m, cv) ------------------------------------------------------------------------------
def lp_gamma(y, m, cv):
    k = 1 / cv**2; th = m * cv**2
    with np.errstate(divide="ignore", invalid="ignore"): return np.where(y > 0, (k - 1) * np.log(np.where(y > 0, y, 1)) - y / th - gammaln(k) - k * np.log(th), -np.inf)
def lp_invgauss(y, m, cv):
    lam = m / cv**2
    with np.errstate(divide="ignore", invalid="ignore"): yy = np.where(y > 0, y, 1); return np.where(y > 0, 0.5 * np.log(lam / (2 * np.pi * yy**3)) - lam * (yy - m) ** 2 / (2 * m * m * yy), -np.inf)
def lp_lognorm(y, m, cv):
    s2 = np.log1p(cv**2); mu = np.log(m) - s2 / 2
    with np.errstate(divide="ignore", invalid="ignore"): yy = np.where(y > 0, y, 1); return np.where(y > 0, -np.log(yy) - 0.5 * np.log(2 * np.pi * s2) - (np.log(yy) - mu) ** 2 / (2 * s2), -np.inf)
LP = {"gamma": lp_gamma, "invgauss": lp_invgauss, "lognorm": lp_lognorm,
      "fpt_jacobi": lambda y, m, cv: fpt("jacobi").logpdf(y, m, cv), "fpt_cle": lambda y, m, cv: fpt("cle").logpdf(y, m, cv)}
def cv_bounds(fam):
    if fam.startswith("fpt_"): f = fpt(fam[4:]); return f.cv[0] * 1.001, f.cv[-1] * 0.999
    return 0.005, 5.0

def loglik(fam, x, q, delta, m, cv, eps=0.0, span=None):
    """sum over intervals of log[(1-eps) <f>_q (x) + eps/span]; <.>_q = triangular average over the frame step"""
    if q > 0: lp = logsumexp(LP[fam](x[:, None] + q * NODES[None, :] - delta, m, cv) + LOGW[None, :], axis=1)
    else: lp = LP[fam](x - delta, m, cv)
    if eps > 0: lp = np.logaddexp(np.log1p(-eps) + lp, np.log(eps / span))
    return float(np.sum(lp))

def fit(fam, x, q, delta, eps=0.0):
    """maximum likelihood over (m, cv); returns (logL, m, cv)"""
    y = x - delta; span = 6 * x.mean()
    good = y[y > 0]; m0 = good.mean() if len(good) else x.mean(); c0 = max(good.std() / m0 if len(good) > 2 else 0.2, 0.03)
    lo, hi = cv_bounds(fam)
    def nll(p):
        m, cv = np.exp(p[0]), np.exp(p[1])
        if not (lo < cv < hi) or m <= 0: return 1e12
        v = -loglik(fam, x, q, delta, m, cv, eps, span); return v if np.isfinite(v) else 1e12
    best = None
    for c in (c0, 0.6 * c0, 1.6 * c0):
        c = float(np.clip(c, lo * 1.05, hi * 0.95)); r = minimize(nll, [np.log(m0), np.log(c)], method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-7, maxiter=600))
        if best is None or r.fun < best.fun: best = r
    return -best.fun, float(np.exp(best.x[0])), float(np.exp(best.x[1]))

def sample(fam, rng, n, m, cv):
    if fam == "gamma": return rng.gamma(1 / cv**2, m * cv**2, n)
    if fam == "invgauss": return rng.wald(m, m / cv**2, n)
    if fam == "lognorm": s2 = np.log1p(cv**2); return rng.lognormal(np.log(m) - s2 / 2, np.sqrt(s2), n)
    return fpt(fam[4:]).sample(rng, n, m, cv)

def detrend(x, t):
    if len(x) >= 4: b = np.polyfit(t, x, 1); return x - np.polyval(b, t) + x.mean()
    return x

def boot_summary(R, obs, names, fail=1e9):
    """Parametric-bootstrap summary.  R: one row per replicate with the total log-likelihood of every fitted family and
    the generating family in column 'truth'; obs: observed totals.  A total beyond `fail` in absolute value marks a
    replicate in which the fit of that family was infeasible in at least one segment (an interval at or below the
    offset after frame rounding and detrending); such replicates are excluded pairwise and counted in n_failed."""
    import pandas as pd
    rows = []
    for tr in names:
        g = R[R.truth == tr]
        for alt in names:
            if alt == tr: continue
            ok = (g[tr].abs() < fail) & (g[alt].abs() < fail); dsim = (g[tr] - g[alt])[ok].values; dobs = obs[tr] - obs[alt]
            rows.append(dict(truth=tr, alt=alt, dlogL_obs=dobs, n_valid=int(ok.sum()), n_failed=int((~ok).sum()),
                             sim_mean=dsim.mean() if len(dsim) else np.nan, sim_lo=np.quantile(dsim, .025) if len(dsim) else np.nan,
                             sim_hi=np.quantile(dsim, .975) if len(dsim) else np.nan, p_sim_le_obs=float(np.mean(dsim <= dobs)) if len(dsim) else np.nan))
    return pd.DataFrame(rows)
