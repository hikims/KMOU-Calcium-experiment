"""
c1_hierarchical_slope.py — the sigma-T_av law separated into record, cell and dose components (pooled slope 1.33,
errors-in-variables and pseudo-replication).

Units: segment (authors' file) < cell-at-one-dose (train) < cell (record, ST) < record.  y = log sd (detrended), x = log T_av.
A  pooled slope with segment- and record-bootstrap intervals; the same with the spike offset removed, x' = log(T_av - 8.64 s)
B  hybrid (between/within) regression: x is split into
     xr  record mean                                   (between records: dose, day, dish)
     xc  cell mean - record mean                       (between cells of one record: the cell axis)
     xs  cell-at-dose mean - cell mean                 (same cell at different doses, step records: the dose axis within a cell)
     xw  segment - cell-at-dose mean                   (same cell, same dose, later segment: drift)
   OLS with record-bootstrap intervals, and a mixed model with random intercepts for record and cell.
C  what the model predicts for the same design: sd = c (T_av - T_sp)^g with g = 1 (common nu) or 3/2 (common sigma_H)
   gives which log-log slope over the observed T_av values?
D  calibration of the estimators by simulation: synthetic segments with the data design (same T_av, lengths, frame
   steps, detrending), inverse-Gaussian intervals with offset, true exponent s and a cell-level log-normal scatter of
   the CV chosen to reproduce the residual spread; bias and spread of every slope estimator.
Output: results/c1_*.csv
"""
import os, sys, numpy as np, pandas as pd, warnings
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import isifam as F
OUT = os.path.join(HERE, "results"); CODE = os.path.join(HERE, ".."); DELTA = 16.4 * 0.527
S = pd.read_csv(os.path.join(CODE, "data_dose_units", "results", "segment_table.csv")); S = S[S.n >= 10].reset_index(drop=True)
S["cell"] = S.record + ":" + S.ST.astype(str); S["y"] = np.log(S.sd_detr); S["x"] = np.log(S.Tav); S["xo"] = np.log(S.Tav - DELTA)
def split(d, col):
    xr = d.groupby("record")[col].transform("mean"); xcl = d.groupby("cell")[col].transform("mean"); xcd = d.groupby("train")[col].transform("mean")
    return pd.DataFrame(dict(xr=xr, xc=xcl - xr, xs=xcd - xcl, xw=d[col] - xcd))
COMP = ["xr", "xc", "xs", "xw"]
def ols(X, y): X1 = np.column_stack([np.ones(len(y)), X]); return np.linalg.lstsq(X1, y, rcond=None)[0][1:]
def estimates(d, col="x"):
    Z = split(d, col); b = ols(Z[COMP].values, d.y.values); return dict(pooled=float(np.polyfit(d[col], d.y, 1)[0]), **dict(zip(COMP, b)))
def rec_boot(d, col, nb=2000, seed=0):
    rng = np.random.default_rng(seed); recs = d.record.unique(); g = {r: d[d.record == r] for r in recs}; out = []
    for _ in range(nb):
        pick = rng.choice(recs, len(recs)); dd = pd.concat([g[r].assign(record=f"{r}#{i}", cell=g[r].cell + f"#{i}", train=g[r].train + f"#{i}") for i, r in enumerate(pick)]); out.append(estimates(dd, col))
    return pd.DataFrame(out)

G = {}
def _init(meta, xbar, a0): G.update(meta=meta, xbar=xbar, a0=a0)
def simulate(a):
    seed, s_true, tau_cell = a; rng = np.random.default_rng(seed); meta, XBAR, A0 = G["meta"], G["xbar"], G["a0"]
    lcv_cell = {c: rng.normal(0, tau_cell) for c in S.cell.unique()}; rows = []
    for i, r in S.iterrows():
        m = r.Tav - DELTA; sd = np.exp(A0 + s_true * (np.log(r.Tav) - XBAR) + lcv_cell[r.cell]); cv = float(np.clip(sd / m, 0.03, 3.0)); n = int(r.n) - 1; q = meta[i]
        x = F.sample("invgauss", rng, n, m, cv) + DELTA; t = np.cumsum(x)
        if q > 0: tt = np.round(np.concatenate([[0.0], t]) / q) * q; x = np.diff(tt); t = tt[1:]
        xd = F.detrend(x, t); rows.append((xd.mean(), xd.std(ddof=1)))
    d = S[["record", "cell", "train"]].copy(); d["Tav"], d["sd"] = zip(*rows); d = d[d.sd > 0]; d["y"] = np.log(d.sd); d["x"] = np.log(d.Tav)
    return dict(s_true=s_true, **estimates(d, "x"), resid_sd=float(np.std(d.y - np.polyval(np.polyfit(d.x, d.y, 1), d.x))))

if __name__ == "__main__":
    pd.set_option("display.width", 220)
    print(f"{len(S)} segments, {S.train.nunique()} cell-at-dose trains, {S.cell.nunique()} cells, {S.record.nunique()} records; cells seen at >= 2 doses: {(S.groupby('cell').dose.nunique() >= 2).sum()}; trains with >= 2 segments: {(S.groupby('train').size() >= 2).sum()}")
    Z = split(S, "x"); print("variance of the components of x = log T_av:", Z[COMP].var().round(4).to_dict(), " total", round(S.x.var(), 4))
    rows = []
    for col, lab in (("x", "log T_av"), ("xo", "log (T_av - 8.64 s)")):
        e = estimates(S, col); B = rec_boot(S, col); rng = np.random.default_rng(1)
        segb = np.array([np.polyfit(S[col].values[i], S.y.values[i], 1)[0] for i in (rng.integers(0, len(S), len(S)) for _ in range(2000))])
        rows.append(dict(regressor=lab, component="pooled", slope=e["pooled"], rec_lo=B.pooled.quantile(.025), rec_hi=B.pooled.quantile(.975), seg_lo=np.quantile(segb, .025), seg_hi=np.quantile(segb, .975)))
        for c in COMP: rows.append(dict(regressor=lab, component=c, slope=e[c], rec_lo=B[c].quantile(.025), rec_hi=B[c].quantile(.975)))
    R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "c1_slopes.csv"), index=False); print("\nA/B. slopes (95 % intervals: record bootstrap; pooled also segment bootstrap)"); print(R.round(3).to_string(index=False))
    # mixed model
    try:
        import statsmodels.formula.api as smf
        d = pd.concat([S[["y", "record", "cell", "dose"]], Z], axis=1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore"); mm = smf.mixedlm("y ~ xr + xc + xs + xw", d, groups=d["record"], re_formula="1", vc_formula={"cell": "0 + C(cell)"}).fit(reml=True)
        ci = mm.conf_int(); M = pd.DataFrame(dict(slope=mm.params[COMP], lo=ci.loc[COMP, 0], hi=ci.loc[COMP, 1])); M.to_csv(os.path.join(OUT, "c1_mixed_model.csv"))
        print("\nmixed model (random intercepts: record, cell within record)"); print(M.round(3).to_string()); print("variance components: record %.3f, cell %.3f, residual %.3f" % (float(mm.cov_re.iloc[0, 0]), float(mm.vcomp[0]), mm.scale))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore"); m2 = smf.mixedlm("y ~ x + C(dose)", pd.concat([S[["y", "x", "record", "cell", "dose"]]], axis=1), groups=S["record"], vc_formula={"cell": "0 + C(cell)"}).fit(reml=True)
        print("mixed model with dose as a fixed factor: slope of x %.3f [%.3f, %.3f]" % (m2.params["x"], *m2.conf_int().loc["x"]))
    except Exception as ex: print("mixed model failed:", ex)
    # C: model prediction over the observed T_av values
    T = S.Tav.values; rows = []
    for g, lab in ((1.0, "common nu: sd ∝ (T_av - T_sp)"), (1.5, "common sigma_H: sd ∝ (T_av - T_sp)^1.5")):
        for tsp in (0.0, DELTA):
            rows.append(dict(law=lab, T_sp=tsp, loglog_slope_over_data_Tav=float(np.polyfit(np.log(T), g * np.log(T - tsp), 1)[0])))
    Cm = pd.DataFrame(rows); Cm.to_csv(os.path.join(OUT, "c1_model_prediction.csv"), index=False); print("\nC. log-log slope implied over the observed T_av values"); print(Cm.round(3).to_string(index=False))
    lin = np.polyfit(S.Tav, S.sd_detr, 1); print("linear fit sd = %.3f (T_av - %.1f s)" % (lin[0], -lin[1] / lin[0]))
    # D: calibration
    from b2_family_comparison import load_segments
    segs = load_segments(trim=0); qmap = {s["id"]: s["q"] for s in segs}; SIMMETA = [qmap.get(f"{r.train}|{r.segment}", 0.0) for r in S.itertuples()]
    XBAR = S.x.mean(); fit = np.polyfit(S.x, S.y, 1); A0 = np.polyval(fit, XBAR); resid = float(np.std(S.y - np.polyval(fit, S.x))); print(f"\nD. residual sd of log sd about the pooled line: {resid:.3f}")
    with Pool(26, initializer=_init, initargs=(SIMMETA, XBAR, A0)) as pool:
        tune = pd.DataFrame(pool.map(simulate, [(s, 1.33, tau) for tau in (0.0, 0.2, 0.3, 0.4, 0.5) for s in range(12)]))
        tune["tau"] = np.repeat([0.0, 0.2, 0.3, 0.4, 0.5], 12); tt = tune.groupby("tau").resid_sd.mean(); tau = float(np.interp(resid, tt.values, tt.index.values)); print("cell-level scatter needed:", tt.round(3).to_dict(), "-> tau =", round(tau, 3))
        sim = pd.DataFrame(pool.map(simulate, [(1000 + i, s, tau) for s in (1.0, 1.33, 1.5) for i in range(300)]))
    sim.to_csv(os.path.join(OUT, "c1_calibration_raw.csv"), index=False)
    cal = sim.groupby("s_true")[["pooled"] + COMP].agg(["mean", "std"]).round(3); cal.to_csv(os.path.join(OUT, "c1_calibration.csv")); print(cal.to_string())
