"""
b2_family_comparison.py — which interval density describes the HEK293 segments?  (gamma / inverse
Gaussian / log-normal / reduced first-passage density)

Data: the 196 segments with >= 10 spikes; intervals detrended inside the segment (as in the manuscript); intervals
below 0.3 of the segment mean (doublets) dropped in the primary analysis.  Each family has two free parameters per
segment (mean and CV of the shifted interval) and a fixed common offset delta (0, or the model's spike-and-early-recovery
offset 16.4 model units = 8.64 s).  Likelihood accounts for the frame step of each record (isifam.loglik).
Parts:  A  total log-likelihood per family, differences with segment- and record-bootstrap intervals, by dose
        B  profile over a common offset delta
        C  sensitivity: all intervals with a 1 % uniform contamination component; no detrending; two burst segments removed
        D  model-as-data: full model (Jacobi, chemical Langevin) cut into the data segment lengths
        E  parametric bootstrap: data regenerated from each fitted family, same lengths and frame steps -> distribution of
           the log-likelihood differences if that family were true
Output: results/b2_*.csv
"""
import os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "round2_checks"))
import isifam as F
OUT = os.path.join(HERE, "results"); CODE = os.path.join(HERE, "..")
T0 = 0.527; DELTA = 16.4 * T0
FAMS = [("gamma", 0.0), ("invgauss", 0.0), ("lognorm", 0.0), ("gamma", DELTA), ("invgauss", DELTA), ("lognorm", DELTA), ("fpt_jacobi", DELTA), ("fpt_cle", DELTA)]
name = lambda fam, d: fam + ("+off" if d > 0 else "")
NAMES = [name(*f) for f in FAMS]

def quantum(x):
    for q in (4.0, 2.0, 1.0, 0.5, 0.25):
        if np.mean(np.abs(x / q - np.round(x / q)) * q < 0.02) > 0.97: return q
    return 0.0
def load_segments(detr=True, trim=0.3):
    E = pd.read_csv(os.path.join(CODE, "data_dose_units", "results", "events_canonical.csv")); E["segment"] = E.file.str.replace(r"\.csv$", "", regex=True)
    S = pd.read_csv(os.path.join(CODE, "data_dose_units", "results", "segment_table.csv")); keep = set(zip(S[S.n >= 10].train, S[S.n >= 10].segment)); segs = []
    for (tr, sg), d in E.groupby(["train", "segment"], sort=False):
        if (tr, sg) not in keep: continue
        d = d.sort_values("time"); x = d.time.diff().values[1:]; t = d.time.values[1:]; q = quantum(x)
        xd = F.detrend(x, t) if detr else x; r = xd / xd.mean(); ok = r >= trim if trim else np.ones(len(x), bool)
        segs.append(dict(id=f"{tr}|{sg}", record=d.record.iloc[0], dose=int(d.dose.iloc[0]), q=q, x=xd[ok], n=int(ok.sum()), ndrop=int((~ok).sum()), burst=bool((r < 0.3).sum() >= 5)))
    return segs
def fit_segment(a):
    x, q, eps, fams = a
    return [F.fit(fam, x, q, d, eps) for fam, d in fams]
def fit_all(pool, segs, eps=0.0, fams=FAMS):
    res = pool.map(fit_segment, [(s["x"], s["q"], eps, fams) for s in segs], chunksize=4)
    L = np.array([[r[0] for r in rs] for rs in res]); M = np.array([[r[1] for r in rs] for rs in res]); C = np.array([[r[2] for r in rs] for rs in res])
    return L, M, C
def summary(L, segs, names, ref, nb=4000, seed=0):
    rng = np.random.default_rng(seed); n = np.array([s["n"] for s in segs]); rec = np.array([s["record"] for s in segs]); urec = np.unique(rec)
    tot = L.sum(0); iref = names.index(ref); D = L - L[:, [iref]]; rows = []
    idx = rng.integers(0, len(segs), (nb, len(segs))); bs = D[idx].sum(1)                                   # segment bootstrap
    byrec = np.array([D[rec == r].sum(0) for r in urec]); br = byrec[rng.integers(0, len(urec), (nb, len(urec)))].sum(1)   # record (cluster) bootstrap
    best = L.argmax(1)
    for j, nm in enumerate(names):
        rows.append(dict(family=nm, logL=tot[j], dlogL_vs_best=tot[j] - tot.max(), per_interval=tot[j] / n.sum(), dlogL_vs_ref=D[:, j].sum(),
                         seg_lo=np.quantile(bs[:, j], .025), seg_hi=np.quantile(bs[:, j], .975), rec_lo=np.quantile(br[:, j], .025), rec_hi=np.quantile(br[:, j], .975),
                         n_best=int((best == j).sum()), frac_seg_better_than_ref=float(np.mean(D[:, j] > 0))))
    return pd.DataFrame(rows)

def synth_segment(rng, fam, d, m, cv, n, q):
    """n intervals from the fitted family, spike times rounded to the frame grid, detrended like the data"""
    x = F.sample(fam, rng, n, m, cv) + d; t = np.cumsum(x)
    if q > 0: tt = np.round(np.concatenate([[0.0], t]) / q) * q; x = np.diff(tt); t = tt[1:]
    x = F.detrend(x, t); return x[x / x.mean() >= 0.3]
def boot_job(a):
    seed, fam_i, pars, meta = a; rng = np.random.default_rng(seed); fam, d = FAMS[fam_i]; tot = np.zeros(len(FAMS))
    for (m, cv), (n, q) in zip(pars, meta):
        x = synth_segment(rng, fam, d, m, cv, n, q)
        tot += np.array([F.fit(f2, x, q, d2)[0] for f2, d2 in FAMS])
    return fam_i, tot

if __name__ == "__main__":
    pd.set_option("display.width", 250); fmt = lambda v: f"{v:.2f}"
    with Pool(26) as pool:
        segs = load_segments(); ntot = sum(s["n"] for s in segs)
        print(f"segments {len(segs)}, intervals {ntot}, dropped doublets {sum(s['ndrop'] for s in segs)}; frame step counts {pd.Series([s['q'] for s in segs]).value_counts().to_dict()}")
        # ---- A
        L, M, C = fit_all(pool, segs); np.save(os.path.join(OUT, "b2_logL_primary.npy"), L)
        A = summary(L, segs, NAMES, "fpt_jacobi+off"); A.to_csv(os.path.join(OUT, "b2_A_primary.csv"), index=False); print("\nA. primary (detrended, doublets dropped); reference = fpt_jacobi+off"); print(A.to_string(index=False, float_format=fmt))
        P = pd.DataFrame(dict(id=[s["id"] for s in segs], record=[s["record"] for s in segs], dose=[s["dose"] for s in segs], n=[s["n"] for s in segs], q=[s["q"] for s in segs]))
        for j, nm in enumerate(NAMES): P["logL_" + nm] = L[:, j]; P["m_" + nm] = M[:, j]; P["cv_" + nm] = C[:, j]
        P["nu_jacobi"] = [F.fpt("jacobi").nu_of_cv(c) for c in P["cv_fpt_jacobi+off"]]; P["nu_cle"] = [F.fpt("cle").nu_of_cv(c) for c in P["cv_fpt_cle+off"]]
        P.to_csv(os.path.join(OUT, "b2_A_segments.csv"), index=False)
        bd = P.groupby("dose").apply(lambda g: pd.Series({nm: g["logL_" + nm].sum() - g["logL_fpt_jacobi+off"].sum() for nm in NAMES} | {"intervals": g.n.sum()}), include_groups=False)
        bd.to_csv(os.path.join(OUT, "b2_A_by_dose.csv")); print("\nby dose: logL - logL(fpt_jacobi+off)"); print(bd.round(1).to_string())
        print("\nfitted nu (Jacobi form): median %.3f, IQR %.3f-%.3f;  implied N = 1/nu^2 (chemical Langevin form): median %.0f, IQR %.0f-%.0f" % (
            P.nu_jacobi.median(), P.nu_jacobi.quantile(.25), P.nu_jacobi.quantile(.75), (1 / P.nu_cle**2).median(), (1 / P.nu_cle**2).quantile(.25), (1 / P.nu_cle**2).quantile(.75)))
        # ---- B
        rows = []
        for d in (0.0, 2.0, 4.0, 6.0, 8.0, DELTA, 10.0, 11.0, 12.0):
            fams = [(f, d) for f in ("gamma", "invgauss", "lognorm", "fpt_jacobi", "fpt_cle")]; Ld, _, _ = fit_all(pool, segs, fams=fams)
            rows.append(dict(delta=d, **{f: Ld[:, j].sum() for j, (f, _) in enumerate(fams)}))
        B = pd.DataFrame(rows); B.to_csv(os.path.join(OUT, "b2_B_offset_profile.csv"), index=False); print("\nB. total logL against a common offset delta (s)"); print(B.to_string(index=False, float_format=fmt))
        # ---- C
        out = []
        for lab, sg, eps in (("all intervals + 1% contamination", load_segments(trim=0), 0.01), ("no detrending, doublets dropped", load_segments(detr=False), 0.0),
                             ("two burst segments removed", [s for s in segs if not s["burst"]], 0.0)):
            Lc, _, _ = fit_all(pool, sg, eps=eps); c = summary(Lc, sg, NAMES, "fpt_jacobi+off"); c.insert(0, "variant", lab); out.append(c)
        Cc = pd.concat(out); Cc.to_csv(os.path.join(OUT, "b2_C_sensitivity.csv"), index=False); print("\nC. sensitivity"); print(Cc[["variant", "family", "logL", "dlogL_vs_best", "dlogL_vs_ref", "seg_lo", "seg_hi", "rec_lo", "rec_hi"]].to_string(index=False, float_format=fmt))
        # ---- D
        import common as cm
        BIN = os.path.join(HERE, "m1_noiseform_sde"); star = pd.read_csv(os.path.join(OUT, "a3_matched_segmented.csv")).set_index("form").amp
        out = []
        for form, key, extra in (("jacobi", "sigH", dict(hnoise=0)), ("cle", "Nrec", dict(hnoise=2))):
            msegs = []
            for seed in range(1, 7):
                lst = cm.run(BIN, {key: float(star[form]), "seed": seed, **extra}, "isi=1"); x, t = lst.ISI.values * T0, lst.t.values * T0; rng = np.random.default_rng(seed); i = 0
                for n_ in rng.permutation([s["n"] + s["ndrop"] for s in segs]):
                    if i + n_ > len(x): break
                    xd = F.detrend(x[i:i + n_], t[i:i + n_]); xd = xd[xd / xd.mean() >= 0.3]; msegs.append(dict(record=f"seed{seed}", n=len(xd), x=xd, q=0.0)); i += n_
            Lm, _, _ = fit_all(pool, msegs); d = summary(Lm, msegs, NAMES, "fpt_jacobi+off"); d.insert(0, "truth", f"full model, {form} ({len(msegs)} segments, {sum(s['n'] for s in msegs)} intervals)"); out.append(d)
        Dd = pd.concat(out); Dd.to_csv(os.path.join(OUT, "b2_D_model_as_data.csv"), index=False); print("\nD. model-as-data (full model at the CV-matched amplitude, reference geometry, T0 = 0.527 s)")
        print(Dd[["truth", "family", "logL", "dlogL_vs_best", "dlogL_vs_ref", "seg_lo", "seg_hi", "per_interval"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        # ---- E
        NREP = int(os.environ.get("NREP", 60)); meta = [(s["n"] + s["ndrop"], s["q"]) for s in segs]; jobs = []
        for fi in range(len(FAMS)):
            pars = list(zip(M[:, fi], C[:, fi]))
            jobs += [(1000 * fi + r, fi, pars, meta) for r in range(NREP)]
        res = pool.map(boot_job, jobs, chunksize=1)
        R = pd.DataFrame([dict(truth=NAMES[fi], **dict(zip(NAMES, tot))) for fi, tot in res]); R.to_csv(os.path.join(OUT, "b2_E_parametric_bootstrap_raw.csv"), index=False)
        obs = dict(zip(NAMES, L.sum(0))); rows = []
        for tr in NAMES:
            g = R[R.truth == tr]
            for alt in NAMES:
                if alt == tr: continue
                dsim = (g[tr] - g[alt]).values; dobs = obs[tr] - obs[alt]
                rows.append(dict(truth=tr, alt=alt, dlogL_obs=dobs, sim_mean=dsim.mean(), sim_lo=np.quantile(dsim, .025), sim_hi=np.quantile(dsim, .975), p_sim_le_obs=float(np.mean(dsim <= dobs))))
        Eb = pd.DataFrame(rows); Eb.to_csv(os.path.join(OUT, "b2_E_parametric_bootstrap.csv"), index=False)
        print(f"\nE. parametric bootstrap ({NREP} replicates per truth): logL(truth) - logL(alt); observed value vs its distribution if 'truth' generated the data")
        print(Eb.to_string(index=False, float_format=fmt))
