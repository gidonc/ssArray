"""Simulated designs with one feature changed at a time: is each feature reliably harder or easier for the scheme it targets?

A design is (margins, prior centre, prior sds); the model sees only these.  Knobs and the mechanism each one targets:
  N        table total (eps = 1): margin tightness.                         Target: ILR and log cells.
  sigma_b  prior width.                                                     Target: position logit, odds-ratio logit.
  tie      1 = row and column margins have the same shape (remaining row and column amounts tie, so the binding
           bound flips); 0 = staggered so they never tie.                   Target: position logit (kinks).
  kappa    log odds ratio on the diagonal of the prior centre (centre = exp(kappa * I) raked to the margins, so its
           margins agree with the data).                                    Target: position logit with a wide prior;
                                                                            log cells should gain on ILR.
  tiny     share of the smallest row margin.                                Target: everything under an interactions-only
                                                                            prior; adjusted table's solve.
  conflict centre's row and column effects moved this many prior sds away from the margins.   Target: everything, past ~15.
  prior    'equal' or 'inter' (sd 50 on row and column effects).            Target: adjusted table gains under 'inter'.
  size     R = C.                                                           Target: cell-by-cell schemes (and adjusted table's cost).
Baseline: 5x5, N 30,000, sigma_b 0.5, staggered margins with ratio 1.6, kappa 0, equal sds, no conflict.
Each arm changes one or two knobs.  Three replicates per arm: the margins are jittered (log-normal, sd 0.15) so a result
has to hold across margins of the same kind, and the sampler seed changes too.
Run:  python3 experiments/17_sim_features.py     Writes experiments/output/17_sim_features.csv and .md
Environment: MARGIN_COORDS = last (default) | largest | whitened  (see margincoords.py; output gets the suffix _<coords>),
             SEQ_ONLY = 1 to run only the four sequential schemes,
             MODEL = current for the model as it now stands (largest-column margin coordinates, centred interior,
             split basis without its matrix; output suffix _current).
"""
import os, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2
import sys; sys.path.insert(0, os.path.join(ROOT, "experiments")); import margincoords as mc
CURRENT = os.environ.get("MODEL") == "current"      # MODEL=current: largest-column margins, centred interior, split basis without its matrix
COORDS = "largest" if CURRENT else os.environ.get("MARGIN_COORDS", "last"); PARAMS = (1, 2, 3, 4) if os.environ.get("SEQ_ONLY") == "1" else (0, 5, 1, 2, 3, 4)
BASE = dict(size=5, N=3e4, sigma_b=0.5, tie=0, kappa=0.0, tiny=None, conflict=0.0, prior="equal", ratio=1.6)
ARMS = [("baseline", {}),
        ("low volume", dict(N=30.0)), ("high volume", dict(N=3e6)),
        ("wide prior", dict(sigma_b=2.0)),
        ("diagonal centre", dict(kappa=3.0)), ("diagonal centre, wide prior", dict(kappa=3.0, sigma_b=2.0)),
        ("tied margins, diagonal, wide", dict(tie=1, kappa=3.0, sigma_b=2.0)), ("tied margins, wide", dict(tie=1, sigma_b=2.0)),
        ("tiny margin", dict(tiny=1e-4)), ("tiny margin, interactions-only", dict(tiny=1e-4, prior="inter")),
        ("interactions-only, wide", dict(prior="inter", sigma_b=2.0)),
        ("conflict 5 sd", dict(conflict=5.0)), ("conflict 20 sd", dict(conflict=20.0)),
        ("3x3", dict(size=3)), ("8x8", dict(size=8))]


def ipf(x, w, m, it=3000):
    x = x.copy()
    for _ in range(it):
        x *= (w / x.sum(1))[:, None]; x *= (m / x.sum(0))[None, :]
    return x


def split_basis(R, C):
    Hr, Hc = e2.helmert(R), e2.helmert(C); one = lambda n: np.ones((n, 1)) / np.sqrt(n)
    return np.c_[np.kron(Hr, one(C)), np.kron(one(R), Hc), np.kron(Hr, Hc)]


def design(rng, size, N, sigma_b, tie, kappa, tiny, conflict, prior, ratio):
    """returns observed margins w, m; prior centre (table); prior sds in the split basis; the basis"""
    R = C = size; D = R * C; nm = R + C - 2
    w = ratio ** -np.arange(R, dtype=float); m = w.copy() if tie else ratio ** -(np.arange(C) + 0.5)
    jw = np.exp(rng.normal(0, 0.15, R)); jm = jw.copy() if tie else np.exp(rng.normal(0, 0.15, C))   # tied margins stay tied
    w, m = w * jw, m * jm
    if tiny is not None: w[-1] = tiny * w[:-1].sum() / (1 - tiny)
    w, m = N * w / w.sum(), N * m / m.sum()
    V = split_basis(R, C)
    cen = ipf(np.exp(kappa * np.eye(R, C)), w, m)
    sig = np.r_[np.full(nm, sigma_b if prior == "equal" else 50.0), np.full(D - 1 - nm, sigma_b)]
    if conflict > 0:
        u = rng.standard_normal(nm); u /= np.linalg.norm(u)
        cen = np.exp(np.log(cen) + (V[:, :nm] @ (conflict * sigma_b * u)).reshape(R, C)); cen *= N / cen.sum()
    return w, m, cen, sig, V


def flips(Tt, order, rem_col):
    """share of visited cells whose binding upper bound (row or column remainder) changes across the posterior"""
    n, R, C = Tt.shape; sr = Tt.sum(2).copy(); sc = Tt.sum(1).copy(); f = []
    for k in range(R - 1):
        r = order[k] - 1; last = rem_col[r] - 1
        for c in range(C):
            if c == last: continue
            b = (sc[:, c] < sr[:, r]).mean(); f.append(2 * min(b, 1 - b))
            x = Tt[:, r, c]; sc[:, c] -= x; sr[:, r] -= x
        sc[:, last] -= sr[:, r]
    return float(np.mean(f))


def main():
    rows = []; out = os.path.join(ROOT, "experiments", "output", "17_sim_features" + ("_current" if CURRENT else "" if COORDS == "last" else "_" + COORDS))
    for arm, ch in ARMS:
        k = dict(BASE); k.update(ch)
        for rep in (1, 2, 3):
            rng = np.random.default_rng(100 * rep + 7)
            w, m, cen, sig, V = design(rng, **k); R = C = k["size"]; D = R * C
            T = np.outer(w, m) / w.sum()
            Km, scm = mc.margin_K(w, m, 1.0, COORDS, 1)
            for p in PARAMS:
                dd = e2.base_data(T, 1.0, 0 if p in (0, 5) else scm); dd["V"] = V; dd["K_margin"] = Km
                if CURRENT: dd["centred_interior"] = 1; dd["split_basis"] = 1; dd["V"] = np.zeros((0, 0))
                dd.update(param=p, mu_b=V.T @ np.log(cen).ravel(), sigma_b=sig, mu_logv=float(np.log(k["N"])), sigma_logv=1.0)
                row = dict(arm=arm, rep=rep, param=r4.NAMES[p], min_margin=float(min(w.min(), m.min())), **{a: b for a, b in k.items() if a != "ratio"}); t0 = time.time()
                try:
                    f = e2.model.sample(data=dd, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=rep,
                                        inits=[{"theta": e2.start(p, dd)}] * 2, show_progress=False, show_console=False, timeout=200)
                    sm = f.summary(); c = sm[sm.index.str.startswith("log_T[")]; dr = f.draws_pd()
                    ess = float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()); lf = float(dr["n_leapfrog__"].mean())
                    row.update(seconds=time.time() - t0, max_rhat=float(c["R_hat"].max()), min_ess=ess, div=int(dr["divergent__"].sum()), leapfrogs=lf,
                               ess_per_1000_grad=ess / lf, ess_per_s=ess / (time.time() - t0))
                    if p == 3: row["flip"] = flips(np.exp(f.stan_variable("log_T")), dd["row_order"], dd["rem_col"])
                except Exception as ex:
                    row.update(seconds=time.time() - t0, max_rhat=np.inf, min_ess=0.0, error=type(ex).__name__)
                rows.append(row); pd.DataFrame(rows).to_csv(out + ".csv", index=False)
            print(arm, rep, "done", flush=True)


if __name__ == "__main__":
    main()
