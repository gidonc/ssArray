"""Where does each parameterisation have the advantage?  Designed margins: table size x margin family x volume x prior width.

The model sees only the margins, so a design is (R, C, row margins, column margins, total N); eps = 1 (count-like noise),
so volume sets the tightness tau = n_max * sigma_b^2 / eps^2.  Prior centred on the independence table.

Predictions written before the run:
  P1  Direct (ILR, log cells) converges at low volume and fails at high volume, at every size and family.
  P2  At low volume direct is at least as good as sequential per gradient, more so in larger tables (margins are a
      shrinking share, (R+C-1)/RC, of the directions).
  P3  The nonlinearity of the map to the ILR coordinates grows with table size for the cell-by-cell schemes (position
      logit, odds-ratio logit: a cell depends on every earlier cell), more slowly for adjusted row (one solve per row),
      least for adjusted table (one joint solve, no ordering).
  P4  Position logit loses most against odds-ratio logit where the active bound flips across the posterior
      (kinks at sc = sr and sr = rest); they are close where the active bounds are stable.
  P5  Adjusted table has the best ESS per gradient among sequential schemes in large tables but the worst cost per gradient.

Per fit: rhat, min ESS of log cells, divergences, leapfrogs, seconds, ESS per gradient, and (sequential only)
  r_ilr      nonlinearity of theta -> ILR over one posterior sd (as in 07), median and p90
and from the posterior tables, by replaying the between-bounds fill:
  flip_up    mean over visited cells of 2*min(f, 1-f), f = share of draws where the upper bound is the column (sc < sr)
  flip_lo    same for the lower bound being active (sr > rest)
  lo_active  share of (cell, draw) with an active lower bound
  rel_width  median (upper - lower) / upper
  near       share of (cell, draw) within 2% of either bound
Run:  BRIDGESTAN=/path/to/bridgestan python3 experiments/08_regimes.py
Writes experiments/output/08_regimes.csv
"""
import os, sys, time, argparse, importlib.util, numpy as np, pandas as pd
import bridgestan as bs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("c6", os.path.join(ROOT, "experiments", "06_curvature.py"))
c6 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c6)
r4, e2, jsonable, SO = c6.r4, c6.e2, c6.jsonable, c6.SO


def margins(family, R, C):
    g = lambda n, q: q ** -np.arange(n, dtype=float)
    if family == "uniform":   w, m = np.ones(R), np.ones(C)
    elif family == "geometric": w, m = g(R, 2), g(C, 2)             # both margins spread, same shape
    elif family == "dominant":  w, m = np.r_[1.5 * (R - 1), np.ones(R - 1)], np.r_[1.5 * (C - 1), np.ones(C - 1)]   # 60% in one row and one column
    elif family == "cross":     w, m = np.ones(R), g(C, 3)           # equal rows against very unequal columns
    else: raise ValueError(family)
    return w / w.sum(), m / m.sum()


def replay(Tt, order, rem_col):
    """between-bounds fill replayed on posterior tables (draws x R x C); order, rem_col are 1-based as in the Stan data"""
    n, R, C = Tt.shape
    sr = Tt.sum(2).copy(); sc = Tt.sum(1).copy()
    upcol, loact, relw, pos = [], [], [], []
    for k in range(R - 1):
        r = order[k] - 1; last = rem_col[r] - 1
        rest = sc.sum(1)
        for c in range(C):
            if c == last: continue
            rest = rest - sc[:, c]
            lo = np.maximum(0, sr[:, r] - rest); up = np.minimum(sc[:, c], sr[:, r])
            x = Tt[:, r, c]
            upcol.append(sc[:, c] < sr[:, r]); loact.append(sr[:, r] > rest)
            relw.append((up - lo) / up); pos.append((x - lo) / (up - lo))
            sc[:, c] -= x; sr[:, r] -= x
        sc[:, last] -= sr[:, r]
    upcol, loact, relw, pos = map(np.array, (upcol, loact, relw, pos))
    fl = lambda b: float(np.mean(2 * np.minimum(b.mean(1), 1 - b.mean(1))))
    return dict(flip_up=fl(upcol), flip_lo=fl(loact), lo_active=float(loact.mean()), rel_width=float(np.median(relw)),
                near=float(((pos < 0.02) | (pos > 0.98)).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="2x2,4x4,3x8,8x8"); ap.add_argument("--families", default="uniform,geometric,dominant,cross")
    ap.add_argument("--N", default="1e3,1e6"); ap.add_argument("--sigma", default="1,2")
    ap.add_argument("--params", default="0,5,1,2,3,4"); ap.add_argument("--label", default="regimes")
    ap.add_argument("--timeout", type=float, default=150)
    a = ap.parse_args()
    out = os.path.join(ROOT, "experiments", "output", f"08_{a.label}.csv")
    rows = []; rng = np.random.default_rng(1)
    for size in a.sizes.split(","):
        R, C = [int(x) for x in size.split("x")]; D = R * C
        for fam in a.families.split(","):
            w, m = margins(fam, R, C)
            for N in [float(x) for x in a.N.split(",")]:
                T = N * np.outer(w, m)
                for sb in [float(x) for x in a.sigma.split(",")]:
                    for p in [int(x) for x in a.params.split(",")]:
                        d = e2.base_data(T, 1.0, 0 if p in (0, 5) else 1)
                        V = d["V"]
                        d.update(param=p, mu_b=V.T @ np.log(T).ravel(), sigma_b=np.full(D - 1, sb), mu_logv=float(np.log(N)), sigma_logv=1.0)
                        row = dict(size=size, R=R, C=C, family=fam, N=N, sigma_b=sb, param=r4.NAMES[p],
                                   tau=N * max(w.max(), m.max()) * sb ** 2, stiff_share=(R + C - 1) / D)
                        t0 = time.time()
                        try:
                            f = e2.model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=1,
                                                inits=[{"theta": e2.start(p, d)}] * 2, show_progress=False, show_console=False, timeout=a.timeout)
                        except Exception as ex:
                            row.update(seconds=time.time() - t0, max_rhat=np.inf, min_ess=0.0, error=type(ex).__name__)
                            rows.append(row); print(row, flush=True); pd.DataFrame(rows).to_csv(out, index=False); continue
                        row["seconds"] = time.time() - t0
                        s = f.summary(); c = s[s.index.str.startswith("log_T[")]; dr = f.draws_pd()
                        ess = float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()); lf = float(dr["n_leapfrog__"].mean())
                        row.update(max_rhat=float(c["R_hat"].max()), min_ess=ess, div=int(dr["divergent__"].sum()), leapfrogs=lf,
                                   ess_per_1000_grad=ess / lf, ess_per_s=ess / row["seconds"])
                        row.update(replay(np.exp(f.stan_variable("log_T")), d["row_order"], d["rem_col"]))
                        if p in (1, 2, 3, 4):
                            th = f.stan_variable("theta"); th = th.reshape(-1, th.shape[-1])
                            model = bs.StanModel(SO, jsonable(d), seed=1)
                            def z(t):
                                full = model.param_constrain(np.ascontiguousarray(t, dtype=np.float64), include_tp=True)
                                return V.T @ np.log(full[D:2 * D].reshape((R, C), order="F")).ravel()
                            M = np.atleast_2d(np.cov(th.T)); L = np.linalg.cholesky(M + 1e-10 * np.trace(M) / D * np.eye(D))
                            rs = []
                            for i in rng.choice(len(th), 30, replace=False):
                                z0 = z(th[i])
                                for _ in range(4):
                                    u = rng.standard_normal(D); dv = L @ (u / np.linalg.norm(u))
                                    lin = (z(th[i] + 1e-4 * dv) - z(th[i] - 1e-4 * dv)) / 2e-4
                                    zs = z(th[i] + dv)
                                    if np.all(np.isfinite(zs)) and np.linalg.norm(lin) > 1e-12:
                                        rs.append(np.linalg.norm(zs - z0 - lin) / np.linalg.norm(lin))
                            row.update(r_ilr=float(np.median(rs)), r_ilr_p90=float(np.quantile(rs, 0.9)))
                            del model
                        rows.append(row)
                        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
                        pd.DataFrame(rows).to_csv(out, index=False)


if __name__ == "__main__":
    main()
