"""Nonlinearity of the map from each parameterisation's coordinates to the quantities the density is built from.

The prior is Gaussian in the ILR coordinates z = V' log(cells) and the penalty is on the margins, so what matters is how
far the map theta -> (z, log row sums, log column sums) departs from linear over a posterior-sized step.
At posterior draws theta, take a random direction u in the posterior-whitened space (theta = theta0 + s*L*u, |u| = 1,
cov(theta) = L L'), and compare the actual change in x with its linearisation:
    r = |x(theta0 + s L u) - x(theta0) - s J L u| / |s J L u|       (J L u by a central difference)
r = 0 means the map is linear over that step.  Reported: median and 90th percentile over draws and directions, for s = 1, 2.
Run from the repository root:  BRIDGESTAN=/path/to/bridgestan python3 experiments/07_map_curvature.py [--eps 1,0.1] [--sigma 1,2]
Writes experiments/output/07_map_curvature.csv and .md.
"""
import os, sys, argparse, importlib.util, numpy as np, pandas as pd
import bridgestan as bs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("c6", os.path.join(ROOT, "experiments", "06_curvature.py"))
c6 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c6)
r4, e2, jsonable, SO = c6.r4, c6.e2, c6.jsonable, c6.SO


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--path", default="data/scot_2007_all_areas.csv")
    ap.add_argument("--eps", default="1,0.1"); ap.add_argument("--sigma", default="1,2")
    ap.add_argument("--params", default="0,5,1,2,3,4")
    ap.add_argument("--draws", type=int, default=200); ap.add_argument("--dirs", type=int, default=20)
    ap.add_argument("--label", default="map_curvature")
    a = ap.parse_args()
    T = r4.load(os.path.join(ROOT, a.path)); R, C = T.shape; D = R * C
    rng = np.random.default_rng(1)
    rows = []
    for eps in [float(x) for x in a.eps.split(",")]:
        for sb in [float(x) for x in a.sigma.split(",")]:
            for p in [int(x) for x in a.params.split(",")]:
                sc = 0 if p in (0, 5) else 1
                d = e2.base_data(T, eps, sc)
                V = d["V"]
                d.update(param=p, mu_b=V.T @ r4.centre_table(T, "indep").ravel(), sigma_b=np.full(D - 1, sb), mu_logv=float(np.log(T.sum())), sigma_logv=1.0)
                f = e2.model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=1,
                                    inits=[{"theta": e2.start(p, d)}] * 2, show_progress=False, show_console=False)
                s = f.summary(); c = s[s.index.str.startswith("log_T[")]
                rhat = float(c["R_hat"].max())
                th = f.stan_variable("theta"); th = th.reshape(-1, th.shape[-1])
                model = bs.StanModel(SO, jsonable(d), seed=1)

                def xmap(t):
                    t = np.ascontiguousarray(t, dtype=np.float64)
                    full = model.param_constrain(t, include_tp=True)
                    Tm = full[D:2 * D].reshape((R, C), order="F")             # Stan matrices are column-major
                    lt = np.log(Tm)
                    z = V.T @ lt.ravel()                                       # cells row by row, as in the prior
                    return dict(ilr=z, margins=np.r_[np.log(Tm.sum(1)), np.log(Tm.sum(0))])

                M = np.cov(th.T); L = np.linalg.cholesky(M + 1e-12 * np.trace(M) / D * np.eye(D))
                idx = rng.choice(len(th), a.draws, replace=False)
                res = {(q, sstep): [] for q in ("ilr", "margins") for sstep in (1.0, 2.0)}
                for i in idx:
                    t0 = th[i]; x0 = xmap(t0)
                    for _ in range(a.dirs):
                        u = rng.standard_normal(D); u /= np.linalg.norm(u); dv = L @ u
                        h = 1e-4
                        xp, xm = xmap(t0 + h * dv), xmap(t0 - h * dv)
                        for sstep in (1.0, 2.0):
                            xs = xmap(t0 + sstep * dv)
                            for q in ("ilr", "margins"):
                                lin = sstep * (xp[q] - xm[q]) / (2 * h)
                                den = np.linalg.norm(lin)
                                if den > 1e-12 and np.all(np.isfinite(xs[q])):
                                    res[(q, sstep)].append(np.linalg.norm(xs[q] - x0[q] - lin) / den)
                for (q, sstep), v in res.items():
                    v = np.array(v)
                    rows.append(dict(eps=eps, sigma_b=sb, param=r4.NAMES[p], rhat=rhat, quantity=q, step=sstep,
                                     median=float(np.median(v)), p90=float(np.quantile(v, 0.9)), n=len(v)))
                print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
                pd.DataFrame(rows).to_csv(os.path.join(ROOT, "experiments", "output", f"07_{a.label}.csv"), index=False)
    df = pd.DataFrame(rows)
    wide = df[df.step == 1.0].pivot_table(index=["eps", "sigma_b", "param", "rhat"], columns="quantity", values=["median", "p90"]).round(3)
    with open(os.path.join(ROOT, "experiments", "output", f"07_{a.label}.md"), "w") as fh:
        fh.write("# Nonlinearity of the coordinate maps over one posterior-sd step (random whitened directions)\n\n"
                 "r = |actual change - linear prediction| / |linear prediction|; 0 = linear. Rows with rhat >= 1.05 have an unreliable "
                 "posterior covariance.\nilr = coordinates the prior is Gaussian in; margins = log row and column sums the penalty acts on.\n\n"
                 + wide.to_markdown() + "\n\nStep 2 sd: see the csv.\n")
    print(wide.to_string())


if __name__ == "__main__":
    main()
