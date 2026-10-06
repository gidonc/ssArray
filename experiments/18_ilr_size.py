"""Why do the direct parameterisations look better in larger tables?

Two candidate reasons, separated here:
  (1) per-margin volume: at a fixed table total, a larger table has smaller margins, so each margin is looser.
      Check: hold the total fixed, or hold each margin's count fixed (total grows with R).
  (2) dilution: a log cell is mostly made of loose, prior-only directions when the table is large (the margins are
      R+C-1 of RC directions), so the cells mix well even if the margin directions do not.
      Check: ESS of the log cells against ESS of the log margins.
Equal margins (jittered), independence centre, sigma_b 0.5, eps 1, sizes 3, 5, 8, 10.  ILR, log cells, odds-ratio logit
(reference).  Three replicates.  ESS by arviz (bulk), minimum over cells and over margins.
Run:  python3 experiments/18_ilr_size.py     Writes experiments/output/18_ilr_size.csv and .md
"""
import os, time, importlib.util, numpy as np, pandas as pd, arviz as az
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("s17", os.path.join(ROOT, "experiments", "17_sim_features.py"))
s17 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s17); e2, r4 = s17.e2, s17.r4


def ess_min(x):                                     # x: chains x draws x k
    return float(np.min([az.ess(x[:, :, j], method="bulk") for j in range(x.shape[2])]))


def main():
    rows = []; out = os.path.join(ROOT, "experiments", "output", "18_ilr_size")
    for scheme, tot in (("total fixed at 30,000", lambda R: 3e4), ("each margin fixed at 6,000", lambda R: 6e3 * R)):
        for R in (3, 5, 8, 10):
            for rep in (1, 2, 3):
                rng = np.random.default_rng(10 * rep + R); N = tot(R); D = R * R
                w = np.exp(rng.normal(0, 0.15, R)); m = np.exp(rng.normal(0, 0.15, R)); w, m = N * w / w.sum(), N * m / m.sum()
                T = np.outer(w, m) / N; V = s17.split_basis(R, R)
                for p in (0, 5, 2):
                    dd = e2.base_data(T, 1.0, 0 if p in (0, 5) else 1); dd["V"] = V
                    dd.update(param=p, mu_b=V.T @ np.log(T).ravel(), sigma_b=np.full(D - 1, 0.5), mu_logv=float(np.log(N)), sigma_logv=1.0)
                    t0 = time.time()
                    f = e2.model.sample(data=dd, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=rep,
                                        inits=[{"theta": e2.start(p, dd)}] * 2, show_progress=False, show_console=False)
                    sec = time.time() - t0
                    lt = f.stan_variable("log_T").reshape(2, 500, R, R) if False else f.draws_pd()
                    cols = [c for c in lt.columns if c.startswith("log_T[")]
                    a = lt[cols].values.reshape(2, 500, R, R, order="F") if False else None
                    x = np.stack([lt[lt["chain__"] == ch][cols].values for ch in (1.0, 2.0)])            # chains x draws x cells (column-major names)
                    Tt = np.exp(x).reshape(2, 500, R, R, order="F")                                       # log_T[r,c] columns run r fastest
                    marg = np.log(np.concatenate([Tt.sum(3), Tt.sum(2)], axis=2))
                    sd_cell = float(np.median(x.reshape(-1, D).std(0))); sd_marg = float(np.median(marg.reshape(-1, 2 * R).std(0)))
                    lf = float(lt["n_leapfrog__"].mean())
                    rows.append(dict(scheme=scheme, size=R, rep=rep, param=r4.NAMES[p], N=N, margin=N / R, ess_cells=ess_min(x), ess_margins=ess_min(marg),
                                     leapfrogs=lf, step_size=float(np.mean(f.step_size)), sd_log_cell=sd_cell, sd_log_margin=sd_marg, seconds=sec,
                                     treedepth_hits=float((lt["treedepth__"] >= 10).mean())))
                    pd.DataFrame(rows).to_csv(out + ".csv", index=False)
                print(scheme, R, rep, "done", flush=True)
    df = pd.DataFrame(rows)
    for v in ("ess_cells", "ess_margins"): df[v + "_per_1000_grad"] = df[v] / df.leapfrogs
    g = df.groupby(["scheme", "param", "size"], sort=False)[["ess_cells", "ess_margins", "leapfrogs", "step_size", "ess_cells_per_1000_grad", "ess_margins_per_1000_grad", "sd_log_cell", "sd_log_margin"]].median()
    open(out + ".md", "w").write("# Direct parameterisations and table size (medians of three replicates)\n\n" + g.round(4).to_markdown() + "\n")
    pd.set_option("display.width", 250); print(g.round(4).to_string())


if __name__ == "__main__":
    main()
