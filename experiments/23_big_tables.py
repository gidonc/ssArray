"""Do the parameterisations run in reasonable time on large tables?

Square tables 10, 16, 24, 32 and one wide table 5 x 40.  6,000 per row margin on average (so the total grows with the
table), margins jittered (sd 0.15), independence centre, sigma_b 0.5, eps 1, largest-column margin coordinates.
All six parameterisations, 2 chains, 500 warm-up + 500 draws, time limit 900 s per fit.
Reported: wall time, leapfrogs per draw, time per gradient (wall time / gradients per chain, warm-up included),
min ESS of the log cells, ESS per second.
Note: every parameterisation evaluates the prior as V' log(cells) with a dense RC x (RC-1) basis, which costs RC^2 per
gradient and is part of all the times here; a structured basis would remove it.
Run:  python3 experiments/23_big_tables.py     Writes experiments/output/23_big_tables.csv and .md
"""
import os, sys, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, "experiments"))
import simdesign as sd, margincoords as mc
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2
SIZES = [(10, 10), (5, 40), (16, 16), (24, 24), (32, 32)]; LIMIT = 900


def main():
    rows = []; out = os.path.join(ROOT, "experiments", "output", "23_big_tables")
    for R, C in SIZES:
        d = sd.make_design(np.random.default_rng(R * 100 + C), R=R, C=C, volume=("margin", 6000), jitter=0.15); w, m = d["w"], d["m"]; D = R * C
        T = np.outer(w, m) / d["N"]; K, sc = mc.margin_K(w, m, 1.0, "largest", 1)
        for p in (2, 1, 3, 0, 5, 4):
            dd = e2.base_data(T, 1.0, 0 if p in (0, 5) else sc); dd["V"] = d["V"]; dd["K_margin"] = K
            dd.update(param=p, mu_b=d["V"].T @ np.log(T).ravel(), sigma_b=d["sigma"], mu_logv=float(np.log(d["N"])), sigma_logv=1.0)
            row = dict(R=R, C=C, cells=D, total=d["N"], param=r4.NAMES[p]); t0 = time.time()
            try:
                fit = e2.model.sample(data=dd, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=1, save_warmup=True,
                                      inits=[{"theta": e2.start(p, dd)}] * 2, show_progress=False, show_console=False, timeout=LIMIT)
                wall = time.time() - t0
                dr = fit.draws_pd(inc_warmup=True); post = dr[dr["iter__"] > 500] if "iter__" in dr else dr.groupby("chain__").tail(500)
                cols = [c for c in dr.columns if c.startswith("log_T[")]
                import arviz as az
                x = np.stack([post[post["chain__"] == ch][cols].values for ch in sorted(post["chain__"].unique())])
                ess = float(np.min([az.ess(x[:, :, j], method="bulk") for j in range(0, x.shape[2], max(1, x.shape[2] // 200))]))      # up to ~200 cells
                rh = float(np.max([az.rhat(x[:, :, j]) for j in range(0, x.shape[2], max(1, x.shape[2] // 200))]))
                grads = dr.groupby("chain__")["n_leapfrog__"].sum().mean()
                row.update(seconds=wall, leapfrogs=float(post["n_leapfrog__"].mean()), ms_per_gradient=1000 * wall / grads, min_ess=ess, max_rhat=rh,
                           ess_per_s=ess / wall, div=int(post["divergent__"].sum()), treedepth_hits=float((post["treedepth__"] >= 10).mean()))
            except Exception as ex:
                row.update(seconds=time.time() - t0, error=type(ex).__name__)
            rows.append(row); pd.DataFrame(rows).to_csv(out + ".csv", index=False)
            print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    df = pd.DataFrame(rows)
    open(out + ".md", "w").write("# Large tables: 2 chains, 500 + 500 draws\n\n" + df.round(2).to_markdown(index=False) + "\n")


if __name__ == "__main__":
    main()
