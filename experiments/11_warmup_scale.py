"""Is the failure of the unscaled sequential versions a warm-up problem?  (Side issue for the main ssEI models.)

The tight margin directions are isolated in the margin parameters (10_isolation), so a diagonal metric can in principle
fix their scale.  Scottish grand total, eps 0.01, sigma_b 1, margins as plain log deviations (scale_margins = 0):
  default      500 warm-up iterations, unit initial metric
  long         2000 warm-up iterations
  init metric  500 warm-up, initial diagonal metric set to the known penalty scale of each margin parameter
  scaled       scale_margins = 1 (the reparameterisation), 500 warm-up
Also reports how far the adapted metric ends from the known scale for the margin parameters (ratio of sds, worst case).
Run:  python3 experiments/11_warmup_scale.py     Writes experiments/output/11_warmup_scale.csv and .md
"""
import os, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2
EPS, SB = 0.01, 1.0


def main():
    T = r4.load(os.path.join(ROOT, "data", "scot_2007_all_areas.csv")); R, C = T.shape; D = R * C; k = R + C - 1
    w, m = T.sum(1), T.sum(0)
    s = np.r_[EPS / np.sqrt(w), EPS * np.sqrt(1 / m[:C - 1] + 1 / m[C - 1])]            # penalty scale of each margin parameter
    rows = []
    for p in (1, 2):
        for name, sc, warm, init in (("default", 0, 500, False), ("long warm-up", 0, 2000, False), ("init metric", 0, 500, True), ("scaled", 1, 500, False)):
            for seed in (1, 2):
                d = e2.base_data(T, EPS, sc)
                d.update(param=p, mu_b=d["V"].T @ r4.centre_table(T, "indep").ravel(), sigma_b=np.full(D - 1, SB), mu_logv=float(np.log(T.sum())), sigma_logv=1.0)
                kw = dict(inv_metric=np.r_[s ** 2, np.ones(D - k)]) if init else {}
                t0 = time.time()
                f = e2.model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=warm, iter_sampling=500, seed=seed, metric="diag_e",
                                    inits=[{"theta": e2.start(p, d)}] * 2, show_progress=False, show_console=False, **kw)
                sm = f.summary(); c = sm[sm.index.str.startswith("log_T[")]; dr = f.draws_pd()
                im = np.asarray(f.metric)                                                # chains x D adapted inverse metric (variances)
                target = (s if sc == 0 else np.ones(k))
                off = np.exp(np.abs(np.log(np.sqrt(im[:, :k]) / target)).max())          # worst factor between adapted sd and penalty scale
                rows.append(dict(param=r4.NAMES[p], setup=name, seed=seed, max_rhat=float(c["R_hat"].max()),
                                 min_ess=float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()), div=int(dr["divergent__"].sum()),
                                 leapfrogs=float(dr["n_leapfrog__"].mean()), seconds=time.time() - t0, metric_off_by=float(off)))
                print({a: (round(b, 3) if isinstance(b, float) else b) for a, b in rows[-1].items()}, flush=True)
    df = pd.DataFrame(rows); out = os.path.join(ROOT, "experiments", "output", "11_warmup_scale")
    df.to_csv(out + ".csv", index=False)
    g = df.assign(fail=df.max_rhat > 1.05).groupby(["param", "setup"], sort=False).agg(failed=("fail", "sum"), max_rhat=("max_rhat", "max"), min_ess=("min_ess", "median"),
            leapfrogs=("leapfrogs", "median"), seconds=("seconds", "median"), metric_off_by=("metric_off_by", "max")).round(2)
    open(out + ".md", "w").write(f"# Warm-up and margin scale, Scottish grand total (eps {EPS}, sigma_b {SB})\n\n" + g.to_markdown() + "\n")
    print(g.to_string())


if __name__ == "__main__":
    main()
