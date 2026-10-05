"""Stage A of the evaluation: tables simulated from the model's own prior, so there is no prior-data conflict.

For each simulated table
  1. draw ILR coordinates ~ normal(0, sigma_b) and the log total ~ normal(log TOTAL, 1), and form the table;
  2. observe its row and column totals with the model's own margin noise,  obs = total + eps * sqrt(total) * N(0, 1);
  3. fit the same density in each parameterisation, with the same data and seed.
Because the data are generated from the model, the rank of the true log cell among the posterior draws
should be uniform (simulation-based calibration). Sampling cost is recorded for every run.

Run from the repository root:  python3 experiments/03_simulated_tables.py [n_tables] [sigma_b ...]
Writes experiments/output/03_simulated_tables.csv (one row per run) and a summary .md.
"""
import os, sys, time, importlib.util, logging, numpy as np, pandas as pd

logging.getLogger("cmdstanpy").disabled = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("e2", os.path.join(ROOT, "experiments", "02_single_area_models.py"))
e2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(e2)

R = C = 7
D = R * C
TOTAL = 20000.0
V = e2.helmert(D)
PARAMS = {0: "ILR + log volume", 5: "log cells", 1: "position logit", 2: "odds-ratio logit", 3: "adjusted row", 4: "adjusted table"}


def simulate(rng, sigma_b, eps):
    z = rng.normal(0, sigma_b, D - 1)
    lt = V @ z
    lt += rng.normal(np.log(TOTAL), 1.0) - np.log(np.exp(lt).sum())
    T = np.exp(lt).reshape(R, C)
    w, m = T.sum(1), T.sum(0)
    ow = np.maximum(w + eps * np.sqrt(w) * rng.normal(size=R), 1.0)
    om = np.maximum(m + eps * np.sqrt(m) * rng.normal(size=C), 1.0)
    return T, ow, om


def fit_one(T, ow, om, eps, sigma_b, param, seed, warm=500, draws=500):
    d = e2.base_data(T, eps, 0)
    d["obs_w"], d["obs_m"] = ow, om
    d["row_order"] = [int(i) + 1 for i in np.argsort(ow)]
    d.update(param=param, mu_b=np.zeros(D - 1), sigma_b=np.full(D - 1, sigma_b), mu_logv=float(np.log(TOTAL)), sigma_logv=1.0)
    t0 = time.time()
    f = e2.model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=warm, iter_sampling=draws, seed=seed,
                        inits=[{"theta": e2.start(param, d)}] * 2, show_progress=False, show_console=False)
    wall = time.time() - t0
    s = f.summary(); c = s[s.index.str.startswith("log_T[")]
    ess = float(c["ESS_bulk" if "ESS_bulk" in s.columns else "N_Eff"].min())
    dr = f.draws_pd()
    lT = f.stan_variable("log_T").reshape(-1, D)
    true = np.log(T).ravel()
    ranks = (lT < true).mean(0)
    return dict(min_ess=ess, max_rhat=float(c["R_hat"].max()), div=int(dr["divergent__"].sum()),
                leap=float(dr["n_leapfrog__"].mean()), seconds=wall, rank_mean=float(ranks.mean()),
                rank_sd=float(ranks.std()), _ranks=ranks)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    sigmas = [float(x) for x in sys.argv[2:]] or [1.0, 2.0]
    eps_list = [1.0, 0.1]
    rows, allranks = [], []
    for sb in sigmas:
        for eps in eps_list:
            rng = np.random.default_rng(int(sb * 1000 + eps * 100))
            for k in range(n):
                T, ow, om = simulate(rng, sb, eps)
                for p in (0, 5, 1):
                    r = fit_one(T, ow, om, eps, sb, p, seed=k + 1)
                    ranks = r.pop("_ranks")
                    r.update(sigma_b=sb, eps=eps, table=k, param=PARAMS[p])
                    rows.append(r); allranks.append((PARAMS[p], sb, eps, ranks))
                    print({k_: (round(v, 3) if isinstance(v, float) else v) for k_, v in r.items()}, flush=True)
                    pd.DataFrame(rows).to_csv(os.path.join(ROOT, "experiments", "output", "03_simulated_tables.csv"), index=False)
    df = pd.DataFrame(rows)
    df["fail"] = (df.max_rhat > 1.05) | (df["div"] > 0)
    df["ess_per_s"] = df.min_ess / df.seconds
    g = df.groupby(["sigma_b", "eps", "param"]).agg(runs=("fail", "size"), failed=("fail", "sum"),
            median_min_ess=("min_ess", "median"), median_seconds=("seconds", "median"),
            median_ess_per_s=("ess_per_s", "median"), median_rhat=("max_rhat", "median"))
    out = os.path.join(ROOT, "experiments", "output", "03_simulated_tables.md")
    with open(out, "w") as fh:
        fh.write("# Simulated tables, prior known\n\n" + g.round(2).to_markdown() + "\n")
    print(g.round(2).to_string())


if __name__ == "__main__":
    main()
