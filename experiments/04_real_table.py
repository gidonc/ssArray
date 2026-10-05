"""Sampling performance on a real table, across margin noise (eps), prior width, volume and parameterisation.

The table's own margins are the data. The prior on the ILR coordinates is centred on the INDEPENDENCE table implied by
the margins (so the prior does not fight the margins), with spread sigma_b. The log total has prior sd 1.
All parameterisations target the same density; only the coordinates differ.

Input: a CSV with columns row_no, col_no, votes (an `area` column, if present, is summed over: this gives the grand total).
Run from the repository root:
  python3 experiments/04_real_table.py data/scot_area60.csv area60 [--volumes 1,100] [--eps 1,0.1,0.01]
         [--sigma 1,2] [--params 0,5,1,2,3,4] [--seeds 1,2]
Writes experiments/output/04_<label>.csv and .md
"""
import os, sys, argparse, importlib.util, logging, time, numpy as np, pandas as pd

logging.getLogger("cmdstanpy").disabled = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("e2", os.path.join(ROOT, "experiments", "02_single_area_models.py"))
e2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(e2)
NAMES = {0: "ILR + log volume", 5: "log cells", 1: "position logit", 2: "odds-ratio logit", 3: "adjusted row", 4: "adjusted table"}


def load(path):
    df = pd.read_csv(path)
    g = df.groupby(["row_no", "col_no"])["votes"].sum().reset_index()
    R, C = int(g.row_no.max()), int(g.col_no.max())
    T = np.zeros((R, C))
    for r, c, v in g.itertuples(index=False):
        T[int(r) - 1, int(c) - 1] = v
    return T


def difficulty(w, m):
    """mean over cells of the Frechet range as a share of min(w_r, m_c): small = nearly forced cells"""
    T = w.sum()
    lo = np.maximum(0, w[:, None] + m[None, :] - T)
    hi = np.minimum(w[:, None], m[None, :])
    return float(np.mean((hi - lo) / hi))


def centre_table(T, centre):
    w, m = T.sum(1), T.sum(0); R, C = T.shape
    if centre == 'indep':
        return np.log(np.outer(w, m))
    if centre == 'reverse':                                         # independence table from reversed margins: conflicts with both margins
        l = np.log(np.outer(w[::-1], m[::-1]))
        return l + (np.log(T.sum()) - np.log(np.exp(l).sum()))
    if centre == 'diag':                                            # mass pushed on the diagonal, which the margins cannot support
        l = np.log(np.outer(w, m)); l = l + 3 * np.eye(R, C)
        return l + (np.log(T.sum()) - np.log(np.exp(l).sum()))
    raise ValueError(centre)


def run(T, eps, sigma_b, param, seed, warm=500, draws=500, scale=0, centre='indep'):
    R, C = T.shape; D = R * C
    d = e2.base_data(T, eps, scale)
    V = d["V"]
    lt0 = centre_table(T, centre).ravel()                       # prior centre, cells row by row
    d.update(param=param, mu_b=V.T @ lt0, sigma_b=np.full(D - 1, sigma_b), mu_logv=float(np.log(T.sum())), sigma_logv=1.0)
    t0 = time.time()
    f = e2.model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=warm, iter_sampling=draws, seed=seed,
                        inits=[{"theta": e2.start(param, d)}] * 2, show_progress=False, show_console=False)
    wall = time.time() - t0
    s = f.summary(); c = s[s.index.str.startswith("log_T[")]
    dr = f.draws_pd()
    return dict(min_ess=float(c["ESS_bulk" if "ESS_bulk" in s.columns else "N_Eff"].min()), max_rhat=float(c["R_hat"].max()),
                div=int(dr["divergent__"].sum()), leapfrogs=float(dr["n_leapfrog__"].mean()), seconds=wall)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path"); ap.add_argument("label")
    ap.add_argument("--volumes", default="1"); ap.add_argument("--eps", default="1,0.1,0.01")
    ap.add_argument("--sigma", default="1,2"); ap.add_argument("--params", default="0,5,1,2,3,4"); ap.add_argument("--seeds", default="1,2")
    ap.add_argument("--centre", default="indep", help="prior centre: indep, reverse, diag (comma list ok)")
    ap.add_argument("--scale", default="0", help="scale_margins for the sequential versions: 0, 1 or 0,1")
    a = ap.parse_args()
    f = lambda s: [float(x) for x in s.split(",")]
    T0 = load(os.path.join(ROOT, a.path) if not os.path.isabs(a.path) else a.path)
    print(f"table {T0.shape}, total {T0.sum():.0f}, row margins {T0.sum(1).round()}, col margins {T0.sum(0).round()}", flush=True)
    out_csv = os.path.join(ROOT, "experiments", "output", f"04_{a.label}.csv")
    rows = []
    for vol in f(a.volumes):
        T = T0 * vol
        w, m = T.sum(1), T.sum(0)
        diff = difficulty(w, m)
        for eps in f(a.eps):
            for sb in f(a.sigma):
                tau = max(w.max(), m.max()) * sb ** 2 / eps ** 2
                for cen in a.centre.split(","):
                 for p in [int(x) for x in a.params.split(",")]:
                  for sc in [int(x) for x in a.scale.split(",")]:
                    if sc == 1 and p in (0, 5):
                        continue                                    # the switch does nothing for the direct versions
                    for seed in [int(x) for x in a.seeds.split(",")]:
                        r = run(T, eps, sb, p, seed, scale=sc, centre=cen)
                        r.update(label=a.label, volume=vol, total=T.sum(), eps=eps, sigma_b=sb, tau=tau, frechet_width=diff,
                                 param=NAMES[p], scale_margins=sc, seed=seed, centre=cen)
                        rows.append(r)
                        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k in
                               ("volume", "eps", "sigma_b", "param", "scale_margins", "centre", "seed", "min_ess", "max_rhat", "div", "seconds")}, flush=True)
                        pd.DataFrame(rows).to_csv(out_csv, index=False)
    df = pd.DataFrame(rows)
    df["fail"] = df.max_rhat > 1.05
    g = df.groupby(["centre", "volume", "eps", "sigma_b", "param", "scale_margins"]).agg(tau=("tau", "first"), runs=("fail", "size"), failed=("fail", "sum"),
            median_min_ess=("min_ess", "median"), median_seconds=("seconds", "median"), max_rhat=("max_rhat", "max"))
    g["ess_per_s_if_ok"] = df[~df.fail].groupby(["centre", "volume", "eps", "sigma_b", "param", "scale_margins"]).apply(lambda x: (x.min_ess / x.seconds).median())
    with open(os.path.join(ROOT, "experiments", "output", f"04_{a.label}.md"), "w") as fh:
        fh.write(f"# Real table: {a.label}\n\nfailed = max rhat > 1.05 (of runs); times are wall-clock seconds for 2 chains.\n\n"
                 + g.round(2).to_markdown() + "\n")
    print(g.round(2).to_string())


if __name__ == "__main__":
    main()
