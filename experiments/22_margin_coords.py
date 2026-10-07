"""The three margin coordinate systems (last-column reference, largest-column reference, whitened) on the designs where
the reference column mattered (21) and on the size and balance designs (19).  Sequential schemes only: the direct
parameterisations have no margin parameters.
Part A (as 21): 5x5, both margins geometric at balance 0.5 and 0.35, the column of size rank k moved to the last position.
Part B (as 19): 3x3, 5x5, 8x8 with total / per-margin / per-cell count fixed; unbalanced rows; both sides unbalanced.
Total 30,000 at 5x5, independence centre, sigma_b 0.5, eps 1.  Replicates: margin jitter 0.1 and a new seed.
Run:  python3 experiments/22_margin_coords.py     Writes experiments/output/22_margin_coords.csv and .md
"""
import os, sys, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, "experiments"))
import simdesign as sd, margincoords as mc
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2
spec = importlib.util.spec_from_file_location("e19", os.path.join(ROOT, "experiments", "19_size_margin_dials.py"))
e19 = importlib.util.module_from_spec(spec); spec.loader.exec_module(e19)
rows = []; OUT = os.path.join(ROOT, "experiments", "output", "22_margin_coords")


def fit_all(part, design, rep, w, m, d, extra):
    R, C = len(w), len(m); T = np.outer(w, m) / w.sum(); N = w.sum()
    for coords in ("last", "largest", "whitened"):
        K, sc = mc.margin_K(w, m, 1.0, coords, 1)
        for p in (1, 2, 3, 4):
            dd = e2.base_data(T, 1.0, sc); dd["V"] = d["V"]; dd["K_margin"] = K
            dd.update(param=p, mu_b=d["V"].T @ np.log(T).ravel(), sigma_b=d["sigma"], mu_logv=float(np.log(N)), sigma_logv=1.0)
            t0 = time.time()
            fit = e2.model.sample(data=dd, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=rep,
                                  inits=[{"theta": e2.start(p, dd)}] * 2, show_progress=False, show_console=False)
            sm = fit.summary(); c = sm[sm.index.str.startswith("log_T[")]; dr = fit.draws_pd()
            ess = float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()); lf = float(dr["n_leapfrog__"].mean())
            th = fit.stan_variable("theta"); th = th.reshape(-1, th.shape[-1])[:, :R + C - 1]; cc = np.corrcoef(th.T)
            rows.append(dict(part=part, design=design, rep=rep, coords=coords, param=r4.NAMES[p], max_rhat=float(c["R_hat"].max()), min_ess=ess, leapfrogs=lf,
                             ess_per_1000_grad=ess / lf, ess_per_s=ess / (time.time() - t0), div=int(dr["divergent__"].sum()),
                             max_corr_margin=float(np.abs(cc[np.triu_indices(R + C - 1, 1)]).max()), last_col_share=float(m[-1] / N), **extra))
            pd.DataFrame(rows).to_csv(OUT + ".csv", index=False)


def main():
    for bal in (0.5, 0.35):
        for k in (1, 2, 3, 4, 5):
            for rep in (1, 2, 3):
                d = sd.make_design(np.random.default_rng(rep), jitter=0.1, balance_row=bal, balance_col=bal); C = 5
                m = np.sort(d["m"])[::-1]; m = m[[j for j in range(C) if j != k - 1] + [k - 1]]
                fit_all("A", f"balance {bal}, rank {k} last", rep, d["w"], m, d, dict(balance=bal, last_rank=k))
            print("A", bal, k, "done", flush=True)
    for lab, kw in e19.DESIGNS:
        for rep in (1, 2):
            d = sd.make_design(np.random.default_rng(rep), jitter=0.1, **kw)
            fit_all("B", lab, rep, d["w"], d["m"], d, dict(balance=np.nan, last_rank=np.nan))
        print("B", lab, "done", flush=True)


if __name__ == "__main__":
    main()
