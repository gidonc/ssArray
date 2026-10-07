"""Follow-up to 20: the 'alignment' effect is the size of the reference column in the margin parameterisation.

In single_area.stan the column margins are a softmax of log-ratios against the LAST column.  If that column is small,
its noise enters every column log-ratio, so the margin parameters are strongly correlated and a diagonal metric cannot
separate them.  (Aligned and reversed margins are the same posterior with the columns relabelled: the prior and the
penalty do not depend on column order.  Only the coordinates do.)
Test: one table (5x5, total 30,000, independence centre, sigma_b 0.5, both margins geometric, balance 0.5 and 0.35).
Rows and columns in descending order, except that the column of size rank k is moved to the last position, k = 1 (largest)
to 5 (smallest).  Nothing else changes.  All four sequential schemes, three replicates.
Also reported: the largest correlation between column-margin parameters in the posterior, and the condition number of the
correlation matrix of the margin parameters.
Run:  python3 experiments/21_reference_column.py     Writes experiments/output/21_reference_column.csv and .md
"""
import os, sys, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, "experiments"))
import simdesign as sd
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2


def main():
    rows = []; out = os.path.join(ROOT, "experiments", "output", "21_reference_column")
    for bal in (0.5, 0.35):
        for k in (1, 2, 3, 4, 5):
            for rep in (1, 2, 3):
                d = sd.make_design(np.random.default_rng(rep), jitter=0.1, balance_row=bal, balance_col=bal); R, C = 5, 5; D = R * C
                w = d["w"]; m = np.sort(d["m"])[::-1]; perm = [j for j in range(C) if j != k - 1] + [k - 1]; m = m[perm]      # rank-k column last
                N = d["N"]; T = np.outer(w, m) / N
                for p in (1, 2, 3, 4):
                    dd = e2.base_data(T, 1.0, 1); dd["V"] = d["V"]
                    dd.update(param=p, mu_b=d["V"].T @ np.log(T).ravel(), sigma_b=d["sigma"], mu_logv=float(np.log(N)), sigma_logv=1.0)
                    t0 = time.time()
                    fit = e2.model.sample(data=dd, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=rep,
                                          inits=[{"theta": e2.start(p, dd)}] * 2, show_progress=False, show_console=False)
                    sm = fit.summary(); c = sm[sm.index.str.startswith("log_T[")]; dr = fit.draws_pd()
                    ess = float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()); lf = float(dr["n_leapfrog__"].mean())
                    th = fit.stan_variable("theta"); th = th.reshape(-1, th.shape[-1])
                    cm = np.corrcoef(th[:, :R + C - 1].T); ev = np.linalg.eigvalsh(cm)
                    cc = np.corrcoef(th[:, R:R + C - 1].T); mx = float(np.abs(cc[np.triu_indices(C - 1, 1)]).max())
                    rows.append(dict(balance=bal, last_column_rank=k, last_column_share=float(m[-1] / N), rep=rep, param=r4.NAMES[p], max_rhat=float(c["R_hat"].max()),
                                     min_ess=ess, leapfrogs=lf, ess_per_1000_grad=ess / lf, max_corr_col_params=mx, cond_margin_corr=float(ev[-1] / ev[0]), seconds=time.time() - t0))
                    pd.DataFrame(rows).to_csv(out + ".csv", index=False)
            print(bal, k, "done", flush=True)
    df = pd.DataFrame(rows)
    g = df.pivot_table(index=["balance", "last_column_rank"], columns="param", values="ess_per_1000_grad", aggfunc="median")
    s = df.groupby(["balance", "last_column_rank"])[["last_column_share", "max_corr_col_params", "cond_margin_corr"]].median()
    lf = df.pivot_table(index=["balance", "last_column_rank"], columns="param", values="leapfrogs", aggfunc="median")
    open(out + ".md", "w").write("# Size of the reference (last) column: median ESS per 1000 gradients (three replicates)\n\n" + s.round(3).join(g.round(1)).to_markdown() + "\n\n## Leapfrogs per draw\n\n" + lf.round(0).to_markdown() + "\n")
    pd.set_option("display.width", 250); print(s.round(3).join(g.round(1)).to_string()); print(lf.round(0).to_string())


if __name__ == "__main__":
    main()
