"""Fundamentals check: HMC cost on a Gaussian is set by the ratio of the largest to the smallest scale left after the metric.

10-dimensional Gaussian, nine directions with sd 1 and one tight direction with sd s.
  aligned   the tight direction is one parameter
  rotated   the tight direction is spread equally over all ten parameters (same distribution, rotated)
  diagonal metric  Stan's default adaptation: should remove the aligned case, not the rotated one, where the
                   step size should fall like s and the leapfrogs per draw rise like 1/s (up to the tree-depth cap, 1023)
  dense metric     should remove both
Reports step size, leapfrogs per draw, min ESS, and ESS per 1000 gradients.
Run:  python3 experiments/12_stiffness_basics.py     Writes experiments/output/12_stiffness_basics.csv and .md
"""
import os, numpy as np, pandas as pd, logging
from cmdstanpy import CmdStanModel
logging.getLogger("cmdstanpy").setLevel(logging.ERROR)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
m = CmdStanModel(stan_file=os.path.join(ROOT, "experiments", "stan", "gauss.stan"))
D = 10; rows = []
u = np.ones(D) / np.sqrt(D)
for s in (1, 0.1, 0.01, 0.001):
    for layout in ("aligned", "rotated"):
        if layout == "aligned": S = np.diag(np.r_[s ** 2, np.ones(D - 1)])
        else: S = np.eye(D) - (1 - s ** 2) * np.outer(u, u)
        L = np.linalg.cholesky(S)
        for metric, kw in (("diagonal", dict(metric="diag_e")), ("dense", dict(metric="dense_e"))):
            f = m.sample(data=dict(D=D, L=L), chains=2, parallel_chains=2, iter_warmup=1000, iter_sampling=1000, seed=1,
                         show_progress=False, show_console=False, **kw)
            sm = f.summary(); c = sm[sm.index.str.startswith("theta")]; dr = f.draws_pd()
            ess = float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()); lf = float(dr["n_leapfrog__"].mean())
            rows.append(dict(s=s, layout=layout, metric=metric, step_size=float(np.mean(f.step_size)), leapfrogs=lf, min_ess=ess,
                             max_rhat=float(c["R_hat"].max()), ess_per_1000_grad=ess / (lf * 2000) * 1000))
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
df = pd.DataFrame(rows); out = os.path.join(ROOT, "experiments", "output", "12_stiffness_basics")
df.to_csv(out + ".csv", index=False)
open(out + ".md", "w").write("# Gaussian with one tight direction of sd s (nine others sd 1)\n\n" + df.round(4).to_markdown(index=False) + "\n")
