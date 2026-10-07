"""Does a split basis remove the linearity problem?  Prior on the interactions alone against the equal-sd prior.

Basis: row effects, column effects (12 coordinates), then interactions (36).  With equal sds this is a rotation of the
Helmert basis and changes nothing.  With the main-effect sds made very wide (50) the prior acts on the interactions only,
and the main effects are left to the margins.  Prediction: adjusted table, whose parameters are the interaction
coordinates, gains most; position logit and odds-ratio logit, whose parameters are collapsed quantities, gain least.
Scottish grand total, eps 0.1, scaled margins, two seeds.
Run:  python3 experiments/14_split_prior.py     Writes experiments/output/14_split_prior.csv and .md
"""
import os, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2


def main():
    T = r4.load(os.path.join(ROOT, "data", "scot_2007_all_areas.csv")); R, C = T.shape; D = R * C; nm = R + C - 2
    Hr, Hc = e2.helmert(R), e2.helmert(C); one = lambda n: np.ones((n, 1)) / np.sqrt(n)
    V = np.c_[np.kron(Hr, one(C)), np.kron(one(R), Hc), np.kron(Hr, Hc)]
    assert np.allclose(V.T @ V, np.eye(D - 1))
    rows = []
    for centre in ("indep", "actual"):
        for sb in (1.0, 2.0):
            for prior, s_main in (("equal sd", sb), ("interactions only", 50.0)):
                for p in (1, 2, 3, 4):
                    for seed in (1, 2):
                        d = e2.base_data(T, 0.1, 1); d["V"] = V
                        d.update(param=p, mu_b=V.T @ r4.centre_table(T, centre).ravel(), sigma_b=np.r_[np.full(nm, s_main), np.full(D - 1 - nm, sb)],
                                 mu_logv=float(np.log(T.sum())), sigma_logv=1.0)
                        t0 = time.time()
                        f = e2.model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=seed,
                                            inits=[{"theta": e2.start(p, d)}] * 2, show_progress=False, show_console=False)
                        sec = time.time() - t0
                        sm = f.summary(); c = sm[sm.index.str.startswith("log_T[")]; dr = f.draws_pd()
                        ess = float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()); lf = float(dr["n_leapfrog__"].mean())
                        rows.append(dict(centre=centre, sigma_b=sb, prior=prior, param=r4.NAMES[p], seed=seed, max_rhat=float(c["R_hat"].max()), min_ess=ess,
                                         div=int(dr["divergent__"].sum()), leapfrogs=lf, ess_per_1000_grad=ess / lf, seconds=sec))
                        print({a: (round(b, 3) if isinstance(b, float) else b) for a, b in rows[-1].items()}, flush=True)
    df = pd.DataFrame(rows); out = os.path.join(ROOT, "experiments", "output", "14_split_prior")
    df.to_csv(out + ".csv", index=False)
    g = df.assign(fail=df.max_rhat > 1.05).groupby(["centre", "sigma_b", "param", "prior"], sort=False).agg(failed=("fail", "sum"), max_rhat=("max_rhat", "max"),
            min_ess=("min_ess", "median"), leapfrogs=("leapfrogs", "median"), ess_per_1000_grad=("ess_per_1000_grad", "median"), div=("div", "sum"), seconds=("seconds", "median")).round(2)
    open(out + ".md", "w").write("# Prior on interactions alone against equal sds, split basis, Scottish grand total (eps 0.1)\n\n" + g.to_markdown() + "\n")
    print(g.to_string())


if __name__ == "__main__":
    main()
