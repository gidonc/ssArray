"""What is each allocation scheme linear in?  Split the map's nonlinearity into the interaction and main-effect parts.

Log cells (centred) split into an interaction part (log odds ratios, (R-1)(C-1) dimensions) and a main-effect part (row and
column effects, R+C-2 dimensions).  By construction
  adjusted table    lambda are the interaction parameters themselves: exactly linear in the interaction part
  adjusted row      lambda differences are log odds ratios with only the remaining rows collapsed
  odds-ratio logit  lambda is the log odds ratio of a 2x2 with remaining rows and columns collapsed
  position logit    lambda is a conditional logit within whichever margin binds, and that switches at the kinks
so the nonlinearity in the interaction part should rise in that order, and all of them are nonlinear in the main effects.
Same measure as 07 (one posterior-sd step in a random whitened direction, r = |actual - linear| / |linear|), computed
separately on the two parts.  Scottish grand total, eps 0.1.
Run:  BRIDGESTAN=/path/to/bridgestan python3 experiments/13_linear_in_what.py   Writes experiments/output/13_linear_in_what.csv/.md
"""
import os, importlib.util, numpy as np, pandas as pd
import bridgestan as bs
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("c6", os.path.join(ROOT, "experiments", "06_curvature.py"))
c6 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c6)
r4, e2, jsonable, SO = c6.r4, c6.e2, c6.jsonable, c6.SO


def main():
    T = r4.load(os.path.join(ROOT, "data", "scot_2007_all_areas.csv")); R, C = T.shape; D = R * C; k = R + C - 1
    Hr, Hc = e2.helmert(R), e2.helmert(C); one = lambda n: np.ones((n, 1)) / np.sqrt(n)
    B = {"interaction": np.kron(Hr, Hc), "main effects": np.c_[np.kron(Hr, one(C)), np.kron(one(R), Hc)]}
    rng = np.random.default_rng(1); rows = []
    for centre in ("indep", "actual"):
        for sb in (1.0, 2.0):
            for p in (1, 2, 3, 4):
                d = e2.base_data(T, 0.1, 1); V = d["V"]
                d.update(param=p, mu_b=V.T @ r4.centre_table(T, centre).ravel(), sigma_b=np.full(D - 1, sb), mu_logv=float(np.log(T.sum())), sigma_logv=1.0)
                f = e2.model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=1,
                                    inits=[{"theta": e2.start(p, d)}] * 2, show_progress=False, show_console=False)
                sm = f.summary(); c = sm[sm.index.str.startswith("log_T[")]
                th = f.stan_variable("theta"); th = th.reshape(-1, th.shape[-1])
                model = bs.StanModel(SO, jsonable(d), seed=1)
                lt = lambda t: np.log(model.param_constrain(np.ascontiguousarray(t, dtype=np.float64), include_tp=True)[D:2 * D].reshape((R, C), order="F")).ravel()
                for scope in ("all parameters", "interior only"):                 # interior only: margins held fixed, step in lambda alone
                    idx = np.arange(D) if scope == "all parameters" else np.arange(k, D)
                    L = np.linalg.cholesky(np.cov(th[:, idx].T) + 1e-12 * np.eye(len(idx)))
                    res = {q: [] for q in B}
                    for i in rng.choice(len(th), 100, replace=False):
                        x0 = lt(th[i])
                        for _ in range(8):
                            u = rng.standard_normal(len(idx)); dv = np.zeros(D); dv[idx] = L @ (u / np.linalg.norm(u))
                            lin = (lt(th[i] + 1e-4 * dv) - lt(th[i] - 1e-4 * dv)) / 2e-4; xs = lt(th[i] + dv)
                            if not np.all(np.isfinite(xs)): continue
                            for q, Bq in B.items():
                                a, b = Bq.T @ (xs - x0 - lin), Bq.T @ lin
                                if np.linalg.norm(b) > 1e-12: res[q].append(np.linalg.norm(a) / np.linalg.norm(b))
                    for q in B:
                        rows.append(dict(centre=centre, sigma_b=sb, param=r4.NAMES[p], rhat=float(c["R_hat"].max()), step=scope, part=q,
                                         median=float(np.median(res[q])), p90=float(np.quantile(res[q], 0.9))))
                print(r4.NAMES[p], centre, sb, {(r["step"], r["part"]): round(r["median"], 4) for r in rows[-4:]}, flush=True); del model
    df = pd.DataFrame(rows); out = os.path.join(ROOT, "experiments", "output", "13_linear_in_what")
    df.to_csv(out + ".csv", index=False)
    w = df.pivot_table(index=["centre", "sigma_b", "param"], columns=["step", "part"], values="median", sort=False).round(4)
    open(out + ".md", "w").write("# Nonlinearity of each scheme, split into interaction and main-effect parts (median r, one posterior sd)\n\n" + w.to_markdown() + "\n")
    print(w.to_string())


if __name__ == "__main__":
    main()
