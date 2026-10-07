"""First use of the separate dials (simdesign.py): which aspect of size, and which aspect of the margins, moves which scheme?

Size: 3x3, 5x5, 8x8 with the total, the count per margin, or the count per cell held fixed (they coincide at 5x5).
Margins (5x5, total 30,000): rows unbalanced three ways at the same effective number where possible, and both sides unbalanced
with the large margins aligned or reversed.  Independence centre, sigma_b 0.5, eps 1.  Two replicates (margin jitter 0.1).
ILR, log cells, odds-ratio logit, adjusted row.
Run:  python3 experiments/19_size_margin_dials.py     Writes experiments/output/19_size_margin_dials.csv and .md
"""
import os, sys, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, "experiments"))
import simdesign as sd
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2
DESIGNS = [("5x5 baseline", {})]
for R in (3, 8):
    DESIGNS += [(f"{R}x{R}, total fixed", dict(R=R, C=R)), (f"{R}x{R}, per-margin fixed", dict(R=R, C=R, volume=("margin", 6000))),
                (f"{R}x{R}, per-cell fixed", dict(R=R, C=R, volume=("cell", 1200)))]
DESIGNS += [("rows geometric 0.5", dict(balance_row=0.5)), ("rows one dominant 0.5", dict(balance_row=0.5, family_row="dominant")),
            ("rows one tiny 0.85", dict(balance_row=0.85, family_row="tiny")), ("both geometric 0.5, aligned", dict(balance_row=0.5, balance_col=0.5)),
            ("both geometric 0.5, reversed", dict(balance_row=0.5, balance_col=0.5, align="reversed"))]


def main():
    rows = []; out = os.path.join(ROOT, "experiments", "output", "19_size_margin_dials")
    for lab, kw in DESIGNS:
        for rep in (1, 2):
            d = sd.make_design(np.random.default_rng(rep), jitter=0.1, **kw); f = d["features"]; R, C = f["R"], f["C"]; D = R * C
            T = np.outer(d["w"], d["m"]) / d["N"]
            for p in (0, 5, 2, 3):
                dd = e2.base_data(T, d["eps"], 0 if p in (0, 5) else 1); dd["V"] = d["V"]
                dd.update(param=p, mu_b=d["V"].T @ np.log(d["centre"]).ravel(), sigma_b=d["sigma"], mu_logv=float(np.log(d["N"])), sigma_logv=1.0)
                t0 = time.time()
                fit = e2.model.sample(data=dd, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=rep,
                                      inits=[{"theta": e2.start(p, dd)}] * 2, show_progress=False, show_console=False)
                sm = fit.summary(); c = sm[sm.index.str.startswith("log_T[")]; dr = fit.draws_pd()
                ess = float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()); lf = float(dr["n_leapfrog__"].mean())
                rows.append(dict(design=lab, rep=rep, param=r4.NAMES[p], max_rhat=float(c["R_hat"].max()), min_ess=ess, leapfrogs=lf, step_size=float(np.mean(fit.step_size)),
                                 ess_per_1000_grad=ess / lf, div=int(dr["divergent__"].sum()), seconds=time.time() - t0,
                                 **{k: f[k] for k in ("total", "row_margin_mean", "cell_mean", "margin_min", "margin_max", "balance_row", "balance_col", "tight_index_max", "tight_index_min")}))
                pd.DataFrame(rows).to_csv(out + ".csv", index=False)
        print(lab, "done", flush=True)
    df = pd.DataFrame(rows); labs = [l for l, _ in DESIGNS]
    g = df.pivot_table(index="design", columns="param", values="ess_per_1000_grad", aggfunc="median").reindex(labs)
    h = df.groupby("design")[["total", "row_margin_mean", "cell_mean", "margin_min", "margin_max"]].median().reindex(labs)
    lf = df.pivot_table(index="design", columns="param", values="leapfrogs", aggfunc="median").reindex(labs)
    open(out + ".md", "w").write("# Size and margin dials: median ESS per 1000 gradients (two replicates)\n\n" + h.round(0).join(g.round(1)).to_markdown() + "\n\n## Leapfrogs per draw\n\n" + lf.round(0).to_markdown() + "\n")
    pd.set_option("display.width", 250); print(h.round(0).join(g.round(1)).to_string()); print(lf.round(0).to_string()); print(df.groupby("design").max_rhat.max().reindex(labs).round(2).to_string())


if __name__ == "__main__":
    main()
