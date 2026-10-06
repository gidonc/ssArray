"""The alignment effect: when the large row meets the large column, the sequential schemes slow down.  Why?

Hypothesis: it is the choice of the cell found by subtraction, not alignment as such.  By default the left-over column
in row r is column r (the diagonal), and rows are filled smallest first.  With aligned unbalanced margins the diagonal
cell of a small row is the product of two small shares, so the left-over cell is a small difference of large numbers.
Test, all at 5x5, total 30,000, independence centre, sigma_b 0.5, both margins geometric:
  designs      balance 1, 0.75, 0.5, 0.35 aligned; balance 0.5 reversed and shifted
  left-over    'diagonal' (default), 'largest column' for every row, 'smallest column' for every row
  row order    smallest row first (default; the largest row is found by subtraction), or largest first (aligned 0.5 only)
Adjusted table takes no ordering, so it is the control.  Three replicates (margin jitter 0.1).
Reported with each run: the smallest share of its row that a left-over cell has at the centre.
Run:  python3 experiments/20_alignment.py     Writes experiments/output/20_alignment.csv and .md
"""
import os, sys, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, "experiments"))
import simdesign as sd
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2
DESIGNS = [(f"aligned, balance {b}", dict(balance_row=b, balance_col=b)) for b in (1.0, 0.75, 0.5, 0.35)] + \
          [("reversed, balance 0.5", dict(balance_row=0.5, balance_col=0.5, align="reversed")), ("shifted, balance 0.5", dict(balance_row=0.5, balance_col=0.5, align="shifted"))]


def main():
    rows = []; out = os.path.join(ROOT, "experiments", "output", "20_alignment")
    for lab, kw in DESIGNS:
        for order in (("smallest first", "largest first") if lab == "aligned, balance 0.5" else ("smallest first",)):
            for rem in ("diagonal", "largest column", "smallest column"):
                for rep in (1, 2, 3):
                    d = sd.make_design(np.random.default_rng(rep), jitter=0.1, **kw); f = d["features"]; R, C = f["R"], f["C"]; D = R * C
                    w, m = d["w"], d["m"]; T = np.outer(w, m) / d["N"]
                    ro = [int(i) + 1 for i in (np.argsort(w) if order == "smallest first" else np.argsort(-w))]
                    rc = {"diagonal": [r + 1 if r < C else C for r in range(R)], "largest column": [int(np.argmax(m)) + 1] * R, "smallest column": [int(np.argmin(m)) + 1] * R}[rem]
                    visited = ro[:-1]
                    share = min(d["centre"][r - 1, rc[r - 1] - 1] / w[r - 1] for r in visited)        # left-over cell as a share of its row
                    last_share = float((d["centre"][ro[-1] - 1] / m).min())                           # last row's cells as a share of their columns
                    for p in (1, 2, 3, 4):
                        if p == 4 and not (rem == "diagonal" and order == "smallest first"): continue  # no ordering in adjusted table
                        dd = e2.base_data(T, d["eps"], 1); dd["V"] = d["V"]; dd["row_order"] = ro; dd["rem_col"] = rc
                        dd.update(param=p, mu_b=d["V"].T @ np.log(d["centre"]).ravel(), sigma_b=d["sigma"], mu_logv=float(np.log(d["N"])), sigma_logv=1.0)
                        t0 = time.time()
                        fit = e2.model.sample(data=dd, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=rep,
                                              inits=[{"theta": e2.start(p, dd)}] * 2, show_progress=False, show_console=False)
                        sm = fit.summary(); c = sm[sm.index.str.startswith("log_T[")]; dr = fit.draws_pd()
                        ess = float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()); lf = float(dr["n_leapfrog__"].mean())
                        rows.append(dict(design=lab, row_order=order, left_over=rem, rep=rep, param=r4.NAMES[p], max_rhat=float(c["R_hat"].max()), min_ess=ess, leapfrogs=lf,
                                         ess_per_1000_grad=ess / lf, div=int(dr["divergent__"].sum()), left_over_share=float(share), last_row_share=last_share, seconds=time.time() - t0))
                        pd.DataFrame(rows).to_csv(out + ".csv", index=False)
                print(lab, order, rem, "done", flush=True)
    df = pd.DataFrame(rows)
    g = df.pivot_table(index=["design", "row_order", "left_over"], columns="param", values="ess_per_1000_grad", aggfunc="median", sort=False)
    s = df.groupby(["design", "row_order", "left_over"], sort=False)[["left_over_share", "last_row_share"]].median()
    lf = df.pivot_table(index=["design", "row_order", "left_over"], columns="param", values="leapfrogs", aggfunc="median", sort=False)
    open(out + ".md", "w").write("# Alignment and the left-over cell: median ESS per 1000 gradients (three replicates)\n\n" + s.round(3).join(g.round(1)).to_markdown() + "\n\n## Leapfrogs per draw\n\n" + lf.round(0).to_markdown() + "\n")
    pd.set_option("display.width", 250); print(s.round(3).join(g.round(1)).to_string()); print(lf.round(0).to_string()); print("max rhat", df.max_rhat.max().round(3), "div", df["div"].sum())


if __name__ == "__main__":
    main()
