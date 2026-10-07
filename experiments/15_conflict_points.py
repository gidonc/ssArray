"""A few points on the main expectations, with the cross-area kind of conflict: a common centre, an area's own margins.

Scottish 2007 constituencies collapsed to 5x5 (Con, Lab, LD, SNP, rest) on both votes, so every area has positive margins.
For an area, the observed margins are its own; the prior centre stands in for a cross-area mean:
  indep       the area's own independence table (no conflict, no interaction)
  national    the national table's composition (its margins differ from the area's: the conflict)
  raked       the national table raked to the area's margins (national odds ratios, no margin conflict)
Prior, in a basis split into row/column effects and interactions:
  equal sd            sigma_b on every coordinate
  interactions only   sigma_b on the interactions, sd 50 on the row and column effects (run with the national centre)
eps = 1 (count noise).  All six parameterisations; sequential ones with scaled margins.  One seed, two chains.
Part A: three areas (small, median, large mismatch with the national margins) x sigma_b 0.5, 2 x the four prior set-ups.
Part B: the median area, raked centre, sigma_b 1, with the volume multiplied by 0.01, 1, 100 (margin tightness).
Run:  python3 experiments/15_conflict_points.py     Writes experiments/output/15_conflict_points.csv and .md
"""
import os, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2
RG = [[0], [1], [2], [3], [4, 5, 6]]; CG = [[0], [1], [2], [5], [3, 4, 6]]        # rows: ..., SNP, spoilt, Other, Green; cols: ..., Other, Green, SNP, spoilt


def ipf(x, w, m, it=2000):
    x = x.copy()
    for _ in range(it):
        x *= (w / x.sum(1))[:, None]; x *= (m / x.sum(0))[None, :]
    return x


def main():
    d = pd.read_csv(os.path.join(ROOT, "data", "scot_2007_all_areas.csv")); R0, C0 = d.row_no.max(), d.col_no.max()
    col = lambda t: np.array([[t[np.ix_(r, c)].sum() for c in CG] for r in RG])
    tabs = {a: col(g.sort_values(["row_no", "col_no"]).votes.values.reshape(R0, C0).astype(float)) for a, g in d.groupby("area")}
    name = d.groupby("area").district.first()
    G = sum(tabs.values()); R, C = G.shape; D = R * C; nm = R + C - 2
    Hr, Hc = e2.helmert(R), e2.helmert(C); one = lambda n: np.ones((n, 1)) / np.sqrt(n)
    V = np.c_[np.kron(Hr, one(C)), np.kron(one(R), Hc), np.kron(Hr, Hc)]
    gw, gm = G.sum(1) / G.sum(), G.sum(0) / G.sum()
    mis = {a: float(np.sqrt(((np.log(t.sum(1) / t.sum()) - np.log(gw)) ** 2).sum() + ((np.log(t.sum(0) / t.sum()) - np.log(gm)) ** 2).sum())) for a, t in tabs.items()}
    order = sorted(mis, key=mis.get); areas = [order[3], order[len(order) // 2], order[-4]]
    print("areas:", [(a, name[a], round(mis[a], 2), tabs[a].sum()) for a in areas], flush=True)
    rows = []; out = os.path.join(ROOT, "experiments", "output", "15_conflict_points")

    def fit(part, a, vol, sb, cname, prior):
        t = tabs[a] * vol; w, m = t.sum(1), t.sum(0); N = t.sum()
        cen = {"indep": np.outer(w, m) / N, "national": G / G.sum() * N, "raked": ipf(G / G.sum() * N, w, m)}[cname]
        sig = np.r_[np.full(nm, sb if prior == "equal sd" else 50.0), np.full(D - 1 - nm, sb)]
        # size of the margin conflict: distance of the centre's row/column effects from the raked centre's, in prior sds
        conflict = float(np.linalg.norm(V[:, :nm].T @ (np.log(cen) - np.log(ipf(cen, w, m))).ravel() / sig[:nm]))
        for p in (0, 5, 1, 2, 3, 4):
            dd = e2.base_data(t, 1.0, 0 if p in (0, 5) else 1); dd["V"] = V
            dd.update(param=p, mu_b=V.T @ np.log(cen).ravel(), sigma_b=sig, mu_logv=float(np.log(N)), sigma_logv=1.0)
            row = dict(part=part, area=name[a], mismatch=mis[a], N=N, sigma_b=sb, centre=cname, prior=prior, conflict_sd=conflict, param=r4.NAMES[p]); t0 = time.time()
            try:
                f = e2.model.sample(data=dd, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=1,
                                    inits=[{"theta": e2.start(p, dd)}] * 2, show_progress=False, show_console=False, timeout=150)
                sm = f.summary(); c = sm[sm.index.str.startswith("log_T[")]; dr = f.draws_pd()
                ess = float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()); lf = float(dr["n_leapfrog__"].mean())
                row.update(seconds=time.time() - t0, max_rhat=float(c["R_hat"].max()), min_ess=ess, div=int(dr["divergent__"].sum()), leapfrogs=lf,
                           ess_per_1000_grad=ess / lf, ess_per_s=ess / (time.time() - t0))
            except Exception as ex:
                row.update(seconds=time.time() - t0, max_rhat=np.inf, min_ess=0.0, error=type(ex).__name__)
            rows.append(row); print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
            pd.DataFrame(rows).to_csv(out + ".csv", index=False)

    for a in areas:
        for sb in (0.5, 2.0):
            for cname, prior in (("indep", "equal sd"), ("raked", "equal sd"), ("national", "equal sd"), ("national", "interactions only")):
                fit("A", a, 1.0, sb, cname, prior)
    for vol in (0.01, 100.0):
        fit("B", areas[1], vol, 1.0, "raked", "equal sd")
    fit("B", areas[1], 1.0, 1.0, "raked", "equal sd")


if __name__ == "__main__":
    main()
