"""Check and extend the conflict points (15) on two elections: Scotland 2007 and New Zealand 2017, both collapsed to 5x5.

Same design as 15, with more areas and two seeds.  For an area the observed margins are its own; the prior centre is
  indep     its own independence table                      (no conflict, no interaction)
  raked     the national table raked to the area's margins  (national odds ratios, no margin conflict)
  national  the national table's composition                (margins differ from the area's)
and the prior is equal sd on every coordinate, or (national centre only) on the interactions alone (sd 50 on row/column effects).
eps = 1.  Six parameterisations, sequential ones with scaled margins.  Areas: five per dataset, spread over the mismatch
between area and national margins, among areas with all margins positive.
Part A: areas x sigma_b (0.5, 2) x four set-ups x two seeds.   Part B: the middle area, raked centre, sigma_b 1, volume x 0.01, 1, 100.
Run:  python3 experiments/16_conflict_datasets.py     Writes experiments/output/16_conflict_datasets.csv
"""
import os, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2
DATA = {"Scotland 2007": "scot_2007_5x5.csv", "New Zealand 2017": "nz_2017_5x5.csv"}


def ipf(x, w, m, it=2000):
    x = x.copy()
    for _ in range(it):
        x *= (w / x.sum(1))[:, None]; x *= (m / x.sum(0))[None, :]
    return x


def main():
    rows = []; out = os.path.join(ROOT, "experiments", "output", "16_conflict_datasets.csv")
    for dname, fn in DATA.items():
        d = pd.read_csv(os.path.join(ROOT, "data", fn)); R, C = d.row_no.max(), d.col_no.max(); D = R * C; nm = R + C - 2
        tabs = {a: g.sort_values(["row_no", "col_no"]).votes.values.reshape(R, C).astype(float) for a, g in d.groupby("area")}
        tabs = {a: t for a, t in tabs.items() if t.sum(1).min() > 0 and t.sum(0).min() > 0}
        name = d.groupby("area").district.first()
        G = sum(tabs.values())
        Hr, Hc = e2.helmert(R), e2.helmert(C); one = lambda n: np.ones((n, 1)) / np.sqrt(n)
        V = np.c_[np.kron(Hr, one(C)), np.kron(one(R), Hc), np.kron(Hr, Hc)]
        gw, gm = G.sum(1) / G.sum(), G.sum(0) / G.sum()
        mis = {a: float(np.sqrt(((np.log(t.sum(1) / t.sum()) - np.log(gw)) ** 2).sum() + ((np.log(t.sum(0) / t.sum()) - np.log(gm)) ** 2).sum())) for a, t in tabs.items()}
        order = sorted(mis, key=mis.get); n = len(order)
        areas = [order[int(q * (n - 1))] for q in (0.05, 0.275, 0.5, 0.725, 0.95)]
        print(dname, len(tabs), "areas;", [(str(name[a]), round(mis[a], 2), int(tabs[a].sum())) for a in areas], flush=True)

        def fit(part, a, vol, sb, cname, prior, seed):
            t = tabs[a] * vol; w, m = t.sum(1), t.sum(0); N = t.sum()
            cen = {"indep": np.outer(w, m) / N, "national": G / G.sum() * N, "raked": ipf(G / G.sum() * N, w, m)}[cname]
            sig = np.r_[np.full(nm, sb if prior == "equal sd" else 50.0), np.full(D - 1 - nm, sb)]
            conflict = float(np.linalg.norm(V[:, :nm].T @ (np.log(cen) - np.log(ipf(cen, w, m))).ravel() / sig[:nm]))
            for p in (0, 5, 1, 2, 3, 4):
                dd = e2.base_data(t, 1.0, 0 if p in (0, 5) else 1); dd["V"] = V
                dd.update(param=p, mu_b=V.T @ np.log(cen).ravel(), sigma_b=sig, mu_logv=float(np.log(N)), sigma_logv=1.0)
                row = dict(dataset=dname, part=part, area=str(name[a]), mismatch=mis[a], N=N, sigma_b=sb, centre=cname, prior=prior, conflict_sd=conflict,
                           param=r4.NAMES[p], seed=seed); t0 = time.time()
                try:
                    f = e2.model.sample(data=dd, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=seed,
                                        inits=[{"theta": e2.start(p, dd)}] * 2, show_progress=False, show_console=False, timeout=150)
                    sm = f.summary(); c = sm[sm.index.str.startswith("log_T[")]; dr = f.draws_pd()
                    ess = float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()); lf = float(dr["n_leapfrog__"].mean())
                    row.update(seconds=time.time() - t0, max_rhat=float(c["R_hat"].max()), min_ess=ess, div=int(dr["divergent__"].sum()), leapfrogs=lf,
                               ess_per_1000_grad=ess / lf, ess_per_s=ess / (time.time() - t0))
                except Exception as ex:
                    row.update(seconds=time.time() - t0, max_rhat=np.inf, min_ess=0.0, error=type(ex).__name__)
                rows.append(row); pd.DataFrame(rows).to_csv(out, index=False)
            print(dname, part, name[a], vol, sb, cname, prior, seed, "done", len(rows), flush=True)

        for seed in (1, 2):
            for a in areas:
                for sb in (0.5, 2.0):
                    for cname, prior in (("indep", "equal sd"), ("raked", "equal sd"), ("national", "equal sd"), ("national", "interactions only")):
                        fit("A", a, 1.0, sb, cname, prior, seed)
            for vol in (0.01, 1.0, 100.0):
                fit("B", areas[2], vol, 1.0, "raked", "equal sd", seed)


if __name__ == "__main__":
    main()
