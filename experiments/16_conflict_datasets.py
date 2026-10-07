"""Check and extend the conflict points (15) on two elections: Scotland 2007 and New Zealand 2017, both collapsed to 5x5.

Same design as 15, with more areas and two seeds.  For an area the observed margins are its own; the prior centre is
  indep     its own independence table                      (no conflict, no interaction)
  raked     the national table raked to the area's margins  (national odds ratios, no margin conflict)
  national  the national table's composition                (margins differ from the area's)
and the prior is equal sd on every coordinate, or (national centre only) on the interactions alone (sd 50 on row/column effects).
eps = 1.  Six parameterisations, sequential ones with scaled margins.  Areas: five per dataset, spread over the mismatch
between area and national margins, among areas with all margins positive.
Part A: areas x sigma_b (0.5, 2) x four set-ups x two seeds.   Part B: the middle area, raked centre, sigma_b 1, volume x 0.01, 1, 100.
Run:  python3 experiments/16_conflict_datasets.py [datasets, comma separated] [output label]
Default: the two elections, written to experiments/output/16_conflict_datasets.csv
Environment: MARGIN_COORDS = last (default) | largest | whitened (see margincoords.py), SEQ_ONLY = 1 for the sequential
schemes only, SEEDS = e.g. "1" (default "1,2"), MODEL = current for the model as it now stands (largest-column margin
coordinates, centred interior, split basis without its matrix).
"""
import os, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2
import sys; sys.path.insert(0, os.path.join(ROOT, "experiments")); import margincoords as mc
CURRENT = os.environ.get("MODEL") == "current"      # MODEL=current: largest-column margins, centred interior, split basis without its matrix
COORDS = "largest" if CURRENT else os.environ.get("MARGIN_COORDS", "last"); PARAMS = (1, 2, 3, 4) if os.environ.get("SEQ_ONLY") == "1" else (0, 5, 1, 2, 3, 4)
SEEDS = [int(x) for x in os.environ.get("SEEDS", "1,2").split(",")]
DATA = {"Scotland 2007": "scot_2007_5x5.csv", "New Zealand 2017": "nz_2017_5x5.csv", "senc": "senc_3x3.csv", "redistrict": "redistrict_margins.csv"}
# redistrict has margins only (no known interior): its common centre is the independence table of the pooled margins, and there is no raked set-up


def ipf(x, w, m, it=2000):
    x = x.copy()
    for _ in range(it):
        x *= (w / x.sum(1))[:, None]; x *= (m / x.sum(0))[None, :]
    return x


def main():
    import sys
    which = sys.argv[1].split(",") if len(sys.argv) > 1 else ["Scotland 2007", "New Zealand 2017"]
    label = sys.argv[2] if len(sys.argv) > 2 else "16_conflict_datasets"
    rows = []; out = os.path.join(ROOT, "experiments", "output", label + ".csv")
    for dname in which:
        fn = DATA[dname]
        d = pd.read_csv(os.path.join(ROOT, "data", fn)); interior = "votes" in d.columns
        if interior:
            R, C = d.row_no.max(), d.col_no.max()
            tabs = {a: g.sort_values(["row_no", "col_no"]).votes.values.reshape(R, C).astype(float) for a, g in d.groupby("area")}
        else:                                                       # margins only: carry each area as its independence table
            tabs = {}
            for a, g in d.groupby("area"):
                w = g[g.margin == "row"].sort_values("no")["count"].values.astype(float); m = g[g.margin == "col"].sort_values("no")["count"].values.astype(float)
                tabs[a] = np.outer(w, m) / w.sum()
            R, C = next(iter(tabs.values())).shape
        D = R * C; nm = R + C - 2
        tabs = {a: t for a, t in tabs.items() if t.sum(1).min() > 0 and t.sum(0).min() > 0}
        name = d.groupby("area").district.first()
        G = sum(tabs.values())
        if not interior: G = np.outer(G.sum(1), G.sum(0)) / G.sum()
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
            Km, scm = mc.margin_K(w, m, 1.0, COORDS, 1)
            for p in PARAMS:
                dd = e2.base_data(t, 1.0, 0 if p in (0, 5) else scm); dd["V"] = V; dd["K_margin"] = Km
                if CURRENT: dd["centred_interior"] = 1; dd["split_basis"] = 1; dd["V"] = np.zeros((0, 0))
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

        for seed in SEEDS:
            for a in areas:
                for sb in (0.5, 2.0):
                    for cname, prior in (("indep", "equal sd"), ("raked", "equal sd"), ("national", "equal sd"), ("national", "interactions only")):
                        if cname == "raked" and not interior: continue
                        fit("A", a, 1.0, sb, cname, prior, seed)
            for vol in (0.01, 1.0, 100.0):
                fit("B", areas[2], vol, 1.0, "raked" if interior else "indep", "equal sd", seed)


if __name__ == "__main__":
    main()
