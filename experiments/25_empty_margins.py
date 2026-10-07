"""Areas with an empty row or column, left out of experiment 16: drop the empty rows and columns and fit the smaller table.

The common centre for such an area is the national table closed over the area's active cells (its subcomposition: every
log-ratio between active cells is kept).  Otherwise as experiment 16 with the current model (largest-column margin
coordinates, centred interior, split basis): equal-sd prior, centres indep / raked / national, sigma_b 0.5 and 2, eps 1,
four sequential schemes, up to six areas per dataset spread over the mismatch with the national margins, two seeds.
Run:  python3 experiments/25_empty_margins.py     Writes experiments/output/25_empty_margins.csv and .md
"""
import os, sys, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, "experiments"))
import simdesign as sd, margincoords as mc
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2
DATA = {"New Zealand 2017": "nz_2017_5x5.csv", "senc": "senc_3x3.csv", "redistrict": "redistrict_margins.csv"}


def main():
    rows = []; out = os.path.join(ROOT, "experiments", "output", "25_empty_margins")
    for dname, fn in DATA.items():
        d = pd.read_csv(os.path.join(ROOT, "data", fn)); interior = "votes" in d.columns; name = d.groupby("area").district.first()
        if interior:
            R0, C0 = d.row_no.max(), d.col_no.max()
            tabs = {a: g.sort_values(["row_no", "col_no"]).votes.values.reshape(R0, C0).astype(float) for a, g in d.groupby("area")}
        else:
            tabs = {}
            for a, g in d.groupby("area"):
                w = g[g.margin == "row"].sort_values("no")["count"].values.astype(float); m = g[g.margin == "col"].sort_values("no")["count"].values.astype(float)
                tabs[a] = np.outer(w, m) / w.sum()
        G = sum(tabs.values())
        if not interior: G = np.outer(G.sum(1), G.sum(0)) / G.sum()
        cand = {}
        for a, t in tabs.items():
            rk, ck = t.sum(1) > 0, t.sum(0) > 0
            if (rk.all() and ck.all()) or rk.sum() < 2 or ck.sum() < 2: continue
            tt = t[np.ix_(rk, ck)]; Ga = G[np.ix_(rk, ck)]; Ga = Ga / Ga.sum()                      # subcomposition of the national table
            mis = float(np.sqrt(((np.log(tt.sum(1) / tt.sum()) - np.log(Ga.sum(1))) ** 2).sum() + ((np.log(tt.sum(0) / tt.sum()) - np.log(Ga.sum(0))) ** 2).sum()))
            cand[a] = (tt, Ga, mis)
        order = sorted(cand, key=lambda a: cand[a][2]); n = len(order)
        areas = [order[int(q * (n - 1))] for q in np.linspace(0.05, 0.95, min(6, n))]
        print(dname, n, "areas with an empty row or column;", [(str(name[a]), cand[a][0].shape, round(cand[a][2], 2), int(cand[a][0].sum()), int(min(cand[a][0].sum(1).min(), cand[a][0].sum(0).min()))) for a in areas], flush=True)
        for seed in (1, 2):
            for a in areas:
                t, Ga, mis = cand[a]; R, C = t.shape; D = R * C; nm = R + C - 2; w, m = t.sum(1), t.sum(0); N = t.sum()
                V = sd.split_basis(R, C); Km, scm = mc.margin_K(w, m, 1.0, "largest", 1)
                for sb in (0.5, 2.0):
                    for cname in (("indep", "raked", "national") if interior else ("indep", "national")):
                        cen = {"indep": np.outer(w, m) / N, "national": Ga * N, "raked": sd.ipf(Ga * N, w, m)}[cname]
                        conflict = float(np.linalg.norm(V[:, :nm].T @ (np.log(cen) - np.log(sd.ipf(cen, w, m))).ravel() / sb))
                        for p in (1, 2, 3, 4):
                            dd = e2.base_data(t, 1.0, scm); dd["K_margin"] = Km; dd["centred_interior"] = 1; dd["split_basis"] = 1; dd["V"] = np.zeros((0, 0))
                            dd.update(param=p, mu_b=V.T @ np.log(cen).ravel(), sigma_b=np.full(D - 1, sb), mu_logv=float(np.log(N)), sigma_logv=1.0)
                            row = dict(dataset=dname, area=str(name[a]), R=R, C=C, mismatch=mis, N=N, min_margin=float(min(w.min(), m.min())), sigma_b=sb, centre=cname,
                                       conflict_sd=conflict, param=r4.NAMES[p], seed=seed); t0 = time.time()
                            try:
                                f = e2.model.sample(data=dd, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=seed,
                                                    inits=[{"theta": e2.start(p, dd)}] * 2, show_progress=False, show_console=False, timeout=150)
                                sm = f.summary(); c = sm[sm.index.str.startswith("log_T[")]; dr = f.draws_pd()
                                ess = float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()); lf = float(dr["n_leapfrog__"].mean())
                                row.update(seconds=time.time() - t0, max_rhat=float(c["R_hat"].max()), min_ess=ess, div=int(dr["divergent__"].sum()), leapfrogs=lf,
                                           ess_per_1000_grad=ess / lf, ess_per_s=ess / (time.time() - t0))
                            except Exception as ex:
                                row.update(seconds=time.time() - t0, max_rhat=np.inf, min_ess=0.0, error=type(ex).__name__)
                            rows.append(row); pd.DataFrame(rows).to_csv(out + ".csv", index=False)
            print(dname, "seed", seed, "done", len(rows), flush=True)


if __name__ == "__main__":
    main()
