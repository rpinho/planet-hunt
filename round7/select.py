"""Round 7, stage 1: pick Eumemic's weak / unpromoted signals worth a second look.

    .venv/bin/python round7/select.py stage1     # table-only cuts -> runs/round7_eumemic_weak/stage1.csv
    .venv/bin/python round7/select.py products   # MAST product inventory for stage-1 stars (cached)
    .venv/bin/python round7/select.py rank       # expected power with all sectors -> targets.csv

Source: github.com/eumemic/exoplanets results_public/all_signals*.csv.gz (every BLS
signal of his earlier searches: K/M dwarfs T<=11.5, TOI hosts, faint M dwarfs
11.5<T<=13.5, SPOC-TCE re-check, 20-100 d). His promotion path starts with a
pre-filter at SNR >= 9, so the SNR 6-9 band was never vetted (TOI-6000 c sat
there at SNR 7.05).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from planet_hunt import config  # noqa: E402

EUM = Path.home() / "dev/planet-hunt-data/r7/eumemic_src/repo"
OUT = config.RUNS / "round7_eumemic_weak"
OUT.mkdir(parents=True, exist_ok=True)
R5 = config.CACHE / "r5"
TABLES = {"km": "all_signals", "toih": "all_signals_toi_hosts", "f1": "all_signals_faint_t1",
          "f2": "all_signals_faint_t2", "tce": "all_signals_tce", "long": "all_signals_long",
          "longtoi": "all_signals_long_toi_hosts"}
ALIASES = [1, 2, 3, 0.5, 1 / 3, 1.5, 2 / 3, 4, 0.25]


def load_signals() -> pd.DataFrame:
    out = []
    for key, f in TABLES.items():
        d = pd.read_csv(EUM / "results_public" / f"{f}.csv.gz")
        d["src"] = key
        out.append(d)
    d = pd.concat(out, ignore_index=True)
    for c in ("toi_alias", "repeat", "rot_alias"):
        d[c] = d[c].fillna(False).astype(bool)
    return d


def his_candidates() -> pd.DataFrame:
    """Every signal he has listed in any tier (counted or not)."""
    rp = EUM / "results_public"
    rows = []
    for f in ["candidates.csv", "candidates_2026-10-08.csv", "candidates_2026-10-08b.csv",
              "candidates_2026-10-08c.csv", "candidates_2026-10-08d.csv"]:
        d = pd.read_csv(rp / f)
        rows.append(pd.DataFrame({"tic": d["TIC"], "P": d["P_d"], "his_list": f}))
    d = pd.read_csv(rp / "candidates_allsky.csv")
    rows.append(pd.DataFrame({"tic": d["TIC"], "P": d["P_d"], "his_list": "allsky:" + d["tier"].astype(str)}))
    d = pd.read_csv(rp / "tiers_earlier_candidates.csv")
    rows.append(pd.DataFrame({"tic": d["TIC"], "P": d["P_d"], "his_list": "tiers:" + d["tier"].astype(str)}))
    for f in sorted(rp.glob("followup_*.csv")):
        d = pd.read_csv(f)
        rows.append(pd.DataFrame({"tic": d["tic"], "P": d["period"], "his_list": f.stem + ":" + d["flags"].fillna("clean").astype(str)}))
    c = pd.concat(rows, ignore_index=True)
    c["tic"] = pd.to_numeric(c["tic"], errors="coerce")
    return c.dropna(subset=["tic"]).astype({"tic": "int64"})


def catalogs():
    toi = pd.read_csv(R5 / "toi_r5.csv", comment="#", low_memory=False)
    toi = pd.DataFrame({"tic": toi["TIC ID"], "P": pd.to_numeric(toi["Period (days)"], errors="coerce"),
                        "name": "TOI-" + toi["TOI"].astype(str), "disp": toi["TFOPWG Disposition"].astype(str)})
    ctoi = pd.read_csv(R5 / "ctoi_r5.csv", comment="#", low_memory=False)
    ccol = [c for c in ctoi.columns if c.lower().startswith("period")][0]
    tcol = [c for c in ctoi.columns if c.upper().startswith("TIC")][0]
    ncol = [c for c in ctoi.columns if c.upper().startswith("CTOI")][0]
    ctoi = pd.DataFrame({"tic": ctoi[tcol], "P": pd.to_numeric(ctoi[ccol], errors="coerce"),
                         "name": "CTOI-" + ctoi[ncol].astype(str), "disp": "CTOI"})
    ps = pd.read_csv(R5 / "pscomp_r5.csv")
    ps = pd.DataFrame({"tic": pd.to_numeric(ps["tic_id"].astype(str).str.replace("TIC", "").str.strip(), errors="coerce"),
                       "P": ps["pl_orbper"], "name": ps["pl_name"], "disp": "confirmed"})
    cat = pd.concat([toi, ctoi, ps], ignore_index=True).dropna(subset=["tic"])
    cat["tic"] = cat["tic"].astype("int64")
    return cat


def alias_match(P, Ps, tol=0.003):
    Ps = np.asarray([p for p in Ps if np.isfinite(p) and p > 0])
    if not len(Ps):
        return None
    for a in ALIASES:
        r = np.abs(P / (Ps * a) - 1)
        if (r < tol).any():
            return float(Ps[np.argmin(r)]), a
    return None


def cmd_stage1(a):
    d = load_signals()
    n0 = len(d)
    f = d
    steps = [("all signals", len(f))]
    f = f[(f.snr >= a.snr_min) & (f.snr < a.snr_max)]; steps.append((f"{a.snr_min} <= snr < {a.snr_max}", len(f)))
    f = f[f.n_good_events >= 3]; steps.append((">= 3 good transits", len(f)))
    f = f[(f.period >= 0.5)]; steps.append(("P >= 0.5 d", len(f)))
    f = f[~f.repeat & ~f.rot_alias & ~f.toi_alias]; steps.append(("not repeat / rotation alias / TOI alias", len(f)))
    f = f[f.rho > a.rho_min]; steps.append((f"stellar density > {a.rho_min} g/cc (small star proxy)", len(f)))
    f = f[(f.odd_even_sigma.abs() < 3) & ((f.secondary / f.secondary_err) < 3)]; steps.append(("odd/even < 3 sigma, secondary < 3 sigma", len(f)))
    lim = np.minimum(0.8, 1 / np.sqrt(f.n_good_events) + 0.35)
    f = f[f.max_ses_frac < lim]; steps.append(("no single dominant event", len(f)))
    r = f.duration_h / f.exp_duration_h
    f = f[(r > 0.25) & (r < 1.6)]; steps.append(("duration consistent with the star", len(f)))
    f = f[f.sde >= a.sde_min]; steps.append((f"SDE >= {a.sde_min}", len(f)))
    # one row per (star, period): the same star can sit in several of his tables
    f = f.copy()
    f["pkey"] = (f.tic.astype(str) + "_" + np.round(np.log(f.period) / 0.002).astype(int).astype(str))
    f = f.sort_values("snr", ascending=False).drop_duplicates("pkey")
    steps.append(("deduplicated across his tables", len(f)))
    # his own lists (any tier) and the catalogs
    his = his_candidates()
    cat = catalogs()
    hisg = his.groupby("tic")["P"].apply(list).to_dict()
    hisl = his.groupby("tic")["his_list"].apply(lambda s: ";".join(sorted(set(s)))).to_dict()
    catg = cat.groupby("tic")
    keep, why = [], []
    for row in f.itertuples():
        m = alias_match(row.period, hisg.get(row.tic, []))
        if m:
            keep.append(False); why.append(f"his list {hisl[row.tic]}"); continue
        if row.tic in catg.groups:
            g = catg.get_group(row.tic)
            m = alias_match(row.period, g["P"].values)
            if m:
                nm = g.loc[(g["P"] - m[0]).abs().idxmin(), "name"]
                keep.append(False); why.append(f"catalog {nm} x{m[1]:.3g}"); continue
        keep.append(True); why.append("")
    f["excluded"] = why
    ex = f[~np.array(keep)]
    f = f[np.array(keep)]
    steps.append(("not in his lists (any tier) nor TOI/CTOI/confirmed at an alias", len(f)))
    f.to_csv(OUT / "stage1.csv", index=False)
    ex.to_csv(OUT / "stage1_excluded.csv", index=False)
    pd.DataFrame(steps, columns=["step", "n"]).to_csv(OUT / "stage1_funnel.csv", index=False)
    for s, n in steps:
        print(f"{n:>8}  {s}")
    print(f.groupby("src").size())
    print("stars:", f.tic.nunique())


def prerank(s1: pd.DataFrame, n_faint: int) -> pd.DataFrame:
    """Which stage-1 signals get a MAST product inventory: every TOI-host,
    long-period and bright K/M signal, and the best n_faint of each faint table
    by his SNR + SDE (the faint tables hold 33k signals, nearly all noise)."""
    s1 = s1.copy()
    s1["score0"] = s1.snr + 0.1 * s1.sde.clip(upper=25)
    keep = s1[s1.src.isin(["toih", "longtoi", "long", "km", "tce"])]
    faint = [s1[s1.src == k].nlargest(n_faint, "score0") for k in ("f1", "f2")]
    return pd.concat([keep] + faint, ignore_index=True)


def cmd_products(a):
    from planet_hunt import ffi
    s1 = prerank(pd.read_csv(OUT / "stage1.csv"), a.n_faint)
    s1.to_csv(OUT / "stage1_inventoried.csv", index=False)
    tics = s1.tic.unique().tolist()
    print(len(s1), "signals on", len(tics), "stars", flush=True)
    prod = ffi.lightcurve_products(tics, chunk=150)
    print(len(prod), "product rows for", prod.tic.nunique(), "stars")


HIS_CUT = {"km": 105, "toih": 105, "longtoi": 105, "long": 105, "tce": 105, "f1": 105, "f2": 105}


def coverage(tics, prod, two):
    """Per star: chosen product per sector (2-min > TESS-SPOC > QLP)."""
    rows = []
    two = two[two.tic.isin(tics)].drop_duplicates(["tic", "sector"])
    p = prod[prod.tic.isin(tics)].dropna(subset=["dataURL"])
    t2 = two.groupby("tic")["sector"].apply(set).to_dict()
    ts = p[p.provenance_name == "TESS-SPOC"].groupby("tic")["sequence_number"].apply(lambda x: set(x.astype(int))).to_dict()
    tq = p[p.provenance_name == "QLP"].groupby("tic")["sequence_number"].apply(lambda x: set(x.astype(int))).to_dict()
    for t in tics:
        a, b, c = t2.get(t, set()), ts.get(t, set()), tq.get(t, set())
        allsec = {s for s in (a | b | c) if 0 < s < 200}; b = {s for s in b if s < 200}; c = {s for s in c if s < 200}
        qlp_only = c - a - b
        rows.append({"tic": t, "sectors": sorted(allsec), "n_all": len(allsec), "n_2min": len(a),
                     "max_sector": max(allsec) if allsec else 0,
                     "qlp_only_dropped_range": len([s for s in qlp_only if 94 <= s <= 97 or 101 <= s <= 104])})
    return pd.DataFrame(rows)


def cmd_rank(a):
    from planet_hunt import data as pdata
    s1 = pd.read_csv(OUT / "stage1_inventoried.csv")
    prod = pd.read_csv(config.TABLE_DIR / "lcprod.csv")
    two = pd.read_csv(R5 / "twomin_index.csv.gz")
    cov = coverage(s1.tic.unique().tolist(), prod, two).set_index("tic")
    s1 = s1.join(cov, on="tic")
    # his cutoff: which (104 or 105) reproduces his n_sectors best, per table
    for src, g in s1.groupby("src"):
        res = {}
        for cut in (104, 105):
            n_le = g.sectors.apply(lambda ss: sum(x <= cut for x in ss))
            res[cut] = float(np.mean(n_le == g.n_sectors))
        print(f"{src}: fraction with n(sectors<=cut)==his n_sectors: {res}  max sector seen {g.max_sector.max()}")
    s1["his_cut"] = s1.src.map(HIS_CUT)
    s1["n_le_cut"] = [sum(x <= c for x in ss) for ss, c in zip(s1.sectors, s1.his_cut)]
    s1["n_after"] = [sum(x > c for x in ss) for ss, c in zip(s1.sectors, s1.his_cut)]
    s1["n_missing_le_cut"] = (s1.n_le_cut - s1.n_sectors).clip(lower=0)
    s1["n_new"] = s1.n_missing_le_cut + s1.n_after
    s1["n_total"] = s1.n_sectors + s1.n_new
    s1["E_hold"] = s1.snr * np.sqrt(s1.n_new / s1.n_sectors)
    s1["E_full"] = s1.snr * np.sqrt(s1.n_total / s1.n_sectors)
    print(s1[["n_sectors", "n_le_cut", "n_after", "n_new", "E_hold", "E_full"]].describe().round(2))
    s1.to_csv(OUT / "stage2_power.csv", index=False)
    # stellar parameters for the plausible ones
    pool = s1[(s1.n_new >= 1) & (s1.E_hold >= a.min_ehold)].copy()
    tic = pdata.tic_info(pool.tic.unique().tolist())
    tic = tic.rename(columns={"ID": "tic", "Teff": "teff", "rad": "rstar", "mass": "mstar", "Tmag": "tmag"})
    pool = pool.merge(tic[["tic", "teff", "rstar", "mstar", "tmag", "ra", "dec", "d", "contratio"]], on="tic", how="left")
    small = (pool.teff < 4200) | (pool.rstar < 0.7)
    pool = pool[small & pool.rstar.notna() & pool.mstar.notna()].copy()
    pool["teff"] = pool.teff.fillna(3500.0)
    pool["rp_his"] = np.sqrt(pool.depth.clip(lower=0)) * pool.rstar * 109.1
    pool["score"] = pool.E_full + 0.5 * pool.E_hold
    pool = pool.sort_values("score", ascending=False)
    print(len(pool), "small-star signals with E_hold >=", a.min_ehold)
    tg = pool.head(a.n).copy()
    tg.to_csv(OUT / "targets.csv", index=False)
    pool.to_csv(OUT / "stage2_pool.csv", index=False)
    print(tg.src.value_counts())
    print(tg[["tic", "src", "period", "snr", "n_sectors", "n_new", "E_hold", "E_full", "teff", "rstar", "tmag"]].head(40).to_string())


def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("stage1")
    s.add_argument("--snr-min", type=float, default=6.0)
    s.add_argument("--snr-max", type=float, default=9.0)
    s.add_argument("--rho-min", type=float, default=2.5)
    s.add_argument("--sde-min", type=float, default=8.0)
    s = sp.add_parser("products")
    s.add_argument("--n-faint", type=int, default=900)
    s = sp.add_parser("rank")
    s.add_argument("--min-ehold", type=float, default=2.0)
    s.add_argument("--n", type=int, default=120)
    a = p.parse_args()
    {"stage1": cmd_stage1, "products": cmd_products, "rank": cmd_rank}[a.cmd](a)


if __name__ == "__main__":
    main()
