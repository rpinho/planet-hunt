"""Round 8, stage 2: merge every source, exclude the known/claimed, rank by detection power with all data now on MAST.

    .venv/bin/python round8/rank.py merge          # signals_all.csv, signals_excluded.csv, funnel
    .venv/bin/python round8/rank.py tic            # TIC rows for the merged pool (own cache file)
    .venv/bin/python round8/rank.py rank --n 900   # power ranking -> targets.csv

Expected power (as round 7): their SNR scaled by sqrt(data):
    E_hold = SNR_src * sqrt(N_new / N_theirs),  E_full = SNR_src * sqrt((N_theirs + N_new) / N_theirs).
N_new comes from the MAST inventory with the per-source labelling rule (holdout.py); for Eumemic's tables N_theirs is his
own n_sectors, so QLP files his reader dropped show up as N_le_cut - n_sectors (round 7's estimator, which matched his
counts for 113/120).  Single-sector sources are penalised for the look-elsewhere window (x0.85).
SNR is capped at 15 for ranking. Score (the 2x2): E_full + 0.5 E_hold, x1.15 for R* < 0.45 (mid/late M), x1.1 for P < 3 d (many transits);
Eumemic-listed (E3) rows x0.8 (someone already lists them).
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd

from common import DATA, OUT, R8C, config

KNOWN = DATA / "novelty" / "known.csv"
ALIASES = [1, 2, 3, 0.5, 1 / 3, 1.5, 2 / 3, 4, 0.25]
TIC_CACHE = R8C / "tic_rows_r8.csv"
TIC_COLS = ["ID", "ra", "dec", "Tmag", "Teff", "rad", "mass", "d", "contratio", "logg", "GAIA", "plx", "lumclass",
            "disposition"]
PRIORITY = {"E1": 0, "E2": 1, "C1": 2, "C2": 3, "C3": 4, "C4": 5, "C5": 6, "C6": 7, "C7": 8, "E3": 9}


def alias_match(P, Ps, tol=0.003):
    Ps = np.asarray([p for p in Ps if np.isfinite(p) and p > 0])
    if not np.isfinite(P) or not len(Ps):
        return None
    for a in ALIASES:
        r = np.abs(P / (Ps * a) - 1)
        if (r < tol).any():
            return int(np.argmin(r)), a
    return None


def cmd_merge(a):
    e = pd.read_csv(OUT / "src_eumemic.csv")
    c = pd.read_csv(OUT / "src_community.csv")
    d = pd.concat([e, c], ignore_index=True)
    d["mode"] = d["mode"].fillna("fixed") if "mode" in d else "fixed"
    d["src_sectors"] = d["src_sectors"].fillna("")
    steps = [("all source rows", len(d), d.tic.nunique())]
    for s, g in d.groupby("stream"):
        steps.append((f"  {s}", len(g), g.tic.nunique()))
    d = d[np.isfinite(d.period) & (d.period >= 0.5) & (d.period <= 30)]
    steps.append(("0.5 <= P <= 30 d", len(d), d.tic.nunique()))
    k = pd.read_csv(KNOWN)
    kg = k.groupby("tic")
    why = []
    for r in d.itertuples():
        w = ""
        if r.tic in kg.groups:
            g = kg.get_group(r.tic)
            if r.stream == "E3":
                g = g[g.kind != "eumemic_list"]
            g_nan = g[g.P.isna() & (g.kind == "community_claim")]
            if len(g_nan):
                w = g_nan.name.iloc[0]
            else:
                m = alias_match(r.period, g.P.values)
                if m:
                    row = g.iloc[[i for i, p in enumerate(g.P.values) if np.isfinite(p) and p > 0][m[0]]]
                    w = f"{row.kind}: {row['name']} x{m[1]:.3g}"
        why.append(w)
    d["excluded"] = why
    ex = d[d.excluded != ""]
    d = d[d.excluded == ""].copy()
    for kind in ("TOI", "CTOI", "confirmed", "TCE", "eumemic_list", "community_claim"):
        steps.append((f"  excluded as {kind}", int(ex.excluded.str.startswith(kind).sum()), None))
    steps.append(("not TOI/CTOI/confirmed/TCE/claimed at an alias", len(d), d.tic.nunique()))
    # dedupe: same star, same period (or alias) across streams -> keep the higher-priority stream, credit all
    d["prio"] = d.stream.map(PRIORITY)
    d = d.sort_values(["tic", "prio", "snr_src"], ascending=[True, True, False])
    keep, also = [], []
    for tic, g in d.groupby("tic", sort=False):
        kept = []
        for i, r in g.iterrows():
            m = alias_match(r.period, [g.loc[j, "period"] for j in kept], tol=0.002)
            if m:
                also.append((kept[m[0]], f"{r.stream}:{r.author.split(' (')[0]}"))
            else:
                kept.append(i)
        keep += kept
    d2 = d.loc[keep].copy()
    extra = pd.DataFrame(also, columns=["idx", "who"]).groupby("idx").who.apply(lambda s: "; ".join(sorted(set(s))))
    d2["also_in"] = d2.index.map(extra).fillna("")
    steps.append(("deduplicated across sources (same star, period or alias)", len(d2), d2.tic.nunique()))
    d2.to_csv(OUT / "signals_all.csv", index=False)
    ex.to_csv(OUT / "signals_excluded.csv", index=False)
    f = pd.DataFrame(steps, columns=["step", "signals", "stars"])
    f.to_csv(OUT / "funnel_merge.csv", index=False)
    print(f.to_string(index=False))
    print(d2.groupby("stream").size().to_string())


def tic_rows(tics) -> pd.DataFrame:
    from astroquery.mast import Catalogs
    parts = []
    for p in (config.CATALOG_DIR / "tic_rows.csv", TIC_CACHE):
        if p.exists():
            x = pd.read_csv(p, low_memory=False)
            parts.append(x[[c for c in TIC_COLS if c in x.columns]])
    have = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=TIC_COLS)
    have["ID"] = have["ID"].astype("int64")
    have = have.drop_duplicates("ID")
    todo = sorted({int(t) for t in tics} - set(have.ID))
    print(f"[tic] {len(todo)} to query", flush=True)
    for i in range(0, len(todo), 400):
        ids = [str(t) for t in todo[i:i + 400]]
        t0, k = time.time(), 0
        while True:
            try:
                tab = Catalogs.query_criteria(catalog="Tic", ID=ids)
                break
            except Exception as exc:
                k += 1
                if time.time() - t0 > 900:
                    print("[tic] give up chunk", exc, flush=True)
                    tab = None
                    break
                time.sleep(min(60, 5 * k))
        if tab is None:
            continue
        df = tab.to_pandas()
        df = df[[c for c in TIC_COLS if c in df.columns]]
        df.to_csv(TIC_CACHE, mode="a", header=not TIC_CACHE.exists(), index=False)
        if (i // 400) % 10 == 0:
            print(f"[tic] {min(i + 400, len(todo))}/{len(todo)}", flush=True)
    x = pd.concat([have, pd.read_csv(TIC_CACHE, low_memory=False)] if TIC_CACHE.exists() else [have], ignore_index=True)
    x["ID"] = x["ID"].astype("int64")
    return x.drop_duplicates("ID")[lambda z: z.ID.isin({int(t) for t in tics})]


def cmd_tic(a):
    d = pd.read_csv(OUT / "signals_all.csv")
    t = tic_rows(d.tic.unique())
    print(len(t), "TIC rows for", d.tic.nunique(), "stars")


def counts(row, sectors: dict) -> tuple[int, int, int]:
    """(n_avail, n_theirs_est, n_new_est) from the inventory, by the source's labelling rule."""
    allsec = sorted(sectors)
    rule = row["rule"]
    if rule == "eum_early":
        n_le = sum(s <= 105 for s in allsec)
        nsrc = int(row["n_sec_src"]) if np.isfinite(row["n_sec_src"]) else n_le
        n_new = max(n_le - nsrc, 0) + sum(s > 105 for s in allsec)
        return len(allsec), nsrc, n_new
    if rule == "eum_allsky":
        n_th = sum(s <= 104 and not (sectors[s] == {"QLP"} and 74 <= s <= 79) for s in allsec)
        n_new = len(allsec) - n_th
        # dropped empty-error QLP files are unknown until downloaded: count half of the QLP-only S94-97/S101-104
        q = sum(s in (set(range(94, 98)) | set(range(101, 105))) and sectors[s] == {"QLP"} for s in allsec)
        return len(allsec), max(n_th - q // 2, 1), n_new + q // 2
    if rule == "lk_ffi":
        n_th_max = sum(s <= 104 and bool(sectors[s] & {"TESS-SPOC", "QLP"}) for s in allsec)
        nsrc = int(row["n_sec_src"]) if np.isfinite(row["n_sec_src"]) else n_th_max
        n_new = max(n_th_max - nsrc, 0) + sum(s > 104 or not (sectors[s] & {"TESS-SPOC", "QLP"}) for s in allsec)
        return len(allsec), nsrc, n_new
    src = set(json.loads(row["src_sectors"])) if str(row["src_sectors"]).startswith("[") else set()
    n_th = len(src & set(allsec)) or len(src)
    return len(allsec), n_th, len(set(allsec) - src)


def cmd_rank(a):
    from inventory import load_all
    d = pd.read_csv(OUT / "signals_all.csv")
    prod = load_all()
    prod = prod[prod.tic.isin(set(d.tic))].dropna(subset=["dataURL"])
    two = pd.read_csv(config.CACHE / "r5" / "twomin_index.csv.gz")
    two = two[two.tic.isin(set(d.tic))].drop_duplicates(["tic", "sector"])
    secmap = {}
    for r in prod.itertuples():
        s = int(r.sequence_number)
        if 0 < s < 200:
            secmap.setdefault(int(r.tic), {}).setdefault(s, set()).add(r.provenance_name)
    for r in two.itertuples():
        secmap.setdefault(int(r.tic), {}).setdefault(int(r.sector), set()).add("SPOC-2min")
    inv_done = set(load_all().tic)
    c = [counts(r, secmap.get(int(r["tic"]), {})) for r in d.to_dict("records")]
    d["n_avail"], d["n_theirs"], d["n_new"] = zip(*c)
    d["inventoried"] = d.tic.isin(inv_done)
    # C7 has no SNR column (crude stand-in); capped at 15: the round is about weak signals, and a huge signal
    # someone rejected is a variability/EB suspect - it still enters if it has new data, it just does not crowd the top
    snr = d.snr_src.fillna(d.depth_ppm / 300.0).clip(upper=15)
    d["snr_rank"] = snr
    d["E_hold"] = snr * np.sqrt(d.n_new / d.n_theirs.clip(lower=1))
    d["E_full"] = snr * np.sqrt((d.n_theirs + d.n_new) / d.n_theirs.clip(lower=1))
    d.loc[d["mode"] == "window", ["E_hold", "E_full"]] *= 0.85
    d0 = d.copy()
    d0.to_csv(OUT / "power_all.csv", index=False)
    d = d[d.tic.isin(inv_done) & (d.n_new >= 1) & (d.E_hold >= a.min_ehold)]   # TIC only for the plausible pool
    tic = tic_rows(d.tic.unique()).rename(columns={"ID": "tic", "Teff": "teff", "rad": "rstar", "mass": "mstar",
                                                    "Tmag": "tmag"})
    d = d.merge(tic[["tic", "teff", "rstar", "mstar", "tmag", "ra", "dec", "d", "contratio", "logg", "GAIA"]],
                on="tic", how="left")
    steps = [("merged, not known/claimed", len(d0), d0.tic.nunique()),
             ("MAST inventory done", int(d0.inventoried.sum()), d0[d0.inventoried].tic.nunique()),
             ("any independent new data (N_new >= 1)", int((d0.inventoried & (d0.n_new >= 1)).sum()), None),
             (f"E_hold >= {a.min_ehold}", len(d), d.tic.nunique())]
    small = ((d.teff < 4200) | (d.rstar < 0.7)) & d.rstar.notna() & d.mstar.notna() & (d.rstar < 1.0)
    d = d[small].copy()
    steps.append(("small star (TIC Teff < 4200 K or R* < 0.7, R* < 1)", len(d), d.tic.nunique()))
    d["teff"] = d.teff.fillna(3500.0)
    d["score"] = (d.E_full + 0.5 * d.E_hold) * np.where(d.rstar < 0.45, 1.15, 1.0) * np.where(d.period < 3, 1.1, 1.0) \
        * np.where(d.stream == "E3", 0.8, 1.0)
    d = d.sort_values("score", ascending=False)
    d.to_csv(OUT / "rank_pool.csv", index=False)
    # one star may carry several signals: all of them go in when the star is selected
    stars = list(dict.fromkeys(d.tic))[:a.n]
    tg = d[d.tic.isin(set(stars))].copy()
    tg["star_rank"] = tg.tic.map({t: i + 1 for i, t in enumerate(stars)})
    tg = tg.sort_values(["star_rank", "score"], ascending=[True, False])
    steps.append((f"top {a.n} stars by score (all their signals)", len(tg), tg.tic.nunique()))
    tg.to_csv(OUT / "targets.csv", index=False)
    f = pd.DataFrame(steps, columns=["step", "signals", "stars"])
    f.to_csv(OUT / "funnel_rank.csv", index=False)
    print(f.to_string(index=False))
    print(tg.groupby("stream").agg(n=("tic", "size"), E_hold_med=("E_hold", "median")).to_string())
    print(tg[["tic", "stream", "period", "snr_src", "n_theirs", "n_new", "E_hold", "E_full", "teff", "rstar", "tmag",
              "score"]].head(30).round(2).to_string())


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sp = p.add_subparsers(dest="cmd", required=True)
    sp.add_parser("merge")
    sp.add_parser("tic")
    s = sp.add_parser("rank")
    s.add_argument("--n", type=int, default=900)
    s.add_argument("--min-ehold", type=float, default=2.0)
    a = p.parse_args()
    {"merge": cmd_merge, "tic": cmd_tic, "rank": cmd_rank}[a.cmd](a)
