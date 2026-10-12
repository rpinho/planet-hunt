"""Method A: faint red dwarfs in TESS full-frame images, many sectors stacked.

    .venv/bin/python mdwarf_hunt.py select --n 300      # pick stars, list their FFI light curves
    .venv/bin/python mdwarf_hunt.py fetch               # download + shrink them (resumable)
    .venv/bin/python mdwarf_hunt.py search              # stitch all sectors, BLS, vet, cross-match
    .venv/bin/python mdwarf_hunt.py vet                 # centroid, crowding, neighbours, SPOC, diff image
    .venv/bin/python mdwarf_hunt.py vet --eb --no-diffimage   # method C: eclipsing binaries
    .venv/bin/python mimic_report.py                    # one report for methods A and B

Idea credited to Mehran Jalali's public thread
(https://x.com/mehran__jalali/status/2108339123669111033): a red dwarf 1/6
the Sun's width, 3% deep transits every 11.7 d, found by averaging dozens of
TESS transits from 2018-2025. Small planets give deep transits on small
stars, so faint M dwarfs that were never 2-minute targets are a gap SPOC's
planet search does not cover.

Selection: TIC dwarfs with Teff 2700-3900 K, R 0.08-0.6 Rsun, Tmag 11-15,
d < 100 pc, within 12 deg of an ecliptic pole (TESS's continuous viewing
zones: 30-70 sectors per star). Stars WITHOUT any SPOC 2-minute light curve
are preferred (never searched by SPOC); known TOI hosts among the M dwarfs
are added as positive controls.
Outputs runs/mdwarf/{targets.csv, chosen_lcs.csv, results.csv, lc/, vet/}.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import time

import numpy as np
import pandas as pd

from planet_hunt import config
from planet_hunt import ffi

RUN = config.RUNS / "mdwarf"


def cmd_select(args):
    RUN.mkdir(parents=True, exist_ok=True)
    stars = pd.concat([ffi.mdwarfs_near_pole(p, radius_deg=args.radius) for p in ("south", "north")],
                      ignore_index=True).drop_duplicates("ID")
    print(f"[select] {len(stars)} M dwarfs (Tmag 11-15, d<100 pc) within {args.radius} deg of the poles",
          flush=True)
    prod = ffi.lightcurve_products(list(stars["ID"]))
    cov = ffi.coverage_table(prod)
    stars = stars.rename(columns={"ID": "tic"}).astype({"tic": "int64"})
    cov["tic"] = cov["tic"].astype("int64")
    stars = stars.merge(cov, on="tic", how="left").fillna(
        {"n_spoc": 0, "n_tessspoc": 0, "n_qlp": 0, "n_ffi": 0})
    stars.to_csv(RUN / "all_stars.csv", index=False)

    from planet_hunt.catalogs import load_toi
    toi = load_toi()
    hosts = set(toi["tic"])
    stars["toi_host"] = stars["tic"].isin(hosts)
    pool = stars[stars["n_ffi"] >= args.min_sectors]
    ctrl = pool[pool["toi_host"]].sample(min(args.controls, int(pool["toi_host"].sum())),
                                         random_state=args.seed)
    ctrl = ctrl.assign(control=True)
    fresh = pool[~pool["toi_host"] & (pool["n_spoc"] == 0)]
    n_rand = args.n - len(ctrl)
    rand = fresh.sample(min(n_rand, len(fresh)), random_state=args.seed).assign(control=False)
    if len(rand) < n_rand:   # top up with stars that do have some 2-min data
        rest = pool[~pool["toi_host"] & (pool["n_spoc"] > 0)]
        rand = pd.concat([rand, rest.sample(min(n_rand - len(rand), len(rest)),
                                            random_state=args.seed).assign(control=False)])
    sel = pd.concat([ctrl, rand], ignore_index=True)
    sel.to_csv(RUN / "targets.csv", index=False)
    chosen = ffi.choose_sectors(prod[prod["tic"].isin(sel["tic"])])
    chosen.to_csv(RUN / "chosen_lcs.csv", index=False)
    print(f"[select] pool with >= {args.min_sectors} FFI sectors: {len(pool)} "
          f"({int((pool['n_spoc'] == 0).sum())} never 2-min). Picked {len(sel)}: "
          f"{len(ctrl)} TOI-host controls + {len(rand)} others; {len(chosen)} light curves "
          f"(median {sel['n_ffi'].median():.0f} sectors/star)", flush=True)


def cmd_fetch(args):
    chosen = pd.read_csv(RUN / "chosen_lcs.csv")
    n = ffi.fetch_many(chosen, workers=args.workers, max_cache_gb=args.max_gb)
    print(f"[fetch] {n} light curves extracted; cache {ffi.cache_gb():.2f} GB", flush=True)


# --------------------------------------------------------------------------
# search
# --------------------------------------------------------------------------
def search_star(job):
    """Stitch all sectors of one star, search, mask, search again."""
    from planet_hunt import multisector as ms
    tic, star, max_signals, pmax = job
    rows = []
    try:
        sectors = ffi.load_star(tic)
        if len(sectors) < 2:
            return [{"tic": tic, "signal": 1, "error": f"only {len(sectors)} sectors"}]
        rstar = float(star.get("rad") or 0.3)
        mstar = float(star.get("mass") or rstar)
        teff = float(star.get("Teff") or 3300)
        masks = []

        def mask_fn(t):
            m = np.zeros(len(t), bool)
            for (P, t0, d) in masks:
                m |= ms.in_tr(t, P, t0, d, 1.5)
            return m
        for k in range(max_signals):
            lc = ms.stitch(sectors, transit_mask=mask_fn if masks else None)
            if not lc or len(lc["t"]) < 1000:
                break
            # cut out transits of earlier signals before searching again
            keep = ~mask_fn(lc["t"]) if masks else np.ones(len(lc["t"]), bool)
            lc_s = {kk: (v[keep] if isinstance(v, np.ndarray) and len(v) == len(keep) else v)
                    for kk, v in lc.items()}
            b = ms.search(lc_s, rstar, mstar, pmax=pmax)
            # re-detrend with this signal masked too, then measure
            def measure(P, t0, dur):
                masks_try = masks + [(P, t0, dur)]
                lc2 = ms.stitch(sectors, transit_mask=lambda t: np.any(
                    [ms.in_tr(t, P_, t0_, d_, 1.5) for (P_, t0_, d_) in masks_try], axis=0))
                keep2 = np.ones(len(lc2["t"]), bool)
                for (P_, t0_, d_) in masks:
                    keep2 &= ~ms.in_tr(lc2["t"], P_, t0_, d_, 1.5)
                lc2 = {kk: (v[keep2] if isinstance(v, np.ndarray) and len(v) == len(keep2) else v)
                       for kk, v in lc2.items()}
                return lc2, ms.transit_stats(lc2, P, t0, dur, rstar, mstar, teff)
            lc2, st = measure(b["period"], b["t0"], b["duration"])
            # BLS often locks onto TWICE the true period of a short-period
            # planet: then the "secondary eclipse" at phase 0.5 is just as deep
            # as the primary. Halve the period in that case.
            e_sec = st["depth_sec_ppm"] / st["sec_sigma"] if st["sec_sigma"] and np.isfinite(st["sec_sigma"]) else np.inf
            if (np.isfinite(st["sec_sigma"]) and st["sec_sigma"] > 3 and abs(st["depth_sec_ppm"] - st["depth_ppm"])
                    < 3 * np.hypot(e_sec, st["depth_err_ppm"]) and b["period"] / 2 >= 0.25):
                b["period"] /= 2
                b["halved"] = True
                lc2, st = measure(b["period"], b["t0"], b["duration"])
            eb, caution = ms.vet(st, b["sde"])
            tce = bool(st["snr"] >= ms.SNR_MIN and st["n_transits"] >= 3 and b["sde"] >= ms.SDE_MIN)
            np.savez_compressed(RUN / "lc" / f"{tic}_{k + 1}.npz", t=lc2["t"],
                                flat=lc2["flat"].astype(np.float32), dy=lc2["dy"].astype(np.float32),
                                cad=lc2["cad"].astype(np.float32), sector=lc2["sector"].astype(np.int16),
                                cx=lc2["cx"].astype(np.float32), cy=lc2["cy"].astype(np.float32),
                                periods=b["periods"].astype(np.float32),
                                power=b["power"].astype(np.float32))
            row = {"tic": tic, "signal": k + 1, **{kk: v for kk, v in st.items()
                                                    if kk not in ("season_depths", "per_transit")},
                   "sde": b["sde"], "power_coh": b["power_coh"], "power_semi": b["power_semi"],
                   "n_seasons": b["n_seasons"], "T_all": b["T_all"], "n_periods": b["n_periods"],
                   "n_sectors": len(lc2["sector_info"]), "n_points": len(lc2["t"]),
                   "eb_flags": "; ".join(eb), "caution_flags": "; ".join(caution), "tce": tce,
                   "season_depths": json.dumps(st["season_depths"]),
                   "rstar": rstar, "mstar": mstar, "teff": teff, "halved": bool(b.get("halved", False))}
            rows.append(row)
            if not tce:
                break
            masks.append((b["period"], b["t0"], b["duration"]))
    except Exception as exc:
        import traceback
        rows.append({"tic": tic, "signal": len(rows) + 1,
                     "error": (repr(exc) + traceback.format_exc()[-300:])[:500]})
    return rows


def cmd_search(args):
    from planet_hunt.catalogs import crossmatch, load_confirmed, load_ctoi, load_toi
    (RUN / "lc").mkdir(parents=True, exist_ok=True)
    tg = pd.read_csv(RUN / "targets.csv")
    if args.tic:
        tg = tg[tg["tic"].isin(args.tic)]
    jobs = [(int(r.tic), r._asdict(), args.max_signals, args.pmax) for r in tg.itertuples()]
    t0 = time.time()
    rows = []
    with cf.ProcessPoolExecutor(args.workers) as ex:
        for k, r in enumerate(ex.map(search_star, jobs, chunksize=1), 1):
            rows.extend(r)
            if k % 10 == 0 or k == len(jobs):
                print(f"[search] {k}/{len(jobs)} stars ({time.time() - t0:.0f} s)", flush=True)
    res = pd.DataFrame(rows)
    for c in ("error", "period", "snr", "tce", "eb_flags", "caution_flags"):
        if c not in res:
            res[c] = np.nan
    toi, conf, ctoi = load_toi(), load_confirmed(), load_ctoi()
    labs = [crossmatch(int(r.tic), float(r.period) if pd.notna(r.period) else np.nan, toi, conf, ctoi)
            for r in res.itertuples()]
    res["catalog_status"] = [x["status"] for x in labs]
    res["catalog_label"] = [x["label"] for x in labs]
    res["tce"] = res["tce"].fillna(False).astype(bool)
    res = res.merge(tg[["tic", "Tmag", "d", "ra", "dec", "pole", "n_spoc", "n_ffi", "control",
                        "contratio"]], on="tic", how="left")
    res = flag_shared_periods(res)
    res = flag_repeats(res)
    out = RUN / ("results.csv" if not args.tic else "results_subset.csv")
    res.sort_values(["tce", "snr"], ascending=False).to_csv(out, index=False)
    t = res[res["tce"]]
    print(f"[search] {len(tg)} stars, {len(t)} TCEs; catalog: "
          f"{t['catalog_status'].value_counts().to_dict()}; clean: "
          f"{int((t['eb_flags'].fillna('') == '').sum())}; {time.time() - t0:.0f} s", flush=True)


def flag_repeats(res: pd.DataFrame) -> pd.DataFrame:
    """A period that comes back after its transits were cut out was never a
    transit: rotational modulation, pulsation or an eclipsing binary
    leaking through. Flag the later signal and the original."""
    from planet_hunt.catalogs import period_match

    def add(i, msg):
        old = res.at[i, "eb_flags"] if isinstance(res.at[i, "eb_flags"], str) else ""
        if msg not in old:
            res.at[i, "eb_flags"] = "; ".join(filter(None, [old, msg]))
    for tic, g in res[res["tce"]].groupby("tic"):
        g = g.sort_values("signal")
        for i, r in g.iterrows():
            for j, q in g[g["signal"] < r["signal"]].iterrows():
                if period_match(r["period"], q["period"], tol=0.02) is not None:
                    add(i, f"repeat of signal {int(q['signal'])} (persists after masking)")
                    add(j, "persists after masking: variability or EB, not a transit")
    return res


def flag_shared_periods(res: pd.DataFrame, tol: float = 5e-4) -> pd.DataFrame:
    """The same period (within 0.05%) on >= 2 other stars = a systematic."""
    ok = res[res["tce"]].dropna(subset=["period"])
    counts = {}
    for i, r in ok.iterrows():
        others = ok[(ok["tic"] != r["tic"]) & (np.abs(ok["period"] / r["period"] - 1) < tol)]
        counts[i] = others["tic"].nunique()
    res["shared_period"] = pd.Series(counts).reindex(res.index).fillna(0).astype(int)
    for i in res.index[res["shared_period"] >= 2]:
        old = res.at[i, "eb_flags"] if isinstance(res.at[i, "eb_flags"], str) else ""
        res.at[i, "eb_flags"] = "; ".join(filter(None, [old, f"same period on {res.at[i, 'shared_period']} other stars (systematic)"]))
    return res


# --------------------------------------------------------------------------
# vet
# --------------------------------------------------------------------------
def spoc_tces_any(tic: int, period: float) -> dict:
    """SPOC TCEs for this star from ANY DV run (only 2-min stars have them)."""
    import urllib.request
    from planet_hunt import spoc
    from planet_hunt.catalogs import period_match
    spoc.DV_DIR.mkdir(parents=True, exist_ok=True)
    runs = spoc._dv_products(tic)
    found = []
    for s0, s1, uri, fn in runs:
        path = spoc.DV_DIR / fn
        if not path.exists():
            with urllib.request.urlopen(config.MAST_DOWNLOAD + uri, timeout=120) as r:
                path.write_bytes(r.read())
        for t in spoc._parse_dvr(path.read_text(errors="ignore")):
            t["run"] = f"S{s0}" if s0 == s1 else f"S{s0}-{s1}"
            t["match"] = period_match(period, t["period"], tol=0.01)
            found.append(t)
    hits = [t for t in found if t["match"] is not None]
    if hits:
        h = hits[0]
        v = f"SPOC TCE matches: P={h['period']:.4f} d, {h['depth_ppm']:.0f} ppm ({h['run']}{h['match']})"
    elif found:
        v = "SPOC TCE(s) at other periods: " + "; ".join(f"P={t['period']:.3f} d ({t['run']})" for t in found[:4])
    elif runs:
        v = "SPOC searched this star (2-min data) and raised no TCE at this period"
    else:
        v = "no SPOC DV products (not a 2-min target, or no TCE ever)"
    return {"n_runs": len(runs), "tces": found, "match": bool(hits), "verdict": v}


def gaia_context(gaia_id, ra: float, dec: float) -> dict:
    """Gaia DR3 RUWE (>1.4 suggests an unresolved binary), and whether Gaia's
    eclipsing-binary table or the AAVSO VSX lists anything within 30\"."""
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    out = {}
    v = Vizier(columns=["Source", "RUWE", "Gmag", "BP-RP", "Plx", "IPDfmp", "Dup"], row_limit=5)
    c = SkyCoord(ra * u.deg, dec * u.deg)
    try:
        r = v.query_region(c, radius=5 * u.arcsec, catalog="I/355/gaiadr3")
        if len(r):
            t = r[0].to_pandas()
            if pd.notna(gaia_id):
                t = t[t["Source"] == int(gaia_id)] if (t["Source"] == int(gaia_id)).any() else t.head(1)
            out["ruwe"] = float(t["RUWE"].iloc[0])
            out["n_gaia_5as"] = int(len(r[0]))
    except Exception as exc:
        out["gaia_error"] = repr(exc)[:100]
    # also: published FFI planet-candidate lists not in ExoFOP (DTARPS-S,
    # Melton et al. 2024; neural-network FFI candidates, Olmschenk et al. 2021)
    for cat, key in (("I/358/veb", "gaia_eb"), ("B/vsx/vsx", "vsx"), ("J/AJ/167/203", "dtarps"),
                     ("J/AJ/161/273", "olmschenk")):
        try:
            r = Vizier(row_limit=5, columns=["**"]).query_region(c, radius=30 * u.arcsec, catalog=cat)
            out[key] = [{k: (str(row[k]) if k in row.colnames else "") for k in
                         ("Name", "Type", "Period", "Source", "Freq")} for row in r[0]] if len(r) else []
        except Exception as exc:
            out[key + "_error"] = repr(exc)[:100]
    return out


def crowding_test(tic: int, lc: dict, P: float, t0: float, dur: float) -> dict:
    """Does the transit depth track TESS-SPOC's crowding correction?

    PDCSAP flux is corrected for the light of neighbours in the aperture
    (CROWDSAP = fraction of the aperture flux that is the target, which
    changes from sector to sector with the aperture and the pointing). For
    a dip on the target the corrected depth is the same in every sector. For
    a dip on a neighbour the "corrected" depth is larger in sectors where
    more neighbour light is in the aperture, so it anti-correlates with
    CROWDSAP. A blend signature."""
    from scipy import stats
    from planet_hunt import multisector as ms
    meta = {d["meta"]["sector"]: d["meta"] for d in ffi.load_star(tic)}
    t, y, dy, sec = lc["t"], lc["flat"].astype(float), lc["dy"].astype(float), lc["sector"]
    intr = ms.in_tr(t, P, t0, dur)
    oot = ~ms.in_tr(t, P, t0, dur, 2.0)
    rows = []
    for s in np.unique(sec):
        m = meta.get(int(s), {})
        if m.get("provenance") != "TESS-SPOC" or not np.isfinite(m.get("crowdsap", np.nan)):
            continue
        a, b = (sec == s) & intr, (sec == s) & oot
        if a.sum() < 3 or b.sum() < 50:
            continue
        lev = np.average(y[b], weights=1 / dy[b] ** 2)
        rows.append((int(s), float(m["crowdsap"]), float(lev - np.average(y[a], weights=1 / dy[a] ** 2)) * 1e6,
                     float(1 / np.sqrt(np.sum(1 / dy[a] ** 2))) * 1e6))
    if len(rows) < 6:
        return {"ok": False, "verdict": "too few TESS-SPOC sectors"}
    c = np.array([r[1] for r in rows]); d = np.array([r[2] for r in rows])
    if np.ptp(c) < 0.03:
        return {"ok": False, "verdict": f"crowding nearly constant (CROWDSAP {c.min():.2f}-{c.max():.2f})",
                "per_sector": rows}
    rho, p = stats.spearmanr(c, d)
    bad = p < 1e-3 and rho < -0.4
    v = (f"BLEND SIGNATURE: depth anti-correlates with the crowding correction (Spearman {rho:+.2f}, p={p:.0e}; "
         f"CROWDSAP {c.min():.2f}-{c.max():.2f}): the dip is probably on a neighbour" if bad else
         f"depth does not track crowding (Spearman {rho:+.2f}, p={p:.2f}; CROWDSAP {c.min():.2f}-{c.max():.2f})")
    return {"ok": True, "rho": float(rho), "p": float(p), "blend": bool(bad), "verdict": v, "per_sector": rows}


def to_vet(res: pd.DataFrame, borderline: bool = False) -> pd.DataFrame:
    """TCEs that are not known planets/TOIs and carry no red flag."""
    flags = res["eb_flags"].fillna("")
    ok = res["tce"] & res["catalog_status"].isin(["new", "host_other"])
    return res[ok & ((flags == "") if not borderline else True)]


def cmd_vet(args):
    from planet_hunt import multisector as ms
    from planet_hunt.centroid import neighbour_check
    res = pd.read_csv(RUN / "results.csv")
    tg = pd.read_csv(RUN / "targets.csv").set_index("tic")
    (RUN / "vet").mkdir(exist_ok=True)
    if args.tic:
        todo = res[res["tic"] == args.tic]
    elif args.eb:     # method C: eclipsing-binary-like signals only (catalog context)
        # deep, eclipse-shaped (not sinusoidal), not a handful of outliers,
        # with a binary signature: second eclipse, odd/even, too big, V shape
        f = res["eb_flags"].fillna("")
        todo = res[res["tce"] & (res["depth_ppm"] > 3000) & (res["sin_ratio"] < 0.3)
                   & ~f.str.contains("dominated")
                   & f.str.contains("second eclipse|secondary|odd/even|too big|V-shaped")]
        todo = todo.sort_values("snr", ascending=False).drop_duplicates("tic")
    else:
        todo = to_vet(res, borderline=args.all_flags)
    # also the positive controls' known signals, to calibrate the centroid test
    if not args.tic and not args.eb:
        ctrl = res[res["tce"] & res["catalog_status"].isin(["known_planet", "known_toi"])]
        todo = pd.concat([todo, ctrl]).drop_duplicates(["tic", "signal"])
    print(f"[vet] {len(todo)} signals", flush=True)
    for r in todo.itertuples():
        t = tg.loc[r.tic]
        z = np.load(RUN / "lc" / f"{r.tic}_{r.signal}.npz")
        lc = {k: z[k] for k in ("t", "flat", "dy", "cad", "sector", "cx", "cy")}
        out = {"tic": int(r.tic), "signal": int(r.signal)}
        try:
            out["centroid"] = ms.centroid_multi(lc, r.period, r.t0, r.duration_h / 24, r.depth_ppm)
        except Exception as exc:
            out["centroid"] = {"ok": False, "verdict": f"failed: {exc!r}"[:150]}
        try:
            out["crowding"] = crowding_test(int(r.tic), lc, r.period, r.t0, r.duration_h / 24)
        except Exception as exc:
            out["crowding"] = {"ok": False, "verdict": f"failed: {exc!r}"[:150]}
        try:
            out["neighbours"] = neighbour_check(int(r.tic), float(t["ra"]), float(t["dec"]),
                                                float(t["Tmag"]), float(r.depth_ppm))
        except Exception as exc:
            out["neighbours"] = {"error": repr(exc)[:150]}
        try:
            out["spoc"] = spoc_tces_any(int(r.tic), float(r.period)) if t["n_spoc"] > 0 else \
                {"verdict": "no 2-min data: SPOC never searched this star"}
        except Exception as exc:
            out["spoc"] = {"verdict": f"SPOC lookup failed: {exc!r}"[:150]}
        if args.diffimage and r.catalog_status in ("new", "host_other"):
            from planet_hunt.diffimage import difference_image
            try:   # the 5 most recent sectors: 200-s cadence, sharpest in-transit images
                secs = sorted(set(lc["sector"].tolist()))[-5:]
                out["diffimage"] = difference_image(float(t["ra"]), float(t["dec"]), secs, r.period,
                                                    r.t0, r.duration_h / 24)
            except Exception as exc:
                out["diffimage"] = {"ok": False, "verdict": f"failed: {exc!r}"[:150]}
        try:
            out["gaia"] = gaia_context(t.get("GAIA"), float(t["ra"]), float(t["dec"]))
        except Exception as exc:
            out["gaia"] = {"error": repr(exc)[:150]}
        nb = out["neighbours"]
        print(f"TIC {r.tic} #{r.signal} P={r.period:.4f} {r.catalog_status}: {out['centroid']['verdict']}; "
              f"crowding: {out['crowding']['verdict']}; "
              f"neighbours able: {nb.get('n_able', '?')}/{nb.get('n_within', '?')}; {out['spoc']['verdict']}; "
              f"RUWE {out['gaia'].get('ruwe', '?')}; EB cats: {len(out['gaia'].get('gaia_eb', []))} Gaia, "
              f"{len(out['gaia'].get('vsx', []))} VSX; diff image: "
              f"{out.get('diffimage', {}).get('verdict', 'not run')}", flush=True)
        (RUN / "vet" / f"{r.tic}_{r.signal}.json").write_text(json.dumps(out, indent=1, default=float))


def cmd_reflag(args):
    """Re-apply the cross-signal flags (repeats, shared periods) to results.csv."""
    res = pd.read_csv(RUN / "results.csv")
    res = flag_repeats(res)
    res.to_csv(RUN / "results.csv", index=False)
    t = res[res["tce"] & res["catalog_status"].isin(["new", "host_other"])]
    print(f"[reflag] clean uncatalogued TCEs: {int((t['eb_flags'].fillna('') == '').sum())}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--n", type=int, default=300)
    s.add_argument("--controls", type=int, default=25)
    s.add_argument("--radius", type=float, default=12.0)
    s.add_argument("--min-sectors", type=int, default=8)
    s.add_argument("--seed", type=int, default=42)
    f = sub.add_parser("fetch")
    f.add_argument("--workers", type=int, default=8)
    f.add_argument("--max-gb", type=float, default=8.0)
    r = sub.add_parser("search")
    r.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    r.add_argument("--max-signals", type=int, default=3)
    r.add_argument("--pmax", type=float, default=30.0)
    r.add_argument("--tic", type=int, nargs="*")
    v = sub.add_parser("vet")
    v.add_argument("--tic", type=int)
    v.add_argument("--all-flags", action="store_true", help="also vet red-flagged uncatalogued TCEs")
    v.add_argument("--eb", action="store_true", help="eclipsing-binary-like signals (method C)")
    v.add_argument("--no-diffimage", dest="diffimage", action="store_false",
                   help="skip the TESScut difference-image test (slow: ~5 cutouts per signal)")
    sub.add_parser("reflag")
    args = ap.parse_args()
    {"select": cmd_select, "fetch": cmd_fetch, "reflag": cmd_reflag,
     "search": lambda a: cmd_search(a), "vet": lambda a: cmd_vet(a)}[args.cmd](args)


if __name__ == "__main__":
    main()
