"""Round 5, track A: hidden extra planets around small-star planet hosts.

    .venv/bin/python round5_hosts.py select           # targets.csv
    .venv/bin/python round5_hosts.py fetch            # SPOC 2-min (+TESS-SPOC FFI gap fill), shrunk to npz
    .venv/bin/python round5_hosts.py search           # one search per star, runs/round5_hosts/search/
    .venv/bin/python round5_hosts.py inject --n 300   # injection-recovery
    .venv/bin/python round5_hosts.py collect          # candidates.csv

Targets: every TOI host (TFOPWG not FP/FA) and every confirmed-planet host with
Teff < 4200 K or R* < 0.7 Rsun, Tmag < 14, and >= 8 sectors of SPOC 2-minute
light curves. The 2-minute sector list comes from MAST's per-sector download
scripts (cache/r5/twomin_index.csv.gz), so no MAST query is needed per star.

Per star (planet_hunt/hostsearch.py): every known TOI / planet is masked
(1.5 x duration); each sector is detrended with a transit-masked running
median (window chosen from the star's rotation); a semi-coherent BLS sums the
season periodograms over 0.5-30 d; the top peaks are refined coherently on all
data; statistics on the un-binned light curve with per-sector red-noise
inflation. Up to three signals per star (find, mask, search again).
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

from planet_hunt import config

R5 = config.CACHE / "r5"
LC_DIR = R5 / "lc"             # shrunk light curves, one npz per star
OUT = config.RUNS / "round5_hosts"
ARCHIVE = "https://archive.stsci.edu/missions/tess/tid"
MIN_FREE_GB = 20.0


def free_gb() -> float:
    return shutil.disk_usage("/").free / 1e9


# --------------------------------------------------------------------------
# Target selection
# --------------------------------------------------------------------------
def _num(x):
    return pd.to_numeric(x, errors="coerce")


def known_ephemerides(tic: int, toi: pd.DataFrame, ps: pd.DataFrame) -> list[dict]:
    """Every TOI (any disposition: an FP TOI on the star is still a signal to
    mask) and every transiting confirmed planet with an ephemeris."""
    out = []
    for r in toi[toi["tic"] == tic].itertuples():
        P, t0, d = _num(r.period), _num(r.epoch), _num(r.duration_h)
        out.append({"name": f"TOI-{r.TOI}", "period": P, "t0": t0 - 2457000 if t0 > 2e6 else t0,
                    "duration_h": d, "depth_ppm": _num(r.depth_ppm), "disp": str(r.tfop),
                    "src": "TOI"})
    for r in ps[ps["tic"] == tic].itertuples():
        P, t0, d = _num(r.pl_orbper), _num(r.pl_tranmid), _num(r.pl_trandur)
        if not np.isfinite(P):
            continue
        # skip if a TOI already carries it
        if any(np.isfinite(k["period"]) and abs(k["period"] / P - 1) < 0.002 for k in out):
            continue
        out.append({"name": r.pl_name, "period": P,
                    "t0": t0 - 2457000 if np.isfinite(t0) and t0 > 2e6 else t0,
                    "duration_h": d, "depth_ppm": np.nan, "disp": "confirmed",
                    "transiting": int(r.tran_flag) if r.tran_flag == r.tran_flag else 0,
                    "src": "pscomppars"})
    return out


def load_catalogs():
    toi = pd.read_csv(R5 / "toi_r5.csv", comment="#", low_memory=False)
    toi = toi.rename(columns={"TIC ID": "tic", "Period (days)": "period", "Epoch (BJD)": "epoch",
                              "Duration (hours)": "duration_h", "Depth (ppm)": "depth_ppm",
                              "TFOPWG Disposition": "tfop", "TESS Disposition": "tess_disp",
                              "TESS Mag": "tmag", "Stellar Eff Temp (K)": "teff",
                              "Stellar Radius (R_Sun)": "rstar", "Stellar Mass (M_Sun)": "mstar",
                              "RA": "ra", "Dec": "dec", "Planet Radius (R_Earth)": "rp"})
    toi["tic"] = toi["tic"].astype("int64")
    ps = pd.read_csv(R5 / "pscomp_r5.csv")
    ps["tic"] = _num(ps["tic_id"].astype(str).str.replace("TIC", "").str.strip()).astype("Int64")
    return toi, ps


def _deg(s):
    """ExoFOP RA/Dec are sexagesimal strings."""
    from astropy.coordinates import Angle
    import astropy.units as u
    return s


def cmd_select(a):
    from astropy.coordinates import SkyCoord, BarycentricTrueEcliptic
    import astropy.units as u
    toi, ps = load_catalogs()
    idx = pd.read_csv(R5 / "twomin_index.csv.gz")
    nsec = idx.groupby("tic")["sector"].nunique()
    # TOI hosts: drop hosts whose every TOI is FP/FA
    bad = toi["tfop"].isin(["FP", "FA"]) | (toi["tfop"].isna() & toi["tess_disp"].isin(["FP", "FA", "EB"]))
    good_tics = set(toi.loc[~bad, "tic"])
    rows = {}
    for r in toi[toi["tic"].isin(good_tics)].itertuples():
        if r.tic in rows:
            continue
        rows[r.tic] = {"tic": int(r.tic), "src": "TOI", "tmag": _num(r.tmag), "teff": _num(r.teff),
                       "rstar": _num(r.rstar), "mstar": _num(r.mstar), "ra": r.ra, "dec": r.dec,
                       "host_toi": int(np.floor(r.TOI))}
    for r in ps.dropna(subset=["tic"]).itertuples():
        t = int(r.tic)
        if t in rows:
            rows[t]["confirmed_host"] = r.hostname
            for k, v in (("teff", r.st_teff), ("rstar", r.st_rad), ("mstar", r.st_mass)):
                if not np.isfinite(rows[t][k]) and np.isfinite(_num(v)):
                    rows[t][k] = _num(v)
            continue
        rows[t] = {"tic": t, "src": "confirmed", "tmag": _num(r.sy_tmag), "teff": _num(r.st_teff),
                   "rstar": _num(r.st_rad), "mstar": _num(r.st_mass), "ra": r.ra, "dec": r.dec,
                   "host_toi": None, "confirmed_host": r.hostname}
    df = pd.DataFrame(rows.values())
    df["n2min"] = df["tic"].map(nsec).fillna(0).astype(int)
    small = (df["teff"] < 4200) | (df["rstar"] < 0.7)
    sel = df[small & (df["tmag"] < 14) & (df["n2min"] >= 8)].copy()
    # coordinates (ExoFOP gives sexagesimal strings, the archive degrees)
    def coord(r):
        try:
            return SkyCoord(float(r.ra) * u.deg, float(r.dec) * u.deg)
        except (TypeError, ValueError):
            return SkyCoord(str(r.ra), str(r.dec), unit=(u.hourangle, u.deg))
    cs = [coord(r) for r in sel.itertuples()]
    sel["ra"] = [c.ra.deg for c in cs]
    sel["dec"] = [c.dec.deg for c in cs]
    sel["elat"] = [c.transform_to(BarycentricTrueEcliptic()).lat.deg for c in cs]
    # fill missing masses with a main-sequence radius relation (M dwarfs: M ~ R)
    m = ~np.isfinite(sel["mstar"])
    sel.loc[m, "mstar"] = sel.loc[m, "rstar"].clip(0.1, 1.0)
    m = ~np.isfinite(sel["rstar"])
    sel.loc[m, "rstar"] = 0.6
    sel.loc[~np.isfinite(sel["mstar"]), "mstar"] = 0.6
    sel["known"] = [json.dumps(known_ephemerides(t, toi, ps), default=float) for t in sel["tic"]]
    sel["n_known"] = [len(json.loads(k)) for k in sel["known"]]
    sel = sel.sort_values(["n2min", "tmag"], ascending=[False, True]).reset_index(drop=True)
    OUT.mkdir(parents=True, exist_ok=True)
    sel.to_csv(OUT / "targets.csv", index=False)
    print(f"{len(df)} hosts with catalog entries; {small.sum()} small; "
          f"{len(sel)} selected (Tmag<14, >=8 2-min sectors)")
    print(sel["src"].value_counts().to_dict())
    print("sectors:", sel["n2min"].describe().round(1).to_dict())
    print("2-min files to fetch:", int(sel["n2min"].sum()))


# --------------------------------------------------------------------------
# Fetch: SPOC 2-min light curves from the static archive tree, shrunk
# --------------------------------------------------------------------------
def _download(url: str, path: Path, retries: int = 4) -> bool:
    """curl, not urllib: urllib reads the archive ~7x slower on this machine."""
    import subprocess
    for k in range(retries):
        r = subprocess.run(["curl", "-s", "-f", "-m", "300", "-o", str(path), "-w", "%{http_code}", url],
                           capture_output=True, text=True)
        if r.returncode == 0 and path.exists() and path.stat().st_size > 10000:
            return True
        path.unlink(missing_ok=True)
        if r.stdout.strip() == "404":
            return False
        time.sleep(5 * (k + 1))
    return False


def shrink_spoc(path: Path) -> dict | None:
    from astropy.io import fits
    with fits.open(path) as hd:
        d = hd[1].data
        h0, h1 = hd[0].header, hd[1].header
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d["PDCSAP_FLUX"], float)
        fe = np.asarray(d["PDCSAP_FLUX_ERR"], float)
        q = np.asarray(d["QUALITY"], np.int32)
        cx = np.asarray(d["MOM_CENTR1"], float)
        cy = np.asarray(d["MOM_CENTR2"], float)
        try:
            ap = hd["APERTURE"].data
            c0, r0 = int(hd["APERTURE"].header["CRVAL1P"]), int(hd["APERTURE"].header["CRVAL2P"])
            yy, xx = np.nonzero(ap & 2)
            aperture = np.column_stack([c0 + xx, r0 + yy]).astype(np.int16)
        except Exception:
            aperture = np.zeros((0, 2), np.int16)
        meta = {"sector": int(h0["SECTOR"]), "camera": int(h0.get("CAMERA", 0)), "ccd": int(h0.get("CCD", 0)),
                "crowdsap": float(h1.get("CROWDSAP", np.nan) or np.nan),
                "flfrcsap": float(h1.get("FLFRCSAP", np.nan) or np.nan),
                "cdpp1": float(h1.get("CDPP1_0", np.nan) or np.nan)}
    dump_t = t[((q & 32) > 0) | ((q & 8) > 0)]
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
    if ok.sum() < 500:
        return None
    med = np.nanmedian(f[ok])
    return {"time": t[ok], "flux": (f[ok] / med).astype(np.float32), "ferr": (fe[ok] / med).astype(np.float32),
            "cx": cx[ok].astype(np.float32), "cy": cy[ok].astype(np.float32), "dumps": dump_t[np.isfinite(dump_t)],
            "aperture": aperture, "meta": meta}


def fetch_star(tic: int, files: list[str]) -> dict:
    out_path = LC_DIR / f"{tic}.npz"
    if out_path.exists():
        return {"tic": tic, "status": "cached"}
    tmp = R5 / "tmp"
    tmp.mkdir(exist_ok=True)
    t = f"{tic:016d}"
    secs = []
    for fn in files:
        s = int(fn.split("-s")[1][:4])
        url = f"{ARCHIVE}/s{s:04d}/{t[0:4]}/{t[4:8]}/{t[8:12]}/{t[12:16]}/{fn}"
        cached = config.MAST_DIR / fn
        p = cached if cached.exists() else tmp / fn
        if not p.exists() and not _download(url, p):
            continue
        try:
            d = shrink_spoc(p)
        except Exception:
            d = None
        if p.parent == tmp:
            p.unlink(missing_ok=True)
        if d is not None:
            secs.append(d)
    if not secs:
        return {"tic": tic, "status": "empty"}
    blob = {}
    for i, d in enumerate(sorted(secs, key=lambda d: d["meta"]["sector"])):
        for k in ("time", "flux", "ferr", "cx", "cy", "dumps", "aperture"):
            blob[f"s{i}_{k}"] = d[k]
        blob[f"s{i}_meta"] = np.array(json.dumps(d["meta"]))
    blob["n"] = np.array(len(secs))
    part = out_path.with_suffix(".part.npz")
    np.savez_compressed(part, **blob)
    part.rename(out_path)          # atomic: a search never reads a half-written file
    return {"tic": tic, "status": "ok", "n": len(secs)}


def load_lc(tic: int) -> list[dict]:
    z = np.load(LC_DIR / f"{tic}.npz")
    out = []
    for i in range(int(z["n"])):
        d = {k: z[f"s{i}_{k}"] for k in ("time", "flux", "ferr", "cx", "cy", "dumps", "aperture")}
        d["meta"] = json.loads(str(z[f"s{i}_meta"]))
        d["meta"]["provenance"] = "SPOC-2min"
        out.append(d)
    return out


def cmd_fetch(a):
    LC_DIR.mkdir(parents=True, exist_ok=True)
    tg = pd.read_csv(OUT / "targets.csv")
    if a.limit:
        tg = tg.head(a.limit)
    idx = pd.read_csv(R5 / "twomin_index.csv.gz")
    idx = idx[idx["tic"].isin(set(tg["tic"]))]
    files = {t: sorted(g.drop_duplicates("sector")["fname"]) for t, g in idx.groupby("tic")}
    todo = [t for t in tg["tic"] if not (LC_DIR / f"{t}.npz").exists()]
    print(f"{len(todo)} stars to fetch, {sum(len(files[t]) for t in todo)} files; free {free_gb():.1f} GB", flush=True)
    done = 0
    with cf.ThreadPoolExecutor(a.workers) as ex:
        futs = {}
        it = iter(todo)
        for t in it:
            futs[ex.submit(fetch_star, t, files[t])] = t
            if len(futs) >= a.workers * 2:
                break
        while futs:
            fin, _ = cf.wait(futs, return_when=cf.FIRST_COMPLETED)
            for f in fin:
                t = futs.pop(f)
                try:
                    r = f.result()
                except Exception as exc:
                    r = {"tic": t, "status": f"error {exc!r}"[:120]}
                done += 1
                if done % 10 == 0 or r["status"] not in ("ok", "cached"):
                    print(f"[{done}/{len(todo)}] {r}  free {free_gb():.1f} GB", flush=True)
            if free_gb() < MIN_FREE_GB:
                print("DISK GUARD: free space below 20 GB, stopping fetch", flush=True)
                for f in futs:
                    f.cancel()
                break
            for t in it:
                futs[ex.submit(fetch_star, t, files[t])] = t
                break


# --------------------------------------------------------------------------
# Search
# --------------------------------------------------------------------------
def _target_row(tic: int) -> dict:
    tg = pd.read_csv(OUT / "targets.csv")
    return tg[tg["tic"] == tic].iloc[0].to_dict()


def search_one(row: dict, outdir: Path, inject: dict | None = None, max_signals: int = 3) -> dict:
    from planet_hunt import hostsearch as hs
    tic = int(row["tic"])
    secs = load_lc(tic)
    known = json.loads(row["known"])
    res = hs.search_star(secs, known, float(row["rstar"]), float(row["mstar"]), float(row["teff"])
                         if np.isfinite(row["teff"]) else 3500.0, inject=inject, max_signals=max_signals)
    res["tic"] = tic
    res["row"] = {k: row[k] for k in ("tic", "tmag", "teff", "rstar", "mstar", "ra", "dec", "elat",
                                     "n2min", "host_toi", "src")}
    return res


def _search_worker(args):
    row, outdir = args
    tic = int(row["tic"])
    path = outdir / f"{tic}.json"
    if path.exists():
        return tic, "cached"
    t0 = time.time()
    try:
        res = search_one(row, outdir)
        res["elapsed_s"] = time.time() - t0
        path.write_text(json.dumps(res, default=_js))
        return tic, f"ok {len(res.get('signals', []))} signals {time.time() - t0:.0f}s"
    except Exception as exc:
        (outdir / f"{tic}.err").write_text(traceback.format_exc())
        return tic, f"error {exc!r}"[:150]


def _js(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def cmd_search(a):
    outdir = OUT / "search"
    outdir.mkdir(parents=True, exist_ok=True)
    tg = pd.read_csv(OUT / "targets.csv")
    tg = tg[[(LC_DIR / f"{t}.npz").exists() for t in tg["tic"]]]
    if a.only:
        tg = tg[tg["tic"].isin(a.only)]
    if a.limit:
        tg = tg.head(a.limit)
    jobs = [(r, outdir) for r in tg.to_dict("records")]
    print(f"{len(jobs)} stars to search with {a.workers} processes", flush=True)
    with cf.ProcessPoolExecutor(a.workers) as ex:
        for k, (tic, msg) in enumerate(ex.map(_search_worker, jobs, chunksize=1)):
            print(f"[{k + 1}/{len(jobs)}] TIC {tic}: {msg}", flush=True)


# --------------------------------------------------------------------------
# Injection-recovery
# --------------------------------------------------------------------------
def _inject_worker(args):
    row, inj, outdir = args
    path = outdir / f"{inj['id']}.json"
    if path.exists():
        return inj["id"], "cached"
    t0 = time.time()
    try:
        res = search_one(row, outdir, inject=inj, max_signals=2)
        rec = res.get("injection_result", {})
        out = {"inj": inj, "tic": int(row["tic"]), "tmag": row["tmag"], "rstar": row["rstar"],
               "n2min": row["n2min"], "result": rec, "elapsed_s": time.time() - t0,
               "rotation": res.get("rotation"), "saturated": res.get("rotation_saturated")}
        path.write_text(json.dumps(out, default=_js))
        return inj["id"], f"{'REC' if rec.get('recovered') else 'miss'} {time.time() - t0:.0f}s"
    except Exception as exc:
        (outdir / f"{inj['id']}.err").write_text(traceback.format_exc())
        return inj["id"], f"error {exc!r}"[:150]


def cmd_inject(a):
    outdir = OUT / "inject"
    outdir.mkdir(parents=True, exist_ok=True)
    tg = pd.read_csv(OUT / "targets.csv")
    tg = tg[[(LC_DIR / f"{t}.npz").exists() for t in tg["tic"]]].reset_index(drop=True)
    rng = np.random.default_rng(a.seed)
    jobs = []
    for i in range(a.n):
        row = tg.iloc[rng.integers(len(tg))].to_dict()
        P = float(np.exp(rng.uniform(np.log(1.0), np.log(30.0))))
        rp = float(rng.uniform(0.5, 2.0))
        inj = {"id": f"inj{a.seed}_{i:04d}", "period": P, "rp_earth": rp, "b": float(rng.uniform(0, 0.8)),
               "phase": float(rng.uniform(0, 1))}
        jobs.append((row, inj, outdir))
    print(f"{len(jobs)} injections with {a.workers} processes", flush=True)
    with cf.ProcessPoolExecutor(a.workers) as ex:
        for k, (iid, msg) in enumerate(ex.map(_inject_worker, jobs, chunksize=1)):
            print(f"[{k + 1}/{len(jobs)}] {iid}: {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("select")
    p = sub.add_parser("fetch"); p.add_argument("--workers", type=int, default=6); p.add_argument("--limit", type=int)
    p = sub.add_parser("search"); p.add_argument("--workers", type=int, default=8)
    p.add_argument("--limit", type=int); p.add_argument("--only", type=int, nargs="*")
    p = sub.add_parser("inject"); p.add_argument("--workers", type=int, default=8)
    p.add_argument("--n", type=int, default=300); p.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    {"select": cmd_select, "fetch": cmd_fetch, "search": cmd_search, "inject": cmd_inject}[a.cmd](a)


if __name__ == "__main__":
    main()
