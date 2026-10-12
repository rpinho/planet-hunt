"""Faint red dwarfs from TESS full-frame images: target selection and data.

Most nearby M dwarfs fainter than Tmag ~12 were never 2-minute targets, so
NASA's SPOC planet search (which only runs on 2-minute stars) never looked at
them. Their light curves exist anyway, extracted from the full-frame images
(FFIs) by two public high-level science products on MAST:

  * TESS-SPOC  (Caldwell et al. 2020): SPOC photometry of up to ~160,000 FFI
    stars per sector, PDCSAP flux, 30/10/3.3-min cadence. No public planet
    search is run on these files.
  * QLP        (Huang et al. 2020; Kunimoto et al. 2021): the MIT Quick-Look
    Pipeline, KSPSAP (spline-detrended) flux. QLP's own search feeds the TOI
    list, mostly for Tmag < 13.5.

For each sector this module takes TESS-SPOC when it exists and QLP otherwise.
The FITS files are reduced to small .npz files (time, flux, error, centroid)
and the FITS is then deleted, so hundreds of stars x 30+ sectors fit in a
few GB of ./cache.

Near the ecliptic poles TESS stares for ~1 year at a time (the continuous
viewing zones, CVZ), so a CVZ star has 30-70 sectors: many transits to stack.
"""
from __future__ import annotations

import concurrent.futures as cf
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from . import config

FFI_DIR = config.CACHE / "ffi"          # extracted per-sector .npz files
TMP_DIR = config.CACHE / "ffi_tmp"      # FITS while being read (deleted after)
POLES = {"south": (90.0, -66.5607), "north": (270.0, 66.5607)}   # ecliptic poles (RA, Dec)
TIC_KEEP = ["ID", "ra", "dec", "pmRA", "pmDEC", "plx", "Tmag", "Teff", "logg", "rad",
            "mass", "lum", "d", "eclat", "contratio", "numcont", "disposition", "GAIA",
            "gaiabp", "gaiarp", "GAIAmag", "Kmag", "lumclass"]


# --------------------------------------------------------------------------
# 1. M dwarfs near the ecliptic poles (TIC v8.2 via MAST)
# --------------------------------------------------------------------------
def mdwarfs_near_pole(pole: str, radius_deg: float = 12.0, teff=(2700, 3900),
                      rad=(0.08, 0.6), tmag=(11.0, 15.0), dmax: float = 100.0) -> pd.DataFrame:
    """TIC dwarfs with M-dwarf Teff/radius within `radius_deg` of a pole.

    The TIC's Teff and radius for cool dwarfs come from the Mann et al. (2015,
    2019) colour/absolute-K relations on Gaia distances, good to ~5-10%.
    The query takes ~2.5 min per pole and is cached."""
    path = config.CATALOG_DIR / (f"mdwarfs_{pole}_{radius_deg:g}deg_t{tmag[0]:g}-{tmag[1]:g}"
                                 f"_d{dmax:g}.csv")
    if path.exists():
        return pd.read_csv(path)
    from astroquery.mast import Catalogs, conf
    conf.timeout = 600
    ra, dec = POLES[pole]
    print(f"[tic] M dwarfs within {radius_deg} deg of the {pole} ecliptic pole ...", flush=True)
    tab = Catalogs.query_criteria(catalog="Tic", coordinates=f"{ra} {dec}",
                                  radius=radius_deg, Teff=list(teff), rad=list(rad),
                                  Tmag=list(tmag), d=[0.1, dmax], objType="STAR")
    df = tab.to_pandas()
    df = df[[c for c in TIC_KEEP if c in df.columns]]
    df = df[df["disposition"].isna() | ~df["disposition"].astype(str).isin(["SPLIT", "DUPLICATE", "ARTIFACT"])]
    df["pole"] = pole
    df.to_csv(path, index=False)
    return df


# --------------------------------------------------------------------------
# 2. Which light-curve products exist for each star? (one MAST query / 150 stars)
# --------------------------------------------------------------------------
def lightcurve_products(tics: list[int], chunk: int = 150) -> pd.DataFrame:
    """All TESS light-curve observations (SPOC 2-min, TESS-SPOC, QLP) for these
    TIC IDs. Cached per TIC in cache/mast_tables/lcprod.csv."""
    from astroquery.mast import Observations, conf
    conf.timeout = 300
    path = config.TABLE_DIR / "lcprod.csv"
    have = pd.read_csv(path) if path.exists() else pd.DataFrame(columns=["tic"])
    done = set(have["tic"].astype("int64")) if len(have) else set()
    todo = [int(t) for t in tics if int(t) not in done]
    new = []
    for i in range(0, len(todo), chunk):
        ids = [str(t) for t in todo[i:i + chunk]]
        for attempt in range(3):
            try:
                obs = Observations.query_criteria(
                    obs_collection=["TESS", "HLSP"], dataproduct_type="timeseries",
                    target_name=ids, provenance_name=["SPOC", "TESS-SPOC", "QLP"])
                break
            except Exception as exc:
                print(f"[mast] retry {attempt + 1}: {exc}", flush=True)
                time.sleep(10)
        else:
            continue
        o = obs["target_name", "provenance_name", "sequence_number", "t_exptime",
                "t_min", "t_max", "dataURL"].to_pandas()
        o["tic"] = pd.to_numeric(o["target_name"], errors="coerce").astype("Int64")
        o = o.dropna(subset=["tic"])
        # QLP and SPOC sometimes list a .txt or fast-cadence twin; keep LC FITS
        url = o["dataURL"].astype(str)
        o = o[url.str.endswith("_lc.fits") | url.str.endswith("_llc.fits")]
        o = o[~((o["provenance_name"] == "SPOC") & (o["t_exptime"] < 100))]
        # stars with no products at all still get a row, so they are not re-queried
        missing = set(int(t) for t in ids) - set(o["tic"].astype("int64"))
        if missing:
            o = pd.concat([o, pd.DataFrame({"tic": sorted(missing)})], ignore_index=True)
        new.append(o.drop(columns=["target_name"]))
        print(f"[mast] products for {min(i + chunk, len(todo))}/{len(todo)} stars", flush=True)
    if new:
        have = pd.concat([have] + new, ignore_index=True)
        have.to_csv(path, index=False)
    have["tic"] = have["tic"].astype("int64")
    return have[have["tic"].isin([int(t) for t in tics])]


def choose_sectors(prod: pd.DataFrame) -> pd.DataFrame:
    """One FFI light curve per (star, sector): TESS-SPOC if present, else QLP."""
    p = prod.dropna(subset=["dataURL"]).copy()
    p = p[p["provenance_name"].isin(["TESS-SPOC", "QLP"])]
    p["rank"] = p["provenance_name"].map({"TESS-SPOC": 0, "QLP": 1})
    p = p.sort_values(["tic", "sequence_number", "rank"])
    return p.drop_duplicates(["tic", "sequence_number"]).drop(columns="rank")


def coverage_table(prod: pd.DataFrame) -> pd.DataFrame:
    """Per star: number of FFI sectors (chosen), SPOC 2-min sectors, TESS-SPOC, QLP."""
    p = prod.dropna(subset=["dataURL"])
    out = pd.DataFrame({"tic": sorted(prod["tic"].astype("int64").unique())}).set_index("tic")
    for name in ("SPOC", "TESS-SPOC", "QLP"):
        out[f"n_{name.lower().replace('-', '')}"] = (
            p[p["provenance_name"] == name].groupby("tic")["sequence_number"].nunique())
    out["n_ffi"] = choose_sectors(prod).groupby("tic")["sequence_number"].nunique()
    return out.fillna(0).astype(int).reset_index()


# --------------------------------------------------------------------------
# 3. Download -> extract -> delete FITS
# --------------------------------------------------------------------------
def _npz_path(tic: int, sector: int) -> Path:
    return FFI_DIR / f"{tic}" / f"s{sector:04d}.npz"


def extract_fits(path: Path, provenance: str) -> dict:
    """Read one TESS-SPOC or QLP light curve into plain arrays.

    TESS-SPOC: PDCSAP_FLUX, QUALITY == 0 (strict: FFI data have few flags),
               centroid MOM_CENTR1/2.
    QLP:       KSPSAP_FLUX / DET_FLUX (QLP's detrended flux), QUALITY == 0, centroid
               SAP_X/SAP_Y.
    """
    from astropy.io import fits
    with fits.open(path, memmap=False) as hd:
        h0, h1 = hd[0].header, hd[1].header
        d = hd[1].data
        names = [n.upper() for n in d.columns.names]

        def col(*cands):
            for c in cands:
                if c in names:
                    return np.asarray(d[c], dtype=float)
            return np.full(len(d), np.nan)
        t = col("TIME")
        q = col("QUALITY")
        if provenance == "TESS-SPOC":
            f, e = col("PDCSAP_FLUX"), col("PDCSAP_FLUX_ERR")
            cx, cy = col("MOM_CENTR1"), col("MOM_CENTR2")
        else:
            # older QLP releases call it KSPSAP_FLUX, newer ones DET_FLUX
            f, e = col("KSPSAP_FLUX", "DET_FLUX"), col("KSPSAP_FLUX_ERR", "DET_FLUX_ERR")
            cx, cy = col("SAP_X"), col("SAP_Y")
            if np.all(~np.isfinite(e)):
                e = np.full_like(f, np.nan)
        ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
        med = np.nanmedian(f[ok]) if ok.any() else 1.0
        meta = {"sector": int(h0.get("SECTOR", 0) or 0), "camera": int(h0.get("CAMERA", 0) or 0),
                "ccd": int(h0.get("CCD", 0) or 0), "provenance": provenance,
                "tmag": float(h0.get("TESSMAG", np.nan) or np.nan),
                "crowdsap": float(h1.get("CROWDSAP", np.nan) or np.nan),
                "cadence_min": float(np.nanmedian(np.diff(t[ok])) * 1440) if ok.sum() > 2 else np.nan}
    err = e[ok] / med
    if not np.isfinite(err).any():
        err = np.full(ok.sum(), np.nan)
    return {"time": t[ok], "flux": (f[ok] / med).astype(np.float32),
            "flux_err": err.astype(np.float32), "cx": cx[ok].astype(np.float32),
            "cy": cy[ok].astype(np.float32), "meta": meta}


def fetch_sector(tic: int, sector: int, url: str, provenance: str, retries: int = 3) -> Path | None:
    """Download one FFI light curve, extract it to .npz and delete the FITS."""
    out = _npz_path(tic, sector)
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    tmp = TMP_DIR / f"{tic}_{sector}_{provenance}.fits"
    for attempt in range(retries):
        try:
            urllib.request.urlretrieve(config.MAST_DOWNLOAD + url, tmp)
            z = extract_fits(tmp, provenance)
            if len(z["time"]) < 50:
                tmp.unlink(missing_ok=True)
                return None
            meta = z.pop("meta")
            np.savez_compressed(out, **z, **{f"meta_{k}": v for k, v in meta.items()})
            tmp.unlink(missing_ok=True)   # scratch file this run created
            return out
        except Exception as exc:
            print(f"[dl] retry {attempt + 1} TIC {tic} S{sector}: {exc!r}"[:160], flush=True)
            time.sleep(3)
    tmp.unlink(missing_ok=True)
    return None


def cache_gb() -> float:
    return sum(f.stat().st_size for f in config.CACHE.rglob("*") if f.is_file()) / 1e9


def fetch_many(chosen: pd.DataFrame, workers: int = 8, max_cache_gb: float = 8.0) -> int:
    """Fetch every chosen (tic, sector) light curve; stops if ./cache > max_cache_gb."""
    jobs = [(int(r.tic), int(r.sequence_number), str(r.dataURL), str(r.provenance_name))
            for r in chosen.itertuples() if not _npz_path(int(r.tic), int(r.sequence_number)).exists()]
    print(f"[dl] {len(jobs)} FFI light curves to fetch ({len(chosen) - len(jobs)} cached)", flush=True)
    n = 0
    for i in range(0, len(jobs), 400):
        gb = cache_gb()
        if gb > max_cache_gb:
            print(f"[dl] cache is {gb:.1f} GB > {max_cache_gb} GB: stopping downloads", flush=True)
            break
        with cf.ThreadPoolExecutor(workers) as ex:
            for r in ex.map(lambda j: fetch_sector(*j), jobs[i:i + 400]):
                n += r is not None
        print(f"[dl] {min(i + 400, len(jobs))}/{len(jobs)} fetched (cache {cache_gb():.2f} GB)", flush=True)
    return n


def load_star(tic: int) -> list[dict]:
    """All extracted sectors of one star, sorted by sector."""
    out = []
    for p in sorted((FFI_DIR / str(tic)).glob("s*.npz")):
        z = np.load(p)
        d = {k: z[k] for k in ("time", "flux", "flux_err", "cx", "cy")}
        d["meta"] = {k[5:]: z[k].item() for k in z.files if k.startswith("meta_")}
        out.append(d)
    return out
