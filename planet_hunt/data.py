"""Target selection, downloading and cleaning TESS SPOC 2-minute light curves.

Everything here talks only to NASA MAST (public data). Results of slow queries
are cached under ./cache so re-runs are fast.
"""
from __future__ import annotations

import concurrent.futures as cf
import time
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import config  # noqa: F401  (must precede astropy imports)

import lightkurve as lk  # noqa: E402

lk.conf.cache_dir = str(config.CACHE / "lightkurve")


# --------------------------------------------------------------------------
# 1. Which stars did TESS observe at 2-minute cadence in this sector?
# --------------------------------------------------------------------------
def sector_observations(sector: int) -> pd.DataFrame:
    """All SPOC 2-min light-curve observations in a sector (one row per star).

    The MAST query takes a few minutes the first time, so the table is cached
    as CSV. Columns of interest: target_name (= TIC ID), dataURL (FITS URI).
    """
    path = config.TABLE_DIR / f"sector_{sector:04d}_spoc120.csv"
    if path.exists():
        df = pd.read_csv(path)
    else:
        from astroquery.mast import Observations

        print(f"[mast] querying SPOC 2-min observations for sector {sector} "
              "(first time can take ~4 min)...", flush=True)
        obs = Observations.query_criteria(
            obs_collection="TESS", dataproduct_type="timeseries",
            sequence_number=sector, provenance_name="SPOC", t_exptime=[119, 121])
        if len(obs) == 0:
            raise SystemExit(f"MAST has no SPOC 2-min light curves for sector {sector} (yet?)")
        df = obs["target_name", "obsid", "obs_id", "s_ra", "s_dec", "t_min",
                 "t_max", "dataURL", "proposal_id"].to_pandas()
        df.to_csv(path, index=False)
    df = df[df["dataURL"].astype(str).str.endswith("_lc.fits")]
    df["tic"] = pd.to_numeric(df["target_name"], errors="coerce").astype("Int64")
    return df.dropna(subset=["tic"]).drop_duplicates("tic").reset_index(drop=True)


# --------------------------------------------------------------------------
# 2. Stellar parameters from the TESS Input Catalog (TIC v8.2, via MAST)
# --------------------------------------------------------------------------
TIC_COLS = ["ID", "Tmag", "Teff", "logg", "rad", "mass", "lumclass",
            "contratio", "ra", "dec", "d"]


def tic_info(tics: list[int], chunk: int = 400) -> pd.DataFrame:
    """TIC rows for a list of TIC IDs, cached in cache/catalogs/tic_rows.csv."""
    from astroquery.mast import Catalogs

    cache = config.CATALOG_DIR / "tic_rows.csv"
    have = pd.read_csv(cache) if cache.exists() else pd.DataFrame(columns=TIC_COLS)
    have["ID"] = have["ID"].astype("int64")
    todo = sorted(set(int(t) for t in tics) - set(have["ID"]))
    new = []
    for i in range(0, len(todo), chunk):
        ids = [str(t) for t in todo[i:i + chunk]]
        for attempt in range(3):
            try:
                tab = Catalogs.query_criteria(catalog="Tic", ID=ids)
                break
            except Exception as exc:  # MAST hiccups happen; retry politely
                print(f"[tic] retry {attempt + 1}: {exc}", flush=True)
                time.sleep(5)
        else:
            continue
        df = tab.to_pandas()
        df = df[[c for c in TIC_COLS if c in df.columns]]
        df["ID"] = df["ID"].astype("int64")
        new.append(df)
        print(f"[tic] {min(i + chunk, len(todo))}/{len(todo)} stars", flush=True)
    if new:
        have = pd.concat([have] + new, ignore_index=True).drop_duplicates("ID")
        have.to_csv(cache, index=False)
    return have[have["ID"].isin([int(t) for t in tics])].reset_index(drop=True)


def select_targets(sector: int, n: int, tmag_max: float = 11.0,
                   n_controls: int = 15, toi: pd.DataFrame | None = None,
                   seed: int = 42) -> pd.DataFrame:
    """Pick ~n bright, quiet-ish main-sequence dwarfs observed in `sector`.

    Selection (all from the TIC):
      Tmag < tmag_max, 4.1 < logg, 3800 K < Teff < 6500 K, contratio < 0.1
      (little light from neighbours diluting the transit).
    "Quiet" cannot be judged before downloading; FGK dwarfs are the best proxy.

    Positive controls: up to `n_controls` of the targets are deliberately
    chosen among stars that already host a TOI with 0.5 d < P < 12.5 d and
    TOI SNR >= 10 (i.e. findable in one sector), so you can check the hunt
    re-discovers known signals. The rest are a random (seeded) draw, so the
    "new candidate" rate is not biased by the controls (they are labelled).
    """
    obs = sector_observations(sector)
    rng = np.random.default_rng(seed)
    all_tics = obs["tic"].astype("int64").to_numpy()

    def passes(df: pd.DataFrame) -> pd.Series:
        return ((df["Tmag"] < tmag_max) & (df["logg"] > 4.1)
                & df["Teff"].between(3800, 6500)
                & (df["contratio"].fillna(0) < 0.1))

    picked = []
    # positive controls: TOI hosts observed this sector
    if toi is not None and n_controls > 0:
        # only TOIs a one-sector search CAN find: >= 2 transits in ~27 d,
        # and not too shallow (TOI-reported SNR >= 10)
        findable = toi[toi["period"].between(0.5, 12.5)
                       & (toi["Planet SNR"].fillna(0) >= 10)]
        hosts = np.intersect1d(all_tics, findable["tic"].unique())
        info = tic_info(list(hosts))
        info = info[passes(info)]
        ctrl = info.sample(min(n_controls, len(info)), random_state=seed)
        ctrl["control"] = True
        picked.append(ctrl)
    # random draw from the rest, querying the TIC in shuffled chunks until full
    n_needed = n - sum(len(p) for p in picked)
    pool = rng.permutation(np.setdiff1d(all_tics, np.concatenate(
        [p["ID"].to_numpy() for p in picked]) if picked else []))
    found = []
    for i in range(0, len(pool), 800):
        info = tic_info(list(pool[i:i + 800]))
        found.append(info[passes(info)])
        if sum(len(f) for f in found) >= n_needed:
            break
    rand = pd.concat(found).sample(frac=1, random_state=seed).head(n_needed)
    rand["control"] = False
    picked.append(rand)
    sel = pd.concat(picked, ignore_index=True).rename(columns={"ID": "tic"})
    sel = sel.merge(obs[["tic", "dataURL"]].astype({"tic": "int64"}), on="tic")
    return sel.sort_values("Tmag").reset_index(drop=True)


# --------------------------------------------------------------------------
# 3. Download light-curve FITS files (parallel, resumable)
# --------------------------------------------------------------------------
def download(data_uri: str, retries: int = 3) -> Path:
    """Download one MAST data URI into cache/mast (skips if present)."""
    out = config.MAST_DIR / data_uri.rsplit("/", 1)[-1]
    if out.exists() and out.stat().st_size > 10_000:
        return out
    tmp = out.with_suffix(".part")
    for attempt in range(retries):
        try:
            urllib.request.urlretrieve(config.MAST_DOWNLOAD + data_uri, tmp)
            tmp.rename(out)
            return out
        except Exception as exc:
            print(f"[dl] retry {attempt + 1} {out.name}: {exc}", flush=True)
            time.sleep(3)
    raise RuntimeError(f"failed to download {data_uri}")


def download_many(uris: list[str], workers: int = 8) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    with cf.ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(download, u): u for u in uris}
        for k, fut in enumerate(cf.as_completed(futs), 1):
            try:
                paths[futs[fut]] = fut.result()
            except Exception as exc:
                print(f"[dl] giving up: {exc}", flush=True)
            if k % 25 == 0 or k == len(uris):
                print(f"[dl] {k}/{len(uris)} light curves", flush=True)
    return paths


# --------------------------------------------------------------------------
# 4. Clean + detrend
# --------------------------------------------------------------------------
@dataclass
class CleanLC:
    """A cleaned light curve as plain numpy arrays (easy to save / plot)."""
    tic: int
    sector: int
    time: np.ndarray        # BTJD = BJD - 2457000
    flux: np.ndarray        # PDCSAP flux, normalised to median 1
    flux_err: np.ndarray
    trend: np.ndarray       # Savitzky-Golay trend removed by flatten()
    flat: np.ndarray        # flux / trend: what the transit search sees
    star: dict = field(default_factory=dict)   # Tmag, Teff, logg, R*, M*

    @property
    def cadence(self) -> float:
        return float(np.median(np.diff(self.time)))


def star_from_header(meta: dict) -> dict:
    """Stellar parameters SPOC copies from the TIC into the FITS header.
    Mass is derived from logg and radius (M = g R^2 / G, solar units)."""
    def num(key):
        try:
            v = float(meta.get(key))
            return v if np.isfinite(v) else None
        except (TypeError, ValueError):
            return None
    star = {"Tmag": num("TESSMAG"), "Teff": num("TEFF"),
            "logg": num("LOGG"), "rad": num("RADIUS")}
    if star["logg"] and star["rad"]:
        star["mass"] = 10 ** (star["logg"] - 4.438) * star["rad"] ** 2
    else:
        star["mass"] = None
    return star


def load_lightcurve(path: str | Path, gap_days: float = 0.5,
                    edge_days: float = 0.25):
    """Read a SPOC light curve: steps 1-3 of the cleaning recipe.

    1. PDCSAP flux, default quality mask (drops flagged cadences), drop NaNs.
    2. Normalise to median 1.
    3. Mask big gaps: drop `edge_days` after the start of the data and after
       every gap longer than `gap_days` (thermal ramps after data downlinks).
    Returns (lightkurve.LightCurve, header-dict).
    """
    lc = lk.read(str(path), quality_bitmask="default").remove_nans().normalize()
    t = lc.time.value
    keep = np.ones(len(t), bool)
    starts = np.r_[t[0], t[1:][np.diff(t) > gap_days]]
    for s in starts:
        keep &= ~((t >= s) & (t < s + edge_days))
    return lc[keep], dict(lc.meta)


def in_transit(t: np.ndarray, period: float, t0: float, duration: float,
               widen: float = 1.0) -> np.ndarray:
    """Boolean mask of points within widen*duration/2 of a transit centre."""
    phase = (t - t0 + 0.5 * period) % period - 0.5 * period
    return np.abs(phase) < 0.5 * duration * widen


def detrend(lc, meta: dict, window_days: float = 1.0, gap_days: float = 0.5,
            transit: tuple[float, float, float] | None = None) -> CleanLC:
    """Steps 4-5 of the cleaning recipe.

    4. Flatten with a Savitzky-Golay filter `window_days` wide, iteratively
       ignoring >3-sigma outliers. A window much longer than a transit (1 d vs
       a few hours) removes stellar variability but keeps transits.
       `transit=(period, t0, duration)` excludes those points (widened 1.5x)
       from the trend fit -- used for a second, cleaner pass once BLS has
       found a signal, so the filter does not eat into the transit depth.
    5. Clip upward outliers (flares, cosmic rays) above 4 sigma. Downward
       points are kept: they could be transits.
    """
    t = lc.time.value
    cadence = float(np.median(np.diff(t)))
    wl = int(round(window_days / cadence)) | 1          # must be odd
    mask = in_transit(t, *transit, widen=1.5) if transit else None
    flat, trend = lc.flatten(window_length=wl, polyorder=2, sigma=3, niters=5,
                             break_tolerance=max(5, int(gap_days / cadence)),
                             mask=mask, return_trend=True)
    f = flat.flux.value.astype(float)
    resid = f - np.nanmedian(f)
    sig = 1.4826 * np.nanmedian(np.abs(resid))
    ok = np.isfinite(f) & (resid < 4 * sig)
    return CleanLC(
        tic=int(meta.get("TICID", 0) or 0), sector=int(meta.get("SECTOR", 0) or 0),
        time=t[ok], flux=lc.flux.value[ok].astype(float),
        flux_err=lc.flux_err.value[ok].astype(float),
        trend=trend.flux.value[ok].astype(float), flat=f[ok],
        star=star_from_header(meta))


def clean_lightcurve(path: str | Path, **kw) -> CleanLC:
    """Convenience: load_lightcurve + detrend in one call."""
    lc, meta = load_lightcurve(path)
    return detrend(lc, meta, **kw)
