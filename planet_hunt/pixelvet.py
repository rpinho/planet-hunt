"""Pixel-level vetting of a transit candidate from TESS full-frame-image cutouts.

`diffimage.py` asks one coarse question (does the light centre of the
difference image sit near the star?). This module does the full job, using
every sector the star was observed in:

1. Scene. Every Gaia DR3 source in and around an 11 x 11 pixel TESScut cutout,
   with an estimated TESS magnitude (Stassun et al. 2019 colour relation) and
   proper motion carried to each sector's epoch. A pixel-integrated
   two-Gaussian PSF (core + halo) is fitted to the out-of-transit image of
   each sector, with the relative fluxes of all stars fixed by Gaia and a free
   pointing offset, PSF width, flux scale and sky level. That gives, for every
   pixel, how much light each star puts there.
2. Difference image. For every transit: per-pixel local linear baseline from
   the frames just before and after, minus the in-transit frames. Averaged
   over a sector. Its noise comes from the same operation at ~30 shifted
   (null) ephemerides, so red noise and scattered light are included.
3. Localisation. One point source at a free sky offset (xi, eta) from the
   target, same PSF as the scene, per-sector amplitude and constant, fitted to
   all sectors at once (sectors have different roll angles, so PSF errors
   average out). Gives an offset with an error ellipse (rescaled by the
   reduced chi^2 and checked by a sector jackknife), and the chi^2 of
   "source = star j" for every nearby Gaia star.
4. Flux budget. For each neighbour: the eclipse depth it would need,
   observed depth x (total flux in the aperture) / (its flux in the aperture).
   More than 100% is impossible.
5. Depth vs aperture. Depth in rings of pixels ranked by the target's PSF
   (1 pixel, the next 3, ...), compared with what each hypothesis predicts:
   on-target dips have a constant dilution-corrected depth, a neighbour's dip
   grows where the neighbour's light is.

Calibrate on known planets and known off-target false positives before
trusting a verdict (deep_vet.py does).
"""
from __future__ import annotations

import concurrent.futures as cf
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special

from . import config
from .diffimage import CUT_DIR

PIX_ARCSEC = 21.0
FIT_RADIUS_PX = 3.5


# --------------------------------------------------------------------------
# Cutouts
# --------------------------------------------------------------------------
def cut_dir(ra: float, dec: float) -> Path:
    return CUT_DIR / f"{ra:.5f}_{dec:+.5f}"


def _readable(path: Path) -> bool:
    from astropy.io import fits
    try:
        with fits.open(path) as hd:
            return hd[1].data["FLUX"].shape[0] > 100
    except Exception:
        return False


def available_sectors(ra: float, dec: float) -> list[int]:
    from astroquery.mast import Tesscut
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import time
    for attempt in range(5):
        try:
            t = Tesscut.get_sectors(coordinates=SkyCoord(ra * u.deg, dec * u.deg))
            return sorted(int(s) for s in t["sector"])
        except Exception as exc:          # MAST drops connections now and then
            print(f"[tesscut] sector list attempt {attempt + 1}: {exc!r}"[:140], flush=True)
            time.sleep(10 * (attempt + 1))
    raise RuntimeError("TESScut sector list unavailable")


def fetch_cutouts(ra: float, dec: float, sectors: list[int] | None = None, size: int = 11,
                  workers: int = 3) -> list[tuple[int, Path]]:
    """Download (or reuse) TESScut cutouts. Each sector downloads into its own
    temporary folder: astroquery names the zip after the current second, so
    parallel downloads into one folder overwrite each other."""
    from astroquery.mast import Tesscut
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    d = cut_dir(ra, dec)
    d.mkdir(parents=True, exist_ok=True)
    c = SkyCoord(ra * u.deg, dec * u.deg)
    if sectors is None:
        sectors = available_sectors(ra, dec)

    def one(s):
        for p in sorted(d.glob(f"*s{s:04d}-*_{size}x{size}_astrocut.fits")):
            if _readable(p):
                return s, p
            p.unlink()        # truncated download from an earlier run
        for attempt in range(3):
            tmp = Path(tempfile.mkdtemp(dir=d, prefix=f".s{s:04d}_"))
            try:
                Tesscut.download_cutouts(coordinates=c, size=size, sector=s, path=str(tmp))
                got = sorted(tmp.glob(f"*s{s:04d}-*.fits"))
                if got and _readable(got[0]):
                    dest = d / got[0].name
                    shutil.move(str(got[0]), dest)
                    return s, dest
            except Exception as exc:
                print(f"[tesscut] S{s} attempt {attempt + 1}: {exc!r}"[:140], flush=True)
            finally:
                shutil.rmtree(tmp, ignore_errors=True)   # temp folder created above
        return s, None
    out = []
    with cf.ThreadPoolExecutor(workers) as ex:
        for s, p in ex.map(one, sectors):
            if p is not None:
                out.append((s, p))
            else:
                print(f"[tesscut] S{s}: no cutout", flush=True)
    return sorted(out)


# --------------------------------------------------------------------------
# Gaia DR3 scene
# --------------------------------------------------------------------------
def gaia_tmag(g, bp_rp):
    """TESS magnitude from Gaia G and BP-RP (Stassun et al. 2019, eq. 1);
    G - 0.43 when there is no colour."""
    c = np.clip(np.asarray(bp_rp, float), -0.2, 3.5)
    t = g - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(c), t, np.asarray(g, float) - 0.43)


def gaia_scene(ra: float, dec: float, radius_arcsec: float = 200.0, gaia_id=None,
               tmag: float | None = None) -> pd.DataFrame:
    """Gaia DR3 sources around a target (cached). Adds Tmag, separation and
    position angle from the target (epoch 2016), and flags the target row."""
    path = config.CATALOG_DIR / "gaia_scene" / f"{ra:.5f}_{dec:+.5f}_{int(radius_arcsec)}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        df = pd.read_csv(path)
    else:
        from astroquery.vizier import Vizier
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag", "BPmag", "RPmag", "BP-RP", "Plx",
                            "e_Plx", "pmRA", "pmDE", "RUWE", "Teff", "Dist"], row_limit=-1)
        r = v.query_region(SkyCoord(ra * u.deg, dec * u.deg), radius=radius_arcsec * u.arcsec,
                           catalog="I/355/gaiadr3")
        df = r[0].to_pandas().rename(columns={"BP-RP": "bp_rp", "RA_ICRS": "ra", "DE_ICRS": "dec",
                                              "Source": "source_id"})
        df.to_csv(path, index=False)
    df = df.copy()
    df["tmag"] = gaia_tmag(df["Gmag"], df["bp_rp"])
    for c in ("pmRA", "pmDE"):
        df[c] = df[c].fillna(0.0)
    if gaia_id is not None and pd.notna(gaia_id) and (df["source_id"] == int(gaia_id)).any():
        it = int(np.flatnonzero(df["source_id"].to_numpy() == int(gaia_id))[0])
    else:   # nearest source, penalising a magnitude mismatch
        sep = np.hypot((df["ra"] - ra) * np.cos(np.radians(dec)), df["dec"] - dec) * 3600
        dm = (df["tmag"] - tmag) if tmag is not None else 0
        it = int(np.argmin((sep / 2.0) ** 2 + (dm / 0.5) ** 2))
    df["is_target"] = False
    df.iloc[it, df.columns.get_loc("is_target")] = True
    r0 = df.iloc[it]
    df["xi"] = (df["ra"] - r0["ra"]) * np.cos(np.radians(r0["dec"])) * 3600     # arcsec east
    df["eta"] = (df["dec"] - r0["dec"]) * 3600                                  # arcsec north
    df["sep"] = np.hypot(df["xi"], df["eta"])
    df["pa_deg"] = np.degrees(np.arctan2(df["xi"], df["eta"])) % 360
    df["rel_flux"] = 10 ** (-0.4 * (df["tmag"] - r0["tmag"]))
    return df.sort_values("sep").reset_index(drop=True)


def _at_epoch(sc: pd.DataFrame, year: float):
    dt = year - 2016.0
    dec = sc["dec"] + sc["pmDE"] * dt / 3.6e6
    ra = sc["ra"] + sc["pmRA"] * dt / 3.6e6 / np.cos(np.radians(sc["dec"]))
    return ra.to_numpy(), dec.to_numpy()


# --------------------------------------------------------------------------
# Cutout -> arrays
# --------------------------------------------------------------------------
def load_cutout(path: Path, scattered_light_mad: float = 6.0) -> dict:
    """Good frames of one cutout, with a per-frame sky level removed (the 20th
    percentile of the cutout; the scene fit adds a constant on top)."""
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(path, memmap=False) as hd:
        d = hd[1].data
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d["FLUX"], np.float32)
        q = np.asarray(d["QUALITY"])
        w = WCS(hd[2].header)
        h0 = hd[0].header
        sector, cam, ccd = int(h0.get("SECTOR", 0)), int(h0.get("CAMERA", 0)), int(h0.get("CCD", 0))
    ok = np.isfinite(t) & (q == 0) & np.isfinite(f).all(axis=(1, 2))
    t, f = t[ok], f[ok]
    sky = np.percentile(f.reshape(len(t), -1), 20, axis=1)
    med = np.median(sky)
    mad = 1.4826 * np.median(np.abs(sky - med)) + 1e-6
    keep = sky < med + scattered_light_mad * mad
    t, f, sky = t[keep], f[keep] - sky[keep, None, None], sky[keep]
    cad = float(np.median(np.diff(t)))
    year = 2000.0 + (np.mean(t) + 2457000.0 - 2451545.0) / 365.25
    return {"t": t, "f": f, "sky": sky, "wcs": w, "sector": sector, "camera": cam, "ccd": ccd,
            "cad": cad, "year": year, "path": str(path)}


def smear_factor(dur: float, cad: float, n: int = 400) -> float:
    """Mean fractional coverage of a box transit of length `dur` by exposures
    of length `cad` whose centres fall inside the transit window: the factor
    by which finite exposures dilute an in-window mean depth."""
    u = (np.arange(n) + 0.5) / n * dur - dur / 2
    ov = np.clip(np.minimum(u + cad / 2, dur / 2) - np.maximum(u - cad / 2, -dur / 2), 0, None) / cad
    return float(ov.mean())


# --------------------------------------------------------------------------
# Difference images
# --------------------------------------------------------------------------
def _windows(t, P, t0, dur, cad):
    ph = (t - t0 + 0.5 * P) % P - 0.5 * P
    ep = np.round((t - t0) / P).astype(int)
    half = 0.5 * dur
    o_lo = half + max(cad, 0.25 * dur)
    o_hi = half + max(2.0 * dur, 5 * cad)
    return ph, ep, np.abs(ph) < half, (np.abs(ph) > o_lo) & (np.abs(ph) < o_hi), o_hi


def transit_diff(t, f, P, t0, dur, cad, exclude=None):
    """Mean out-minus-in image over all transits in t (local linear baseline
    per pixel from the frames either side). Positive = light lost in transit."""
    ph, ep, win_in, win_out, _ = _windows(t, P, t0, dur, cad)
    if exclude is not None:
        win_in &= ~exclude
        win_out &= ~exclude
    ny, nx = f.shape[1:]
    D, W = [], []
    for e in np.unique(ep[win_in]):
        mi, mo = win_in & (ep == e), win_out & (ep == e)
        if mi.sum() < 1 or (mo & (ph < 0)).sum() < 2 or (mo & (ph > 0)).sum() < 2:
            continue
        X = np.column_stack([np.ones(mo.sum()), ph[mo]])
        coef = np.linalg.lstsq(X, f[mo].reshape(mo.sum(), -1), rcond=None)[0]
        pred = np.array([1.0, ph[mi].mean()]) @ coef
        D.append(pred - f[mi].reshape(mi.sum(), -1).mean(0))
        W.append(mi.sum())
    if not D:
        return None, 0
    W = np.asarray(W, float)
    return ((W[:, None] * np.asarray(D)).sum(0) / W.sum()).reshape(ny, nx), len(D)


def sector_diff(cut: dict, P: float, t0: float, dur: float, n_null: int = 30,
                exclude_windows=()) -> dict | None:
    """Difference image, its noise (from shifted null ephemerides) and the
    out-of-transit image for one sector. `exclude_windows`: (P, t0, dur) of
    other signals whose frames are left out."""
    t, f, cad = cut["t"], cut["f"], cut["cad"]
    excl = np.zeros(len(t), bool)
    for (P2, t02, d2) in exclude_windows:
        excl |= np.abs((t - t02 + 0.5 * P2) % P2 - 0.5 * P2) < 0.75 * d2 + cad
    diff, n_tr = transit_diff(t, f, P, t0, dur, cad, excl)
    if diff is None or n_tr < 1:
        return None
    _, _, win_in, _, o_hi = _windows(t, P, t0, dur, cad)
    # null ephemerides: same P, shifted away from the transit and from phase 0.5
    span = P - 2 * (o_hi + dur)
    nulls = []
    if span > 0:
        for s in np.linspace(o_hi + dur, P - o_hi - dur, n_null + 2)[1:-1]:
            if abs(s - 0.5 * P) < o_hi + dur and P > 4 * (o_hi + dur):
                continue
            d_, n_ = transit_diff(t, f, P, t0 + s, dur, cad, excl | win_in)
            if d_ is not None and n_ >= max(1, n_tr // 2):
                nulls.append(d_)
    nulls = np.asarray(nulls)
    if len(nulls) >= 8:
        err = 1.4826 * np.median(np.abs(nulls - np.median(nulls, 0)), 0)
        err = np.maximum(err, 0.5 * np.std(nulls, 0))
    else:   # fall back on the frame-to-frame scatter
        r = np.diff(f, axis=0) / np.sqrt(2)
        sig = 1.4826 * np.median(np.abs(r), 0)
        err = sig / np.sqrt(max(win_in.sum(), 1))
    oot = np.median(f[~win_in & ~excl], axis=0)
    # pixels off the detector (constant) or without noise get no weight
    bad = (np.std(f, axis=0) <= 0) | ~np.isfinite(err) | (err <= 0)
    err = np.where(bad, np.inf, err)
    diff = np.where(bad, 0.0, diff)
    return {"bad": bad, "sector": cut["sector"], "camera": cut["camera"], "ccd": cut["ccd"], "cad": cad,
            "year": cut["year"], "diff": diff, "err": err, "oot": oot, "n_tr": n_tr,
            "n_in": int(win_in.sum()), "nulls": nulls, "smear": smear_factor(dur, cad)}


# --------------------------------------------------------------------------
# PSF and scene model
# --------------------------------------------------------------------------
_SUB = 4   # sub-pixels per pixel side


def _subgrid(shape):
    ny, nx = shape
    o = (np.arange(_SUB) + 0.5) / _SUB - 0.5
    yy, xx = np.mgrid[0:ny, 0:nx]
    full = (ny, nx, _SUB, _SUB)
    sx = np.broadcast_to(xx[..., None, None] + o[None, None, None, :], full).ravel()
    sy = np.broadcast_to(yy[..., None, None] + o[None, None, :, None], full).ravel()
    return sx, sy


def psf_images(x, y, shape, p) -> np.ndarray:
    """Unit-flux pixel-integrated PSF images (n, ny, nx) for sources at
    pixel positions x, y. p = (sx, sy, theta, fh): rotated elliptical Gaussian
    core plus a circular halo 3x wider holding fraction fh of the light."""
    sx_, sy_, th, fh = p
    gx, gy = _subgrid(shape)
    x = np.atleast_1d(np.asarray(x, float))[:, None]
    y = np.atleast_1d(np.asarray(y, float))[:, None]
    dx, dy = gx[None, :] - x, gy[None, :] - y
    c, s = np.cos(th), np.sin(th)
    u, v = c * dx + s * dy, -s * dx + c * dy
    core = np.exp(-0.5 * ((u / sx_) ** 2 + (v / sy_) ** 2)) / (2 * np.pi * sx_ * sy_)
    sh = 3.0 * np.sqrt(sx_ * sy_)
    halo = np.exp(-0.5 * (dx ** 2 + dy ** 2) / sh ** 2) / (2 * np.pi * sh ** 2)
    img = ((1 - fh) * core + fh * halo) / _SUB ** 2
    ny, nx = shape
    return img.reshape(len(x), ny, nx, _SUB * _SUB).sum(-1)


def _unpack(q):
    dx, dy, lsx, lsy, th, lfh, lsc, c0 = q
    return dx, dy, (np.exp(lsx), np.exp(lsy), th, 1 / (1 + np.exp(-lfh))), np.exp(lsc), c0


def fit_scene(oot: np.ndarray, err_img: np.ndarray, px: np.ndarray, py: np.ndarray,
              rel: np.ndarray) -> dict:
    """Fit pointing offset, PSF shape, flux scale and sky to the
    out-of-transit image with all Gaia stars at fixed relative fluxes."""
    shape = oot.shape
    inside = (px > -4) & (px < shape[1] + 3) & (py > -4) & (py < shape[0] + 3) & (rel > 10 ** (-0.4 * 7))
    px, py, rel = px[inside], py[inside], rel[inside]
    c0 = float(np.percentile(oot, 10))
    sc0 = max(float(oot.sum() - c0 * oot.size), 1.0) / max(float(rel.sum()), 1e-3)
    sig = np.sqrt(err_img ** 2 + 1.0)

    def model(q):
        dx, dy, p, sc, c = _unpack(q)
        return c + sc * np.tensordot(rel, psf_images(px + dx, py + dy, shape, p), axes=1)

    def resid(q):
        m = model(q)
        r = (oot - m) / np.sqrt(sig ** 2 + (0.03 * np.abs(m)) ** 2)
        return np.where(np.isfinite(sig), r, 0.0).ravel()
    best = None
    for th0 in (0.0, 0.8, -0.8):
        q0 = [0, 0, np.log(0.7), np.log(0.7), th0, -1.5, np.log(sc0), c0]
        r = optimize.least_squares(resid, q0, loss="soft_l1", f_scale=3.0,
                                   bounds=([-1.5, -1.5, np.log(0.3), np.log(0.3), -np.pi, -6, -20, -1e5],
                                           [1.5, 1.5, np.log(2.5), np.log(2.5), np.pi, 2, 30, 1e5]))
        if best is None or r.cost < best.cost:
            best = r
    dx, dy, p, sc, c = _unpack(best.x)
    m = model(best.x)
    g = np.isfinite(sig)
    return {"dx": float(dx), "dy": float(dy), "psf": tuple(float(v) for v in p), "scale": float(sc),
            "c0": float(c), "model": m,
            "frac_resid": float(np.sum(np.abs(oot - m)[g]) / np.sum(np.abs(oot - c)[g]))}


# --------------------------------------------------------------------------
# Per-sector preparation
# --------------------------------------------------------------------------
def prepare_sector(cut: dict, scene: pd.DataFrame, P: float, t0: float, dur: float,
                   exclude_windows=()) -> dict | None:
    sd = sector_diff(cut, P, t0, dur, exclude_windows=exclude_windows)
    if sd is None:
        return None
    ra_e, dec_e = _at_epoch(scene, cut["year"])
    px, py = cut["wcs"].world_to_pixel_values(ra_e, dec_e)
    sd["px"], sd["py"] = np.asarray(px, float), np.asarray(py, float)
    n_oot = max(len(cut["t"]) - sd["n_in"], 1)
    oot_err = 1.4826 * np.median(np.abs(np.diff(cut["f"], axis=0)), 0) / np.sqrt(2 * n_oot)
    oot_err = np.where(sd["bad"], np.inf, oot_err)
    sd["fit"] = fit_scene(sd["oot"], oot_err, sd["px"], sd["py"], scene["rel_flux"].to_numpy())
    # linear map from sky offsets (arcsec E, N of the target) to pixel offsets
    it = int(np.flatnonzero(scene["is_target"].to_numpy())[0])
    r0, d0 = ra_e[it], dec_e[it]
    e = 10.0 / 3600
    x0, y0 = cut["wcs"].world_to_pixel_values(r0, d0)
    xe, ye = cut["wcs"].world_to_pixel_values(r0 + e / np.cos(np.radians(d0)), d0)
    xn, yn = cut["wcs"].world_to_pixel_values(r0, d0 + e)
    sd["J"] = np.array([[xe - x0, xn - x0], [ye - y0, yn - y0]]) / 10.0   # pixels per arcsec
    sd["tpix"] = np.array([float(x0) + sd["fit"]["dx"], float(y0) + sd["fit"]["dy"]])
    sd["it"] = it
    sd["star_flux"] = sd["fit"]["scale"] * scene["rel_flux"].to_numpy()   # counts per star
    # localisation uses pixels within FIT_RADIUS of the target only: variable
    # stars at the cutout edge must not steer a fit meant to separate stars
    # within ~1 pixel
    yy, xx = np.mgrid[0:sd["diff"].shape[0], 0:sd["diff"].shape[1]]
    sd["fitmask"] = (np.hypot(xx - sd["tpix"][0], yy - sd["tpix"][1]) < FIT_RADIUS_PX) & ~sd["bad"]
    # quality: target fraction in its top-4 pixels, distance to off-detector pixels
    fit = sd["fit"]
    timg = fit["scale"] * psf_images(sd["px"][[it]] + fit["dx"], sd["py"][[it]] + fit["dy"],
                                     sd["diff"].shape, fit["psf"])[0]
    order = rank_pixels(timg, sd["bad"])[:4]
    tot = np.clip(sd["oot"] - fit["c0"], 1e-3, None).ravel()[order].sum()
    sd["target_frac"] = float(timg.ravel()[order].sum() / tot)
    badd = np.hypot(xx - sd["tpix"][0], yy - sd["tpix"][1])[sd["bad"]]
    sd["bad_dist_px"] = float(badd.min()) if badd.size else np.inf
    sd["quality_ok"] = bool(0.2 <= sd["target_frac"] <= 1.5 and sd["bad_dist_px"] > 2.5
                            and fit["frac_resid"] < 0.5)
    return sd


# --------------------------------------------------------------------------
# Localisation
# --------------------------------------------------------------------------
def _chi2_at(sectors: list[dict], xi: np.ndarray, eta: np.ndarray, per_sector: bool = False):
    """Profiled chi^2 (amplitude + constant per sector) of a point source at
    sky offsets (xi, eta) [arcsec] from the target, summed over sectors."""
    xi, eta = np.atleast_1d(xi), np.atleast_1d(eta)
    tot = np.zeros(len(xi))
    amps = []
    for sd in sectors:
        x = sd["tpix"][0] + sd["J"][0, 0] * xi + sd["J"][0, 1] * eta
        y = sd["tpix"][1] + sd["J"][1, 0] * xi + sd["J"][1, 1] * eta
        P = psf_images(x, y, sd["diff"].shape, sd["fit"]["psf"]).reshape(len(xi), -1)
        d = sd["diff"].ravel()
        w = np.where(sd["fitmask"].ravel(), 1 / sd["err"].ravel() ** 2, 0.0)
        Sw, Sp, Spp = w.sum(), P @ w, (P ** 2) @ w
        Sd, Spd, Sdd = (w * d).sum(), P @ (w * d), (w * d * d).sum()
        det = Spp * Sw - Sp ** 2
        A = (Spd * Sw - Sp * Sd) / det
        c = (Spp * Sd - Sp * Spd) / det
        chi = Sdd - A * Spd - c * Sd
        # a dip source loses light: no negative amplitudes (then constant only)
        neg = A < 0
        chi = np.where(neg, Sdd - Sd ** 2 / Sw, chi)
        A = np.where(neg, 0.0, A)
        tot += chi
        amps.append((A, np.sqrt(Sw / det)))
    return (tot, amps) if per_sector else tot


def localize(sectors: list[dict], scene: pd.DataFrame, radius: float = 60.0, step: float = 1.5,
             hyp_radius: float = 63.0) -> dict:
    """Joint multi-sector point-source fit to the difference images."""
    g = np.arange(-radius, radius + step / 2, step)
    XI, ETA = np.meshgrid(g, g)
    chi = np.concatenate([_chi2_at(sectors, a, b) for a, b in
                          zip(np.array_split(XI.ravel(), 40), np.array_split(ETA.ravel(), 40))])
    k = int(np.argmin(chi))
    fun = lambda q: float(_chi2_at(sectors, np.array([q[0]]), np.array([q[1]]))[0])
    r = optimize.minimize(fun, [XI.ravel()[k], ETA.ravel()[k]], method="Nelder-Mead",
                          options={"xatol": 0.02, "fatol": 1e-4})
    best = r.x
    chi_min = r.fun
    npix = sum(int(sd["fitmask"].sum()) for sd in sectors)
    dof = npix - 2 * len(sectors) - 2
    red = max(1.0, chi_min / dof)
    # curvature -> covariance (chi^2 rescaled by the reduced chi^2)
    h = 1.0
    H = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            e_i, e_j = np.eye(2)[i] * h, np.eye(2)[j] * h
            H[i, j] = (fun(best + e_i + e_j) - fun(best + e_i - e_j) - fun(best - e_i + e_j)
                       + fun(best - e_i - e_j)) / (4 * h * h)
    try:
        cov = 2 * np.linalg.inv(H) * red
    except np.linalg.LinAlgError:
        cov = np.full((2, 2), np.nan)
    # jackknife over sectors (catches sector-to-sector systematics)
    jk = []
    if len(sectors) >= 4:
        for i in range(len(sectors)):
            sub = sectors[:i] + sectors[i + 1:]
            f2 = lambda q: float(_chi2_at(sub, np.array([q[0]]), np.array([q[1]]))[0])
            jk.append(optimize.minimize(f2, best, method="Nelder-Mead",
                                        options={"xatol": 0.05, "fatol": 1e-3}).x)
        jk = np.array(jk)
        n = len(jk)
        cov_jk = (n - 1) / n * np.cov((jk - jk.mean(0)).T, bias=True) * n / max(n - 1, 1) * (n - 1)
        cov_jk = (n - 1) / n * ((jk - jk.mean(0)).T @ (jk - jk.mean(0)))
    else:
        cov_jk = np.full((2, 2), np.nan)
    cov_use = cov if not np.all(np.isfinite(cov_jk)) else np.where(
        np.trace(cov_jk) > np.trace(cov), cov_jk, cov)
    # hypotheses: each Gaia star within hyp_radius
    hyp = []
    sc = scene[scene["sep"] <= hyp_radius]
    chis, amps = _chi2_at(sectors, sc["xi"].to_numpy(), sc["eta"].to_numpy(), per_sector=True)
    for n_, (i, row) in enumerate(sc.iterrows()):
        dchi = float(chis[n_] - chi_min) / red
        # implied eclipse depth of that star = flux deficit / its flux
        num = den = 0.0
        for sd, (A, eA) in zip(sectors, amps):
            F = sd["star_flux"][i] * sd["smear"]
            if F > 0 and np.isfinite(A[n_]):
                wgt = (F / eA[n_]) ** 2
                num += wgt * A[n_] / F
                den += wgt
        hyp.append({"source_id": int(row["source_id"]), "is_target": bool(row["is_target"]),
                    "sep_arcsec": float(row["sep"]), "pa_deg": float(row["pa_deg"]),
                    "tmag": float(row["tmag"]), "gmag": float(row["Gmag"]),
                    "dchi2": dchi, "implied_depth": float(num / den) if den > 0 else np.nan,
                    "implied_depth_err": float(1 / np.sqrt(den)) if den > 0 else np.nan})
    hyp = sorted(hyp, key=lambda h_: h_["dchi2"])
    off = float(np.hypot(*best))
    e_off = float(np.sqrt(max(np.array([best[0], best[1]]) @ cov_use @ np.array([best[0], best[1]]), 0))
                  / off) if off > 0 else float(np.sqrt(np.trace(cov_use) / 2))
    tgt = next(h_ for h_ in hyp if h_["is_target"])
    return {"xi": float(best[0]), "eta": float(best[1]), "offset_arcsec": off, "offset_err_arcsec": e_off,
            "cov": cov.tolist(), "cov_jackknife": np.asarray(cov_jk).tolist(), "red_chi2": float(red),
            "n_sectors": len(sectors), "target_dchi2": tgt["dchi2"], "hypotheses": hyp,
            "grid": {"xi": g.tolist(), "chi": ((chi.reshape(XI.shape) - chi_min) / red).tolist()},
            "jackknife": np.asarray(jk).tolist()}


# --------------------------------------------------------------------------
# Flux budget and depth vs aperture
# --------------------------------------------------------------------------
RINGS = ((0, 1), (1, 4), (4, 9), (9, 16), (16, 25))


def rank_pixels(img: np.ndarray, bad: np.ndarray) -> np.ndarray:
    """Pixel indices by decreasing value, off-detector pixels last."""
    return np.argsort(np.where(bad, -np.inf, img).ravel())[::-1]


def ring_test(sectors: list[dict], scene: pd.DataFrame, candidates: list[int]) -> dict:
    """Depth in rings of pixels ranked by the target's PSF, and the fit of
    "the dip is on star j" for each candidate source row j."""
    obs, err, pred = [], [], {j: [] for j in candidates}
    curves = []
    for sd in sectors:
        fit = sd["fit"]
        shape = sd["diff"].shape
        rel = scene["rel_flux"].to_numpy()
        allp = psf_images(sd["px"] + fit["dx"], sd["py"] + fit["dy"], shape, fit["psf"])
        star_img = fit["scale"] * allp * rel[:, None, None]          # per-star light per pixel
        tot = np.clip(sd["oot"] - fit["c0"], 1e-3, None)
        order = rank_pixels(star_img[sd["it"]], sd["bad"])
        null_sums = sd["nulls"].reshape(len(sd["nulls"]), -1) if len(sd["nulls"]) else None
        row = []
        for a, b in RINGS:
            pix = order[a:b]
            F = tot.ravel()[pix].sum()
            d = sd["diff"].ravel()[pix].sum() / F / sd["smear"]
            e = (np.std(null_sums[:, pix].sum(1)) if null_sums is not None and len(null_sums) >= 8
                 else np.sqrt((sd["err"].ravel()[pix] ** 2).sum())) / F / sd["smear"]
            obs.append(d)
            err.append(e)
            for j in candidates:
                pred[j].append(star_img[j].ravel()[pix].sum() / F)
            row.append((d, e, star_img[sd["it"]].ravel()[pix].sum() / F))
        curves.append({"sector": sd["sector"], "rings": row})
    obs, err = np.array(obs), np.array(err)
    out = {}
    for j in candidates:
        f = np.array(pred[j])
        w = 1 / err ** 2
        delta = float(np.sum(w * f * obs) / np.sum(w * f * f))
        chi = float(np.sum(w * (obs - delta * f) ** 2))
        out[int(j)] = {"intrinsic_depth": delta, "intrinsic_depth_err": float(1 / np.sqrt(np.sum(w * f * f))),
                       "chi2": chi, "n": len(obs)}
    # dilution-corrected depth per ring (as if on the target), all sectors combined
    corr = []
    for k in range(len(RINGS)):
        d = np.array([c["rings"][k][0] / c["rings"][k][2] for c in curves])
        e = np.array([c["rings"][k][1] / c["rings"][k][2] for c in curves])
        w = 1 / e ** 2
        corr.append((float(np.sum(w * d) / w.sum()), float(1 / np.sqrt(w.sum())),
                     float(np.median([c["rings"][k][2] for c in curves]))))
    return {"hypotheses": out, "ring_corrected_depth": corr, "rings": RINGS}


def flux_budget(sectors: list[dict], scene: pd.DataFrame, depth_on_target: float,
                max_sep: float = 90.0) -> list[dict]:
    """For each neighbour: the eclipse depth it would need to produce the dip,
    = observed depth x total flux in the aperture / neighbour flux in the
    aperture (top-4 target pixels; median over sectors)."""
    rows = []
    req = {}
    for sd in sectors:
        fit = sd["fit"]
        rel = scene["rel_flux"].to_numpy()
        allp = psf_images(sd["px"] + fit["dx"], sd["py"] + fit["dy"], sd["diff"].shape, fit["psf"])
        star_img = fit["scale"] * allp * rel[:, None, None]
        order = rank_pixels(star_img[sd["it"]], sd["bad"])[:4]
        tot = np.clip(sd["oot"] - fit["c0"], 1e-3, None).ravel()[order].sum()
        ft = star_img[sd["it"]].ravel()[order].sum()
        obs = depth_on_target * ft / tot        # what the aperture shows if the dip is on target
        for j in range(len(scene)):
            fj = star_img[j].ravel()[order].sum()
            req.setdefault(j, []).append(obs * tot / fj if fj > 0 else np.inf)
    for j, v in req.items():
        r = scene.iloc[j]
        if r["is_target"] or r["sep"] > max_sep:
            continue
        need = float(np.median(v))
        rows.append({"source_id": int(r["source_id"]), "sep_arcsec": float(r["sep"]),
                     "pa_deg": float(r["pa_deg"]), "gmag": float(r["Gmag"]), "tmag": float(r["tmag"]),
                     "bp_rp": float(r["bp_rp"]) if pd.notna(r["bp_rp"]) else None,
                     "plx": float(r["Plx"]) if pd.notna(r["Plx"]) else None,
                     "ruwe": float(r["RUWE"]) if pd.notna(r["RUWE"]) else None,
                     "required_depth": need, "can_mimic": bool(need < 1.0),
                     "all_flux_required_depth": float(depth_on_target / r["rel_flux"])})
    return sorted(rows, key=lambda x: x["required_depth"])


def per_sector_depths(sectors: list[dict]) -> list[dict]:
    """Dilution-corrected depth (as if on target) in the top-4 target pixels,
    per sector, corrected for exposure smearing. Same method for every sector."""
    out = []
    for sd in sectors:
        fit = sd["fit"]
        timg = fit["scale"] * psf_images(sd["px"][[sd["it"]]] + fit["dx"], sd["py"][[sd["it"]]] + fit["dy"],
                                         sd["diff"].shape, fit["psf"])[0]
        order = rank_pixels(timg, sd["bad"])[:4]
        ft = timg.ravel()[order].sum()
        tot = np.clip(sd["oot"] - fit["c0"], 1e-3, None).ravel()[order].sum()
        d = sd["diff"].ravel()[order].sum() / ft / sd["smear"]
        ns = sd["nulls"].reshape(len(sd["nulls"]), -1)[:, order].sum(1) if len(sd["nulls"]) >= 8 else None
        e = (np.std(ns) if ns is not None else np.sqrt((sd["err"].ravel()[order] ** 2).sum())) / ft / sd["smear"]
        out.append({"sector": sd["sector"], "camera": sd["camera"], "ccd": sd["ccd"], "year": sd["year"],
                    "cadence_min": sd["cad"] * 1440, "n_tr": sd["n_tr"], "depth_ppm": d * 1e6,
                    "err_ppm": e * 1e6, "target_frac_in_ap": ft / tot,
                    "psf_sigma_px": float(np.sqrt(fit["psf"][0] * fit["psf"][1])),
                    "scene_resid": fit["frac_resid"]})
    return out
