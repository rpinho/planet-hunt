"""Is the dip on the target star? Two cheap checks that need no pixel files.

1. Centroid shift (`centroid_test`). SPOC light curves carry the
   flux-weighted centroid of the aperture (MOM_CENTR1/2, in pixels). If the
   dip comes from a neighbour inside or near the aperture, the centroid moves
   during transit, towards the target and away from the neighbour. For a
   source offset D from the light centre and a dip of fractional depth d, the
   shift is about d * D, so D ~ shift / depth. A clear shift with D >~ 1 pixel
   (21") means "probably a neighbour". No significant shift with a small error
   on D is evidence the dip is on (or very near) the target. Shallow dips give
   large errors on D: then the test is simply inconclusive.

2. Neighbours bright enough to fake it (`neighbour_check`). A neighbour can
   only produce an observed depth d if, fully eclipsed, its flux share is at
   least d: F_nb / F_total >= d, i.e. delta-mag <= -2.5 log10(d). This lists
   the TIC stars within ~2.5 TESS pixels that pass that bound. Zero such stars
   is strong evidence for the target; many means the centroid test (or
   higher-resolution data) has to decide.

Neither check replaces a difference-image centroid analysis or TRICERATOPS.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config
from .data import in_transit

import lightkurve as lk  # noqa: E402

PIXEL_ARCSEC = 21.0


def _detrend_series(t, y, mask, window_days=1.0):
    """Subtract a 1-day running median computed on out-of-transit points."""
    trend = np.full_like(y, np.nan)
    ok = ~mask & np.isfinite(y)
    half = window_days / 2
    to, yo = t[ok], y[ok]
    lo = np.searchsorted(to, t - half)
    hi = np.searchsorted(to, t + half)
    for i in range(len(t)):
        if hi[i] - lo[i] >= 10:
            trend[i] = np.median(yo[lo[i]:hi[i]])
    return y - trend


def centroid_test(path, period: float, t0: float, duration_h: float,
                  depth_ppm: float, gap_days: float = 0.5, edge_days: float = 0.25) -> dict:
    """In- vs out-of-transit flux-weighted centroid for one sector's LC file."""
    lc = lk.read(str(path), quality_bitmask="default")
    t = lc.time.value
    cx = np.asarray(lc["mom_centr1"].value, float)
    cy = np.asarray(lc["mom_centr2"].value, float)
    f = np.asarray(lc["pdcsap_flux"].value, float)
    good = np.isfinite(t) & np.isfinite(cx) & np.isfinite(cy) & np.isfinite(f)
    t, cx, cy = t[good], cx[good], cy[good]
    starts = np.r_[t[0], t[1:][np.diff(t) > gap_days]]
    keep = np.ones(len(t), bool)
    for s in starts:
        keep &= ~((t >= s) & (t < s + edge_days))
    t, cx, cy = t[keep], cx[keep], cy[keep]
    dur = duration_h / 24
    intr = in_transit(t, period, t0, dur)
    near = in_transit(t, period, t0, dur, widen=2.0)
    if intr.sum() < 10:
        return {"ok": False, "why": "too few in-transit cadences"}
    out = {}
    for name, c in (("x", cx), ("y", cy)):
        r = _detrend_series(t, c, near)
        oot = ~near & np.isfinite(r)
        rin = r[intr & np.isfinite(r)]
        sig = 1.4826 * np.median(np.abs(r[oot] - np.median(r[oot])))
        # red-noise inflation: scatter of transit-duration averages
        cad = np.median(np.diff(t))
        nb = max(1, int(round(dur / cad)))
        ro = r[oot] - np.median(r[oot])
        binned = np.convolve(ro, np.ones(nb) / nb, mode="valid")[::nb]
        beta = max(1.0, np.std(binned) / (sig / np.sqrt(nb))) if len(binned) > 5 else 1.0
        out[f"d{name}"] = float(np.mean(rin) - np.median(r[oot]))
        out[f"e{name}"] = float(beta * sig / np.sqrt(len(rin)))
    shift = float(np.hypot(out["dx"], out["dy"]))
    sigma = float(np.hypot(out["dx"] / out["ex"], out["dy"] / out["ey"]))
    err = float(np.hypot(out["ex"], out["ey"]) / np.sqrt(2))
    depth = max(depth_ppm, 1) * 1e-6
    return {"ok": True, "shift_px": shift, "shift_err_px": err, "shift_sigma": sigma,
            "shift_arcsec": shift * PIXEL_ARCSEC,
            # where the source would have to be, relative to the light centroid
            "offset_arcsec": shift / depth * PIXEL_ARCSEC,
            "offset_err_arcsec": err / depth * PIXEL_ARCSEC,
            "n_in": int(intr.sum()), **out}


def centroid_verdict(c: dict) -> str:
    """Plain-language read. Calibrated on sector 90: known planets show no
    offset; TOI-409.01 (a TFOPWG false positive, "centroid offset on a
    neighbour") shows a 2.4-sigma shift pointing ~50" away -> MARGINAL."""
    if not c.get("ok"):
        return "not measured"
    off, err, sig = c["offset_arcsec"], c["offset_err_arcsec"], c["shift_sigma"]
    if sig >= 3 and off > PIXEL_ARCSEC:
        return f"SHIFT: centroid moves {sig:.1f} sigma in transit (source ~{off:.0f}\" away): likely a neighbour"
    if sig >= 2 and off > PIXEL_ARCSEC:
        return (f"MARGINAL shift ({sig:.1f} sigma) pointing ~{off:.0f}\" away: "
                "possibly a neighbour, needs a difference image")
    if sig >= 2:   # significant but small: the light centre moves towards the target
        return (f"small shift ({sig:.1f} sigma) pointing ~{off:.0f}\" from the light centre: "
                f"consistent with the target, localised to ~{off + 2 * err:.0f}\"")
    if off + 2 * err < 1.5 * PIXEL_ARCSEC:
        return f"no shift ({sig:.1f} sigma): source within ~{off + 2 * err:.0f}\" of the light centre"
    return (f"inconclusive: no significant shift ({sig:.1f} sigma), but the dip is too "
            f"shallow to localise better than ~{off + 2 * err:.0f}\"")


def neighbour_check(tic: int, ra: float, dec: float, tmag: float, depth_ppm: float,
                    radius_arcsec: float = 2.5 * PIXEL_ARCSEC) -> dict:
    """TIC stars within `radius_arcsec` bright enough to produce the depth."""
    cache = config.CATALOG_DIR / "neighbours"
    cache.mkdir(exist_ok=True)
    path = cache / f"{tic}.csv"
    if path.exists():
        nb = pd.read_csv(path)
    else:
        from astroquery.mast import Catalogs, conf
        conf.timeout = 120
        tab = Catalogs.query_region(f"{ra} {dec}", radius=radius_arcsec / 3600, catalog="TIC")
        nb = tab.to_pandas()[["ID", "Tmag", "ra", "dec", "dstArcSec"]]
        nb.to_csv(path, index=False)
    nb = nb[nb["ID"].astype("int64") != int(tic)].copy()
    nb = nb[nb["dstArcSec"] <= radius_arcsec]
    f_t = 10 ** (-0.4 * tmag)
    f_nb = 10 ** (-0.4 * nb["Tmag"].astype(float))
    total = f_t + f_nb.sum()
    nb["max_depth_ppm"] = f_nb / total * 1e6      # if 100% eclipsed
    able = nb[nb["max_depth_ppm"] >= depth_ppm].sort_values("Tmag")
    dm_max = -2.5 * np.log10(max(depth_ppm, 1) * 1e-6)
    return {"n_within": int(len(nb)), "n_able": int(len(able)),
            "dmag_limit": float(dm_max), "radius_arcsec": radius_arcsec,
            "dilution": float(f_t / total),
            "able": [{"tic": int(r.ID), "Tmag": float(r.Tmag), "sep_arcsec": float(r.dstArcSec),
                      "max_depth_ppm": float(r.max_depth_ppm)} for r in able.head(8).itertuples()]}
