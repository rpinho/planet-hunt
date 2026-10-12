"""Difference-image centroid test from TESS full-frame-image cutouts.

The light-curve centroid test (multisector.centroid_multi) only says whether
the light centre moves. The standard test is a difference image: average the
pixels in transit, average the pixels just before and after each transit,
subtract. The pixels that lose light are where the eclipsing source is. If
the dip is on the target, the difference image peaks at the target's
position; if it peaks one or more pixels away, a neighbour is the source.

Cutouts come from MAST TESScut (11 x 11 pixels per sector, cached in
cache/tesscut). Offsets are converted to arcsec on the sky with each
sector's WCS and then averaged over sectors (TESS pixels are 21").
"""
from __future__ import annotations

import numpy as np

from . import config

CUT_DIR = config.CACHE / "tesscut"


def _cutouts(ra: float, dec: float, sectors: list[int], size: int = 11):
    from astroquery.mast import Tesscut
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    d = CUT_DIR / f"{ra:.5f}_{dec:+.5f}"
    d.mkdir(parents=True, exist_ok=True)
    c = SkyCoord(ra * u.deg, dec * u.deg)
    out = []
    for s in sectors:
        hits = sorted(d.glob(f"*s{s:04d}-*.fits"))
        if not hits:
            try:
                Tesscut.download_cutouts(coordinates=c, size=size, sector=s, path=str(d))
            except Exception as exc:
                print(f"[tesscut] S{s}: {exc!r}"[:120], flush=True)
                continue
            hits = sorted(d.glob(f"*s{s:04d}-*.fits"))
        if hits:
            out.append((s, hits[0]))
    return out


def difference_image(ra: float, dec: float, sectors: list[int], P: float, t0: float,
                     dur: float, size: int = 11) -> dict:
    """Per-sector and combined offset of the difference-image centroid from
    the target (arcsec), plus the out-of-transit centroid offset for scale."""
    from astropy.io import fits
    from astropy.wcs import WCS
    rows = []
    for s, path in _cutouts(ra, dec, sectors, size):
        with fits.open(path) as hd:
            d = hd[1].data
            t = np.asarray(d["TIME"], float)
            flux = np.asarray(d["FLUX"], float)
            q = np.asarray(d["QUALITY"])
            w = WCS(hd[2].header)
        ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(flux.reshape(len(t), -1)), axis=1)
        t, flux = t[ok], flux[ok]
        if len(t) < 100:
            continue
        ph = (t - t0 + 0.5 * P) % P - 0.5 * P
        epoch = np.round((t - t0) / P).astype(int)
        bkg = np.median(flux, axis=(1, 2))[:, None, None]
        f = flux - bkg
        # pair each transit with the frames just before and after it, so slow
        # scattered-light changes cancel
        diffs = []
        for e in np.unique(epoch[np.abs(ph) < 0.4 * dur]):
            m_in = (epoch == e) & (np.abs(ph) < 0.4 * dur)
            m_out = (epoch == e) & (np.abs(ph) > 0.75 * dur) & (np.abs(ph) < 2.5 * dur)
            if m_in.sum() >= 2 and m_out.sum() >= 4:
                diffs.append(np.median(f[m_out], axis=0) - np.median(f[m_in], axis=0))
        if len(diffs) < 2:
            continue
        intr = np.abs(ph) < 0.4 * dur
        diff = np.mean(diffs, axis=0)
        img_out = np.median(f[~intr], axis=0)
        yy, xx = np.mgrid[0:diff.shape[0], 0:diff.shape[1]]
        tx, ty = w.world_to_pixel_values(ra, dec)

        def centroid(img, r=1.0):
            m = (np.hypot(xx - tx, yy - ty) < r + 1.5) & np.isfinite(img)
            v = np.clip(img[m], 0, None)
            if v.sum() <= 0:
                return np.nan, np.nan
            return float(np.sum(v * xx[m]) / v.sum()), float(np.sum(v * yy[m]) / v.sum())
        dx, dy = centroid(diff)
        ox, oy = centroid(img_out)
        # noise: scatter of the difference image far from the star
        far = np.hypot(xx - tx, yy - ty) > 4
        snr = float(np.nanmax(diff[np.hypot(xx - tx, yy - ty) < 2]) / (1.4826 * np.nanmedian(np.abs(diff[far] - np.nanmedian(diff[far])))))
        sky = w.pixel_to_world_values([dx, ox, tx], [dy, oy, ty])
        ra_d, dec_d = sky[0], sky[1]

        def sep(i, j):
            return float(np.hypot((ra_d[i] - ra_d[j]) * np.cos(np.radians(dec)), dec_d[i] - dec_d[j]) * 3600)
        rows.append({"sector": int(s), "n_in": int(intr.sum()), "diff_snr": snr,
                     "diff_vs_target_arcsec": sep(0, 2), "oot_vs_target_arcsec": sep(1, 2),
                     "diff_vs_oot_arcsec": sep(0, 1),
                     "dra_arcsec": float((ra_d[0] - ra_d[1]) * np.cos(np.radians(dec)) * 3600),
                     "ddec_arcsec": float((dec_d[0] - dec_d[1]) * 3600),
                     "peak_pixel_offset": float(np.hypot(*np.subtract(np.unravel_index(np.nanargmax(diff), diff.shape)[::-1], (tx, ty))))})
    if not rows:
        return {"ok": False, "verdict": "no usable cutouts"}
    good = [r for r in rows if r["diff_snr"] > 3]
    use = good or rows
    w_ = np.array([max(r["diff_snr"], 0.1) ** 2 for r in use])
    mra = float(np.sum(w_ * [r["dra_arcsec"] for r in use]) / w_.sum())
    mdec = float(np.sum(w_ * [r["ddec_arcsec"] for r in use]) / w_.sum())
    scatter = float(np.std([np.hypot(r["dra_arcsec"], r["ddec_arcsec"]) for r in use]) / np.sqrt(max(len(use) - 1, 1)))
    off = float(np.hypot(mra, mdec))
    if off < max(7.0, 2 * scatter):
        v = (f"difference image centred on the star: offset {off:.1f}\" ± {scatter:.1f}\" from the "
             f"out-of-transit centroid ({len(use)} sectors with a visible dip)")
    elif off < 15:
        v = f"difference image offset {off:.1f}\" ± {scatter:.1f}\": inconclusive at TESS resolution"
    else:
        v = f"difference image offset {off:.0f}\" ± {scatter:.0f}\": the dip is probably on a neighbour"
    return {"ok": True, "n_sectors": len(rows), "n_used": len(use), "offset_arcsec": off,
            "offset_err_arcsec": scatter, "dra_arcsec": mra, "ddec_arcsec": mdec, "verdict": v,
            "per_sector": rows}
