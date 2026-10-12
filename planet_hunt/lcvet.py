"""Light-curve vetting on fixed-aperture TESScut photometry.

Pipeline light curves change product across years (TESS-SPOC where it exists,
QLP otherwise), and the two differ in aperture, dilution correction and
detrending, so a depth that changes between years can be a product change.
Here every sector gets the SAME treatment: the top-N pixels of the target's
PSF (from pixelvet's scene fit), the scene's own sky constant removed, a
running-median detrend with the transits masked, and the dilution corrected
with the scene model (flux of the target in the aperture / total).

On that light curve:
  * per-sector and per-season depths (same pixels-ranking rule every sector),
  * odd/even depths and a secondary eclipse at phase 0.5 and at any phase,
  * a limb-darkened transit model (numerical integration over the planet's
    disc, quadratic limb darkening, exposure-time integration per cadence),
    fitted by MCMC with the stellar density free; the posterior density is
    compared with the TIC's. A grazing eclipsing binary or a dip on a bigger
    background star needs a density (or impact parameter) the target cannot
    give.
"""
from __future__ import annotations

import numpy as np

from . import multisector as ms
from .pixelvet import psf_images

G_CGS = 6.674e-8
DAY_S = 86400.0


# --------------------------------------------------------------------------
# Fixed-aperture light curves from cutouts
# --------------------------------------------------------------------------
def cutout_lightcurve(cut: dict, sd: dict, npix: int = 4, row: int | None = None) -> dict:
    """Dilution-corrected relative flux of the target (or of scene star `row`)
    in its top-N PSF pixels."""
    fit = sd["fit"]
    it = sd["it"] if row is None else row
    timg = fit["scale"] * psf_images(sd["px"][[it]] + fit["dx"], sd["py"][[it]] + fit["dy"],
                                     sd["diff"].shape, fit["psf"])[0]
    from .pixelvet import rank_pixels
    order = rank_pixels(timg, sd["bad"])[:npix]
    f = cut["f"].reshape(len(cut["t"]), -1)[:, order].sum(1) - fit["c0"] * npix
    ft = timg.ravel()[order].sum()
    tot = np.clip(sd["oot"] - fit["c0"], 1e-3, None).ravel()[order].sum()
    frac = float(np.clip(ft / tot, 0.05, 1.0))
    rel = f / np.median(f) if np.median(f) > 0 else np.full(len(f), np.nan)
    return {"time": cut["t"], "flux": 1 + (rel - 1) / frac, "cx": np.zeros_like(rel), "cy": np.zeros_like(rel),
            "meta": {"sector": sd["sector"], "camera": sd["camera"], "ccd": sd["ccd"],
                     "provenance": "TESScut-fixed", "crowdsap": frac, "cadence_min": sd["cad"] * 1440}}


# --------------------------------------------------------------------------
# Transit model
# --------------------------------------------------------------------------
def ld_coeffs(teff: float) -> tuple[float, float]:
    """Rough TESS-band quadratic limb darkening (Claret 2017 trend, logg 4.5)."""
    tt = np.array([3000, 3400, 4000, 5000, 5800, 6500])
    u1 = np.array([0.15, 0.20, 0.35, 0.42, 0.38, 0.32])
    u2 = np.array([0.40, 0.38, 0.28, 0.22, 0.21, 0.22])
    return float(np.interp(teff, tt, u1)), float(np.interp(teff, tt, u2))


def _disc_points(n_r: int = 12, n_t: int = 24):
    r = np.sqrt((np.arange(n_r) + 0.5) / n_r)          # equal-area rings
    th = (np.arange(n_t) + 0.5) / n_t * 2 * np.pi
    R, T = np.meshgrid(r, th)
    return (R * np.cos(T)).ravel(), (R * np.sin(T)).ravel()


_DX, _DY = _disc_points()


def occultation(z: np.ndarray, k: float, u1: float, u2: float) -> np.ndarray:
    """Fractional flux loss for sky separation z (stellar radii), radius ratio k."""
    z = np.asarray(z, float)
    out = np.zeros_like(z)
    m = z < 1 + k
    if not m.any():
        return out
    x = z[m, None] + k * _DX[None, :]
    y = k * _DY[None, :]
    r2 = x ** 2 + y ** 2
    mu = np.sqrt(np.clip(1 - r2, 0, None))
    inten = np.where(r2 < 1, 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2, 0.0)
    norm = 1 - u1 / 3 - u2 / 6
    out[m] = k ** 2 * inten.mean(1) / norm
    return out


def transit_model(ph_days: np.ndarray, cad_days: np.ndarray, P: float, k: float, b: float,
                  rho: float, u1: float, u2: float, nsub: int = 7) -> np.ndarray:
    """Exposure-integrated model at phase offsets (days from mid-transit)."""
    a_rs = (G_CGS * rho * (P * DAY_S) ** 2 / (3 * np.pi)) ** (1 / 3)
    off = ((np.arange(nsub) + 0.5) / nsub - 0.5)
    tt = ph_days[:, None] + cad_days[:, None] * off[None, :]
    phi = 2 * np.pi * tt / P
    xs = a_rs * np.sin(phi)
    ys = b * np.cos(phi)
    z = np.sqrt(xs ** 2 + ys ** 2)
    z = np.where(np.cos(phi) > 0, z, 10.0)
    return 1 - occultation(z.ravel(), k, u1, u2).reshape(z.shape).mean(1)


def _binned(lc, P, t0, dur, bin_min=3.0, half_window=3.0):
    """Phase-binned data near transit, separately per cadence group."""
    ph = ms.phase(lc["t"], P, t0)
    near = np.abs(ph) < half_window * dur
    out = []
    for c in np.unique(np.round(lc["cad"][near] * 1440, 1)):
        m = near & (np.abs(lc["cad"] * 1440 - c) < 0.2)
        if m.sum() < 20:
            continue
        w = 1 / lc["dy"][m] ** 2
        edges = np.arange(-half_window * dur, half_window * dur + 1e-9, bin_min / 1440)
        idx = np.digitize(ph[m], edges)
        ws = np.bincount(idx, w, len(edges) + 1)
        ok = ws > 0
        xb = np.bincount(idx, w * ph[m], len(edges) + 1)[ok] / ws[ok]
        yb = np.bincount(idx, w * lc["flat"][m], len(edges) + 1)[ok] / ws[ok]
        out.append((xb, yb, 1 / np.sqrt(ws[ok]), np.full(ok.sum(), c / 1440)))
    x, y, e, cad = (np.concatenate([o[i] for o in out]) for i in range(4))
    return x, y, e, cad


def fit_transit(lc: dict, P: float, t0: float, dur: float, depth: float, teff: float,
                rho_tic: float, n_steps: int = 12000, seed: int = 1, profile: bool = False) -> dict:
    """MCMC (adaptive Metropolis) over (dt0, log k, b, log rho), rho free.
    Returns posterior summaries and the tension with the TIC density."""
    u1, u2 = ld_coeffs(teff)
    x, y, e, cad = _binned(lc, P, t0, dur)
    # red-noise-ish inflation: scale errors so the out-of-transit chi2/dof = 1
    oot = np.abs(x) > 0.75 * dur
    s = np.sqrt(max(1.0, np.mean(((y[oot] - np.average(y[oot], weights=1 / e[oot] ** 2)) / e[oot]) ** 2)))
    e = e * s

    def loglike(q):
        dt, lk, b, lr, f0 = q
        k = np.exp(lk)
        if not (0 <= b < 1 + k) or not (-4 < lr < 4.5) or k > 0.6 or abs(dt) > 0.5 * dur:
            return -np.inf
        m = f0 * transit_model(x - dt, cad, P, k, b, np.exp(lr), u1, u2)
        return -0.5 * np.sum(((y - m) / e) ** 2)
    rng = np.random.default_rng(seed)
    q = np.array([0.0, np.log(np.sqrt(max(depth, 1e-5))), 0.3, np.log(rho_tic), 1.0])
    # quick optimisation from a few starts
    from scipy import optimize
    best = None
    for b0 in (0.1, 0.5, 0.8):
        for lr0 in (np.log(rho_tic), np.log(rho_tic) - 1.5, np.log(rho_tic) + 1.0):
            q0 = q.copy(); q0[2] = b0; q0[3] = lr0
            r = optimize.minimize(lambda qq: -loglike(qq), q0, method="Nelder-Mead",
                                  options={"maxiter": 3000, "xatol": 1e-5, "fatol": 1e-3})
            if np.isfinite(r.fun) and (best is None or r.fun < best.fun):
                best = r
    q = best.x
    lp = loglike(q)
    step = np.array([0.02 * dur, 0.05, 0.05, 0.15, 2e-5])
    chain = []
    acc = 0
    for i in range(n_steps):
        prop = q + step * rng.standard_normal(5)
        lpp = loglike(prop)
        if np.log(rng.random()) < lpp - lp:
            q, lp = prop, lpp
            acc += 1
        chain.append(q.copy())
        if i < n_steps // 3 and i % 500 == 499:      # adapt during burn-in
            ar = acc / 500
            step *= np.exp(ar - 0.3)
            acc = 0
    ch = np.array(chain[n_steps // 3:])
    k = np.exp(ch[:, 1]); b = ch[:, 2]; rho = np.exp(ch[:, 3])
    pct = lambda v: [float(np.percentile(v, p)) for p in (16, 50, 84)]
    # tension: where does the TIC density fall in the posterior of log rho?
    lr = ch[:, 3]
    frac_below = float(np.mean(lr < np.log(rho_tic)))
    from scipy import stats
    z = float(stats.norm.ppf(np.clip(frac_below, 1e-4, 1 - 1e-4)))
    # likelihood with rho fixed at TIC vs free (best fits)
    r_fix = optimize.minimize(lambda qq: -loglike(np.r_[qq[:3], np.log(rho_tic), qq[3]]),
                              np.r_[best.x[:3], best.x[4]], method="Nelder-Mead",
                              options={"maxiter": 3000})
    grazing = float(np.mean(b > 1 - k))
    # profile likelihood in stellar density: robust where the posterior is broad
    prof = []
    if profile:
        for lr0 in np.linspace(np.log(rho_tic) - 3.0, np.log(rho_tic) + 1.5, 13):
            bestp = np.inf
            for b0 in (0.2, 0.6, 0.9):
                q0 = np.r_[best.x[0], best.x[1], b0, best.x[4]]
                rr = optimize.minimize(lambda qq: -loglike(np.r_[qq[:3], lr0, qq[3]]), q0,
                                       method="Nelder-Mead", options={"maxiter": 600})
                bestp = min(bestp, rr.fun)
            prof.append((float(np.exp(lr0)), float(2 * (bestp - best.fun))))
    return {"k": pct(k), "b": pct(b), "rho_cgs": pct(rho), "rho_tic_cgs": float(rho_tic),
            "rho_tic_z": z, "p_grazing": grazing, "depth_model_ppm": float(np.median(k) ** 2 * 1e6),
            "chi2_best": float(2 * best.fun), "chi2_rho_tic": float(2 * r_fix.fun),
            "dchi2_rho_tic": float(2 * (r_fix.fun - best.fun)), "n_bins": int(len(x)),
            "rho_profile": prof,
            "err_scale": float(s), "u": [u1, u2],
            "best": {"dt0": float(best.x[0]), "k": float(np.exp(best.x[1])), "b": float(best.x[2]),
                     "rho": float(np.exp(best.x[3]))},
            "binned": {"x": x.tolist(), "y": y.tolist(), "e": e.tolist(), "cad": cad.tolist()}}


# --------------------------------------------------------------------------
# Simple box statistics on any stitched light curve
# --------------------------------------------------------------------------
def box_stats(lc: dict, P: float, t0: float, dur: float) -> dict:
    """Odd/even, secondary (0.5 and any phase), per-sector and per-season depths,
    with exposure-smearing corrections per cadence."""
    from .pixelvet import smear_factor
    t, y, dy, cad = lc["t"], lc["flat"], lc["dy"], lc["cad"]
    sm = np.array([smear_factor(dur, c) for c in cad]) if len(np.unique(cad)) < 2 else None
    if sm is None:
        uc = np.unique(np.round(cad, 6))
        tab = {c: smear_factor(dur, c) for c in uc}
        sm = np.array([tab[c] for c in np.round(cad, 6)])
    yc = 1 - (1 - y) / sm          # in-transit points rescaled; OOT points ~1 so barely change
    w = 1 / dy ** 2
    intr = ms.in_tr(t, P, t0, dur)
    oot = ~ms.in_tr(t, P, t0, dur, 2.5) & ~ms.in_tr(t, P, t0 + 0.5 * P, dur, 2.5)
    lev = np.average(y[oot], weights=w[oot])

    def dep(sel, yy=yc):
        if sel.sum() < 3:
            return np.nan, np.nan
        return float(lev - np.average(yy[sel], weights=w[sel])), float(1 / np.sqrt(w[sel].sum()))
    ep = np.round((t - t0) / P).astype(int)
    d, e = dep(intr)
    d_odd, e_odd = dep(intr & (ep % 2 == 1))
    d_even, e_even = dep(intr & (ep % 2 == 0))
    d_sec, e_sec = dep(ms.in_tr(t, P, t0 + 0.5 * P, dur), y)
    best = (-np.inf, np.nan, np.nan)
    ph01 = ((t - t0) / P) % 1
    q = dur / P
    for c in np.arange(0.1, 0.9, max(q / 3, 0.001)):
        m = np.abs(ph01 - c) < 0.5 * q
        dd, ee = dep(m, y)
        if np.isfinite(dd) and dd / ee > best[0]:
            best = (dd / ee, c, dd)
    secs = []
    for s in np.unique(lc["sector"]):
        m = lc["sector"] == s
        dd, ee = dep(intr & m)
        if np.isfinite(dd):
            secs.append({"sector": int(s), "depth_ppm": dd * 1e6, "err_ppm": ee * 1e6,
                         "t_mid": float(np.median(t[m]))})
    seas = []
    for i, idx in enumerate(ms.seasons(t)):
        m = np.zeros(len(t), bool); m[idx] = True
        dd, ee = dep(intr & m)
        if np.isfinite(dd):
            seas.append({"season": i + 1, "t_mid": float(np.median(t[idx])), "depth_ppm": dd * 1e6,
                         "err_ppm": ee * 1e6, "sectors": sorted(set(int(x) for x in lc["sector"][idx]))})
    return {"depth_ppm": d * 1e6, "err_ppm": e * 1e6, "snr": d / e if e else np.nan,
            "odd_ppm": d_odd * 1e6, "even_ppm": d_even * 1e6,
            "oe_sigma": float(abs(d_odd - d_even) / np.hypot(e_odd, e_even)),
            "sec_ppm": d_sec * 1e6, "sec_sigma": float(d_sec / e_sec) if e_sec else np.nan,
            "sec_any_sigma": float(best[0]), "sec_any_phase": float(best[1]),
            "sec_any_ppm": float(best[2] * 1e6), "per_sector": secs, "per_season": seas}


def chi2_const(rows, key="depth_ppm", ekey="err_ppm") -> dict:
    from scipy import stats
    d = np.array([r[key] for r in rows]); e = np.array([r[ekey] for r in rows])
    ok = np.isfinite(d) & np.isfinite(e) & (e > 0)
    d, e = d[ok], e[ok]
    if len(d) < 2:
        return {"n": int(len(d))}
    w = 1 / e ** 2
    mu = float(np.sum(w * d) / w.sum())
    c2 = float(np.sum(((d - mu) / e) ** 2))
    return {"n": int(len(d)), "mean": mu, "err": float(1 / np.sqrt(w.sum())), "chi2": c2,
            "p": float(stats.chi2.sf(c2, len(d) - 1))}
