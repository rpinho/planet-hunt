"""Multi-sector, multi-year transit search for one star (method A).

A red dwarf near an ecliptic pole can have 30-70 TESS sectors spread over
eight years. Folding all of them coherently needs a period grid fine enough
that a transit does not drift by more than a fraction of its duration over
the whole baseline: for 8 years that is ~10^6 trial periods per star, too
many for BLS in Python. So the search is done in two stages:

  1. semi-coherent: split the data into observing seasons (sectors less than
     60 days apart, i.e. up to ~1 year), run BLS on each season on a shared
     period grid, and SUM the BLS log-likelihood powers across seasons. A
     real planet adds power at the same period in every season (with its own
     best phase per season); noise peaks do not line up.
  2. coherent: for the top peaks, refine on ALL the data with a fine local
     grid. A real planet keeps one ephemeris for eight years.

The detection statistics are then measured on the un-binned data with the
transits masked from the detrending (so the filter does not eat the depth),
with per-sector noise weights (30-min, 10-min and 200-s cadences mix).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares

G, MSUN, RSUN, REARTH, AU, DAY = 6.674e-11, 1.989e30, 6.957e8, 6.371e6, 1.496e11, 86400.0
TESS_ORBIT = 13.7


# --------------------------------------------------------------------------
# Cleaning
# --------------------------------------------------------------------------
def running_median(t, y, window: float, mask=None):
    """Centred running median over a time window (days), ignoring masked points."""
    s = pd.Series(np.where(mask, np.nan, y) if mask is not None else y,
                  index=pd.to_datetime(t, unit="D"))
    return s.rolling(f"{int(window * 86400)}s", center=True, min_periods=5).median().to_numpy()


def clean_sector(d: dict, window: float = 0.75, gap: float = 0.5, edge: float = 0.25,
                 transit_mask=None) -> dict | None:
    """Detrend one sector. Returns t, flat, dy (robust per-sector scatter),
    cadence (days) and the trend. `transit_mask(t) -> bool array` excludes
    known transits from the trend fit."""
    t = np.asarray(d["time"], float)
    f = np.asarray(d["flux"], float)
    ok = np.isfinite(t) & np.isfinite(f)
    t, f = t[ok], f[ok]
    cx, cy = np.asarray(d["cx"], float)[ok], np.asarray(d["cy"], float)[ok]
    if len(t) < 100:
        return None
    starts = np.r_[t[0], t[1:][np.diff(t) > gap]]
    keep = np.ones(len(t), bool)
    for s in starts:
        keep &= ~((t >= s) & (t < s + edge))
    t, f, cx, cy = t[keep], f[keep], cx[keep], cy[keep]
    m = transit_mask(t) if transit_mask is not None else np.zeros(len(t), bool)
    trend = running_median(t, f, window, m)
    r = f / trend - 1
    sig = 1.4826 * np.nanmedian(np.abs(r[~m] - np.nanmedian(r[~m])))
    # second pass: also ignore strong outliers (flares, dips) when fitting the trend
    m2 = m | (np.abs(r) > 3 * sig)
    trend = running_median(t, f, window, m2)
    flat = f / trend
    resid = flat - 1
    sig = 1.4826 * np.nanmedian(np.abs(resid[~m2]))
    # clip flares: points > 4 sigma up, plus one cadence either side
    up = resid > 4 * sig
    up = up | np.r_[False, up[:-1]] | np.r_[up[1:], False]
    good = np.isfinite(flat) & ~up
    cad = float(np.median(np.diff(t)))
    return {"t": t[good], "flat": flat[good], "trend": trend[good], "flux": f[good],
            "cx": cx[good], "cy": cy[good], "sigma": float(sig), "cadence": cad,
            "sector": int(d["meta"].get("sector", 0)), "provenance": d["meta"].get("provenance", ""),
            "crowdsap": d["meta"].get("crowdsap", np.nan)}


def stitch(sectors: list[dict], transit_mask=None, window: float = 0.75) -> dict:
    """Clean every sector and concatenate (per-point dy, cadence, sector id)."""
    parts = [p for p in (clean_sector(d, window=window, transit_mask=transit_mask) for d in sectors)
             if p is not None and len(p["t"]) > 50 and np.isfinite(p["sigma"]) and p["sigma"] > 0]
    if not parts:
        return {}
    cat = lambda k: np.concatenate([p[k] for p in parts])
    out = {k: cat(k) for k in ("t", "flat", "trend", "flux", "cx", "cy")}
    out["dy"] = np.concatenate([np.full(len(p["t"]), p["sigma"]) for p in parts])
    out["cad"] = np.concatenate([np.full(len(p["t"]), p["cadence"]) for p in parts])
    out["sector"] = np.concatenate([np.full(len(p["t"]), p["sector"]) for p in parts])
    o = np.argsort(out["t"])
    out = {k: v[o] for k, v in out.items()}
    out["sector_info"] = [{k: p[k] for k in ("sector", "sigma", "cadence", "provenance", "crowdsap")}
                          | {"n": len(p["t"])} for p in parts]
    return out


def seasons(t: np.ndarray, gap_days: float = 60.0) -> list[np.ndarray]:
    """Index arrays of observing seasons (no internal gap > gap_days)."""
    br = np.where(np.diff(t) > gap_days)[0]
    edges = np.r_[0, br + 1, len(t)]
    return [np.arange(a, b) for a, b in zip(edges[:-1], edges[1:]) if b - a > 200]


def bin_weighted(t, y, dy, minutes: float = 20.0):
    w = 1 / dy ** 2
    idx = np.floor((t - t[0]) / (minutes / 1440)).astype(int)
    ws = np.bincount(idx, w)
    good = ws > 0
    tb = np.bincount(idx, w * t)[good] / ws[good]
    yb = np.bincount(idx, w * y)[good] / ws[good]
    return tb, yb, 1 / np.sqrt(ws[good])


# --------------------------------------------------------------------------
# Star-aware durations and period grid
# --------------------------------------------------------------------------
def central_duration(P, rstar: float, mstar: float):
    """Duration (days) of a central transit of a small planet."""
    a = (G * mstar * MSUN * (np.asarray(P) * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    return np.asarray(P) / np.pi * np.arcsin(np.minimum(1.0, rstar * RSUN / a))


def period_grid(pmin, pmax, T, rstar, mstar, oversample: float = 2.0, qfrac: float = 0.5):
    """Frequency grid whose step keeps the phase drift over the baseline T
    below 1/oversample of (qfrac x) the central-transit duration."""
    f, fmin = 1 / pmin, 1 / pmax
    out = []
    while f > fmin:
        P = 1 / f
        dur = qfrac * float(central_duration(P, rstar, mstar))
        out.append(P)
        f -= dur / (P * T * oversample)
    return np.sort(np.array(out))


def durations_for(P, rstar, mstar, cad_min: float = 0.0139):
    d0 = float(central_duration(P, rstar, mstar))
    d = np.unique(np.clip(d0 * np.array([0.35, 0.55, 0.8, 1.1, 1.5]), cad_min, 0.4))
    return d


def _bls_chunks(t, y, dy, periods, rstar, mstar, n_chunks: int = 10):
    """BLS likelihood power at `periods`, with star-aware duration lists per chunk."""
    bls = BoxLeastSquares(t, y, dy)
    edges = np.geomspace(periods.min(), periods.max() * 1.0001, n_chunks + 1)
    P_, POW, T0, DUR, DEP = [], [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        ch = periods[(periods >= lo) & (periods < hi)]
        if len(ch) == 0:
            continue
        durs = np.unique(np.r_[durations_for(lo, rstar, mstar), durations_for(hi, rstar, mstar)])
        durs = durs[durs < 0.5 * lo]
        r = bls.power(ch, durs, objective="likelihood")
        pw = np.where(np.asarray(r.depth) > 0, np.asarray(r.power), 0.0)
        P_.append(ch); POW.append(pw); T0.append(np.asarray(r.transit_time))
        DUR.append(np.asarray(r.duration)); DEP.append(np.asarray(r.depth))
    P_, POW, T0, DUR, DEP = map(np.concatenate, (P_, POW, T0, DUR, DEP))
    o = np.argsort(P_)
    return P_[o], POW[o], T0[o], DUR[o], DEP[o]


def sde(power):
    from scipy.ndimage import median_filter
    d = power - median_filter(power, size=501, mode="nearest")
    return (d - np.nanmean(d)) / np.nanstd(d)


def search(lc: dict, rstar: float, mstar: float, pmin: float = 0.5, pmax: float = 30.0,
           n_peaks: int = 5, bin_minutes: float = 20.0) -> dict:
    """Semi-coherent season-summed BLS + coherent refinement of the top peaks."""
    t, y, dy = lc["t"], lc["flat"], lc["dy"]
    seas = seasons(t)
    T_seas = max(t[s].max() - t[s].min() for s in seas)
    pmax = min(pmax, 0.5 * T_seas)
    periods = period_grid(pmin, pmax, T_seas, rstar, mstar)
    total = np.zeros(len(periods))
    for s in seas:
        tb, yb, eb = bin_weighted(t[s], y[s], dy[s], bin_minutes)
        sel = periods < 0.5 * (tb.max() - tb.min())
        if sel.sum() < 10:
            continue
        _, pw, *_ = _bls_chunks(tb, yb, eb, periods[sel], rstar, mstar)
        total[sel] += pw
    s_de = sde(total)
    # top distinct peaks
    # FFI light curves carry scattered-light systematics on TESS's 13.7-day
    # orbit; peaks within 2% of 6.85, 13.7 or 27.4 d are skipped (real planets
    # there are lost, a known blind spot of FFI searches)
    veto = lambda P: any(abs(P / (TESS_ORBIT * k) - 1) < 0.02 for k in (0.5, 1, 2))
    peaks = []
    for i in np.argsort(total)[::-1]:
        if veto(periods[i]):
            continue
        if all(abs(periods[i] / p - 1) > 0.002 for p in peaks):
            peaks.append(periods[i])
        if len(peaks) >= n_peaks:
            break
    # coherent refinement on all the data
    tb, yb, eb = bin_weighted(t, y, dy, bin_minutes)
    T_all = tb.max() - tb.min()
    best = None
    for P0 in peaks:
        dur = float(central_duration(P0, rstar, mstar))
        df_s = 0.5 * dur / (P0 * T_seas)          # semi-coherent peak half-width (freq)
        df_c = 0.5 * dur / (P0 * T_all * 3)       # coherent step
        fgrid = np.arange(1 / P0 - 3 * df_s, 1 / P0 + 3 * df_s, df_c)
        pg = np.sort(1 / fgrid)
        P_, pw, T0, DUR, DEP = _bls_chunks(tb, yb, eb, pg, rstar, mstar, n_chunks=1)
        j = int(np.argmax(pw))
        cand = {"period": float(P_[j]), "t0": float(T0[j]), "duration": float(DUR[j]),
                "depth_bls": float(DEP[j]), "power_coh": float(pw[j]),
                "power_semi": float(total[np.argmin(np.abs(periods - P0))]),
                "sde": float(s_de[np.argmin(np.abs(periods - P0))])}
        if best is None or cand["power_coh"] > best["power_coh"]:
            best = cand
    best.update(periods=periods, power=total, sde_curve=s_de, n_seasons=len(seas),
                T_season=float(T_seas), T_all=float(T_all), n_periods=len(periods))
    return best


# --------------------------------------------------------------------------
# Statistics on the un-binned, multi-cadence light curve
# --------------------------------------------------------------------------
def phase(t, P, t0):
    return (t - t0 + 0.5 * P) % P - 0.5 * P


def in_tr(t, P, t0, dur, widen=1.0):
    return np.abs(phase(t, P, t0)) < 0.5 * dur * widen


def _wdepth(y, w, sel, level):
    if sel.sum() == 0:
        return np.nan, np.nan
    ws = w[sel].sum()
    return float(level - np.sum(w[sel] * y[sel]) / ws), float(1 / np.sqrt(ws))


def red_noise_beta(lc, dur, oot):
    """Per-sector ratio of the measured scatter of transit-length averages to
    the white-noise expectation (>= 1)."""
    beta = np.ones(len(lc["t"]))
    for s in np.unique(lc["sector"]):
        m = (lc["sector"] == s) & oot
        if m.sum() < 50:
            continue
        r = lc["flat"][m] - np.median(lc["flat"][m])
        nb = max(1, int(round(dur / np.median(lc["cad"][m]))))
        if len(r) < 5 * nb:
            continue
        binned = np.convolve(r, np.ones(nb) / nb, mode="valid")[::nb]
        exp = 1.4826 * np.median(np.abs(r)) / np.sqrt(nb)
        beta[lc["sector"] == s] = max(1.0, float(np.std(binned) / exp))
    return beta


def transit_stats(lc: dict, P: float, t0: float, dur: float, rstar: float, mstar: float,
                  teff: float) -> dict:
    t, y = lc["t"], lc["flat"]
    intr = in_tr(t, P, t0, dur)
    oot = ~in_tr(t, P, t0, dur, 2.0) & ~in_tr(t, P, t0 + 0.5 * P, dur, 2.0)
    beta = red_noise_beta(lc, dur, oot)
    w_white = 1 / lc["dy"] ** 2
    w = 1 / (lc["dy"] * beta) ** 2
    level = float(np.sum(w_white[oot] * y[oot]) / w_white[oot].sum())
    depth, err_white = _wdepth(y, w_white, intr, level)
    _, err = _wdepth(y, w, intr, level)
    snr_white = depth / err_white if err_white > 0 else 0.0
    snr = depth / err if err > 0 else 0.0

    epoch = np.round((t - t0) / P).astype(int)
    cov = {}
    for e in np.unique(epoch[intr]):
        m = intr & (epoch == e)
        cov[int(e)] = float(np.sum(lc["cad"][m]))
    epochs = [e for e, c in cov.items() if c >= 0.5 * dur]
    n_tr = len(epochs)

    odd, even = intr & (epoch % 2 == 1), intr & (epoch % 2 == 0)
    d_odd, e_odd = _wdepth(y, w, odd, level)
    d_even, e_even = _wdepth(y, w, even, level)
    oe = abs(d_odd - d_even) / np.hypot(e_odd, e_even) if np.isfinite(d_odd + d_even) else np.nan
    sec = in_tr(t, P, t0 + 0.5 * P, dur)
    d_sec, e_sec = _wdepth(y, w, sec, level)
    sec_sigma = d_sec / e_sec if np.isfinite(d_sec) and e_sec > 0 else np.nan
    # secondary eclipse at ANY phase (eccentric binaries put it away from 0.5)
    sec_any, sec_any_phase = np.nan, np.nan
    ph01 = ((t - t0) / P) % 1.0
    q = dur / P
    best = -np.inf
    for c in np.arange(0.15, 0.85, max(q / 3, 0.002)):
        m = np.abs(ph01 - c) < 0.5 * q
        if m.sum() < 5:
            continue
        d_, e_ = _wdepth(y, w, m, level)
        if e_ > 0 and d_ / e_ > best:
            best, sec_any_phase = d_ / e_, float(c)
            sec_any = d_
    sec_any_sigma = float(best) if np.isfinite(best) else np.nan
    ph = phase(t, P, t0)
    inner = np.abs(ph) < 0.25 * dur
    d_in, _ = _wdepth(y, w, inner, level)
    d_outer, _ = _wdepth(y, w, intr & ~inner, level)
    shape = d_outer / d_in if d_in and d_in > 0 else np.nan

    # per-transit and per-season depth consistency
    per_tr = []
    for e in epochs:
        m = intr & (epoch == e)
        d, de = _wdepth(y, w, m, level)
        per_tr.append((e, d, de))
    per_tr = np.array(per_tr) if per_tr else np.zeros((0, 3))
    seas_rows = []
    for k, s in enumerate(seasons(t)):
        m = np.zeros(len(t), bool); m[s] = True
        d, de = _wdepth(y, w, intr & m, level)
        if np.isfinite(d):
            seas_rows.append({"season": k + 1, "t_mid": float(np.median(t[s])), "depth_ppm": d * 1e6,
                              "err_ppm": de * 1e6})
    from scipy import stats

    def chi2_test(d, e):
        if len(d) < 2:
            return np.nan, np.nan
        wt = 1 / e ** 2
        mu = np.sum(wt * d) / wt.sum()
        c2 = float(np.sum(((d - mu) / e) ** 2))
        return c2, float(stats.chi2.sf(c2, len(d) - 1))
    tr_chi2, tr_p = chi2_test(per_tr[:, 1], per_tr[:, 2]) if len(per_tr) else (np.nan, np.nan)
    sd = np.array([r["depth_ppm"] for r in seas_rows]); se = np.array([r["err_ppm"] for r in seas_rows])
    se_chi2, se_p = chi2_test(sd, se)
    # max relative deviation of the per-transit depths (only big ones matter)
    if len(per_tr) >= 2 and depth > 0:
        rel = float(np.max(np.abs(per_tr[:, 1] - depth)) / depth)
        top = float(np.max(per_tr[:, 1]) / depth)
    else:
        rel = top = np.nan

    # sinusoid vs box (variability leaking through)
    wv = 2 * np.pi * t / P
    X = np.column_stack([np.ones_like(t), np.sin(wv), np.cos(wv), np.sin(2 * wv), np.cos(2 * wv)])
    sw = np.sqrt(w_white)
    coef, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    chi_sin = float(np.sum(w_white * (y - np.average(y, weights=w_white)) ** 2)
                    - np.sum(w_white * (y - X @ coef) ** 2))
    sin_ratio = chi_sin / snr_white ** 2 if snr_white > 0 else np.nan

    # physics
    k = np.sqrt(max(depth, 0))
    rp = rstar * RSUN * k / REARTH
    a = (G * mstar * MSUN * (P * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    a_rs = a / (rstar * RSUN)
    tmax = P / np.pi * np.arcsin(min(1.0, (1 + k) / a_rs))
    lum = rstar ** 2 * (teff / 5772.0) ** 4
    s_earth = lum / (a / AU) ** 2
    teq = teff * np.sqrt(rstar * RSUN / (2 * a)) * (1 - 0.3) ** 0.25
    return {"period": P, "t0": t0, "duration_h": dur * 24, "depth_ppm": depth * 1e6,
            "depth_err_ppm": err * 1e6, "snr": snr, "snr_white": snr_white,
            "beta_med": float(np.median(beta)), "n_transits": n_tr,
            "depth_odd_ppm": d_odd * 1e6, "depth_even_ppm": d_even * 1e6, "oe_sigma": oe,
            "depth_sec_ppm": d_sec * 1e6, "sec_sigma": sec_sigma, "shape": shape,
            "sec_any_ppm": sec_any * 1e6 if np.isfinite(sec_any) else np.nan,
            "sec_any_sigma": sec_any_sigma, "sec_any_phase": sec_any_phase,
            "sin_ratio": sin_ratio, "tr_chi2_p": tr_p, "tr_max_rel_dev": rel, "tr_max_over_mean": top,
            "season_chi2_p": se_p, "season_depths": seas_rows, "rp_earth": rp, "a_rs": a_rs,
            "dur_ratio": dur / tmax if tmax > 0 else np.nan, "a_au": a / AU,
            "insolation": s_earth, "teq_k": teq, "hz": hz_status(s_earth, teff),
            "per_transit": per_tr.tolist()}


# --------------------------------------------------------------------------
# Habitable zone (Kopparapu et al. 2014, 1 Earth-mass coefficients)
# --------------------------------------------------------------------------
_HZ = {"recent_venus": (1.776, 2.136e-4, 2.533e-8, -1.332e-11, -3.097e-15),
       "runaway_greenhouse": (1.107, 1.332e-4, 1.580e-8, -8.308e-12, -1.931e-15),
       "maximum_greenhouse": (0.356, 6.171e-5, 1.698e-9, -3.198e-12, -5.575e-16),
       "early_mars": (0.320, 5.547e-5, 1.526e-9, -2.874e-12, -5.011e-16)}


def hz_seff(teff: float) -> dict:
    T = np.clip(teff, 2600, 7200) - 5780
    return {k: s + a * T + b * T ** 2 + c * T ** 3 + d * T ** 4 for k, (s, a, b, c, d) in _HZ.items()}


def hz_status(s_earth: float, teff: float) -> str:
    h = hz_seff(teff)
    if h["runaway_greenhouse"] >= s_earth >= h["maximum_greenhouse"]:
        return "conservative HZ"
    if h["recent_venus"] >= s_earth >= h["early_mars"]:
        return "optimistic HZ"
    return "too hot" if s_earth > h["recent_venus"] else "too cold"


# --------------------------------------------------------------------------
# Flags
# --------------------------------------------------------------------------
SNR_MIN, SDE_MIN = 7.1, 7.0


def vet(s: dict, sde_val: float) -> tuple[list[str], list[str]]:
    eb, caution = [], []
    if np.isfinite(s["rp_earth"]) and s["rp_earth"] > 15:
        eb.append(f"Rp~{s['rp_earth']:.0f} Re: too big for a planet (stellar/BD companion?)")
    if np.isfinite(s["oe_sigma"]) and s["oe_sigma"] > 3:
        eb.append(f"odd/even {s['oe_sigma']:.1f}σ")
    if (np.isfinite(s["sec_sigma"]) and s["sec_sigma"] > 3.5
            and s["depth_sec_ppm"] > 0.05 * s["depth_ppm"]):
        eb.append(f"secondary {s['sec_sigma']:.1f}σ")
    if (np.isfinite(s.get("sec_any_sigma", np.nan)) and s["sec_any_sigma"] > 5
            and s["sec_any_ppm"] > 0.1 * s["depth_ppm"]):
        eb.append(f"second eclipse at phase {s['sec_any_phase']:.2f} ({s['sec_any_ppm']:.0f} ppm, "
                  f"{s['sec_any_sigma']:.1f}σ): eccentric eclipsing binary?")
    if np.isfinite(s["sin_ratio"]) and s["sin_ratio"] > 0.5:
        eb.append(f"sinusoidal (variability; sine/box = {s['sin_ratio']:.1f})")
    if np.isfinite(s["shape"]) and s["shape"] < 0.45 and s["snr"] > 15:
        eb.append("V-shaped")
    if (np.isfinite(s["tr_chi2_p"]) and s["tr_chi2_p"] < 1e-6 and s["tr_max_over_mean"] > 3):
        eb.append(f"dominated by a few deep events (largest {s['tr_max_over_mean']:.0f}x mean)")
    if np.isfinite(s["season_chi2_p"]) and s["season_chi2_p"] < 1e-4:
        caution.append(f"depth differs between seasons (p={s['season_chi2_p']:.0e})")
    if s["n_transits"] < 3:
        caution.append(f"only {s['n_transits']} transits")
    if np.isfinite(s["dur_ratio"]):
        if s["dur_ratio"] > 1.6:
            caution.append(f"too long for star ({s['dur_ratio']:.1f}x)")
        elif s["dur_ratio"] < 0.15:
            caution.append(f"very short ({s['dur_ratio']:.2f}x)")
    for k in (1, 2, 0.5):
        if abs(s["period"] / (TESS_ORBIT * k) - 1) < 0.02:
            caution.append("P≈TESS orbit" + {1: "", 2: "x2", 0.5: "/2"}[k])
    if sde_val < SDE_MIN:
        caution.append(f"weak peak (SDE {sde_val:.1f})")
    return eb, caution


# --------------------------------------------------------------------------
# Centroid test across sectors
# --------------------------------------------------------------------------
def centroid_multi(lc: dict, P, t0, dur, depth_ppm) -> dict:
    """Per-sector in- vs out-of-transit centroid shift (pixels), combined as a
    chi^2 over sectors (pixel axes rotate between sectors, so shifts are not
    averaged directly). offset = |shift| / depth x 21"/px, where the source
    would have to sit relative to the light centroid."""
    from scipy import stats
    rows = []
    for s in np.unique(lc["sector"]):
        m = lc["sector"] == s
        t = lc["t"][m]
        near = in_tr(t, P, t0, dur, 2.0)
        intr = in_tr(t, P, t0, dur)
        if intr.sum() < 3 or (~near).sum() < 50:
            continue
        r = {}
        for ax in ("cx", "cy"):
            c = lc[ax][m].astype(float)
            if not np.isfinite(c).any():
                break
            tr = running_median(t, c, 1.0, near)
            res = c - tr
            o = ~near & np.isfinite(res)
            sig = 1.4826 * np.nanmedian(np.abs(res[o] - np.nanmedian(res[o])))
            ii = intr & np.isfinite(res)
            if ii.sum() < 3 or not np.isfinite(sig) or sig == 0:
                break
            r["d" + ax] = float(np.nanmean(res[ii]) - np.nanmedian(res[o]))
            r["e" + ax] = float(sig / np.sqrt(ii.sum()))
        if len(r) == 4:
            r["sector"] = int(s)
            rows.append(r)
    if not rows:
        return {"ok": False, "verdict": "not measured"}
    chi2 = sum((r["dcx"] / r["ecx"]) ** 2 + (r["dcy"] / r["ecy"]) ** 2 for r in rows)
    dof = 2 * len(rows)
    p = float(stats.chi2.sf(chi2, dof))
    depth = max(depth_ppm, 1) * 1e-6
    w = np.array([1 / (r["ecx"] ** 2 + r["ecy"] ** 2) for r in rows])
    shift = np.array([np.hypot(r["dcx"], r["dcy"]) for r in rows])
    noise = np.array([np.sqrt(r["ecx"] ** 2 + r["ecy"] ** 2) for r in rows])
    # bias-corrected mean shift magnitude
    mshift = float(np.sqrt(max(np.sum(w * (shift ** 2 - noise ** 2)) / w.sum(), 0)))
    off = mshift / depth * 21.0
    off_lim = float(np.sqrt(1 / w.sum()) / depth * 21.0 * 2)
    # Calibration (300 faint M dwarfs, 28 re-found known planets/TOIs): on
    # crowded fields the light centre moves in transit even for on-target
    # planets, with implied offsets of 2-60". So this is a weak test here;
    # depth-vs-crowding (mdwarf_hunt.crowding_test) is the one that discriminates.
    if p < 1e-3:
        v = (f"light centre moves in transit (p={p:.0e}); implied source ~{off:.0f}\" from it "
             "(known planets on crowded faint stars show 2-60\": not decisive)")
    else:
        v = f"no significant shift (p={p:.2f}); source within ~{max(off, 0) + off_lim:.0f}\" of the light centre"
    return {"ok": True, "n_sectors": len(rows), "chi2": float(chi2), "dof": dof, "p": p,
            "offset_arcsec": off, "offset_lim_arcsec": off_lim, "verdict": v, "per_sector": rows}
