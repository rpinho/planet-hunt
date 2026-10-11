"""Transit search (Box Least Squares) plus simple, explainable vetting.

The vetting here is the "first-pass triage" every transit survey does. It
catches obvious eclipsing binaries and junk; it is NOT enough to call
anything a planet (see README / report for what real vetting involves).
"""
from __future__ import annotations

import numpy as np
from astropy.timeseries import BoxLeastSquares

from .data import CleanLC, detrend, in_transit

# Trial transit durations in days (0.75 h ... 8 h). BLS tries each one.
DURATIONS = np.array([0.03, 0.045, 0.06, 0.08, 0.10, 0.125, 0.16, 0.20, 0.25, 0.33])

G = 6.674e-11
MSUN, RSUN, REARTH, DAY = 1.989e30, 6.957e8, 6.371e6, 86400.0
TESS_ORBIT = 13.7  # days; scattered light / momentum dumps repeat on this clock


# --------------------------------------------------------------------------
# BLS
# --------------------------------------------------------------------------
def bin_lc(t, f, minutes: float = 10.0):
    """Average into fixed time bins (speeds BLS ~5x; fine for >1 h transits)."""
    width = minutes / 1440.0
    idx = np.floor((t - t[0]) / width).astype(int)
    counts = np.bincount(idx)
    good = counts > 0
    tb = np.bincount(idx, t)[good] / counts[good]
    fb = np.bincount(idx, f)[good] / counts[good]
    return tb, fb


def max_duration(period, density_factor: float = 2.0):
    """Longest plausible transit (days) at this period for a Sun-like star,
    times a generous `density_factor`: T ~ 13 h * (P / 1 yr)^(1/3)."""
    return density_factor * (13 / 24) * (np.asarray(period) / 365.25) ** (1 / 3)


def bls_search(t, flat, pmin=0.5, pmax=None, durations=DURATIONS,
               bin_minutes=10.0, frequency_factor=1.0, n_chunks=12):
    """Run astropy BoxLeastSquares on a flattened light curve.

    Period grid: uniform in frequency, spacing ~ min(duration)/baseline^2
    (astropy's `autoperiod`), so a transit cannot slip between trial periods.
    pmax defaults to half the baseline, guaranteeing >= 2 transits are
    possible.

    Durations: only physically plausible ones are tried at each period
    (<= 2x a central transit across a Sun-like star). Without this, BLS
    happily fits 8-hour "transits" at 0.5-day periods to stellar variability.
    Implemented by running BLS on ~12 period chunks, each with its own
    duration list; the log-likelihood power is comparable across chunks.
    Returns a dict with the peak and the full periodogram.
    """
    tb, fb = bin_lc(t, flat, bin_minutes)
    baseline = t.max() - t.min()
    pmax = min(pmax or 15.0, baseline / 2)
    bls = BoxLeastSquares(tb, fb)
    periods = bls.autoperiod(durations, minimum_period=pmin,
                             maximum_period=pmax,
                             frequency_factor=frequency_factor)
    periods = np.sort(np.asarray(periods))
    edges = np.geomspace(pmin, pmax, n_chunks + 1)
    P, POW, T0, DUR, DEP = [], [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        chunk = periods[(periods >= lo) & (periods <= hi)]
        if len(chunk) == 0:
            continue
        durs = durations[(durations <= max_duration(hi)) & (durations < 0.5 * lo)]
        if len(durs) == 0:
            durs = durations[:1]
        r = bls.power(chunk, durs, objective="likelihood")
        P.append(np.asarray(r.period)); POW.append(np.asarray(r.power))
        T0.append(np.asarray(r.transit_time)); DUR.append(np.asarray(r.duration))
        DEP.append(np.asarray(r.depth))
    P, POW, T0, DUR, DEP = map(np.concatenate, (P, POW, T0, DUR, DEP))
    order = np.argsort(P)
    P, POW, T0, DUR, DEP = P[order], POW[order], T0[order], DUR[order], DEP[order]
    i = int(np.nanargmax(POW))
    # Signal Detection Efficiency: how far the peak stands above the
    # periodogram's typical level (after removing a slow running median).
    from scipy.ndimage import median_filter
    detr = POW - median_filter(POW, size=301, mode="nearest")
    sde = (detr[i] - np.nanmean(detr)) / np.nanstd(detr)
    return {"period": float(P[i]), "t0": float(T0[i]), "duration": float(DUR[i]),
            "depth_bls": float(DEP[i]), "sde": float(sde), "periods": P,
            "power": POW}


# --------------------------------------------------------------------------
# Statistics on the unbinned light curve
# --------------------------------------------------------------------------
def _depth(f, sel, oot_level, sigma):
    n = int(sel.sum())
    if n == 0:
        return np.nan, np.nan, 0
    return oot_level - float(np.mean(f[sel])), sigma / np.sqrt(n), n


def transit_stats(lc: CleanLC, period: float, t0: float, duration: float) -> dict:
    """Depth, SNR and vetting statistics for one (P, t0, duration) signal."""
    t, f = lc.time, lc.flat
    cad = lc.cadence
    intr = in_transit(t, period, t0, duration)
    oot = ~in_transit(t, period, t0, duration, widen=2.0)  # keep clear of edges
    oot &= ~in_transit(t, period, t0 + 0.5 * period, duration, widen=2.0)
    if oot.sum() < 100:            # duty cycle too high to define a baseline
        oot = ~intr
    level = float(np.median(f[oot]))
    sigma = 1.4826 * float(np.median(np.abs(f[oot] - level)))   # robust scatter

    depth, depth_err, n_in = _depth(f, intr, level, sigma)

    # Which individual transits actually have data? (>= half the cadences)
    epoch = np.round((t - t0) / period).astype(int)
    expected = duration / cad
    epochs = sorted({e for e in np.unique(epoch[intr])
                     if (intr & (epoch == e)).sum() >= 0.5 * expected})
    n_tr = len(epochs)

    # Transits sitting at the edge of a data gap (or the start/end of the
    # data) are the classic systematic: thermal ramps and imperfect
    # detrending at segment boundaries make fake dips.
    gaps = np.where(np.diff(t) > 0.1)[0]
    edges = np.r_[t[0], t[gaps], t[gaps + 1], t[-1]]
    n_edge = sum(1 for e in epochs
                 if np.min(np.abs(edges - (t0 + e * period))) < 0.5 * duration + 0.25)

    # White-noise SNR and a red-noise-aware SNR. The second uses the scatter
    # of the light curve averaged over one transit duration, which captures
    # correlated noise (stellar granulation, residual systematics).
    snr_white = depth / depth_err if depth_err > 0 else 0.0
    nbin = max(1, int(round(duration / cad)))
    fo = f[oot] - level
    if len(fo) > 3 * nbin:
        running = np.convolve(fo, np.ones(nbin) / nbin, mode="valid")[::nbin]
        sigma_dur = float(np.std(running))
    else:
        sigma_dur = sigma / np.sqrt(nbin)
    snr_red = depth / (sigma_dur / np.sqrt(max(n_tr, 1))) if sigma_dur > 0 else 0.0
    snr = float(min(snr_white, snr_red))

    # Noise inflation factor for the per-transit comparisons below: if the
    # noise is red, white-noise error bars are too small by ~snr_white/snr_red.
    beta = max(1.0, snr_white / snr_red) if snr_red > 0 else 1.0

    # Odd vs even transits: an eclipsing binary at twice the BLS period shows
    # alternating primary/secondary eclipses of different depth.
    odd = intr & (epoch % 2 == 1)
    even = intr & (epoch % 2 == 0)
    d_odd, e_odd, _ = _depth(f, odd, level, sigma)
    d_even, e_even, _ = _depth(f, even, level, sigma)
    oe_sigma = (abs(d_odd - d_even) / (beta * np.hypot(e_odd, e_even))
                if np.isfinite(d_odd) and np.isfinite(d_even) else np.nan)

    # Secondary eclipse at phase 0.5: planets are (almost always) invisible
    # there in TESS data; a stellar companion usually is not.
    sec = in_transit(t, period, t0 + 0.5 * period, duration)
    d_sec, e_sec, _ = _depth(f, sec, level, sigma)
    sec_sigma = d_sec / (beta * e_sec) if e_sec and np.isfinite(d_sec) else np.nan

    # Shape: depth in the central half vs the outer half of the transit.
    # A flat-bottomed planet transit gives ~0.7-1; a V-shape (grazing EB) << 1.
    phase = (t - t0 + 0.5 * period) % period - 0.5 * period
    inner = np.abs(phase) < 0.25 * duration
    outer = intr & ~inner
    d_in, e_in, _ = _depth(f, inner, level, sigma)
    d_out, _, _ = _depth(f, outer, level, sigma)
    shape = d_out / d_in if d_in and d_in > 0 else np.nan

    # Sinusoid test: how much of the signal does a smooth sine wave at P (and
    # P/2, for ellipsoidal binaries) explain, compared with the box? For a
    # genuine short transit the ratio is ~2x the duty cycle (a few %); for
    # spots/pulsation/contact binaries leaking through the flattening it is ~1.
    w = 2 * np.pi * t / period
    X = np.column_stack([np.ones_like(t), np.sin(w), np.cos(w), np.sin(2 * w), np.cos(2 * w)])
    coef, *_ = np.linalg.lstsq(X, f, rcond=None)
    chi_sin = (np.sum((f - f.mean()) ** 2) - np.sum((f - X @ coef) ** 2)) / sigma ** 2
    sin_ratio = float(chi_sin / snr_white ** 2) if snr_white > 0 else np.nan

    # Physical plausibility from the star's radius and mass.
    star = lc.star
    rp_earth = a_rs = dur_ratio = np.nan
    if star.get("rad"):
        rp_earth = star["rad"] * RSUN * np.sqrt(max(depth, 0)) / REARTH
    if star.get("rad") and star.get("mass"):
        a = (G * star["mass"] * MSUN * (period * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
        a_rs = a / (star["rad"] * RSUN)
        k = np.sqrt(max(depth, 0))
        tmax = period / np.pi * np.arcsin(min(1.0, (1 + k) / a_rs))  # b = 0
        dur_ratio = duration / tmax

    return {"period": period, "t0": t0, "duration_h": duration * 24,
            "depth_ppm": depth * 1e6, "depth_err_ppm": depth_err * 1e6,
            "snr": snr, "snr_white": snr_white, "snr_red": snr_red,
            "n_transits": n_tr, "n_edge": n_edge, "n_in": n_in, "sigma_ppm": sigma * 1e6,
            "depth_odd_ppm": d_odd * 1e6, "depth_even_ppm": d_even * 1e6,
            "oe_sigma": oe_sigma, "depth_sec_ppm": d_sec * 1e6,
            "sec_sigma": sec_sigma, "shape": shape, "sin_ratio": sin_ratio, "rp_earth": rp_earth,
            "a_rs": a_rs, "dur_ratio": dur_ratio}


# --------------------------------------------------------------------------
# Vetting flags
# --------------------------------------------------------------------------
SNR_MIN = 7.1   # the classic Kepler/TESS detection threshold
SDE_MIN = 7.0   # periodogram peak must stand out from the forest of noise peaks


def vet(s: dict, sde: float) -> tuple[list[str], list[str]]:
    """Return (eb_flags, caution_flags). EB flags = probably not a planet.
    Caution flags = could still be a planet but deserves suspicion."""
    eb, caution = [], []
    if s["depth_ppm"] > 30_000:
        eb.append("deep (>3%)")
    if np.isfinite(s["rp_earth"]) and s["rp_earth"] > 22:
        eb.append(f"Rp~{s['rp_earth']:.0f} Re (>2 RJup)")
    if np.isfinite(s["oe_sigma"]) and s["oe_sigma"] > 3:
        eb.append(f"odd/even {s['oe_sigma']:.1f}σ")
    if (np.isfinite(s["sec_sigma"]) and s["sec_sigma"] > 3
            and s["depth_sec_ppm"] > 0.1 * s["depth_ppm"]):
        eb.append(f"secondary {s['sec_sigma']:.1f}σ")
    if np.isfinite(s["sin_ratio"]) and s["sin_ratio"] > 0.5:
        eb.append(f"sinusoidal (variability; sine/box = {s['sin_ratio']:.1f})")
    if (np.isfinite(s["shape"]) and s["shape"] < 0.45 and s["snr"] > 12):
        eb.append("V-shaped")
    if s["n_transits"] >= 2 and s["n_edge"] >= 0.5 * s["n_transits"]:
        eb.append(f"{s['n_edge']}/{s['n_transits']} transits at data-gap edges")
    elif s["n_edge"] > 0:
        caution.append(f"{s['n_edge']}/{s['n_transits']} transit(s) at a gap edge")
    if s["n_transits"] < 2:
        caution.append("<2 transits")
    if np.isfinite(s["dur_ratio"]):
        if s["dur_ratio"] > 1.6:
            caution.append(f"too long for star ({s['dur_ratio']:.1f}x)")
        elif s["dur_ratio"] < 0.2:
            caution.append(f"very short ({s['dur_ratio']:.2f}x)")
    for k in (1, 2):
        if abs(s["period"] - TESS_ORBIT / k) / (TESS_ORBIT / k) < 0.04:
            caution.append("P≈TESS orbit" + ("/2" if k == 2 else ""))
    if sde < SDE_MIN:
        caution.append(f"weak peak (SDE {sde:.1f})")
    elif sde < 9:
        caution.append(f"modest peak (SDE {sde:.1f})")
    return eb, caution


def analyze(lc0: CleanLC, raw_lc, meta, pmin=0.5, pmax=15.0) -> dict:
    """Full per-star pipeline: BLS on the first-pass flattened light curve,
    then re-flatten with the transit masked and measure everything on that.
    Returns a result dict plus the arrays needed for plotting."""
    b = bls_search(lc0.time, lc0.flat, pmin=pmin, pmax=pmax)
    lc = detrend(raw_lc, meta, transit=(b["period"], b["t0"], b["duration"]))
    s = transit_stats(lc, b["period"], b["t0"], b["duration"])
    eb, caution = vet(s, b["sde"])
    s.update(sde=b["sde"], eb_flags=eb, caution_flags=caution,
             tce=bool(s["snr"] >= SNR_MIN and s["n_transits"] >= 2
                      and b["sde"] >= SDE_MIN))
    return {"stats": s, "bls": b, "lc": lc}
