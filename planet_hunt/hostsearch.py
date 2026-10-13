"""Search for extra transiting planets around a star that already hosts one (round 5).

Input: the star's SPOC 2-minute PDCSAP sectors (round5_hosts.load_lc), the
known TOI/planet ephemerides, and TIC radius / mass / Teff.

1. Known planets are REMOVED (every cadence within 0.75 duration of a known
   transit, i.e. a window 1.5x the duration). Planets without a transit
   ephemeris (RV-only) are not masked but their periods are flagged.
2. Rotation: Lomb-Scargle per sector on 30-min bins (0.1-12 d). On active
   stars (semi-amplitude > max(1000 ppm, 5x the 30-min noise)) a harmonic
   series at P_rot (3 harmonics, amplitudes linear in time) is fitted per
   sector and subtracted before the running filter, and the filter window is
   shortened to P_rot/6 (>= 0.25 d).
3. Detrending: wotan biweight per sector (gaps > 0.5 d split segments,
   0.25 d trimmed after each gap, flares > 4 sigma clipped). Two windows:
   0.5 d feeds the 0.5-6 d search, 1.0 d the 6-30 d search (a 3-4 h transit
   at 25 d around an early M dwarf needs the long window).
4. Search: fully coherent BLS on 10-min bins of ALL sectors (astropy, star-aware
   durations), on a frequency grid whose step keeps the phase drift over the
   whole baseline below 1/4 of a central transit duration. The 25 strongest
   distinct peaks per range get a quick red-noise SNR; the best is measured
   on the un-binned light curve re-detrended with its transits masked
   (multisector.transit_stats: per-sector red-noise beta, odd/even,
   secondary, per-transit and per-season consistency, sinusoid test).
5. The signal is masked (and its 1/3-3x aliases dropped from the peak pool)
   and the next best peak is taken, up to max_signals. The BLS itself runs
   once per star (a full re-run per iteration costs 3x for little gain on
   weak signals; injection-recovery uses exactly the same procedure).
"""
from __future__ import annotations

import numpy as np
from astropy.timeseries import BoxLeastSquares, LombScargle

from . import lcvet
from . import multisector as ms

REARTH, RSUN = 6.371e6, 6.957e8
SPLIT_P = 6.0
W_SHORT, W_LONG = 0.5, 1.0


# --------------------------------------------------------------------------
# masks and helpers
# --------------------------------------------------------------------------
def known_list(known, rstar, mstar):
    out = []
    for k in known:
        P, t0, d = k.get("period"), k.get("t0"), k.get("duration_h")
        P = float(P) if P is not None else np.nan
        t0 = float(t0) if t0 is not None else np.nan
        d = float(d) / 24 if d is not None and np.isfinite(float(d)) else np.nan
        if not np.isfinite(P) or P <= 0:
            continue
        if not np.isfinite(d):
            d = float(ms.central_duration(P, rstar, mstar))
        masked = np.isfinite(t0) and not (k.get("src") == "pscomppars" and k.get("transiting", 1) == 0)
        out.append({"name": k.get("name"), "period": P, "t0": t0, "dur": d, "masked": masked,
                    "disp": k.get("disp")})
    return out


def mask_known(t, kl, widen=1.5):
    m = np.zeros(len(t), bool)
    for k in kl:
        if k["masked"]:
            m |= ms.in_tr(t, k["period"], k["t0"], k["dur"], widen)
    return m


def inject_transit(secs, inj, rstar, mstar, teff):
    """Multiply a limb-darkened transit into the raw normalised flux."""
    P = inj["period"]
    k = inj["rp_earth"] * REARTH / (rstar * RSUN)
    rho = 1.41 * mstar / rstar ** 3
    u1, u2 = lcvet.ld_coeffs(teff)
    tmin = min(float(d["time"].min()) for d in secs)
    t0 = tmin + inj["phase"] * P
    inj["t0"] = t0
    a_rs = (6.674e-8 * rho * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    b = inj["b"]
    inj["duration"] = float(P / np.pi * np.arcsin(min(1.0, np.sqrt(max((1 + k) ** 2 - b ** 2, 1e-6)) / a_rs)))
    inj["depth_ppm_nominal"] = float(k ** 2 * 1e6)
    out = []
    for d in secs:
        d = dict(d)
        t = d["time"]
        ph = ms.phase(t, P, t0)
        near = np.abs(ph) < 1.5 * inj["duration"] + 0.01
        f = d["flux"].astype(float).copy()
        if near.any():
            cad = np.full(near.sum(), float(np.median(np.diff(t))))
            f[near] *= lcvet.transit_model(ph[near], cad, P, k, b, rho, u1, u2, nsub=3)
        d["flux"] = f
        out.append(d)
    return out


# --------------------------------------------------------------------------
# rotation
# --------------------------------------------------------------------------
def _bin(t, y, minutes):
    idx = np.floor((t - t[0]) / (minutes / 1440)).astype(int)
    n = np.bincount(idx)
    g = n > 0
    return np.bincount(idx, t)[g] / n[g], np.bincount(idx, y)[g] / n[g], n[g]


def rotation(secs, kmask_fn) -> dict:
    rows = []
    for d in secs:
        t, f = d["time"], d["flux"].astype(float)
        m = ~kmask_fn(t)
        if m.sum() < 1000:
            continue
        tb, fb, n = _bin(t[m], f[m], 30)
        if len(tb) < 200:
            continue
        span = tb.max() - tb.min()
        ls = LombScargle(tb, fb)
        fr = np.linspace(1 / min(12, 0.5 * span), 10.0, 6000)
        pw = ls.power(fr)
        j = int(np.argmax(pw))
        model = ls.model(tb, fr[j])
        amp = float(np.sqrt(2) * np.std(model) * 1e6)
        noise30 = float(1.4826 * np.median(np.abs(np.diff(fb))) / np.sqrt(2) * 1e6)
        rows.append({"sector": d["meta"]["sector"], "P": float(1 / fr[j]), "power": float(pw[j]),
                     "amp_ppm": amp, "noise30_ppm": noise30})
    if not rows:
        return {"P_rot": np.nan, "amp_ppm": 0.0, "active": False, "rows": []}
    P = np.array([r["P"] for r in rows]); pw = np.array([r["power"] for r in rows])
    amp = np.array([r["amp_ppm"] for r in rows]); nz = np.array([r["noise30_ppm"] for r in rows])
    good = pw > 0.1
    if good.sum() >= 2:
        Prot = float(np.median(P[good]))
        consist = float(np.mean(np.abs(P[good] / Prot - 1) < 0.1))
    else:
        Prot, consist = float(np.median(P)), 0.0
    a_med = float(np.median(amp)); n_med = float(np.median(nz))
    active = bool(a_med > max(1000.0, 5 * n_med) and Prot < 10 and consist >= 0.3)
    return {"P_rot": Prot, "amp_ppm": a_med, "noise30_ppm": n_med, "consistency": consist,
            "active": active, "rows": rows}


def harmonic_fit(t, f, Prot, mask, nh=3):
    """Per-sector harmonic series at Prot with linearly varying amplitudes."""
    tm = t - np.median(t)
    cols = [np.ones_like(t), tm]
    for h in range(1, nh + 1):
        w = 2 * np.pi * h * t / Prot
        for base in (np.sin(w), np.cos(w)):
            cols += [base, base * tm]
    X = np.column_stack(cols)
    ok = ~mask
    coef, *_ = np.linalg.lstsq(X[ok], f[ok], rcond=None)
    return X @ coef


# --------------------------------------------------------------------------
# detrending
# --------------------------------------------------------------------------
def detrend(secs, window, mask_fn, rot=None, gap=0.5, edge=0.25):
    """Masked biweight per sector. mask_fn(t) -> cadences excluded from the
    trend (known planets, candidate transits); the trend under them is
    interpolated. Returns the stitched light curve dict for multisector."""
    from wotan import flatten
    parts = []
    for d in secs:
        t = np.asarray(d["time"], float)
        f = np.asarray(d["flux"], float)
        ok = np.isfinite(t) & np.isfinite(f)
        t, f = t[ok], f[ok]
        cx, cy = np.asarray(d["cx"], float)[ok], np.asarray(d["cy"], float)[ok]
        if len(t) < 500:
            continue
        starts = np.r_[t[0], t[1:][np.diff(t) > gap]]
        keep = np.ones(len(t), bool)
        for s in starts:
            keep &= ~((t >= s) & (t < s + edge))
        t, f, cx, cy = t[keep], f[keep], cx[keep], cy[keep]
        m = mask_fn(t)
        w = window
        if rot is not None and rot.get("active"):
            f = f / harmonic_fit(t, f, rot["P_rot"], m)
            w = min(window, max(0.25, rot["P_rot"] / 6))
        fm = np.where(m, np.nan, f)
        # first pass trend on unmasked cadences, then clip outliers and redo
        trend = np.full(len(t), np.nan)
        good = ~m
        if good.sum() < 200:
            continue
        _, tr = flatten(t[good], f[good], method="biweight", window_length=w, edge_cutoff=0,
                        break_tolerance=gap, return_trend=True)
        trend[good] = tr
        ok_tr = np.isfinite(trend)
        if ok_tr.sum() < 200:
            continue
        trend = np.interp(t, t[ok_tr], trend[ok_tr])
        flat = f / trend
        resid = flat - 1
        sig = 1.4826 * np.nanmedian(np.abs(resid[~m] - np.nanmedian(resid[~m])))
        up = resid > 4 * sig
        up = up | np.r_[False, up[:-1]] | np.r_[up[1:], False]
        g = np.isfinite(flat) & ~up
        parts.append({"t": t[g], "flat": flat[g], "cx": cx[g], "cy": cy[g], "sigma": float(sig),
                      "cadence": float(np.median(np.diff(t))), "sector": int(d["meta"]["sector"]),
                      "camera": int(d["meta"].get("camera", 0)), "crowdsap": d["meta"].get("crowdsap", np.nan),
                      "window": w})
    if not parts:
        return {}
    cat = lambda k: np.concatenate([p[k] for p in parts])
    out = {k: cat(k) for k in ("t", "flat", "cx", "cy")}
    out["dy"] = np.concatenate([np.full(len(p["t"]), p["sigma"]) for p in parts])
    out["cad"] = np.concatenate([np.full(len(p["t"]), p["cadence"]) for p in parts])
    out["sector"] = np.concatenate([np.full(len(p["t"]), p["sector"]) for p in parts])
    o = np.argsort(out["t"])
    out = {k: v[o] for k, v in out.items()}
    out["sector_info"] = [{k: p[k] for k in ("sector", "sigma", "cadence", "camera", "crowdsap", "window")}
                          | {"n": len(p["t"])} for p in parts]
    return out


def subset(lc, keep):
    return {k: (v[keep] if isinstance(v, np.ndarray) and len(v) == len(keep) else v) for k, v in lc.items()}


# --------------------------------------------------------------------------
# coherent BLS
# --------------------------------------------------------------------------
def coherent_bls(tb, yb, eb, pmin, pmax, rstar, mstar, oversample=2.0, qfrac=0.5):
    T = tb.max() - tb.min()
    periods = ms.period_grid(pmin, pmax, T, rstar, mstar, oversample=oversample, qfrac=qfrac)
    P, pw, T0, DUR, DEP = ms._bls_chunks(tb, yb, eb, periods, rstar, mstar, n_chunks=12)
    return {"period": P, "power": pw, "t0": T0, "dur": DUR, "depth": DEP}


def top_peaks(res, n=15, sep=0.002):
    order = np.argsort(res["power"])[::-1]
    picked = []
    for i in order:
        P = res["period"][i]
        if all(abs(P / res["period"][j] - 1) > sep for j in picked):
            picked.append(i)
        if len(picked) >= n:
            break
    return picked


def quick_snr(tb, yb, eb, sec_b, P, t0, dur):
    """Depth and red-noise SNR on the binned light curve (for ranking)."""
    intr = ms.in_tr(tb, P, t0, dur)
    if intr.sum() < 3:
        return 0.0, 0.0, 0
    oot = ~ms.in_tr(tb, P, t0, dur, 3.0)
    if oot.sum() < 10:
        return 0.0, 0.0, 0
    w = 1 / eb ** 2
    lev = np.average(yb[oot], weights=w[oot])
    dep = lev - np.average(yb[intr], weights=w[intr])
    # per-transit depths -> scatter-based error (absorbs red noise)
    ep = np.round((tb - t0) / P).astype(int)
    ue = np.unique(ep[intr])
    d_e, w_e = [], []
    for e in ue:
        m = intr & (ep == e)
        d_e.append(lev - np.average(yb[m], weights=w[m])); w_e.append(w[m].sum())
    d_e, w_e = np.array(d_e), np.array(w_e)
    err_w = 1 / np.sqrt(w_e.sum())
    if len(d_e) >= 4:
        chi = np.sum(w_e * (d_e - dep) ** 2) / (len(d_e) - 1)
        err = err_w * max(1.0, np.sqrt(chi))
    else:
        err = err_w
    return float(dep), float(dep / err) if err > 0 else 0.0, int(len(ue))


# --------------------------------------------------------------------------
# signal-level checks
# --------------------------------------------------------------------------
def ratio_flags(P, kl, rot):
    out = []
    for k in kl:
        r = P / k["period"]
        for num in range(1, 6):
            for den in range(1, 6):
                if np.gcd(num, den) != 1 or (num > 3 and den > 1) or (den > 3 and num > 1):
                    continue
                if abs(r / (num / den) - 1) < 0.003:
                    out.append(f"P = {num}/{den} x {k['name']}" if den > 1 else
                               (f"P = {k['name']}" if num == 1 else f"P = {num} x {k['name']}"))
    if rot and np.isfinite(rot.get("P_rot", np.nan)) and rot.get("amp_ppm", 0) > 300:
        for f in (0.25, 1 / 3, 0.5, 1, 2, 3):
            if abs(P / (rot["P_rot"] * f) - 1) < 0.02:
                out.append(f"P near rotation x {f:.2g} (P_rot {rot['P_rot']:.2f} d)")
    return out


def known_overlap(t, P, t0, dur, kl):
    """Fraction of the candidate's in-transit cadences that fall within 1.5
    durations of a known planet's transit (after masking these are removed,
    so this uses the raw times)."""
    intr = ms.in_tr(t, P, t0, dur)
    if intr.sum() == 0:
        return 0.0
    m = np.zeros(intr.sum(), bool)
    for k in kl:
        if k["masked"]:
            m |= ms.in_tr(t[intr], k["period"], k["t0"], k["dur"], 3.0)
    return float(m.mean())


def dump_coincidence(secs, P, t0, dur):
    dumps = np.concatenate([np.asarray(d["dumps"], float) for d in secs if len(d["dumps"])]) \
        if any(len(d["dumps"]) for d in secs) else np.zeros(0)
    if len(dumps) == 0:
        return {"n_dumps": 0}
    ev = dumps[np.r_[True, np.diff(np.sort(dumps)) > 0.02]]
    near = np.abs(ms.phase(ev, P, t0)) < dur
    # expected: fraction of phase within +-dur
    t_all = np.concatenate([d["time"] for d in secs])
    n_tr = len(np.unique(np.round((t_all - t0) / P)))
    exp = len(ev) * min(1.0, 2 * dur / P)
    return {"n_dumps": int(len(ev)), "n_near_transit": int(near.sum()), "expected": float(exp),
            "n_epochs": int(n_tr)}


def robust_clip_snr(stats_, frac=0.05):
    """Drop the deepest 5% of per-transit depths (Kepler-robovetter style) and
    recompute the weighted depth / SNR from the per-transit table."""
    pt = np.array(stats_["per_transit"])
    if len(pt) < 5:
        return np.nan
    d, e = pt[:, 1], pt[:, 2]
    ok = np.isfinite(d) & np.isfinite(e) & (e > 0)
    d, e = d[ok], e[ok]
    nd = max(1, int(np.ceil(frac * len(d))))
    keep = np.argsort(d)[:-nd]
    w = 1 / e[keep] ** 2
    dep = np.sum(w * d[keep]) / w.sum()
    err = 1 / np.sqrt(w.sum())
    # scale by the white/red ratio of the full measurement
    scale = stats_["snr"] / stats_["snr_white"] if stats_["snr_white"] > 0 else 1.0
    return float(dep / err * scale)


# --------------------------------------------------------------------------
# main per-star search
# --------------------------------------------------------------------------
def search_star(secs, known, rstar, mstar, teff, inject=None, max_signals=3, n_peaks=25,
                snr_stop=6.0, bin_min=10.0, pmin=0.5, pmax=30.0):
    kl = known_list(known, rstar, mstar)
    if inject is not None:
        secs = inject_transit(secs, inject, rstar, mstar, teff)
    kmask = lambda t: mask_known(t, kl, 1.5)
    rot = rotation(secs, kmask)
    lcs = {}
    for name, W in (("short", W_SHORT), ("long", W_LONG)):
        lc = detrend(secs, W, kmask, rot)
        lcs[name] = subset(lc, ~kmask(lc["t"]))       # known transits removed
    T_all = float(lcs["long"]["t"].max() - lcs["long"]["t"].min())
    sig_mask = []                                     # found signals (P, t0, dur)

    def sigm(t, widen=2.0):
        m = np.zeros(len(t), bool)
        for (P, t0, d) in sig_mask:
            m |= ms.in_tr(t, P, t0, d, widen)
        return m

    # one coherent BLS per range; the top distinct peaks are then ranked by a
    # quick red-noise SNR with the already-found signals masked (iterations)
    pool, star_peaks, binned = [], [], {}
    for name, (lo, hi) in (("short", (pmin, SPLIT_P)), ("long", (SPLIT_P, min(pmax, 0.5 * T_all)))):
        lc = lcs[name]
        tb, yb, eb = ms.bin_weighted(lc["t"], lc["flat"], lc["dy"], bin_min)
        binned[name] = (tb, yb, eb)
        res = coherent_bls(tb, yb, eb, lo, hi, rstar, mstar)
        sde = ms.sde(res["power"])
        for i in top_peaks(res, n_peaks):
            pool.append({"range": name, "period": float(res["period"][i]), "t0": float(res["t0"][i]),
                         "duration": float(res["dur"][i]), "depth_bls": float(res["depth"][i]),
                         "power": float(res["power"][i]), "sde": float(sde[i])})
    star_peaks = sorted(pool, key=lambda c: -c["power"])[:6]
    signals = []
    t_raw = np.concatenate([d["time"] for d in secs])
    for it in range(max_signals):
        cands = []
        for c in pool:
            if any(abs(c["period"] / P_ * f - 1) < 0.003 for (P_, _, _) in sig_mask
                   for f in (1, 2, 0.5, 3, 1 / 3, 1.5, 2 / 3)):
                continue
            tb, yb, eb = binned[c["range"]]
            keep = ~sigm(tb)
            dep, snr, ntr = quick_snr(tb[keep], yb[keep], eb[keep], None, c["period"], c["t0"], c["duration"])
            cands.append(c | {"quick_depth_ppm": dep * 1e6, "quick_snr": snr, "quick_ntr": ntr})
        cands = [c for c in cands if c["quick_ntr"] >= 2]
        if not cands:
            break
        best = max(cands, key=lambda c: c["quick_snr"])
        P, t0, dur = best["period"], best["t0"], best["duration"]
        # final measurement: re-detrend with this transit masked too
        W = W_SHORT if best["range"] == "short" else W_LONG
        fm = lambda t: kmask(t) | sigm(t) | ms.in_tr(t, P, t0, dur, 1.6)
        lc = detrend(secs, W, fm, rot)
        lc = subset(lc, ~(kmask(lc["t"]) | sigm(lc["t"])))
        s = ms.transit_stats(lc, P, t0, dur, rstar, mstar, teff)
        eb_flags, caution = ms.vet(s, best["sde"])
        rf = ratio_flags(P, kl, rot)
        s.update({k: best[k] for k in ("range", "sde", "power", "quick_snr", "depth_bls")})
        s["snr_clip5"] = robust_clip_snr(s)
        s["eb_flags"], s["cautions"], s["ratio_flags"] = eb_flags, caution, rf
        s["known_overlap"] = known_overlap(t_raw, P, t0, dur, kl)
        s["dumps"] = dump_coincidence(secs, P, t0, dur)
        s["iteration"] = it + 1
        s["per_sector"] = per_sector_depths(lc, P, t0, dur)
        signals.append(s)
        sig_mask.append((P, t0, dur))
        if best["quick_snr"] < snr_stop and s["snr"] < snr_stop:
            break
    # star-level rotation saturation: strongest raw peaks at rotation harmonics,
    # or the strongest signal mostly sinusoidal
    sat_reasons = []
    if rot.get("active") or rot.get("amp_ppm", 0) > 1000:
        harm = [c for c in star_peaks if any("rotation" in f for f in ratio_flags(c["period"], [], rot))]
        if len(harm) >= 2:
            sat_reasons.append(f"{len(harm)} of the 6 strongest peaks at rotation harmonics")
    if signals and np.isfinite(signals[0].get("sin_ratio", np.nan)) and signals[0]["sin_ratio"] > 0.5:
        sat_reasons.append(f"strongest signal sinusoidal (sine/box {signals[0]['sin_ratio']:.1f})")
    if signals and signals[0].get("beta_med", 1) > 2.0:
        sat_reasons.append(f"red noise beta {signals[0]['beta_med']:.1f}")
    for s in signals:
        s.pop("season_depths_raw", None)
    out = {"rotation": {k: v for k, v in rot.items() if k != "rows"}, "rotation_rows": rot.get("rows", []),
           "rotation_saturated": bool(sat_reasons), "saturation_reasons": sat_reasons,
           "known": kl, "n_sectors": len(secs), "sectors": [int(d["meta"]["sector"]) for d in secs],
           "T_all": T_all, "n_points": int(len(lcs["long"]["t"])),
           "noise_ppm_2min": float(np.median(lcs["long"]["dy"]) * 1e6),
           "signals": signals, "star_peaks": star_peaks}
    if inject is not None:
        out["injection"] = inject
        out["injection_result"] = injection_match(inject, signals, lcs, rstar, mstar, teff)
    return out


def per_sector_depths(lc, P, t0, dur):
    rows = []
    intr_all = ms.in_tr(lc["t"], P, t0, dur)
    oot_all = ~ms.in_tr(lc["t"], P, t0, dur, 3.0)
    for s in np.unique(lc["sector"]):
        m = lc["sector"] == s
        intr, oot = intr_all & m, oot_all & m
        if intr.sum() < 3:
            continue
        w = 1 / lc["dy"][m] ** 2
        wi = 1 / lc["dy"][intr] ** 2
        lev = np.average(lc["flat"][oot], weights=1 / lc["dy"][oot] ** 2)
        d = lev - np.average(lc["flat"][intr], weights=wi)
        rows.append({"sector": int(s), "depth_ppm": float(d * 1e6), "err_ppm": float(1e6 / np.sqrt(wi.sum())),
                     "n_in": int(intr.sum())})
    return rows


def injection_match(inj, signals, lcs, rstar, mstar, teff):
    P, t0, dur = inj["period"], inj["t0"], inj["duration"]
    lc = lcs["short"] if P < SPLIT_P else lcs["long"]
    tb, yb, eb = ms.bin_weighted(lc["t"], lc["flat"], lc["dy"], 10.0)
    dep, snr_true, ntr_true = quick_snr(tb, yb, eb, None, P, t0, dur)
    out = {"snr_true_quick": snr_true, "depth_true_ppm": dep * 1e6, "ntr_true": ntr_true,
           "recovered": False, "matched": False}
    for s in signals:
        for fac, lab in ((1.0, ""), (0.5, "half"), (2.0, "double")):
            if abs(s["period"] / (P * fac) - 1) < 0.01:
                dph = abs(ms.phase(np.array([s["t0"]]), P * fac if fac < 1 else P, t0))[0]
                if dph < max(dur, s["duration_h"] / 24):
                    out.update({"matched": True, "alias": lab, "iteration": s["iteration"],
                                "snr": s["snr"], "n_transits": s["n_transits"], "depth_ppm": s["depth_ppm"],
                                "eb_flags": s["eb_flags"], "period_found": s["period"]})
                    out["recovered"] = bool(lab == "" and s["snr"] >= 7.5 and s["n_transits"] >= 3)
                    out["recovered_clean"] = bool(out["recovered"] and not s["eb_flags"])
                    return out
    return out
