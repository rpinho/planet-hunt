"""Deep follow-up of a second transit signal on a star that already hosts a known planet.

    python second_signal.py twomin --tic 259233660 --period 4.213463 --t0 1686.714262 \
        --duration-h 1.68 --known 0.448956:1684.47826:0.734 --tag tic259
    python second_signal.py ffi    ... (same arguments)
    python second_signal.py mcmc   ... --rstar 0.3775:0.011 --mstar 0.366:0.009 --teff 3419
    python second_signal.py null   ... ;  python second_signal.py sectornull ...
    python second_signal.py report --tag tic259 --control-tag tic259_ctl --verdict vetting/<verdict>.json

Stages (each writes runs/<tag>/<stage>.json, plots as PNG):

twomin   SPOC 2-minute PDCSAP of every sector. The known planet is masked. Then:
         per-sector depths with red-noise-inflated errors (beta per sector),
         constancy chi2, per-year and per-camera depths, depth vs CROWDSAP
         (a blend's PDCSAP depth tracks the contamination correction),
         coincidence with momentum dumps (QUALITY bit 32), odd/even, secondary,
         alias scan (P/3 .. 3P), period ratio to the known planet, and a
         local ephemeris scan. Saves the in-window data for `mcmc`.
ffi      Fixed-aperture TESScut photometry (pixelvet scene fit + lcvet) for
         EVERY TESScut sector, same pixel rule in each. Per-sector depth at the
         2-min ephemeris propagated with its error; sectors without 2-min
         data are combined separately (independent of the 2-min photons, but
         check whether the discovery search already used them).
null     Noise calibration: box SNR, sign test and Student-t depth at
         thousands of random ephemerides (same light curve, same duration).
sectornull  Per-sector depth scatter and depth-vs-time trend vs random ephemerides.
report   Collects every stage (candidate, an on-target control run with the
         same stages, the pixel localisation and TRICERATOPS outputs) and a
         hand-written verdict into runs/<tag>.json.
mcmc     batman + emcee. (a) the candidate alone, unbinned in-window 2-min
         data, stellar density free; (b) the known planet on its phase-folded
         2-min light curve, density free; (c) both together with ONE stellar
         density and shared limb darkening (circular orbits). If the two
         signals orbit the same star, (a) and (b) agree and (c) is a good fit.
         Gives the period/epoch posterior and Rp with errors.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from planet_hunt import config
from planet_hunt import multisector as ms
from planet_hunt import lcvet

G = 6.674e-8
RSUN_REARTH = 109.076


def _args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["twomin", "ffi", "mcmc", "null", "sectornull", "report"])
    ap.add_argument("--tic", type=int)
    ap.add_argument("--period", type=float)
    ap.add_argument("--t0", type=float)
    ap.add_argument("--duration-h", type=float)
    ap.add_argument("--depth-ppm", type=float, default=758.0)
    ap.add_argument("--known", default=None, help="known planet P:t0:duration_h (masked)")
    ap.add_argument("--rstar", default="0.3775:0.011", help="R*:err (Rsun)")
    ap.add_argument("--mstar", default="0.366:0.009", help="M*:err (Msun)")
    ap.add_argument("--teff", type=float, default=3419.0)
    ap.add_argument("--steps", type=int, default=6000)
    ap.add_argument("--walkers", type=int, default=32)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--control-tag", default=None, help="report: run folder of the on-target control")
    ap.add_argument("--verdict", default=None, help="report: hand-written verdict JSON")
    ap.add_argument("--localization", nargs="*", default=[], help="report: deep_vet per-target JSONs")
    ap.add_argument("--triceratops", nargs="*", default=[], help="report: run_triceratops output JSONs")
    ap.add_argument("--report-name", default=None, help="report: output stem under runs/")
    a = ap.parse_args()
    a.out = config.RUNS / a.tag
    a.out.mkdir(parents=True, exist_ok=True)
    return a


def tess_cycle(sector: int) -> int:
    """TESS observing cycle of a sector (Cycle 1 = S1-13, ..., Cycle 7 = S84-96, Cycle 8 = S97-...)."""
    edges = [13, 26, 39, 55, 69, 83, 96]
    return int(np.searchsorted(edges, sector) + 1)


def _known(a):
    if not a.known:
        return None
    p, t, d = (float(v) for v in a.known.split(":"))
    return p, t, d / 24


def _dumps(tic: int, sectors) -> np.ndarray:
    """Times of reaction-wheel desaturations (QUALITY bit 32) from the 2-min files."""
    from astropy.io import fits
    out = []
    for s in sectors:
        hit = sorted(config.MAST_DIR.glob(f"tess*-s{int(s):04d}-{tic:016d}-*-s_lc.fits"))
        if not hit:
            continue
        with fits.open(hit[0]) as hd:
            q = np.asarray(hd[1].data["QUALITY"]); t = np.asarray(hd[1].data["TIME"], float)
        m = ((q & 32) > 0) & np.isfinite(t)
        if m.any():
            tt = t[m]
            # one time per event (events span a few cadences)
            out.extend(tt[np.r_[True, np.diff(tt) > 0.02]])
    return np.array(out)


def _depth(t, y, w, P, t0, dur, sel=None):
    intr = ms.in_tr(t, P, t0, dur)
    oot = ~ms.in_tr(t, P, t0, dur, 3.0)
    if sel is not None:
        intr &= sel; oot &= sel
    if intr.sum() < 3 or oot.sum() < 10:
        return np.nan, np.nan
    lev = np.average(y[oot], weights=w[oot])
    return float(lev - np.average(y[intr], weights=w[intr])), float(1 / np.sqrt(w[intr].sum()))


def _with_beta(lc, P, t0, dur, known):
    """Per-sector red-noise factor on transit-length bins, applied to dy."""
    oot = ~ms.in_tr(lc["t"], P, t0, dur, 2.0)
    beta = ms.red_noise_beta(lc, dur, oot)
    lc = dict(lc)
    lc["beta"] = beta
    lc["dy_red"] = lc["dy"] * beta
    return lc


def load_2min(a):
    from crossyear import spoc_2min
    P, t0, dur = a.period, a.t0, a.duration_h / 24
    kn = _known(a)

    def mask_fn(t):
        m = ms.in_tr(t, P, t0, dur, 1.6)
        if kn:
            m |= ms.in_tr(t, kn[0], kn[1], kn[2], 1.6)
        return m
    lcs = spoc_2min(a.tic)
    lc = ms.stitch(list(lcs.values()), transit_mask=lambda t: mask_fn(t))
    lc["camera"] = np.zeros(len(lc["t"]), int)
    return lcs, lc


def stage_twomin(a):
    from scipy import stats
    from astropy.io import fits
    P, t0, dur = a.period, a.t0, a.duration_h / 24
    kn = _known(a)
    lcs, lc_all = load_2min(a)
    keep = np.ones(len(lc_all["t"]), bool)
    if kn:
        keep &= ~ms.in_tr(lc_all["t"], kn[0], kn[1], kn[2], 1.6)
    lc = {k: (v[keep] if isinstance(v, np.ndarray) and len(v) == len(keep) else v) for k, v in lc_all.items()}
    lc = _with_beta(lc, P, t0, dur, kn)
    t, y = lc["t"], lc["flat"]
    w_white = 1 / lc["dy"] ** 2
    w = 1 / lc["dy_red"] ** 2
    res = {"tic": a.tic, "ephemeris_in": {"period": P, "t0": t0, "duration_h": a.duration_h},
           "known_planet_mask": a.known, "sectors": sorted(int(s) for s in lcs)}
    # camera per sector from file headers
    cams = {}
    for s in lcs:
        hit = sorted(config.MAST_DIR.glob(f"tess*-s{int(s):04d}-{a.tic:016d}-*-s_lc.fits"))
        with fits.open(hit[0]) as hd:
            cams[int(s)] = int(hd[0].header["CAMERA"])
    res["camera"] = cams
    d, e_white = _depth(t, y, w_white, P, t0, dur)
    d_r, e_r = _depth(t, y, w, P, t0, dur)
    res["fold"] = {"depth_ppm": d * 1e6, "err_white_ppm": e_white * 1e6, "snr_white": d / e_white,
                   "depth_red_weighted_ppm": d_r * 1e6, "err_red_ppm": e_r * 1e6, "snr_red": d_r / e_r,
                   "median_beta": float(np.median(lc["beta"]))}
    # per sector (red-noise errors)
    rows = []
    for s in np.unique(lc["sector"]):
        m = lc["sector"] == s
        dd, ee = _depth(t, y, w, P, t0, dur, m)
        if np.isfinite(dd):
            rows.append({"sector": int(s), "camera": cams.get(int(s)), "depth_ppm": dd * 1e6, "err_ppm": ee * 1e6,
                         "beta": float(np.median(lc["beta"][m])),
                         "crowdsap": float(lcs[int(s)]["meta"]["crowdsap"]),
                         "n_transits": int(len(np.unique(np.round((t[m & ms.in_tr(t, P, t0, dur)] - t0) / P)))),
                         "t_mid": float(np.median(t[m]))})
    res["per_sector"] = rows
    res["per_sector_chi2"] = lcvet.chi2_const(rows)
    # per TESS year (cycle) and per camera
    cyc = tess_cycle
    years = {}
    for lab, fn in (("cycle", cyc), ("camera", lambda s: cams.get(s))):
        out = []
        for g in sorted(set(fn(r["sector"]) for r in rows)):
            secs = [r["sector"] for r in rows if fn(r["sector"]) == g]
            dd, ee = _depth(t, y, w, P, t0, dur, np.isin(lc["sector"], secs))
            out.append({lab: g, "sectors": secs, "depth_ppm": dd * 1e6, "err_ppm": ee * 1e6})
        years[lab] = out
        years[lab + "_chi2"] = lcvet.chi2_const(out)
    res["groups"] = years
    # depth vs contamination: weighted linear fit depth = a + b * (1 - crowdsap)
    x = np.array([1 - r["crowdsap"] for r in rows]); dd = np.array([r["depth_ppm"] for r in rows])
    ee = np.array([r["err_ppm"] for r in rows])
    A = np.vstack([np.ones_like(x), x]).T / ee[:, None]
    coef, *_ = np.linalg.lstsq(A, dd / ee, rcond=None)
    cov = np.linalg.inv(A.T @ A)
    res["depth_vs_contamination"] = {"intercept_ppm": float(coef[0]), "slope_ppm_per_unit": float(coef[1]),
                                     "slope_err": float(np.sqrt(cov[1, 1])),
                                     "slope_sigma": float(coef[1] / np.sqrt(cov[1, 1])),
                                     "contamination_range": [float(x.min()), float(x.max())]}
    # momentum dumps
    dumps = _dumps(a.tic, lcs)
    ep = np.round((t - t0) / P).astype(int)
    intr = ms.in_tr(t, P, t0, dur)
    eps = np.unique(ep[intr])
    tc = t0 + eps * P
    near = np.array([np.any(np.abs(dumps - c) < dur) for c in tc]) if len(dumps) else np.zeros(len(tc), bool)
    # expectation: random phases
    rng = np.random.default_rng(3)
    rand = []
    for _ in range(200):
        sh = rng.uniform(0.15, 0.85) * P
        rand.append(np.mean([np.any(np.abs(dumps - (c + sh)) < dur) for c in tc]))
    bad_eps = eps[near]
    d_nd, e_nd = _depth(t, y, w, P, t0, dur, ~np.isin(ep, bad_eps))
    res["momentum_dumps"] = {"n_dumps": int(len(dumps)), "n_transits": int(len(tc)),
                             "transits_within_1dur_of_dump": int(near.sum()),
                             "expected_random": float(np.mean(rand) * len(tc)),
                             "depth_without_those_ppm": d_nd * 1e6, "err_ppm": e_nd * 1e6}
    # per-epoch scatter with red-noise errors
    rows_e = []
    for e_ in eps:
        m = intr & (ep == e_)
        if m.sum() >= 0.5 * dur / (2 / 1440):
            dd_, ee_ = _depth(t, y, w, P, t0, dur, (ep == e_))
            if np.isfinite(dd_):
                rows_e.append((int(e_), dd_, ee_))
    rows_e = np.array(rows_e)
    z = (rows_e[:, 1] - np.average(rows_e[:, 1], weights=1 / rows_e[:, 2] ** 2)) / rows_e[:, 2]
    res["per_epoch"] = {"n": int(len(rows_e)), "chi2": float(np.sum(z ** 2)),
                        "p": float(stats.chi2.sf(np.sum(z ** 2), len(z) - 1)), "max_sigma": float(np.abs(z).max()),
                        "depth_median_epoch_ppm": float(np.median(rows_e[:, 1]) * 1e6),
                        "frac_positive": float(np.mean(rows_e[:, 1] > 0))}
    # robust depth: median-of-epochs bootstrap and clipping 5% both tails
    order = np.argsort(rows_e[:, 1] / rows_e[:, 2])
    n5 = max(1, int(round(0.05 * len(rows_e))))
    keep_e = rows_e[order[n5:-n5], 0].astype(int)
    d_c, e_c = _depth(t, y, w, P, t0, dur, np.isin(ep, keep_e) | ~intr)
    res["clip_both_tails_5pct"] = {"depth_ppm": d_c * 1e6, "err_ppm": e_c * 1e6, "snr": d_c / e_c}
    # odd / even / secondary (red-noise weights)
    bs = lcvet.box_stats({"t": t, "flat": y, "dy": lc["dy_red"], "cad": lc["cad"], "sector": lc["sector"]},
                         P, t0, dur)
    res["odd_even_secondary"] = {k: bs[k] for k in ("odd_ppm", "even_ppm", "oe_sigma", "sec_ppm", "sec_sigma",
                                                    "sec_any_sigma", "sec_any_phase", "sec_any_ppm")}
    # aliases: best SNR over t0 for P*f
    alias = []
    for f in (1 / 3, 1 / 2, 2 / 3, 1.0, 3 / 2, 2.0, 3.0):
        p = P * f
        best = (-np.inf, np.nan)
        for sh in np.arange(0, p, dur / 4):
            dd_, ee_ = _depth(t, y, w, p, t0 + sh, dur)
            if np.isfinite(dd_) and dd_ / ee_ > best[0]:
                best = (dd_ / ee_, sh)
        dd0, ee0 = _depth(t, y, w, p, t0, dur)
        alias.append({"factor": f, "period": p, "snr_at_t0": dd0 / ee0, "best_snr_any_phase": best[0],
                      "best_shift_d": best[1]})
    res["aliases"] = alias
    if kn:
        r = P / kn[0]
        fr = []
        for q in range(1, 9):
            pnum = round(r * q)
            fr.append({"p_over_q": f"{pnum}/{q}", "value": pnum / q, "frac_diff": (r - pnum / q) / r})
        fr.sort(key=lambda x: abs(x["frac_diff"]))
        # fraction of candidate in-transit cadences lost to the known-planet mask
        lost = np.mean(ms.in_tr(lc_all["t"][ms.in_tr(lc_all["t"], P, t0, dur)], kn[0], kn[1], kn[2], 1.6))
        res["known_planet_relation"] = {"period_ratio": r, "nearest_fractions": fr[:4],
                                        "candidate_in_transit_fraction_masked": float(lost),
                                        "n_times_known_period": [float(n * kn[0]) for n in (9, 10)]}
    # local ephemeris scan and Delta chi2 1-sigma box for P, t0 (box model)
    base = t.max() - t.min()
    Ps = P + np.linspace(-6, 6, 61) * dur / base * P / 4
    t0s = t0 + np.linspace(-0.6, 0.6, 49) * dur
    g = np.full((len(Ps), len(t0s)), np.nan)
    for i, p in enumerate(Ps):
        for j, tt in enumerate(t0s):
            dd_, ee_ = _depth(t, y, w, p, tt, dur)
            g[i, j] = dd_ / ee_
    i, j = np.unravel_index(np.nanargmax(g), g.shape)
    res["local_scan"] = {"best_period": float(Ps[i]), "best_t0": float(t0s[j]), "best_snr": float(g[i, j])}
    # save in-window data for the MCMC (candidate windows, and the known planet's fold)
    win = np.abs(ms.phase(t, P, t0)) < 2.5 * dur
    np.savez(a.out / "twomin_windows.npz", t=t[win], y=y[win], dy=lc["dy_red"][win], sector=lc["sector"][win])
    if kn:
        tk, yk, dyk = lc_all["t"], lc_all["flat"], lc_all["dy"]
        mk = ~ms.in_tr(tk, P, t0, dur, 1.6)
        np.savez(a.out / "twomin_known.npz", t=tk[mk], y=yk[mk], dy=dyk[mk])
    np.savez(a.out / "twomin_epochs.npz", rows=rows_e)
    (a.out / "twomin.json").write_text(json.dumps(res, indent=1, default=float))
    _plot_twomin(res, lc, a)
    print(json.dumps({k: v for k, v in res.items() if k not in ("per_sector",)}, indent=1, default=float))


def _plot_twomin(res, lc, a):
    import matplotlib.pyplot as plt
    P, t0, dur = a.period, a.t0, a.duration_h / 24
    fig, ax = plt.subplots(2, 2, figsize=(13, 8))
    rows = res["per_sector"]
    s = np.array([r["sector"] for r in rows]); d = np.array([r["depth_ppm"] for r in rows])
    e = np.array([r["err_ppm"] for r in rows]); cam = np.array([r["camera"] for r in rows])
    for c in np.unique(cam):
        m = cam == c
        ax[0, 0].errorbar(s[m], d[m], e[m], fmt="o", ms=4, label=f"camera {c}")
    mu = res["per_sector_chi2"]["mean"]
    ax[0, 0].axhline(mu, color="k", ls=":"); ax[0, 0].axhline(0, color="0.7")
    ax[0, 0].set_xlabel("sector"); ax[0, 0].set_ylabel("depth (ppm), red-noise errors")
    ax[0, 0].set_title(f"2-min per sector: chi2 p = {res['per_sector_chi2']['p']:.2g}", fontsize=10)
    ax[0, 0].legend(fontsize=8)
    x = np.array([1 - r["crowdsap"] for r in rows])
    ax[0, 1].errorbar(x, d, e, fmt="o", ms=4, color="0.3")
    dv = res["depth_vs_contamination"]
    xx = np.linspace(x.min(), x.max(), 10)
    ax[0, 1].plot(xx, dv["intercept_ppm"] + dv["slope_ppm_per_unit"] * xx, "r-")
    ax[0, 1].set_xlabel("1 - CROWDSAP (contamination corrected by PDC)"); ax[0, 1].set_ylabel("depth (ppm)")
    ax[0, 1].set_title(f"depth vs contamination: slope {dv['slope_sigma']:.1f} sigma", fontsize=10)
    ph = ms.phase(lc["t"], P, t0) * 24
    m = np.abs(ph) < 6
    so = np.argsort(ph[m])
    tb, yb, eb = ms.bin_weighted(ph[m][so] / 24, lc["flat"][m][so], lc["dy_red"][m][so], minutes=12)
    o = np.argsort(tb)
    ax[1, 0].errorbar(tb[o] * 24, (yb[o] - 1) * 1e6, eb[o] * 1e6, fmt="k.", ms=3)
    ax[1, 0].axvspan(-a.duration_h / 2, a.duration_h / 2, color="r", alpha=0.08)
    ax[1, 0].set_xlabel("hours from mid-transit"); ax[1, 0].set_ylabel("ppm")
    f = res["fold"]
    ax[1, 0].set_title(f"fold: {f['depth_red_weighted_ppm']:.0f} ± {f['err_red_ppm']:.0f} ppm "
                       f"(SNR {f['snr_red']:.1f} with red noise)", fontsize=10)
    al = res["aliases"]
    ax[1, 1].bar([f"{r['factor']:.2g}P" for r in al], [r["best_snr_any_phase"] for r in al], color="0.6",
                 label="best over phase")
    ax[1, 1].plot([f"{r['factor']:.2g}P" for r in al], [r["snr_at_t0"] for r in al], "ro", label="at t0")
    ax[1, 1].set_ylabel("box SNR"); ax[1, 1].legend(fontsize=8); ax[1, 1].set_title("aliases", fontsize=10)
    fig.tight_layout(); fig.savefig(a.out / "twomin.png", dpi=110); plt.close(fig)


# --------------------------------------------------------------------------
def stage_ffi(a):
    from planet_hunt import pixelvet as pv
    from crossyear import tesscut_lcs
    from astroquery.mast import Catalogs
    P, t0, dur = a.period, a.t0, a.duration_h / 24
    kn = _known(a)
    tw = json.loads((a.out / "twomin.json").read_text())
    post = a.out / "mcmc.json"
    if post.exists():
        mc = json.loads(post.read_text())["candidate_free"]
        P, t0 = mc["period"][1], mc["t0"][1]
        sP, st0 = (mc["period"][2] - mc["period"][0]) / 2, (mc["t0"][2] - mc["t0"][0]) / 2
    elif abs(tw["ephemeris_in"]["period"] - a.period) < 1e-3:
        P, t0 = tw["local_scan"]["best_period"], tw["local_scan"]["best_t0"]
        sP, st0 = 2e-5, 0.004
    else:
        sP, st0 = 2e-6, 0.002
    row = Catalogs.query_criteria(catalog="TIC", ID=a.tic).to_pandas().iloc[0]
    lcs = tesscut_lcs(float(row.ra), float(row.dec), row.GAIA, float(row.Tmag), P, t0, dur, set())
    twomin_secs = set(tw["sectors"])
    rows = []
    allt, ally, alldy, alls = [], [], [], []
    for s, d in sorted(lcs.items()):
        def mask_fn(tt):
            m = ms.in_tr(tt, P, t0, dur, 1.6)
            if kn:
                m |= ms.in_tr(tt, kn[0], kn[1], kn[2], 1.6)
            return m
        c = ms.clean_sector(d, window=0.75, transit_mask=mask_fn)
        if not c:
            continue
        keep = ~ms.in_tr(c["t"], kn[0], kn[1], kn[2], 1.6) if kn else np.ones(len(c["t"]), bool)
        tt, yy = c["t"][keep], c["flat"][keep]
        tmp = {"t": tt, "flat": yy, "cad": np.full(len(tt), c["cadence"]), "sector": np.full(len(tt), s)}
        beta = float(np.median(ms.red_noise_beta(tmp, dur, ~ms.in_tr(tt, P, t0, dur, 2.0))))
        dy = np.full(len(tt), c["sigma"] * beta)
        w = 1 / dy ** 2
        from planet_hunt.pixelvet import smear_factor
        sm = smear_factor(dur, c["cadence"])
        dd, ee = _depth(tt, yy, w, P, t0, dur)
        n = np.round((np.median(tt) - t0) / P)
        win = float(np.hypot(st0, abs(n) * sP))
        if d["meta"]["crowdsap"] <= 0.06 or not np.isfinite(ee) or ee / sm > 0.02:
            print(f"[ffi] S{s}: scene fit failed (dilution {d['meta']['crowdsap']:.2f}), dropped", flush=True)
            continue
        rows.append({"sector": int(s), "cadence_min": c["cadence"] * 1440, "beta": beta,
                     "dilution_frac": float(d["meta"]["crowdsap"]),
                     "depth_ppm": dd / sm * 1e6, "err_ppm": ee / sm * 1e6,
                     "has_2min": int(s) in twomin_secs, "ephemeris_window_min": win * 1440,
                     "camera": int(d["meta"]["camera"])})
        allt.append(tt); ally.append(yy); alldy.append(dy / 1); alls.append(np.full(len(tt), s))
    t = np.concatenate(allt); y = np.concatenate(ally); dy = np.concatenate(alldy); sec = np.concatenate(alls)
    w = 1 / dy ** 2
    res = {"ephemeris": {"period": P, "period_err": sP, "t0": t0, "t0_err": st0}, "per_sector": rows}
    res["per_sector_chi2"] = lcvet.chi2_const(rows)
    indep = [r["sector"] for r in rows if not r["has_2min"]]
    for lab, secs in (("all_ffi", [r["sector"] for r in rows]), ("ffi_only_sectors", indep),
                      ("sectors_with_2min", [r["sector"] for r in rows if r["has_2min"]])):
        m = np.isin(sec, secs)
        dd, ee = _depth(t, y, w, P, t0, dur, m)
        res[lab] = {"sectors": secs, "depth_ppm": dd * 1e6, "err_ppm": ee * 1e6, "snr": dd / ee}
    cyc = tess_cycle
    groups = []
    for g in sorted(set(cyc(r["sector"]) for r in rows)):
        secs = [r["sector"] for r in rows if cyc(r["sector"]) == g]
        dd, ee = _depth(t, y, w, P, t0, dur, np.isin(sec, secs))
        groups.append({"cycle": g, "sectors": secs, "depth_ppm": dd * 1e6, "err_ppm": ee * 1e6})
    res["per_cycle"] = groups
    res["per_cycle_chi2"] = lcvet.chi2_const(groups)
    lcd = {"t": t, "flat": y, "dy": dy, "cad": np.full(len(t), 200 / 86400.0), "sector": sec}
    bs = lcvet.box_stats(lcd, P, t0, dur)
    res["odd_even_secondary"] = {k: bs[k] for k in ("odd_ppm", "even_ppm", "oe_sigma", "sec_ppm", "sec_sigma",
                                                    "sec_any_sigma", "sec_any_phase", "sec_any_ppm")}
    # does the FFI-only fold peak at the predicted phase? scan t0 over the full period
    m = np.isin(sec, indep)
    scan = []
    for sh in np.arange(-0.5 * P, 0.5 * P, dur / 4):
        dd, ee = _depth(t, y, w, P, t0 + sh, dur, m)
        if np.isfinite(dd):
            scan.append((sh, dd / ee))
    scan = np.array(scan)
    res["ffi_only_phase_scan"] = {"snr_at_prediction": float(scan[np.argmin(np.abs(scan[:, 0])), 1]),
                                  "best_snr": float(scan[:, 1].max()),
                                  "best_shift_h": float(scan[np.argmax(scan[:, 1]), 0] * 24),
                                  "rank_of_prediction": int(np.sum(scan[:, 1] >= scan[np.argmin(np.abs(scan[:, 0])), 1])),
                                  "n_phases": int(len(scan)), "p99": float(np.percentile(scan[:, 1], 99))}
    (a.out / "ffi.json").write_text(json.dumps(res, indent=1, default=float))
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    for flag, col, lab in ((True, "C0", "sector has 2-min (same photons)"), (False, "C3", "FFI only (no 2-min)")):
        rr = [r for r in rows if r["has_2min"] == flag]
        ax[0].errorbar([r["sector"] for r in rr], [r["depth_ppm"] for r in rr], [r["err_ppm"] for r in rr],
                       fmt="o", ms=4, color=col, label=lab)
    ax[0].axhline(tw["fold"]["depth_red_weighted_ppm"], color="k", ls=":", label="2-min fold")
    ax[0].axhline(0, color="0.7"); ax[0].legend(fontsize=7); ax[0].set_xlabel("sector"); ax[0].set_ylabel("depth (ppm)")
    ax[0].set_title("fixed-aperture TESScut depth per sector", fontsize=10)
    for k_, (lab, col) in enumerate((("ffi_only_sectors", "C3"), ("all_ffi", "k"))):
        mm = np.isin(sec, res[lab]["sectors"])
        ph = ms.phase(t[mm], P, t0) * 24
        q = np.abs(ph) < 6
        so = np.argsort(ph[q])
        tb, yb, eb = ms.bin_weighted(ph[q][so] / 24, y[mm][q][so], dy[mm][q][so], minutes=20)
        o = np.argsort(tb)
        ax[1 + k_].errorbar(tb[o] * 24, (yb[o] - 1) * 1e6, eb[o] * 1e6, fmt=".", color=col, ms=3)
        ax[1 + k_].axvspan(-a.duration_h / 2, a.duration_h / 2, color="r", alpha=0.08)
        ax[1 + k_].set_title(f"{lab}: {res[lab]['depth_ppm']:.0f} ± {res[lab]['err_ppm']:.0f} ppm", fontsize=10)
        ax[1 + k_].set_xlabel("hours from predicted mid-transit")
    fig.tight_layout(); fig.savefig(a.out / "ffi.png", dpi=110); plt.close(fig)
    print(json.dumps({k: v for k, v in res.items() if k != "per_sector"}, indent=1, default=float))
    for r in rows:
        print(r)


# --------------------------------------------------------------------------
def _batman_setup(t, P, ld, supersample=1, exp=None):
    import batman
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w = 0.0, P, 0.03, 20.0, 89.0, 0.0, 90.0
    p.limb_dark, p.u = "quadratic", list(ld)
    if supersample > 1:
        m = batman.TransitModel(p, t, supersample_factor=supersample, exp_time=exp)
    else:
        m = batman.TransitModel(p, t)
    return p, m


def _a_over_rs(rho, P):
    return (G * rho * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)


def stage_mcmc(a):
    import emcee
    from scipy import optimize
    rng = np.random.default_rng(11)
    P0, t00, dur = a.period, a.t0, a.duration_h / 24
    kn = _known(a)
    tw = json.loads((a.out / "twomin.json").read_text())
    P0, t00 = tw["local_scan"]["best_period"], tw["local_scan"]["best_t0"]
    W = np.load(a.out / "twomin_windows.npz")
    t, y, dy = W["t"], W["y"], W["dy"]
    # reference epoch near the middle of the data (decorrelates P and t0)
    nmid = np.round((np.median(t) - t00) / P0)
    tref = t00 + nmid * P0
    u0 = lcvet.ld_coeffs(a.teff)
    rs, rse = (float(v) for v in a.rstar.split(":"))
    msun, mse = (float(v) for v in a.mstar.split(":"))
    rho_tic = 1.989e33 * msun / (4 / 3 * np.pi * (6.957e10 * rs) ** 3)
    rho_tic_err = rho_tic * np.hypot(mse / msun, 3 * rse / rs)
    # known planet: fold on its own (refined) ephemeris, 1-min phase bins
    K = np.load(a.out / "twomin_known.npz")
    Pk, tk0, durk = kn
    # refine known planet's period/epoch with a box scan
    best = (-np.inf, Pk, tk0)
    wk = 1 / K["dy"] ** 2
    basek = K["t"].max() - K["t"].min()
    for p in Pk + np.linspace(-3, 3, 31) * durk / basek * Pk / 6:
        for sh in np.linspace(-0.3, 0.3, 25) * durk:
            dd_, ee_ = _depth(K["t"], K["y"], wk, p, tk0 + sh, durk)
            if dd_ / ee_ > best[0]:
                best = (dd_ / ee_, p, tk0 + sh)
    Pk, tk0 = best[1], best[2]
    phk = ms.phase(K["t"], Pk, tk0)
    mk = np.abs(phk) < 3 * durk
    edges = np.arange(-3 * durk, 3 * durk + 1e-9, 1 / 1440)
    idx = np.digitize(phk[mk], edges)
    ws = np.bincount(idx, wk[mk], len(edges) + 1)
    ok = ws > 0
    xk = np.bincount(idx, wk[mk] * phk[mk], len(edges) + 1)[ok] / ws[ok]
    yk = np.bincount(idx, wk[mk] * K["y"][mk], len(edges) + 1)[ok] / ws[ok]
    ek = 1 / np.sqrt(ws[ok])
    # inflate the fold errors so the out-of-transit chi2/dof = 1
    oo = np.abs(xk) > durk
    sk = np.sqrt(max(1.0, np.mean(((yk[oo] - np.average(yk[oo], weights=1 / ek[oo] ** 2)) / ek[oo]) ** 2)))
    ek = ek * sk
    pk_, mkb = _batman_setup(xk, Pk, u0, supersample=5, exp=2 / 1440)   # 2-min exposures smear the fold
    pc_, mcb = _batman_setup(t - tref, P0, u0)

    def model_c(q):
        dt, dP, k, b, lrho, u1, u2, f0 = q
        P = P0 + dP
        aR = _a_over_rs(np.exp(lrho), P)
        if b >= aR:
            return None
        pc_.per, pc_.rp, pc_.a = P, k, aR
        pc_.inc = np.degrees(np.arccos(b / aR)); pc_.u = [u1, u2]
        pc_.t0 = dt          # batman repeats the transit every P
        return f0 * mcb.light_curve(pc_)

    def model_k(q):
        dt, k, b, lrho, u1, u2, f0 = q
        aR = _a_over_rs(np.exp(lrho), Pk)
        if b >= aR:
            return None
        pk_.rp, pk_.a, pk_.inc, pk_.u, pk_.t0 = k, aR, np.degrees(np.arccos(b / aR)), [u1, u2], dt
        return f0 * mkb.light_curve(pk_)

    def ld_prior(u1, u2):
        if not (0 < u1 + u2 < 1 and u1 > 0):
            return -np.inf
        return -0.5 * (((u1 - u0[0]) / 0.1) ** 2 + ((u2 - u0[1]) / 0.1) ** 2)

    def lp_c(q, rho_prior=False):
        dt, dP, k, b, lrho, u1, u2, f0 = q
        if not (0.005 < k < 0.2 and 0 <= b < 1 + k and np.log(0.5) < lrho < np.log(60) and abs(dt) < 0.05
                and abs(dP) < 1e-3):
            return -np.inf
        l = ld_prior(u1, u2)
        if not np.isfinite(l):
            return l
        m = model_c(q)
        if m is None:
            return -np.inf
        r = -0.5 * np.sum(((y - m) / dy) ** 2) + l
        if rho_prior:
            r += -0.5 * ((np.exp(lrho) - rho_tic) / rho_tic_err) ** 2
        return r

    def lp_k(q):
        dt, k, b, lrho, u1, u2, f0 = q
        if not (0.005 < k < 0.3 and 0 <= b < 1 + k and np.log(0.5) < lrho < np.log(60) and abs(dt) < 0.01):
            return -np.inf
        l = ld_prior(u1, u2)
        if not np.isfinite(l):
            return l
        m = model_k(q)
        if m is None:
            return -np.inf
        return -0.5 * np.sum(((yk - m) / ek) ** 2) + l

    def lp_joint(q):
        # shared lrho, u1, u2: [dtc, dPc, kc, bc, dtk, kk, bk, lrho, u1, u2, f0c, f0k]
        dtc, dPc, kc, bc, dtk, kk, bk, lrho, u1, u2, f0c, f0k = q
        a1 = lp_c([dtc, dPc, kc, bc, lrho, u1, u2, f0c])
        if not np.isfinite(a1):
            return -np.inf
        a2 = lp_k([dtk, kk, bk, lrho, u1, u2, f0k])
        return a1 + a2 - ld_prior(u1, u2)        # LD prior counted once

    def run(lp, x0, scale, steps, label):
        nd = len(x0)
        # start: optimise, then a small ball
        r = optimize.minimize(lambda q: -lp(q) if np.isfinite(lp(q)) else 1e25, x0, method="Nelder-Mead",
                              options={"maxiter": 8000, "xatol": 1e-7, "fatol": 1e-4})
        xb = r.x if np.isfinite(lp(r.x)) else np.array(x0)
        p0 = []
        while len(p0) < a.walkers:
            c = xb + scale * rng.standard_normal(nd)
            if np.isfinite(lp(c)):
                p0.append(c)
        s = emcee.EnsembleSampler(a.walkers, nd, lp)
        s.run_mcmc(np.array(p0), steps, progress=False)
        try:
            tau = float(np.max(s.get_autocorr_time(tol=0)))
        except Exception:
            tau = np.nan
        burn = steps // 3
        ch = s.get_chain(discard=burn, flat=True, thin=5)
        print(f"[mcmc] {label}: acc {np.mean(s.acceptance_fraction):.2f}, tau {tau:.0f}, n {len(ch)}", flush=True)
        return ch, float(np.mean(s.acceptance_fraction)), tau, r

    pct = lambda v: [float(np.percentile(v, q)) for q in (16, 50, 84)]
    out = {"rho_star_input_cgs": [rho_tic, rho_tic_err], "rstar": [rs, rse], "mstar": [msun, mse],
           "ld_prior_center": list(u0), "tref": tref, "known_refined": {"period": Pk, "t0": tk0},
           "known_fold_err_scale": float(sk)}
    kc0 = np.sqrt(tw["fold"]["depth_red_weighted_ppm"] * 1e-6)
    xc = [0.0, 0.0, kc0, 0.4, np.log(rho_tic), u0[0], u0[1], 1.0]
    scc = np.array([5e-4, 2e-7, 2e-3, 0.03, 0.05, 0.02, 0.02, 5e-6])
    chc, acc, tau, rc = run(lp_c, xc, scc, a.steps, "candidate, density free")
    Pp = P0 + chc[:, 1]
    rp = chc[:, 2] * (rs + rse * rng.standard_normal(len(chc))) * RSUN_REARTH
    aR = _a_over_rs(np.exp(chc[:, 4]), Pp)
    inc = np.arccos(chc[:, 3] / aR)
    k = chc[:, 2]
    T14 = Pp / np.pi * np.arcsin(np.sqrt(np.clip((1 + k) ** 2 - chc[:, 3] ** 2, 0, None)) / (aR * np.sin(inc))) * 24
    lr = chc[:, 4]
    from scipy import stats
    frac = float(np.mean(np.exp(lr) < rho_tic))
    out["candidate_free"] = {"t0": pct(tref + chc[:, 0]), "period": pct(Pp), "k": pct(k),
                             "depth_ppm": pct(k ** 2 * 1e6), "b": pct(chc[:, 3]),
                             "rho_cgs": pct(np.exp(lr)), "T14_h": pct(T14), "rp_earth": pct(rp),
                             "rho_star_z": float(stats.norm.ppf(np.clip(frac, 1e-4, 1 - 1e-4))),
                             "p_grazing": float(np.mean(chc[:, 3] > 1 - k)),
                             "acceptance": acc, "tau": tau, "n_samples": int(len(chc)),
                             "chi2_best": float(2 * rc.fun), "n_points": int(len(t))}
    # with the stellar-density prior: the radius posterior to quote
    chp, acc, tau, rp_ = run(lambda q: lp_c(q, True), xc, scc, a.steps, "candidate, TIC density prior")
    Pp2 = P0 + chp[:, 1]
    rp2 = chp[:, 2] * (rs + rse * rng.standard_normal(len(chp))) * RSUN_REARTH
    out["candidate_rho_prior"] = {"t0": pct(tref + chp[:, 0]), "period": pct(Pp2), "k": pct(chp[:, 2]),
                                  "depth_ppm": pct(chp[:, 2] ** 2 * 1e6), "b": pct(chp[:, 3]),
                                  "rho_cgs": pct(np.exp(chp[:, 4])), "rp_earth": pct(rp2),
                                  "acceptance": acc, "tau": tau}
    # null model: flat line chi2 vs best transit -> delta chi2 / BIC
    chi2_flat = float(np.sum(((y - np.average(y, weights=1 / dy ** 2)) / dy) ** 2))
    out["candidate_free"]["dchi2_vs_flat"] = chi2_flat - out["candidate_free"]["chi2_best"]
    out["candidate_free"]["dBIC_vs_flat"] = out["candidate_free"]["dchi2_vs_flat"] - 6 * np.log(len(t))
    # known planet alone
    xk0 = [0.0, np.sqrt(max(1e-5, 1 - np.min(yk))), 0.4, np.log(rho_tic), u0[0], u0[1], 1.0]
    sck = np.array([2e-4, 2e-3, 0.03, 0.05, 0.02, 0.02, 5e-6])
    chk, acc, tau, rk = run(lp_k, xk0, sck, a.steps, "known planet, density free")
    out["known_free"] = {"t0_shift_d": pct(chk[:, 0]), "k": pct(chk[:, 1]), "depth_ppm": pct(chk[:, 1] ** 2 * 1e6),
                         "b": pct(chk[:, 2]), "rho_cgs": pct(np.exp(chk[:, 3])),
                         "rp_earth": pct(chk[:, 1] * (rs + rse * rng.standard_normal(len(chk))) * RSUN_REARTH),
                         "acceptance": acc, "tau": tau, "n_bins": int(len(xk))}
    # joint, shared density
    xj = [chc[:, 0].mean(), chc[:, 1].mean(), np.median(chc[:, 2]), np.median(chc[:, 3]),
          np.median(chk[:, 0]), np.median(chk[:, 1]), np.median(chk[:, 2]), np.median(chk[:, 3]),
          u0[0], u0[1], 1.0, 1.0]
    scj = np.array([5e-4, 2e-7, 2e-3, 0.03, 2e-4, 2e-3, 0.03, 0.03, 0.02, 0.02, 5e-6, 5e-6])
    chj, acc, tau, rj = run(lp_joint, xj, scj, a.steps, "joint, shared density")
    out["joint_shared_density"] = {"rho_cgs": pct(np.exp(chj[:, 7])), "candidate_b": pct(chj[:, 3]),
                                   "candidate_k": pct(chj[:, 2]), "candidate_depth_ppm": pct(chj[:, 2] ** 2 * 1e6),
                                   "candidate_rp_earth": pct(chj[:, 2] * (rs + rse * rng.standard_normal(len(chj)))
                                                             * RSUN_REARTH),
                                   "candidate_period": pct(P0 + chj[:, 1]), "candidate_t0": pct(tref + chj[:, 0]),
                                   "known_b": pct(chj[:, 6]), "known_k": pct(chj[:, 5]),
                                   "acceptance": acc, "tau": tau,
                                   "chi2_best": float(2 * rj.fun),
                                   "dchi2_vs_separate": float(2 * rj.fun - 2 * rc.fun - 2 * rk.fun)}
    # density agreement between the two independent fits
    lc_ = np.log(out["candidate_free"]["rho_cgs"]); lk_ = np.log(out["known_free"]["rho_cgs"])
    sc_ = (lc_[2] - lc_[0]) / 2; sk_ = (lk_[2] - lk_[0]) / 2
    out["density_agreement_sigma"] = float((lc_[1] - lk_[1]) / np.hypot(sc_, sk_))
    # propagate the ephemeris to the next TESS visits (2027+, sectors 117-132)
    pe = out["candidate_free"]
    sP = (pe["period"][2] - pe["period"][0]) / 2; st0 = (pe["t0"][2] - pe["t0"][0]) / 2
    nxt = 4300.0   # BTJD ~ early 2027 (sector 117)
    n = np.ceil((nxt - pe["t0"][1]) / pe["period"][1])
    out["ephemeris_2027"] = {"btjd": float(pe["t0"][1] + n * pe["period"][1]),
                             "sigma_min": float(np.hypot(st0, n * sP) * 1440), "n_cycles": int(n)}
    (a.out / "mcmc.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez(a.out / "mcmc_chains.npz", candidate=chc, candidate_rho=chp, known=chk, joint=chj)
    _plot_mcmc(a, out, t, y, dy, tref, chc, model_c, xk, yk, ek, chk, model_k)
    print(json.dumps(out, indent=1, default=float))


def _plot_mcmc(a, out, t, y, dy, tref, chc, model_c, xk, yk, ek, chk, model_k):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    q = np.median(chc, axis=0)
    P = a.period + 0 if False else out["candidate_free"]["period"][1]
    ph = ms.phase(t, P, out["candidate_free"]["t0"][1]) * 24
    so = np.argsort(ph)
    tb, yb, eb = ms.bin_weighted(ph[so] / 24, y[so], dy[so], minutes=10)
    m = model_c(q)
    ax[0].errorbar(tb * 24, (yb - 1) * 1e6, eb * 1e6, fmt="k.", ms=3)
    tm, ym, _ = ms.bin_weighted(ph[so] / 24, m[so], np.ones(len(so)), minutes=4)
    ax[0].plot(tm * 24, (ym - 1) * 1e6, "r-")
    c = out["candidate_free"]
    ax[0].set_title(f"candidate: P = {c['period'][1]:.6f} d, Rp = {c['rp_earth'][1]:.2f} "
                    f"(+{c['rp_earth'][2]-c['rp_earth'][1]:.2f}/-{c['rp_earth'][1]-c['rp_earth'][0]:.2f}) Re",
                    fontsize=10)
    ax[0].set_xlabel("hours"); ax[0].set_ylabel("ppm")
    qk = np.median(chk, axis=0)
    ax[1].errorbar(xk * 24, (yk - 1) * 1e6, ek * 1e6, fmt="k.", ms=2)
    ax[1].plot(xk * 24, (model_k(qk) - 1) * 1e6, "r-")
    ax[1].set_title(f"known planet fold: depth {out['known_free']['depth_ppm'][1]:.0f} ppm", fontsize=10)
    ax[1].set_xlabel("hours")
    for lab, d_, col in (("candidate (free)", out["candidate_free"]["rho_cgs"], "C0"),
                         ("known planet (free)", out["known_free"]["rho_cgs"], "C1"),
                         ("joint", out["joint_shared_density"]["rho_cgs"], "k")):
        ax[2].errorbar([d_[1]], [lab], xerr=[[d_[1] - d_[0]], [d_[2] - d_[1]]], fmt="o", color=col)
    r, re = out["rho_star_input_cgs"]
    ax[2].axvspan(r - re, r + re, color="g", alpha=0.2, label="catalogue star")
    ax[2].set_xscale("log"); ax[2].set_xlabel("stellar density (g/cc)"); ax[2].legend(fontsize=8)
    ax[2].set_title(f"density agreement {out['density_agreement_sigma']:.1f} sigma", fontsize=10)
    fig.tight_layout(); fig.savefig(a.out / "mcmc.png", dpi=110); plt.close(fig)


def epoch_depths(t, y, w, P, t0, dur, win=3.0):
    """Per-transit depth against a local straight-line baseline fitted to the
    out-of-transit points within +-win durations (vectorised over epochs)."""
    ph = ms.phase(t, P, t0)
    m = np.abs(ph) < win * dur
    ph, yy, ww = ph[m], y[m], w[m]
    ep = np.round((t[m] - t0) / P).astype(int)
    ue, inv = np.unique(ep, return_inverse=True)
    intr = np.abs(ph) < 0.5 * dur
    oot = np.abs(ph) > 0.75 * dur
    n = len(ue)
    S = lambda v, sel: np.bincount(inv[sel], v[sel], n)
    so, sx, sxx = S(ww, oot), S(ww * ph, oot), S(ww * ph * ph, oot)
    sy, sxy = S(ww * yy, oot), S(ww * ph * yy, oot)
    det = so * sxx - sx ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        b1 = (so * sxy - sx * sy) / det
        b0 = (sy - b1 * sx) / so
        wi = S(ww, intr)
        yi = S(ww * yy, intr) / wi
        xi = S(ww * ph, intr) / wi
        d = b0 + b1 * xi - yi
        e = 1 / np.sqrt(wi)
    nin = np.bincount(inv[intr], minlength=n)
    nout = np.bincount(inv[oot], minlength=n)
    ok = np.isfinite(d) & (nin >= 0.5 * nin.max()) & (nout >= 10)
    return ue[ok], d[ok], e[ok]


def robust_stats(d, e, nu=4.0):
    """Sign test z, and the Student-t (nu) maximum-likelihood common depth with
    its curvature error."""
    from scipy import optimize
    npos = int(np.sum(d > 0)); n = len(d)
    z_sign = (npos - 0.5 * n) / np.sqrt(0.25 * n)
    nll = lambda mu: 0.5 * (nu + 1) * np.sum(np.log1p(((d - mu) / e) ** 2 / nu))
    r = optimize.minimize_scalar(nll, bounds=(-5 * np.abs(d).max(), 5 * np.abs(d).max()), method="bounded")
    h = 1e-6
    curv = (nll(r.x + h) - 2 * nll(r.x) + nll(r.x - h)) / h ** 2
    err = 1 / np.sqrt(curv) if curv > 0 else np.nan
    w = 1 / e ** 2
    mu_w = np.sum(w * d) / w.sum()
    chi2 = float(np.sum(((d - mu_w) / e) ** 2))
    return {"n": n, "n_positive": npos, "sign_z": float(z_sign), "t_depth": float(r.x), "t_err": float(err),
            "t_snr": float(r.x / err), "wmean": float(mu_w), "wmean_err": float(1 / np.sqrt(w.sum())),
            "chi2_per_dof": chi2 / max(1, n - 1)}


def stage_null(a):
    P, t0, dur = a.period, a.t0, a.duration_h / 24
    kn = _known(a)
    lcs, lc_all = load_2min(a)
    keep = np.ones(len(lc_all["t"]), bool)
    if kn:
        keep &= ~ms.in_tr(lc_all["t"], kn[0], kn[1], kn[2], 1.6)
    lc = {k: (v[keep] if isinstance(v, np.ndarray) and len(v) == len(keep) else v) for k, v in lc_all.items()}
    lc = _with_beta(lc, P, t0, dur, kn)
    t, y = lc["t"], lc["flat"]
    w = 1 / lc["dy_red"] ** 2
    rng = np.random.default_rng(7)
    d0, e0 = _depth(t, y, w, P, t0, dur)
    ue, de, ee = epoch_depths(t, y, w, P, t0, dur)
    obs = robust_stats(de, ee)
    obs.update({"box_snr": d0 / e0, "box_depth_ppm": d0 * 1e6})
    # null A: same period, other phases (away from the transit)
    nullA = []
    for sh in np.arange(1.5 * dur, P - 1.5 * dur, dur / 2):
        dd, e_ = _depth(t, y, w, P, t0 + sh, dur)
        if np.isfinite(dd):
            nullA.append(dd / e_)
    nullA = np.array(nullA)
    # null B: random periods 2-8 d away from the candidate's aliases, random phases
    bad = [P * f for f in (1 / 3, 1 / 2, 2 / 3, 1, 3 / 2, 2, 3)]
    if kn:
        bad += [kn[0] * k for k in range(1, 20)]
    nullB, robB = [], []
    while len(nullB) < 3000:
        p = rng.uniform(2.0, 8.0)
        if any(abs(p - b) / b < 0.01 for b in bad):
            continue
        tt0 = t0 + rng.uniform(0, p)
        dd, e_ = _depth(t, y, w, p, tt0, dur)
        if not np.isfinite(dd):
            continue
        nullB.append(dd / e_)
        if len(robB) < 600:
            u_, d_, e2 = epoch_depths(t, y, w, p, tt0, dur)
            if len(d_) > 30:
                robB.append(robust_stats(d_, e2))
    nullB = np.array(nullB)
    summ = lambda v: {"n": int(len(v)), "mean": float(np.mean(v)), "std": float(np.std(v)),
                      "p99": float(np.percentile(v, 99)), "max": float(np.max(v)),
                      "n_ge_observed": int(np.sum(v >= obs["box_snr"]))}
    res = {"observed": obs, "null_same_period_other_phases": summ(nullA), "null_random_ephemerides": summ(nullB)}
    for key in ("sign_z", "t_snr", "chi2_per_dof"):
        v = np.array([r[key] for r in robB])
        res[f"null_{key}"] = {"mean": float(np.mean(v)), "std": float(np.std(v)),
                              "n_ge_observed": int(np.sum(v >= obs[key])), "n": int(len(v))}
    sB = res["null_random_ephemerides"]["std"]
    res["calibrated_snr"] = float((obs["box_snr"] - res["null_random_ephemerides"]["mean"]) / sB)
    res["calibrated_sign_z"] = float((obs["sign_z"] - res["null_sign_z"]["mean"]) / res["null_sign_z"]["std"])
    res["calibrated_t_snr"] = float((obs["t_snr"] - res["null_t_snr"]["mean"]) / res["null_t_snr"]["std"])
    # look-elsewhere: how many independent box trials did the discovery search make?
    base = t.max() - t.min()
    pg = np.geomspace(0.5, 30, 4000)
    dlnp = np.log(pg[1] / pg[0])
    # independent periods per d(ln P): base/(dur/4) ... times independent phases P/(dur/2)
    n_trials = float(np.sum(pg * dlnp / (pg * (dur / 4) / base) * (pg / (dur / 2))))
    from scipy import stats
    res["look_elsewhere"] = {"approx_independent_trials_0p5_30d": float(n_trials),
                             "gaussian_threshold_1pct_FA": float(stats.norm.isf(0.01 / n_trials)),
                             "note": "trial count is a rough estimate (phase steps of dur/2, period steps "
                                     "of dur/4 drift over the baseline); compare with the calibrated SNR"}
    np.savez(a.out / "null.npz", nullA=nullA, nullB=nullB, epochs=ue, d=de, e=ee)
    (a.out / "null.json").write_text(json.dumps(res, indent=1, default=float))
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
    ax[0].hist(nullB, bins=60, color="0.6", density=True, label="random ephemerides")
    ax[0].hist(nullA, bins=30, histtype="step", color="C0", density=True, label="same P, other phases")
    ax[0].axvline(obs["box_snr"], color="r", label=f"candidate {obs['box_snr']:.1f}")
    ax[0].set_xlabel("box SNR (red-noise errors)"); ax[0].legend(fontsize=8)
    ax[0].set_title(f"null std {sB:.2f}: calibrated SNR {res['calibrated_snr']:.1f}", fontsize=10)
    tt_ = t0 + ue * P
    ax[1].errorbar(tt_, de * 1e6, ee * 1e6, fmt="k.", ms=3, alpha=0.6)
    ax[1].axhline(obs["t_depth"] * 1e6, color="r"); ax[1].axhline(0, color="0.6")
    ax[1].set_xlabel("BTJD"); ax[1].set_ylabel("per-transit depth (ppm), local baseline")
    ax[1].set_title(f"{obs['n_positive']}/{obs['n']} transits positive (sign z {obs['sign_z']:.1f})", fontsize=10)
    v = np.array([r["sign_z"] for r in robB])
    ax[2].hist(v, bins=40, color="0.6"); ax[2].axvline(obs["sign_z"], color="r")
    ax[2].set_xlabel("sign-test z (null ephemerides)")
    fig.tight_layout(); fig.savefig(a.out / "null.png", dpi=110); plt.close(fig)
    print(json.dumps(res, indent=1, default=float))


def stage_sectornull(a):
    """Is the per-sector depth scatter (and any trend with time) unusual for a
    box of this duration on this light curve? Same statistics at 300 random
    ephemerides."""
    P, t0, dur = a.period, a.t0, a.duration_h / 24
    kn = _known(a)
    lcs, lc = load_2min(a)
    if kn:
        keep = ~ms.in_tr(lc["t"], kn[0], kn[1], kn[2], 1.6)
        lc = {k: (v[keep] if isinstance(v, np.ndarray) and len(v) == len(keep) else v) for k, v in lc.items()}
    lc = _with_beta(lc, P, t0, dur, kn)
    t, y, w, sec = lc["t"], lc["flat"], 1 / lc["dy_red"] ** 2, lc["sector"]
    secs = np.unique(sec)
    idx = {s: sec == s for s in secs}

    def stats_(p, tt0):
        rows = []
        for s in secs:
            d, e = _depth(t, y, w, p, tt0, dur, idx[s])
            if np.isfinite(d):
                rows.append({"depth_ppm": d * 1e6, "err_ppm": e * 1e6, "t": float(np.median(t[idx[s]]))})
        c = lcvet.chi2_const(rows)
        tm = np.array([r["t"] for r in rows]); dd = np.array([r["depth_ppm"] for r in rows])
        ee = np.array([r["err_ppm"] for r in rows])
        A = np.vstack([np.ones_like(tm), (tm - 2500) / 365.25]).T / ee[:, None]
        co, *_ = np.linalg.lstsq(A, dd / ee, rcond=None)
        cov = np.linalg.inv(A.T @ A)
        return c["chi2"] / (c["n"] - 1), co[1] / np.sqrt(cov[1, 1]), co[1], np.sqrt(cov[1, 1])
    obs = stats_(P, t0)
    rng = np.random.default_rng(5)
    nul = []
    bad = [P * f for f in (1 / 3, .5, 2 / 3, 1, 1.5, 2, 3)]
    while len(nul) < 300:
        p = rng.uniform(2, 8)
        if any(abs(p - b) / b < 0.01 for b in bad):
            continue
        nul.append(stats_(p, t0 + rng.uniform(0, p)))
    nul = np.array(nul)
    out = {"observed": {"sector_chi2_dof": obs[0], "trend_sigma": obs[1], "trend_ppm_per_yr": obs[2],
                        "trend_err": obs[3]},
           "null_sector_chi2_dof": {"mean": float(nul[:, 0].mean()), "std": float(nul[:, 0].std()),
                                    "frac_ge_obs": float(np.mean(nul[:, 0] >= obs[0]))},
           "null_trend_sigma": {"mean": float(nul[:, 1].mean()), "std": float(nul[:, 1].std()),
                                "frac_le_obs": float(np.mean(nul[:, 1] <= obs[1]))}}
    (a.out / "sector_null.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1))


def product_paired(tag_candidate: str, tag_control: str) -> dict:
    """Fixed-aperture FFI minus 2-min depth, sector by sector, for the candidate
    and for an on-target control on the same star (same pixels, same dilution
    correction). Errors from the scatter of the differences (the products
    share photons, so quadrature errors overestimate)."""
    out = {}
    for tag, lab in ((tag_candidate, "candidate"), (tag_control, "control")):
        tw = {r["sector"]: r for r in json.loads((config.RUNS / tag / "twomin.json").read_text())["per_sector"]}
        ff = {r["sector"]: r for r in json.loads((config.RUNS / tag / "ffi.json").read_text())["per_sector"]}
        common = sorted(set(tw) & set(ff))
        d = np.array([ff[s]["depth_ppm"] - tw[s]["depth_ppm"] for s in common])
        w = 1 / np.array([ff[s]["err_ppm"] ** 2 + tw[s]["err_ppm"] ** 2 for s in common])
        mu = float(np.sum(w * d) / w.sum())
        chi = float(np.sum(w * (d - mu) ** 2) / (len(d) - 1))
        err = float(np.sqrt(chi / w.sum()))
        f2 = json.loads((config.RUNS / tag / "twomin.json").read_text())["fold"]["depth_red_weighted_ppm"]
        out[lab] = {"n_sectors": len(common), "mean_diff_ppm": mu, "err_ppm": err, "sigma": mu / err,
                    "relative": mu / f2, "relative_err": err / f2, "chi2_dof_of_differences": chi}
    c, k = out["candidate"], out["control"]
    out["relative_excess_difference_sigma"] = float((c["relative"] - k["relative"]) /
                                                    np.hypot(c["relative_err"], k["relative_err"]))
    return out


def stage_report(a):
    rd = lambda p: json.loads(Path(p).read_text()) if Path(p).exists() else None
    out = {"verdict": rd(a.verdict) if a.verdict else None}
    for st in ("twomin", "ffi", "mcmc", "null", "sector_null"):
        v = rd(a.out / f"{st}.json")
        if v is not None:
            out[st] = v
    if a.control_tag:
        ctl = config.RUNS / a.control_tag
        out["control"] = {st: rd(ctl / f"{st}.json") for st in ("twomin", "ffi") if (ctl / f"{st}.json").exists()}
        out["product_paired"] = product_paired(a.tag, a.control_tag)
    loc = {}
    for p in a.localization:
        v = rd(p)
        if v:
            L = v["localization"]
            loc[v["spec"]["name"]] = {"n_sectors": v["n_sectors"], "offset_arcsec": L["offset_arcsec"],
                                      "offset_err_arcsec": L["offset_err_arcsec"],
                                      "target_dchi2": L["target_dchi2"],
                                      "hypotheses": L["hypotheses"][:12], "rings": v.get("rings"),
                                      "flux_budget": v.get("flux_budget"), "transit_fit": v.get("transit_fit")}
    out["localization"] = loc
    out["triceratops"] = {Path(p).stem: {k: v for k, v in rd(p).items() if k != "runs"} | {"runs": rd(p)["runs"]}
                          for p in a.triceratops if rd(p)}
    name = a.report_name or a.tag
    (config.RUNS / f"{name}.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"wrote runs/{name}.json")


def main():
    a = _args()
    if a.stage == "sectornull":
        return stage_sectornull(a)
    if a.stage == "report":
        return stage_report(a)
    {"twomin": stage_twomin, "ffi": stage_ffi, "mcmc": stage_mcmc, "null": stage_null}.get(a.stage, lambda a: None)(a)


if __name__ == "__main__":
    main()
