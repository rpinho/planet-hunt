"""Independent re-derivation of every TIC 144439305 number quoted in the RNAAS draft.

Deliberately does NOT reuse the round-7 pipeline (no planet_hunt imports, no cached npz):
  * reads the raw MAST HLSP FITS (TESS-SPOC PDCSAP for S2/29/69, QLP DET_FLUX for S96/103/104)
    and the raw TESScut cutouts (S2-S106) itself;
  * detrends with a 0.4 d biweight (the pipeline used other windows), transits masked;
  * blind astropy BLS over 0.5-2 d on all six sectors;
  * box depths with a red-noise (time-averaging) beta, per data set;
  * per-sector mid-times -> linear ephemeris;
  * own background-subtracted aperture photometry on the TESScut cutouts (S105/106 have no HLSP
    light curves yet);
  * batman + emcee limb-darkened fit -> Rp/R*, b, rho*;
  * 300 random ephemerides on the hold-out data -> empirical p of z_hold.

    python round7/rnaas_check.py --raw <dir of HLSP fits> --cut <tesscut dir> --out <json>
"""
from __future__ import annotations

import argparse
import glob
import json
import re

import numpy as np
from astropy.io import fits
from wotan import flatten

P0, T00 = 0.939809448365927, 3155.4570696806513      # BTJD; ephemeris under test (round-7 ephem.py)
DUR = 0.8636701508353961 / 24                          # d, round-7 box duration (re-fit below)
RSTAR, E_RSTAR, MSTAR, TEFF = 0.464827, 0.0137756, 0.462244, 3505.0   # TIC 8.2 (queried live 2026-10-13)
RSUN_REARTH = 109.076
HIS = {2, 29, 69}
NEW = {96, 103, 104}


def load_hlsp(path):
    with fits.open(path) as h:
        d, hdr0 = h[1].data, h[0].header
        sector = int(hdr0["SECTOR"])
        t = np.asarray(d["TIME"], float)
        q = np.asarray(d["QUALITY"])
        if "PDCSAP_FLUX" in d.columns.names:
            f, prov = np.asarray(d["PDCSAP_FLUX"], float), "TESS-SPOC"
        else:
            f, prov = np.asarray(d["DET_FLUX"], float), "QLP"
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
    return sector, prov, t[ok], f[ok] / np.nanmedian(f[ok])


def load_cut(path):
    """3x3 box around the brightest central pixel (pixels > 20% of its flux), after subtracting a
    per-cadence sky = median of pixels >= 4 px from the cutout centre. (S105/106 cutouts carry
    FLUX_BKG = 0, i.e. no background subtraction, so a naive aperture is sky-dominated.)"""
    with fits.open(path) as h:
        d = h[1].data
        t = np.asarray(d["TIME"], float)
        q = np.asarray(d["QUALITY"])
        cube = np.asarray(d["FLUX"], float)
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(cube.reshape(len(t), -1)), 1)
    t, cube = t[ok], cube[ok]
    ny, nx = cube.shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(yy - ny // 2, xx - nx // 2)
    cube = cube - np.nanmedian(cube[:, r >= 4], axis=1)[:, None, None]
    img = np.nanmedian(cube, 0)
    py, px = np.unravel_index(np.argmax(np.where(r <= 1.5, img, -np.inf)), img.shape)
    ap = (np.abs(yy - py) <= 1) & (np.abs(xx - px) <= 1) & (img > 0.2 * img[py, px])
    f = cube[:, ap].sum(1)
    good = f > 0
    sector = int(re.search(r"tess-s(\d{4})", path).group(1))
    return sector, "TESScut", t[good], f[good] / np.median(f[good])


def phase(t, P=P0, T0=T00):
    return ((t - T0) / P + 0.5) % 1 - 0.5


def detrend(t, f, mask_transits=True, window=0.4):
    m = np.abs(phase(t)) * P0 < (1.5 * DUR) if mask_transits else np.zeros(t.size, bool)
    w = np.ones_like(f)
    w[m] = np.nan                      # wotan ignores NaNs when computing the trend
    _, trend = flatten(t, np.where(m, np.nan, f), method="biweight", window_length=window,
                       edge_cutoff=0, break_tolerance=0.5, return_trend=True)
    # interpolate the trend across masked points
    good = np.isfinite(trend)
    trend = np.interp(t, t[good], trend[good])
    r = f / trend - 1
    oot = ~m
    s = 1.4826 * np.nanmedian(np.abs(r[oot] - np.nanmedian(r[oot])))
    keep = np.isfinite(r) & (r < 4 * s) & (r > -10 * s)
    return t[keep], r[keep], s


def box_depth(t, r, P=P0, T0=T00, dur=DUR, sel=None, beta_scale=True):
    """Depth = mean(out) - mean(in); white error from OOT scatter; red-noise beta from
    binning the OOT residuals to the transit duration."""
    ph = phase(t, P, T0) * P
    if sel is not None:
        t, r, ph = t[sel], r[sel], ph[sel]
    inn = np.abs(ph) < dur / 2
    oot = (np.abs(ph) > dur) & (np.abs(ph) < 0.5 * P - dur)
    if inn.sum() < 5:
        return np.nan, np.nan, 0, np.nan
    s = np.std(r[oot])
    depth = np.mean(r[oot]) - np.mean(r[inn])
    err = s * np.sqrt(1 / inn.sum() + 1 / oot.sum())
    beta = 1.0
    if beta_scale:
        cad = np.median(np.diff(np.sort(t)))
        nb = max(int(round(dur / cad)), 1)
        if nb > 1:
            to, ro = t[oot], r[oot]
            idx = np.floor((to - to.min()) / dur).astype(int)
            cnt = np.bincount(idx)
            sums = np.bincount(idx, ro)
            full = cnt >= 0.8 * nb
            if full.sum() > 20:
                bm = sums[full] / cnt[full]
                exp = s / np.sqrt(cnt[full].mean())
                beta = max(1.0, np.std(bm) / exp)
    n_tr = len(np.unique(np.round((t[inn] - T0) / P)))
    return depth, err * beta, n_tr, beta


def tc_fit(t, r, dur=DUR):
    """Per-sector mid-time: chi2 of the best-depth box over a grid of Tc shifts."""
    E = np.round((np.median(t) - T00) / P0)
    tc0 = T00 + E * P0
    shifts = np.linspace(-0.03, 0.03, 241)
    s = np.std(r)
    chi = []
    for dt in shifts:
        ph = phase(t, P0, T00 + dt) * P0
        inn = np.abs(ph) < dur / 2
        d = -np.mean(r[inn])
        model = np.where(inn, -d, 0)
        chi.append(np.sum((r - model) ** 2) / s ** 2)
    chi = np.array(chi)
    i = np.argmin(chi)
    ok = chi <= chi[i] + 1
    lo, hi = shifts[ok].min(), shifts[ok].max()
    return int(E), tc0 + shifts[i], max((hi - lo) / 2, shifts[1] - shifts[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--cut", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mcmc", type=int, default=6000)
    a = ap.parse_args()
    rng = np.random.default_rng(144439305)
    out = {}

    lcs = {}
    for p in sorted(glob.glob(f"{a.raw}/*_lc.fits") + glob.glob(f"{a.raw}/*_llc.fits")):
        s, prov, t, f = load_hlsp(p)
        tt, rr, sig = detrend(t, f)
        lcs[s] = dict(prov=prov, t=tt, r=rr, sig=sig, n=len(tt))
    out["sectors_hlsp"] = {s: dict(prov=v["prov"], n=v["n"], sigma_ppm=v["sig"] * 1e6) for s, v in lcs.items()}

    def cat(keys):
        t = np.concatenate([lcs[k]["t"] for k in keys])
        r = np.concatenate([lcs[k]["r"] for k in keys])
        o = np.argsort(t)
        return t[o], r[o]

    # ---- blind BLS (unmasked detrend so the search is blind) ----
    from astropy.timeseries import BoxLeastSquares
    tb, rb = [], []
    for p in sorted(glob.glob(f"{a.raw}/*_lc.fits") + glob.glob(f"{a.raw}/*_llc.fits")):
        s, prov, t, f = load_hlsp(p)
        tt, rr, _ = detrend(t, f, mask_transits=False)
        tb.append(tt); rb.append(rr)
    for label, keys in (("his", HIS), ("all", HIS | NEW)):
        sel_t = np.concatenate([x for x, p in zip(tb, sorted(lcs)) if p in keys])
        sel_r = np.concatenate([x for x, p in zip(rb, sorted(lcs)) if p in keys])
        # bin to 20 min to keep the grid tractable (transit ~52 min)
        bins = np.floor(sel_t / (20 / 1440)).astype(int)
        u, inv = np.unique(bins, return_inverse=True)
        cnt = np.bincount(inv)
        tbin = np.bincount(inv, sel_t) / cnt
        rbin = np.bincount(inv, sel_r) / cnt
        span = tbin.max() - tbin.min()
        fmin, fmax = 1 / 2.0, 1 / 0.5
        freqs = np.arange(fmin, fmax, 2e-6)          # phase drift over the span << duration at P~1 d
        bls = BoxLeastSquares(tbin, rbin)
        res = bls.power(1 / freqs, [0.025, 0.035, 0.05], objective="snr")
        pw = np.asarray(res.power)
        i = int(np.nanargmax(pw))
        # SDE on the snr periodogram
        sde = (pw[i] - np.nanmean(pw)) / np.nanstd(pw)
        top = np.argsort(pw)[::-1]
        peaks = []
        for j in top:
            Pj = 1 / freqs[j]
            if all(abs(Pj - q) / q > 0.002 for q, _ in peaks):
                peaks.append((Pj, float(pw[j])))
            if len(peaks) == 5:
                break
        out[f"bls_{label}"] = dict(n_freq=len(freqs), best_P=float(1 / freqs[i]), best_snr=float(pw[i]),
                                   sde=float(sde), depth_ppm=float(res.depth[i] * 1e6),
                                   top5=[(round(q, 6), round(v, 2)) for q, v in peaks])

    # ---- depths at the fixed ephemeris ----
    for label, keys in (("his", HIS), ("holdout", NEW), ("all6", HIS | NEW)):
        t, r = cat(keys)
        d, e, ntr, beta = box_depth(t, r)
        out[f"depth_{label}"] = dict(depth_ppm=d * 1e6, err_ppm=e * 1e6, z=d / e, n_tr=ntr, beta=beta)
    for s in sorted(lcs):
        d, e, ntr, beta = box_depth(lcs[s]["t"], lcs[s]["r"])
        out.setdefault("depth_sector", {})[s] = dict(prov=lcs[s]["prov"], depth_ppm=d * 1e6, err_ppm=e * 1e6, n_tr=ntr)
    ta, ra = cat(HIS | NEW)
    # odd / even, phase 0.5
    ep = np.round((ta - T00) / P0).astype(int)
    do, eo, _, _ = box_depth(ta, ra, sel=ep % 2 == 1)
    de, ee, _, _ = box_depth(ta, ra, sel=ep % 2 == 0)
    out["odd_even"] = dict(odd=do * 1e6, even=de * 1e6, sigma=abs(do - de) / np.hypot(eo, ee))
    dh, eh, _, _ = box_depth(ta, ra, T0=T00 + P0 / 2)
    out["phase_half"] = dict(depth_ppm=dh * 1e6, err_ppm=eh * 1e6)

    # hold-out null: 300 random ephemerides in the NEW data
    tn, rn = cat(NEW)
    dhold = out["depth_holdout"]
    zs = []
    while len(zs) < 300:
        Pr = P0 * (1 + rng.uniform(-0.05, 0.05))
        if abs(Pr - P0) / P0 < 0.005:
            continue
        d, e, _, _ = box_depth(tn, rn, P=Pr, T0=T00 + rng.uniform(0, Pr))
        zs.append(d / e)
    zs = np.array(zs)
    out["holdout_null"] = dict(n=300, mean=float(zs.mean()), std=float(zs.std()),
                               max=float(zs.max()), n_ge_obs=int((zs >= dhold["z"]).sum()))

    # ---- per-sector mid-times -> ephemeris ----
    rows = []
    for s in sorted(lcs):
        E, tc, etc = tc_fit(lcs[s]["t"], lcs[s]["r"])
        rows.append((s, E, tc, etc))
    E = np.array([x[1] for x in rows], float)
    tc = np.array([x[2] for x in rows])
    w = 1 / np.array([x[3] for x in rows]) ** 2
    Eref = np.round(np.average(E, weights=w))
    A = np.vstack([np.ones_like(E), E - Eref]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    beta_ = cov @ (A.T @ (w * tc))
    resid = tc - A @ beta_
    chi2 = float(np.sum(w * resid ** 2))
    dof = len(tc) - 2
    sc = max(1.0, chi2 / dof)
    out["ephemeris"] = dict(T0_btjd=float(beta_[0]), T0_bjd=float(beta_[0] + 2457000), T0_err=float(np.sqrt(cov[0, 0] * sc)),
                            Eref=int(Eref), P=float(beta_[1]), P_err=float(np.sqrt(cov[1, 1] * sc)), chi2=chi2, dof=dof,
                            per_sector=[dict(sector=s, epoch=int(e), tc=float(t), err_min=float(x * 1440),
                                             oc_min=float(rr * 1440)) for (s, e, t, x), rr in zip(rows, resid)])
    # propagate to 2027-06-01 (BTJD 4462.5)
    tq = 4462.5
    n = np.round((tq - beta_[0]) / beta_[1])
    v = np.array([1, n])
    out["ephemeris"]["sigma_2027-06-01_min"] = float(np.sqrt(v @ cov @ v * sc) * 1440)
    # compare to the ephemeris under test
    nn = np.round((beta_[0] - T00) / P0)
    out["ephemeris"]["T0_vs_tested_min"] = float((beta_[0] - (T00 + nn * P0)) * 1440)

    # ---- TESScut fixed-aperture photometry, all sectors ----
    cuts = {}
    for p in sorted(glob.glob(f"{a.cut}/*_astrocut.fits")):
        s, prov, t, f = load_cut(p)
        tt, rr, sig = detrend(t, f)
        d, e, ntr, beta = box_depth(tt, rr)
        cuts[s] = dict(t=tt, r=rr)
        out.setdefault("tesscut_sector", {})
        out.setdefault("tesscut_sector", {})[s] = dict(depth_ppm=d * 1e6, err_ppm=e * 1e6, n_tr=ntr, sigma_ppm=sig * 1e6)
    for label, keys in (("S105_106", {105, 106}), ("all8", set(cuts))):
        t = np.concatenate([cuts[k]["t"] for k in keys])
        r = np.concatenate([cuts[k]["r"] for k in keys])
        d, e, ntr, beta = box_depth(t, r)
        out[f"tesscut_{label}"] = dict(depth_ppm=d * 1e6, err_ppm=e * 1e6, z=d / e, n_tr=ntr)
    ds = np.array([out["tesscut_sector"][s]["depth_ppm"] for s in sorted(cuts)])
    es = np.array([out["tesscut_sector"][s]["err_ppm"] for s in sorted(cuts)])
    wm = np.average(ds, weights=1 / es ** 2)
    from scipy.stats import chi2 as C2
    c2 = float(np.sum((ds - wm) ** 2 / es ** 2))
    out["tesscut_sector_consistency"] = dict(mean=wm, err=float(1 / np.sqrt(np.sum(1 / es ** 2))), chi2=c2,
                                            dof=len(ds) - 1, p=float(C2.sf(c2, len(ds) - 1)))

    # ---- limb-darkened fit on the six HLSP sectors (batman + emcee) ----
    import batman
    import emcee
    tf, rf, cadf, sgf = [], [], [], []
    fitsets = [lcs[s] for s in sorted(lcs)] + [cuts[105], cuts[106]]
    for L in fitsets:
        m = np.abs(phase(L["t"])) * P0 < 0.15
        tf.append(L["t"][m]); rf.append(L["r"][m])
        cad = np.median(np.diff(L["t"]))
        cadf.append(np.full(m.sum(), cad))
        oot = np.abs(phase(L["t"][m])) * P0 > 0.05
        sgf.append(np.full(m.sum(), np.std(L["r"][m][oot])))
    tf, rf, cadf, sgf = map(np.concatenate, (tf, rf, cadf, sgf))
    groups = {}
    for c in np.unique(np.round(cadf * 1440, 1)):
        g = np.round(cadf * 1440, 1) == c
        groups[c] = g
    sig = sgf                                                       # per-data-set white noise

    def model(theta):
        dt, k, b, lrho = theta
        rho = 10 ** lrho
        aRs = (6.674e-8 * rho * (P0 * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
        if b >= aRs:
            return None
        inc = np.degrees(np.arccos(b / aRs))
        prm = batman.TransitParams()
        prm.t0, prm.per, prm.rp, prm.a, prm.inc = T00 + dt, P0, k, aRs, inc
        prm.ecc, prm.w, prm.limb_dark, prm.u = 0.0, 90.0, "quadratic", [0.25, 0.35]
        mod = np.empty_like(tf)
        for c, g in groups.items():
            ss = max(1, int(round(c / 2)))
            mm = batman.TransitModel(prm, tf[g], supersample_factor=ss, exp_time=c / 1440)
            mod[g] = mm.light_curve(prm) - 1
        return mod

    def lnp(theta):
        dt, k, b, lrho = theta
        if not (-0.02 < dt < 0.02 and 0.005 < k < 0.08 and 0 <= b < 1 + k and -0.5 < lrho < 1.7):
            return -np.inf
        m = model(theta)
        if m is None:
            return -np.inf
        return -0.5 * np.sum(((rf - m) / sig) ** 2)

    nw = 24
    p0 = np.array([0.0, 0.025, 0.3, np.log10(6.5)]) + 1e-3 * rng.standard_normal((nw, 4))
    p0[:, 2] = np.abs(p0[:, 2])
    sam = emcee.EnsembleSampler(nw, 4, lnp)
    sam.run_mcmc(p0, a.mcmc, progress=False)
    ch = sam.get_chain(discard=a.mcmc // 3, flat=True)
    q = lambda x: [float(v) for v in np.percentile(x, [16, 50, 84])]
    k = ch[:, 1]
    rp = k * RSTAR * RSUN_REARTH
    rp_full = k * (RSTAR + E_RSTAR * rng.standard_normal(k.size)) * RSUN_REARTH
    out["ld_fit"] = dict(k=q(k), b=q(ch[:, 2]), rho_cgs=q(10 ** ch[:, 3]), dt_min=q(ch[:, 0] * 1440),
                         rp_earth_k_only=q(rp), rp_earth_with_Rstar_err=q(rp_full), sigma_ppm=float(np.median(sig) * 1e6), fit_sets='TESS-SPOC S2/29/69 + QLP S96/103/104 + TESScut S105/106',
                         acc=float(np.mean(sam.acceptance_fraction)), n_points=int(tf.size))

    # ---- derived ----
    G, Msun, AU = 6.674e-11, 1.989e30, 1.496e11
    a_m = (G * MSTAR * Msun * (P0 * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    aRs = a_m / (RSTAR * 6.957e8)
    out["derived"] = dict(a_au=a_m / AU, a_over_Rs=aRs, Teq_K_A0=TEFF * np.sqrt(1 / (2 * aRs)), Teq_K_A03=TEFF * np.sqrt(1 / (2 * aRs)) * 0.7 ** 0.25,
                          S_earth=(RSTAR ** 2 * (TEFF / 5772) ** 4) / (a_m / AU) ** 2,
                          rho_tic_cgs=MSTAR * 1.989e33 / (4 / 3 * np.pi * (RSTAR * 6.957e10) ** 3))
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k not in ("sectors_hlsp",)}, indent=1, default=float))


if __name__ == "__main__":
    main()
