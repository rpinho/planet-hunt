"""Round 8A light-curve gate: every TESS product, fixed-ephemeris hold-out test, and the
light-curve false-positive checks, for one P1 candidate.

    .venv/bin/python round8/lcgate.py toi1846 [--mcmc 4000] [--nnull 2000]

Pre-registered (fixed before any new data was loaded, 2026-10-13 19:30 EDT):
  * the ephemeris under test is the DISCOVERY ephemeris in targets.py; nothing is re-fitted
    before the hold-out depth is measured;
  * hold-out data = sectors the discovery search never saw (see targets.py);
  * "reproduced" = hold-out depth z >= 2.5 AND within 3 sigma of the discovery-data depth;
    "refuted" = hold-out depth more than 3 sigma below the discovery-data depth;
    otherwise "not testable / inconclusive" (expected z_hold < 2.5 means the test has no power).

Estimator: per sector and product, the light curve is flattened (wotan biweight, candidate and
known transits masked), then the depth is the weighted linear least-squares amplitude of a box
whose value per cadence is the fraction of the exposure inside the transit (so 2-min, 200-s,
10-min and 30-min data measure the same box depth). Per-sector weights use the robust scatter
times a time-averaging red-noise beta at the transit duration.
"""
from __future__ import annotations

import argparse
import glob
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from planet_hunt import config  # noqa: E402
from planet_hunt import pixelvet as pv  # noqa: E402
from targets import CANDS  # noqa: E402

warnings.filterwarnings("ignore")
OUT = config.RUNS / "round8_p1"
RAW = config.CACHE / "r8" / "raw"
PREF = ["SPOC-2min", "TESS-SPOC", "QLP", "TGLC", "eleanor-lite", "TESScut"]
RSUN_REARTH = 109.076


# ----------------------------------------------------------------------------- loading
def sector_of(fn):
    m = (re.search(r"tess\d{13}-s(\d{4})-", fn) or re.search(r"_s(\d{4})-", fn) or re.search(r"-s(\d{4})_", fn)
         or re.search(r"tess-s(\d{4})-", fn) or re.search(r"-s(\d{4})-cam", fn))
    return int(m.group(1))


def prov_of(fn):
    if fn.startswith("tess") and fn.endswith("-s_lc.fits"):
        return "SPOC-2min"
    if "hlsp_qlp" in fn:
        return "QLP"
    if "hlsp_tess-spoc" in fn:
        return "TESS-SPOC"
    if "hlsp_tglc" in fn:
        return "TGLC"
    if "hlsp_gsfc-eleanor-lite" in fn:
        return "eleanor-lite"
    return None


def read_lc(path):
    from astropy.io import fits
    fn = Path(path).name
    prov = prov_of(fn)
    with fits.open(path, memmap=False) as h:
        d = h[1].data
        cols = d.columns.names
        t = np.asarray(d["TIME"], float)
        if prov in ("SPOC-2min", "TESS-SPOC"):
            f = np.asarray(d["PDCSAP_FLUX"], float)
            q = np.asarray(d["QUALITY"])
        elif prov == "QLP":
            col = "KSPSAP_FLUX" if "KSPSAP_FLUX" in cols else "DET_FLUX"
            f = np.asarray(d[col], float)
            q = np.asarray(d["QUALITY"])
        elif prov == "TGLC":
            col = "cal_psf_flux" if "cal_psf_flux" in cols else [c for c in cols if "flux" in c.lower()][0]
            f = np.asarray(d[col], float)
            q = np.asarray(d["TESS_flags"]) if "TESS_flags" in cols else np.zeros(len(t), int)
            if "TGLC_flags" in cols:
                q = q | np.asarray(d["TGLC_flags"])
        elif prov == "eleanor-lite":
            col = "CORR_FLUX" if "CORR_FLUX" in cols else "RAW_FLUX"
            f = np.asarray(d[col], float)
            q = np.asarray(d["QUALITY"])
        else:
            return None
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0) & (f > 0)
    if ok.sum() < 50:
        return None
    t, f = t[ok], f[ok] / np.nanmedian(f[ok])
    cad = float(np.nanmedian(np.diff(t)))
    return {"sector": sector_of(fn), "prov": prov, "t": t, "f": f, "cad": cad, "file": fn}


def read_cut(path, ra, dec):
    """Fixed aperture on a TESScut cutout: 3x3 box centred on the target's WCS pixel (pixels above
    20% of the brightest), minus a per-cadence sky = median of pixels >= 4 px from the target
    (S105+ cutouts carry no background subtraction)."""
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(path, memmap=False) as h:
        d = h[1].data
        t = np.asarray(d["TIME"], float)
        q = np.asarray(d["QUALITY"])
        cube = np.asarray(d["FLUX"], float)
        w = WCS(h[2].header)
    x, y = w.world_to_pixel_values(ra, dec)
    ny, nx = cube.shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(yy - float(y), xx - float(x))
    core = r <= 2.2                          # the aperture can only use these pixels
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(cube[:, core]), 1)
    t, cube = t[ok], cube[ok]
    if len(t) < 50:
        return None
    cube = cube - np.nanmedian(cube[:, r >= 4], axis=1)[:, None, None]
    img = np.nanmedian(cube, 0)
    py, px = int(round(float(y))), int(round(float(x)))
    ap = (np.abs(yy - py) <= 1) & (np.abs(xx - px) <= 1) & (img > 0.2 * img[py, px]) & np.isfinite(img)
    f = cube[:, ap].sum(1)
    good = f > 0
    if good.sum() < 50:
        return None
    t, f = t[good], f[good] / np.median(f[good])
    return {"sector": sector_of(Path(path).name), "prov": "TESScut", "t": t, "f": f,
            "cad": float(np.median(np.diff(t))), "file": Path(path).name}


def load_all(key, c, inv):
    lcs = []
    files = sorted(glob.glob(str(RAW / str(c["tic"]) / "*.fits")))
    best = {}
    for p in files:
        fn = Path(p).name
        pr = prov_of(fn)
        if pr is None:
            continue
        s = sector_of(fn)
        best.setdefault((pr, s), []).append(p)
    for (pr, s), ps in sorted(best.items()):
        p = sorted(ps)[-1]                      # latest version (v02 > v01)
        try:
            d = read_lc(p)
        except Exception as exc:
            print("  read fail", Path(p).name, repr(exc)[:80])
            d = None
        if d is not None:
            lcs.append(d)
    cdir = pv.cut_dir(inv["ra"], inv["dec"])
    for p in sorted(cdir.glob("*_11x11_astrocut.fits")):
        try:
            d = read_cut(p, inv["ra"], inv["dec"])
        except Exception as exc:
            print("  cut fail", p.name, repr(exc)[:80])
            d = None
        if d is not None:
            lcs.append(d)
    return lcs


# ----------------------------------------------------------------------------- detrend / measure
def phase(t, P, t0):
    return ((t - t0) / P + 0.5) % 1 - 0.5


def in_window(t, P, t0, half):
    return np.abs(phase(t, P, t0) * P) < half


def overlap(t, cad, P, t0, dur):
    """Fraction of each exposure [t-cad/2, t+cad/2] inside the box transit."""
    x = phase(t, P, t0) * P
    lo = np.maximum(x - cad / 2, -dur / 2)
    hi = np.minimum(x + cad / 2, dur / 2)
    return np.clip(hi - lo, 0, None) / cad


def flatten(d, mask, window):
    from wotan import flatten as wf
    f = np.where(mask, np.nan, d["f"])
    _, trend = wf(d["t"], f, method="biweight", window_length=window, edge_cutoff=0.0,
                  break_tolerance=0.5, return_trend=True)
    good = np.isfinite(trend)
    if good.sum() < 20:
        return None
    trend = np.interp(d["t"], d["t"][good], trend[good])
    r = d["f"] / trend - 1
    oot = ~mask & np.isfinite(r)
    s = 1.4826 * np.nanmedian(np.abs(r[oot] - np.nanmedian(r[oot])))
    keep = np.isfinite(r) & (r < 4 * s) & ((r > -10 * s) | mask)
    return {**d, "t": d["t"][keep], "r": r[keep], "sig": s, "m": mask[keep]}


def beta_at(t, r, dur, cad):
    """Time-averaging red-noise factor at the transit duration (>= 1)."""
    nb = max(int(round(dur / cad)), 1)
    if nb < 2:
        return 1.0
    idx = np.floor((t - t.min()) / dur).astype(int)
    cnt = np.bincount(idx)
    sums = np.bincount(idx, r)
    full = cnt >= 0.7 * nb
    if full.sum() < 15:
        return 1.0
    bm = sums[full] / cnt[full]
    exp = np.std(r) / np.sqrt(cnt[full].mean())
    return float(max(1.0, np.std(bm) / exp))


def depth_fit(t, r, e, cad, P, t0, dur, sel=None):
    """Weighted LSQ box amplitude: r = a - depth * overlap. Returns depth, err, n_epochs."""
    if sel is not None:
        t, r, e, cad = t[sel], r[sel], e[sel], cad[sel]
    o = overlap(t, cad, P, t0, dur)
    near = np.abs(phase(t, P, t0) * P) < max(4 * dur, dur + 2 * np.median(cad) if len(cad) else 0)
    if (o > 0.25).sum() < 3 or near.sum() < 20:
        return np.nan, np.nan, 0
    t, r, e, o = t[near], r[near], e[near], o[near]
    w = 1 / e ** 2
    X = np.vstack([np.ones_like(o), -o]).T
    A = X.T @ (X * w[:, None])
    try:
        C = np.linalg.inv(A)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, 0
    b = C @ (X.T @ (w * r))
    n_ep = len(np.unique(np.round((t[o > 0.25] - t0) / P)))
    return float(b[1]), float(np.sqrt(C[1, 1])), int(n_ep)


class Data:
    """Concatenated flattened light curves of one product choice."""

    def __init__(self, parts):
        self.parts = parts
        if not parts:
            self.t = np.array([]); self.r = self.e = self.cad = self.t; self.s = np.array([], int)
            self.prov = np.array([], object); return
        self.t = np.concatenate([p["t"] for p in parts])
        self.r = np.concatenate([p["r"] for p in parts])
        self.e = np.concatenate([np.full(len(p["t"]), p["sig"] * p["beta"]) for p in parts])
        self.cad = np.concatenate([np.full(len(p["t"]), p["cad"]) for p in parts])
        self.s = np.concatenate([np.full(len(p["t"]), p["sector"]) for p in parts])
        self.prov = np.concatenate([np.full(len(p["t"]), p["prov"], dtype=object) for p in parts])
        self.kany = np.concatenate([p["kany"] for p in parts])
        self.kctrl = np.concatenate([p["kctrl"] for p in parts])
        # candidate measurements never use cadences inside known transits
        k = ~self.kany
        self.ctrl = (self.t[~(self.kany & ~self.kctrl)], self.r[~(self.kany & ~self.kctrl)],
                     self.e[~(self.kany & ~self.kctrl)], self.cad[~(self.kany & ~self.kctrl)])
        self.t, self.r, self.e, self.cad, self.s, self.prov = (self.t[k], self.r[k], self.e[k], self.cad[k],
                                                               self.s[k], self.prov[k])

    def depth(self, P, t0, dur, sel=None):
        if len(self.t) == 0:
            return np.nan, np.nan, 0
        return depth_fit(self.t, self.r, self.e, self.cad, P, t0, dur, sel)

    def ctrl_depth(self, P, t0, dur):
        t, r, e, cad = self.ctrl
        if len(t) == 0:
            return np.nan, np.nan, 0
        return depth_fit(t, r, e, cad, P, t0, dur)


def D(d, e, n=None):
    out = {"depth_ppm": d * 1e6, "err_ppm": e * 1e6, "z": d / e if e and np.isfinite(e) and e > 0 else np.nan}
    if n is not None:
        out["n_tr"] = n
    return out


# ----------------------------------------------------------------------------- tests
def null_ephem(dat, P, t0, dur, n, rng, exclude=0.005, span=0.05):
    z = []
    for _ in range(n):
        while True:
            dp = rng.uniform(-span, span)
            if abs(dp) > exclude:
                break
        Pr = P * (1 + dp)
        d, e, _ = dat.depth(Pr, rng.uniform(0, Pr) + t0, dur)
        z.append(d / e if e else np.nan)
    z = np.array(z)
    return z[np.isfinite(z)]


def secondary_scan(dat, P, t0, dur):
    phs = np.arange(0, 1, max(dur / P / 3, 0.002))
    zs = []
    for ph in phs:
        if min(ph, 1 - ph) * P < 1.5 * dur:
            continue
        d, e, _ = dat.depth(P, t0 + ph * P, dur)
        zs.append((ph, d / e if e else np.nan, d))
    a = np.array([(p, z, d) for p, z, d in zs if np.isfinite(z)])
    if not len(a):
        return {}
    i = int(np.argmax(a[:, 1]))
    j = int(np.argmin(a[:, 1]))
    return {"n_phases": int(len(a)), "max_dip_z": float(a[i, 1]), "max_dip_phase": float(a[i, 0]),
            "max_dip_ppm": float(a[i, 2] * 1e6), "min_z": float(a[j, 1]), "scan_std": float(np.std(a[:, 1])),
            "scan_mean": float(np.mean(a[:, 1]))}


def period_scan(dat, P, t0, dur, n=241):
    """SNR at fixed duration over a fine grid of periods (T0 re-optimised per period by
    phase-shifting over one duration) around the candidate period."""
    ncyc = max(np.ptp(dat.t) / P, 1)
    width = dur / (P * ncyc)                 # relative width of the coherent peak
    frac = 40 * width
    Ps = P * (1 + np.linspace(-frac, frac, n))
    tref = np.median(dat.t)
    E = np.round((tref - t0) / P)
    tm = t0 + E * P                          # keep the mid-baseline epoch fixed
    out = []
    for Pp in Ps:
        best = -9
        for dt in np.linspace(-dur, dur, 9):
            d, e, _ = dat.depth(Pp, tm + dt, dur)
            if e and np.isfinite(d / e):
                best = max(best, d / e)
        out.append(best)
    out = np.array(out)
    i = int(np.argmax(out))
    return {"P_best": float(Ps[i]), "z_best": float(out[i]), "z_at_P": float(out[n // 2]),
            "frac_offset": float(Ps[i] / P - 1), "grid": [float(Ps[0]), float(Ps[-1]), n],
            "peak_width_rel": float(width), "offset_in_widths": float((Ps[i] / P - 1) / width),
            "second_peak_z": float(np.max(np.where(np.abs(Ps / Ps[i] - 1) > 3 * width, out, -9))),
            "curve": out.round(2).tolist()}


def per_epoch(dat, P, t0, dur):
    ep = np.round((dat.t - t0) / P).astype(int)
    rows = []
    for k in np.unique(ep):
        sel = ep == k
        o = overlap(dat.t[sel], dat.cad[sel], P, t0, dur)
        if (o > 0.5).sum() < 2:
            continue
        d, e, _ = dat.depth(P, t0, dur, sel)
        if np.isfinite(d):
            rows.append((int(k), float(d), float(e), int(dat.s[sel][0])))
    return rows


def ephemeris(dat, P, t0, dur, groups):
    """Mid-time per group (TESS year) at fixed P and shape: chi^2 of the box LSQ over a grid of
    shifts; 1-sigma from delta-chi2 = 1. Then a weighted linear fit Tc = T0 + E P."""
    rows = []
    for g in np.unique(groups):
        sel = groups == g
        shifts = np.linspace(-1.2 * dur, 1.2 * dur, 121)
        chi = []
        tt, rr, ee, cc = dat.t[sel], dat.r[sel], dat.e[sel], dat.cad[sel]
        near = np.abs(phase(tt, P, t0) * P) < 4 * dur
        tt, rr, ee, cc = tt[near], rr[near], ee[near], cc[near]
        if len(tt) < 30:
            continue
        for dt in shifts:
            o = overlap(tt, cc, P, t0 + dt, dur)
            w = 1 / ee ** 2
            X = np.vstack([np.ones_like(o), -o]).T
            A = X.T @ (X * w[:, None])
            b = np.linalg.solve(A, X.T @ (w * rr))
            res = rr - X @ b
            chi.append(np.sum(w * res ** 2))
        chi = np.array(chi)
        i = int(np.argmin(chi))
        if chi.max() - chi[i] < 4:          # no constraint from this group
            continue
        # box chi^2 surfaces are jagged at the cadence scale: fit a parabola over the region within
        # delta-chi2 < 9 of the minimum (>= 7 points) and take sigma = 1/sqrt(curvature)
        reg = np.flatnonzero(chi <= chi[i] + 9)
        reg = np.arange(max(0, min(reg.min(), i - 3)), min(len(chi), max(reg.max(), i + 3) + 1))
        cf_ = np.polyfit(shifts[reg], chi[reg], 2)
        if cf_[0] <= 0:
            continue
        x0 = -cf_[1] / (2 * cf_[0])
        sig = 1 / np.sqrt(cf_[0])
        if abs(x0) > 1.2 * dur:
            x0 = shifts[i]
        E = int(np.round((np.median(tt) - t0) / P))
        rows.append((int(g), E, t0 + E * P + x0, max(sig, shifts[1] - shifts[0])))
    if len(rows) < 2:
        return None
    E = np.array([r[1] for r in rows], float); tc = np.array([r[2] for r in rows]); w = 1 / np.array([r[3] for r in rows]) ** 2
    Eref = np.round(np.average(E, weights=w))
    A = np.vstack([np.ones_like(E), E - Eref]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    b = cov @ (A.T @ (w * tc))
    resid = tc - A @ b
    chi2 = float(np.sum(w * resid ** 2)); dof = len(tc) - 2
    sc = max(1.0, chi2 / dof) if dof > 0 else 1.0
    return {"T0_btjd": float(b[0]), "T0_bjd": float(b[0] + 2457000), "T0_err": float(np.sqrt(cov[0, 0] * sc)),
            "P": float(b[1]), "P_err": float(np.sqrt(cov[1, 1] * sc)), "chi2": chi2, "dof": dof,
            "P_disc": P, "P_shift_sigma": float((b[1] - P) / np.sqrt(cov[1, 1] * sc)),
            "groups": [dict(group=r[0], epoch=r[1], tc=r[2], err_min=r[3] * 1440, oc_min=float(x * 1440))
                       for r, x in zip(rows, resid)]}


# ----------------------------------------------------------------------------- shape fit
def ld_coeffs(teff):
    if teff < 3400:
        return [0.20, 0.40]
    if teff < 4000:
        return [0.25, 0.35]
    if teff < 5000:
        return [0.35, 0.25]
    return [0.38, 0.22]


def shape_fit(dat, P, t0, dur, c, nsteps, rng, rho_prior=None):
    import batman
    import emcee
    near = np.abs(phase(dat.t, P, t0) * P) < max(3 * dur, 0.1)
    t, r, e, cad = dat.t[near], dat.r[near], dat.e[near], dat.cad[near]
    groups = {cc: np.round(cad * 1440, 1) == cc for cc in np.unique(np.round(cad * 1440, 1))}
    u = ld_coeffs(c["teff"])
    k0 = np.sqrt(max(c["depth_ppm"], 50) * 1e-6)

    def model(theta):
        dt, k, b, lrho = theta
        rho = 10 ** lrho
        aRs = (6.674e-8 * rho * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
        if b >= aRs:
            return None
        prm = batman.TransitParams()
        prm.t0, prm.per, prm.rp, prm.a, prm.inc = t0 + dt, P, k, aRs, np.degrees(np.arccos(b / aRs))
        prm.ecc, prm.w, prm.limb_dark, prm.u = 0.0, 90.0, "quadratic", u
        m = np.empty_like(t)
        for cc, g in groups.items():
            ss = max(1, int(round(cc / 2)))
            m[g] = batman.TransitModel(prm, t[g], supersample_factor=ss, exp_time=cc / 1440).light_curve(prm) - 1
        return m

    def lnp(theta):
        dt, k, b, lrho = theta
        if not (-0.6 * dur < dt < 0.6 * dur and 0.002 < k < 0.3 and 0 <= b < 1 + k and -1.5 < lrho < 2.3):
            return -np.inf
        m = model(theta)
        if m is None:
            return -np.inf
        lp = -0.5 * np.sum(((r - m) / e) ** 2)
        if rho_prior is not None:
            lp += -0.5 * ((10 ** lrho - rho_prior[0]) / rho_prior[1]) ** 2
        return lp

    nw = 24
    rho0 = rho_prior[0] if rho_prior else c["mstar"] / c["rstar"] ** 3 * 1.41
    p0 = np.array([0.0, k0, 0.3, np.log10(rho0)]) + np.array([1e-4, 1e-3, 0.05, 0.02]) * rng.standard_normal((nw, 4))
    p0[:, 2] = np.abs(p0[:, 2])
    sam = emcee.EnsembleSampler(nw, 4, lnp)
    sam.run_mcmc(p0, nsteps, progress=False)
    ch = sam.get_chain(discard=nsteps // 3, flat=True)
    q = lambda x: [float(v) for v in np.percentile(x, [16, 50, 84])]
    k = ch[:, 1]
    er = c.get("e_rstar") or 0.03 * c["rstar"]
    rp = k * (c["rstar"] + er * rng.standard_normal(k.size)) * RSUN_REARTH
    rho = 10 ** ch[:, 3]
    aRs = (6.674e-8 * rho * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    b = ch[:, 2]
    T14 = P / np.pi * np.arcsin(np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / aRs) * 24
    return {"k": q(k), "b": q(b), "rho_cgs": q(rho), "dt_min": q(ch[:, 0] * 1440), "rp_earth": q(rp),
            "T14_h": q(T14), "depth_ppm_k2": q(k ** 2 * 1e6), "acc": float(np.mean(sam.acceptance_fraction)),
            "n_points": int(t.size), "u": u, "rho_prior": rho_prior}


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("key")
    ap.add_argument("--mcmc", type=int, default=4000)
    ap.add_argument("--nnull", type=int, default=2000)
    ap.add_argument("--no-fit", action="store_true")
    a = ap.parse_args()
    c = dict(CANDS[a.key])
    inv = json.load(open(OUT / "inventory.json"))[a.key]
    tic_row = inv["tic_row"]
    if tic_row.get("e_rad"):
        c["e_rstar"] = tic_row["e_rad"]
    rng = np.random.default_rng(c["tic"] % 2 ** 32)
    P, t0, dur = c["P"], c["t0"], c["dur_h"] / 24
    # known transits (TOI table + confirmed planets + hand list), masked everywhere
    import round5_hosts as r5
    toi, ps = r5.load_catalogs()
    known = []
    for k in r5.known_ephemerides(c["tic"], toi, ps):
        if np.isfinite(k["period"]) and np.isfinite(k["t0"]) and abs(k["period"] / P - 1) > 0.002:
            known.append(dict(name=k["name"], P=k["period"], t0=k["t0"],
                              dur_h=k["duration_h"] if np.isfinite(k["duration_h"]) else 2.0))
    for k in c["known"]:
        if not any(abs(k["P"] / q["P"] - 1) < 0.002 for q in known):
            known.append(k)
    print("known masked:", [(k["name"], round(k["P"], 5)) for k in known], flush=True)

    def kmask(t):
        m = np.zeros(len(t), bool)
        for k in known:
            m |= in_window(t, k["P"], k["t0"], 0.75 * k["dur_h"] / 24 + 0.02)
        return m

    def kmask_one(t, k):
        return in_window(t, k["P"], k["t0"], 0.75 * k["dur_h"] / 24 + 0.02)

    lcs = load_all(a.key, c, inv)
    window = 0.5 if P < 6 else 1.0
    window = max(window, 3.5 * dur)
    flat = []
    for d in lcs:
        m = in_window(d["t"], P, t0, 0.75 * dur + d["cad"] / 2) | kmask(d["t"])
        fd = flatten(d, m, window)
        if fd is None or len(fd["t"]) < 50 or not np.isfinite(fd["sig"]) or fd["sig"] > 0.02:
            continue                              # broken product (e.g. an empty/edge TESScut sector)
        fd["kany"] = kmask(fd["t"])
        fd["kctrl"] = kmask_one(fd["t"], known[0]) if known else np.zeros(len(fd["t"]), bool)
        oot = ~in_window(fd["t"], P, t0, 0.75 * dur + d["cad"] / 2) & ~fd["kany"]
        fd["beta"] = beta_at(fd["t"][oot], fd["r"][oot], dur, fd["cad"])
        flat.append(fd)
    # products per sector
    have = {}
    for d in flat:
        have.setdefault(d["sector"], {})[d["prov"]] = d
    sectors = sorted(have)
    disc_secs = set(c.get("disc_sectors") or [])
    if c["disc"] == "r5hosts":
        disc_secs = {s for s in sectors if "SPOC-2min" in have[s] and s <= c["disc_max_sector"]}
    elif c["disc"] == "r5faint":
        disc_secs = {int(p.stem[1:]) for p in (config.CACHE / "ffi" / str(c["tic"])).glob("s*.npz")}
    primary, info = [], []
    for s in sectors:
        provs = have[s]
        if s in disc_secs and c["disc_product"] in provs:
            pick = c["disc_product"]
        else:
            pick = next(p for p in PREF if p in provs)
        d = provs[pick]
        d["role"] = "disc" if s in disc_secs else "hold"
        primary.append(d)
        info.append({"sector": s, "primary": pick, "role": d["role"], "products": sorted(provs),
                     "cad_min": round(d["cad"] * 1440, 2), "sigma_ppm": round(d["sig"] * 1e6), "beta": round(d["beta"], 2),
                     "n": int(len(d["t"]))})
    ALL = Data(primary)
    DISC = Data([d for d in primary if d["role"] == "disc"])
    HOLD = Data([d for d in primary if d["role"] == "hold"])
    res = {"key": a.key, "cand": {k: v for k, v in c.items() if k != "known"}, "known_masked": known,
           "sectors": info, "n_disc": len(disc_secs & set(sectors)), "n_hold": int(sum(d["role"] == "hold" for d in primary)),
           "window_d": window}
    print(f"{len(sectors)} sectors; disc {res['n_disc']}, hold-out {res['n_hold']}: "
          f"{[i['sector'] for i in info if i['role'] == 'hold']}", flush=True)

    # 1. hold-out test at the fixed discovery ephemeris
    dd, de, dn = DISC.depth(P, t0, dur)
    hd, he, hn = HOLD.depth(P, t0, dur)
    ad, ae, an = ALL.depth(P, t0, dur)
    res["disc"] = D(dd, de, dn)
    res["hold"] = D(hd, he, hn)
    res["all"] = D(ad, ae, an)
    if np.isfinite(hd):
        res["hold"]["consistency_sigma"] = float((hd - dd) / np.hypot(he, de))
        res["hold"]["expected_z_if_real"] = float(dd / he)
        zexp = dd / he
        if hd / he >= 2.5 and abs(res["hold"]["consistency_sigma"]) < 3:
            v = "REPRODUCED"
        elif res["hold"]["consistency_sigma"] < -3:
            v = "REFUTED"
        elif zexp < 2.5:
            v = "NO POWER (expected z < 2.5)"
        else:
            v = "INCONCLUSIVE"
        res["hold"]["verdict"] = v
    else:
        res["hold"]["verdict"] = "NO HOLD-OUT DATA"
    print("disc", res["disc"], "\nhold", res["hold"], "\nall", res["all"], flush=True)

    # 2. per sector, per product (all products, every sector)
    rows = []
    for s in sectors:
        for pr, d in have[s].items():
            x = Data([d])
            dp, ep, n = x.depth(P, t0, dur)
            rows.append({"sector": s, "prov": pr, "role": "disc" if s in disc_secs else "hold", **D(dp, ep, n)})
    res["per_sector_product"] = rows
    # product comparison on common sectors (vs SPOC-2min or TESS-SPOC), with the known planet as control
    ctrl = known[0] if known else None
    comp = {}
    base = "SPOC-2min" if any("SPOC-2min" in have[s] for s in sectors) else c["disc_product"]
    for pr in PREF:
        if pr == base:
            continue
        com = [s for s in sectors if pr in have[s] and base in have[s]]
        if len(com) < 2:
            continue
        A_ = Data([have[s][pr] for s in com]); B_ = Data([have[s][base] for s in com])
        da, ea, _ = A_.depth(P, t0, dur); db, eb, _ = B_.depth(P, t0, dur)
        entry = {"sectors": com, "cand_prod_ppm": da * 1e6, "cand_prod_err": ea * 1e6, "cand_base_ppm": db * 1e6,
                 "cand_base_err": eb * 1e6, "ratio": da / db if db else np.nan,
                 "ratio_err": abs(da / db) * np.hypot(ea / da, eb / db) if da and db else np.nan}
        if ctrl is not None:
            # control: the known planet's own transits in the same two products and sectors
            cP, ct0, cdur = ctrl["P"], ctrl["t0"], ctrl["dur_h"] / 24
            ka, kea, _ = A_.ctrl_depth(cP, ct0, cdur); kb, keb, _ = B_.ctrl_depth(cP, ct0, cdur)
            entry["control"] = ctrl["name"]
            entry["ctrl_ratio"] = ka / kb if kb else np.nan
            entry["ctrl_ratio_err"] = abs(ka / kb) * np.hypot(kea / ka, keb / kb) if ka and kb else np.nan
            entry["cand_vs_ctrl_sigma"] = float((entry["ratio"] - entry["ctrl_ratio"]) /
                                                np.hypot(entry["ratio_err"], entry["ctrl_ratio_err"]))
        comp[pr] = entry
    res["product_compare_base"] = base
    res["product_compare"] = comp

    # 3. odd/even, secondary, aliases (all primary data)
    ep = np.round((ALL.t - t0) / P).astype(int)
    do_, eo_, _ = ALL.depth(P, t0, dur, ep % 2 == 1)
    de_, ee_, _ = ALL.depth(P, t0, dur, ep % 2 == 0)
    res["odd_even"] = {"odd": D(do_, eo_), "even": D(de_, ee_), "sigma": float(abs(do_ - de_) / np.hypot(eo_, ee_))}
    h2, h2e, _ = ALL.depth(P, t0 + P / 2, dur)
    res["phase_half"] = D(h2, h2e)
    res["alias_P_over_2"] = {"phase0.5_over_primary": float(h2 / ad), "err": float(h2e / ad),
                             "excluded_sigma": float((ad - h2) / np.hypot(ae, h2e))}
    p13, e13, _ = ALL.depth(P, t0 + P / 3, dur)
    p23, e23, _ = ALL.depth(P, t0 + 2 * P / 3, dur)
    res["alias_P_over_3"] = {"phase1/3": D(p13, e13), "phase2/3": D(p23, e23),
                             "excluded_sigma": float((ad - 0.5 * (p13 + p23)) / np.hypot(ae, 0.5 * np.hypot(e13, e23)))}
    res3 = []
    for k in range(3):
        d3, e3, _ = ALL.depth(P, t0, dur, ep % 3 == k)
        res3.append(D(d3, e3))
    res["alias_3P_residues"] = res3
    res["alias_2P"] = "see odd_even (both parities carry the transit if P is right)"
    res["secondary_scan"] = secondary_scan(ALL, P, t0, dur)
    # nearby periods: is the fixed P the local optimum? (all data)
    ps_ = period_scan(ALL, P, t0, dur)
    res["period_scan"] = {k: v for k, v in ps_.items() if k != "curve"}
    np.save(OUT / f"{a.key}_period_scan.npy", np.array(ps_["curve"]))

    # 4. per-epoch / per-sector / per-season consistency
    from scipy.stats import chi2 as C2
    pe = per_epoch(ALL, P, t0, dur)
    if pe:
        d_ = np.array([x[1] for x in pe]); e_ = np.array([x[2] for x in pe])
        wm = np.sum(d_ / e_ ** 2) / np.sum(1 / e_ ** 2)
        c2 = float(np.sum((d_ - wm) ** 2 / e_ ** 2))
        loo = []
        for i in range(len(pe)):
            m = np.ones(len(pe), bool); m[i] = False
            w_ = 1 / e_[m] ** 2
            loo.append(np.sum(d_[m] * w_) / np.sqrt(np.sum(w_)))
        top = np.argsort(d_ / e_)[::-1]
        drop = max(1, int(round(0.05 * len(pe))))
        m = np.ones(len(pe), bool); m[top[:drop]] = False
        z_drop = float(np.sum(d_[m] / e_[m] ** 2) / np.sqrt(np.sum(1 / e_[m] ** 2)))
        res["per_epoch"] = {"n": len(pe), "chi2": c2, "dof": len(pe) - 1, "p": float(C2.sf(c2, len(pe) - 1)),
                            "frac_positive": float(np.mean(d_ > 0)), "min_loo_z": float(np.min(loo)),
                            "z_drop_top5pct": z_drop, "n_dropped": drop,
                            "max_single_z": float(np.max(d_ / e_)),
                            "epochs": [dict(epoch=x[0], depth_ppm=x[1] * 1e6, err_ppm=x[2] * 1e6, sector=x[3]) for x in pe]
                            if len(pe) <= 60 else "omitted (>60)"}
    ps_rows = [r for r in rows if any(i["sector"] == r["sector"] and i["primary"] == r["prov"] for i in info)]
    d_ = np.array([r["depth_ppm"] for r in ps_rows]); e_ = np.array([r["err_ppm"] for r in ps_rows])
    ok = np.isfinite(d_) & np.isfinite(e_) & (e_ > 0)
    wm = np.sum(d_[ok] / e_[ok] ** 2) / np.sum(1 / e_[ok] ** 2)
    c2 = float(np.sum((d_[ok] - wm) ** 2 / e_[ok] ** 2))
    res["per_sector_chi2"] = {"chi2": c2, "dof": int(ok.sum() - 1), "p": float(C2.sf(c2, ok.sum() - 1)),
                              "n_positive": int((d_[ok] > 0).sum()), "n": int(ok.sum())}
    # seasons = TESS years (~13 sectors)
    yr = (ALL.s - 1) // 13 + 1
    seas = []
    for y in np.unique(yr):
        dy, ey, ny = ALL.depth(P, t0, dur, yr == y)
        if np.isfinite(dy):
            seas.append({"year": int(y), **D(dy, ey, ny)})
    d_ = np.array([x["depth_ppm"] for x in seas]); e_ = np.array([x["err_ppm"] for x in seas])
    wm = np.sum(d_ / e_ ** 2) / np.sum(1 / e_ ** 2)
    c2 = float(np.sum((d_ - wm) ** 2 / e_ ** 2))
    res["seasons"] = {"rows": seas, "chi2": c2, "dof": len(seas) - 1, "p": float(C2.sf(c2, max(len(seas) - 1, 1)))}
    half = ALL.t < np.median(ALL.t)
    h1 = ALL.depth(P, t0, dur, half); h2_ = ALL.depth(P, t0, dur, ~half)
    res["halves"] = {"first": D(*h1), "second": D(*h2_)}

    # 5. random-ephemeris nulls (calibrate the errors; NOT a look-elsewhere correction)
    zf = null_ephem(ALL, P, t0, dur, a.nnull, rng)
    res["null_all"] = {"n": int(zf.size), "mean": float(zf.mean()), "std": float(zf.std()),
                       "p_emp": float(np.mean(zf >= res["all"]["z"])), "max": float(zf.max()),
                       "z_obs_over_null_std": float((res["all"]["z"] - zf.mean()) / zf.std())}
    if res["n_hold"] and np.isfinite(hd):
        zh = null_ephem(HOLD, P, t0, dur, a.nnull, rng)
        res["null_hold"] = {"n": int(zh.size), "mean": float(zh.mean()), "std": float(zh.std()),
                            "p_emp": float(np.mean(zh >= res["hold"]["z"])), "max": float(zh.max())}
    print("odd/even", round(res["odd_even"]["sigma"], 2), "phase0.5", res["phase_half"],
          "\nsec scan", res["secondary_scan"], "\nperiod scan", res["period_scan"],
          "\nper-sector", res["per_sector_chi2"], "\nseasons p", round(res["seasons"]["p"], 3),
          "\nnull_all", res["null_all"], "\nnull_hold", res.get("null_hold"), flush=True)

    # 5b. ephemeris from per-year mid-times (for ground-data predictions; after the hold-out test)
    eph = ephemeris(ALL, P, t0, dur, (ALL.s - 1) // 13 + 1)
    res["ephemeris"] = eph
    if eph:
        json.dump(eph, open(OUT / f"{a.key}_ephem.json", "w"), indent=1)
        print("ephemeris", {k: v for k, v in eph.items() if k != "groups"}, flush=True)

    # 6. shape fits
    if not a.no_fit:
        rho_tic = c["mstar"] / c["rstar"] ** 3 * 1.408
        e_rho = rho_tic * 0.25
        res["fit_free"] = shape_fit(ALL, P, t0, dur, c, a.mcmc, rng)
        res["fit_rhoprior"] = shape_fit(ALL, P, t0, dur, c, a.mcmc, rng, rho_prior=(rho_tic, e_rho))
        res["rho_tic"] = rho_tic
        print("fit free", res["fit_free"]["rp_earth"], res["fit_free"]["rho_cgs"], res["fit_free"]["b"],
              "\nfit rho prior", res["fit_rhoprior"]["rp_earth"], res["fit_rhoprior"]["b"], res["fit_rhoprior"]["T14_h"],
              flush=True)
    # fold for plots / TRICERATOPS
    near = np.abs(phase(ALL.t, P, t0) * P) < max(4 * dur, 0.15)
    np.savez(OUT / f"{a.key}_fold.npz", phase_d=phase(ALL.t[near], P, t0) * P, flux=1 + ALL.r[near],
             err=ALL.e[near], cad=ALL.cad[near], sector=ALL.s[near],
             role=np.array([1 if s in disc_secs else 0 for s in ALL.s[near]]), period=P,
             depth=max(ad, 1e-5), sectors=np.array(sectors))

    def js(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o) if np.isfinite(o) else None
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.bool_):
            return bool(o)
        return str(o)
    (OUT / f"{a.key}_lcgate.json").write_text(json.dumps(res, indent=1, default=js))
    print("wrote", OUT / f"{a.key}_lcgate.json")


if __name__ == "__main__":
    main()
