"""Pixel tests from SPOC 2-min target pixel files.
(1) Per-sector difference images for b, c, d18, d2 with an empirical noise image from fake epochs;
    Gaussian-PSF centroid of the difference image vs the out-of-transit image -> sky offset.
(2) Depth versus aperture size (1 px .. 7x7 + SPOC aperture), relative to planet c / b."""
from common import *
import json, glob
from astropy.wcs import WCS
from scipy.optimize import least_squares
E = json.load(open(OUT / 'data' / 'ephem.json'))
PLS = ['b', 'c', 'd18', 'd2']
rng = np.random.default_rng(42)

def load_tpf(s):
    fn = glob.glob(str(RAW / f'*-s{s:04d}-*-s_tp.fits'))[0]
    with fits.open(fn) as h:
        d = h[1].data
        t = np.array(d['TIME']); F = np.array(d['FLUX'], float); q = np.array(d['QUALITY'])
        ap = np.array(h[2].data)
        w = WCS(h[2].header)
        hdr0 = dict(h[0].header)
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1)
    return dict(t=t[ok], F=F[ok], ap=ap, wcs=w, ra=hdr0['RA_OBJ'], dec=hdr0['DEC_OBJ'])

def all_transit_mask(t, widen=1.5):
    return mask_all(t, E, widen=widen)

def diff_one(tp, tc, dur, nmin=15):
    t, F = tp['t'], tp['F']
    w = np.abs(t - tc) < 0.75 * dur + 0.3
    tin = np.abs(t - tc) < 0.35 * dur
    oot = w & (np.abs(t - tc) > 0.75 * dur)
    # other transits inside the window are excluded from both
    others = np.zeros(len(t), bool)
    for o in PLS:
        e = E[o]
        if abs(((tc - e['t0'] + 0.5 * e['P']) % e['P']) - 0.5 * e['P']) < 1e-3:
            continue
        others |= in_tr(t, e['P'], e['t0'], e['dur'], 1.5)
    tin &= ~others; oot &= ~others
    if tin.sum() < nmin or oot.sum() < 4 * nmin:
        return None
    # need OOT on both sides
    if (t[oot] < tc).sum() < nmin or (t[oot] > tc).sum() < nmin:
        return None
    X = np.vstack([np.ones(oot.sum()), t[oot] - tc]).T
    Y = F[oot].reshape(oot.sum(), -1)
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    Xi = np.vstack([np.ones(tin.sum()), t[tin] - tc]).T
    model_in = Xi @ beta
    diff = (model_in - F[tin].reshape(tin.sum(), -1)).mean(0).reshape(F.shape[1:])
    direct = beta[0].reshape(F.shape[1:])
    return diff, direct

def gauss_fit(img, err, x0, y0, sig=None):
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    def model(p):
        A, x, y, s, c = p if sig is None else (*p[:3], sig, p[3])
        return A * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * s ** 2)) + c
    p0 = [img.max(), x0, y0, 1.0, 0.0] if sig is None else [img[int(round(y0)), int(round(x0))], x0, y0, 0.0]
    r = least_squares(lambda p: ((img - model(p)) / err).ravel(), p0)
    J = r.jac; cov = np.linalg.pinv(J.T @ J) * max(1, np.sum(r.fun ** 2) / (r.fun.size - len(p0)))
    return r.x, np.sqrt(np.diag(cov))

def aperture_phot_depths(tp, centre, epochs, dur):
    ny, nx = tp['F'].shape[1:]
    yy, xx = np.mgrid[0:ny, 0:nx]
    cx, cy = centre
    r = np.hypot(xx - cx, yy - cy)
    aps = {'r<=0.75px(1-2px)': r <= 0.75, 'r<=1.5px': r <= 1.5, 'r<=2.5px': r <= 2.5, 'r<=3.5px': r <= 3.5,
           'annulus 1.5-3.5px': (r > 1.5) & (r <= 3.5), 'SPOC': (tp['ap'] & 2) > 0}
    out = {}
    for name, a in aps.items():
        fl = tp['F'][:, a].sum(1)
        fl = fl / np.median(fl)
        # per-epoch local linear depth (box, central 70%)
        D, W = [], []
        for tc in epochs:
            t = tp['t']; w = np.abs(t - tc) < 0.75 * dur + 0.3
            tin = np.abs(t - tc) < 0.35 * dur; oot = w & (np.abs(t - tc) > 0.75 * dur)
            if tin.sum() < 15 or oot.sum() < 60: continue
            X = np.vstack([np.ones(oot.sum()), t[oot] - tc]).T
            b, *_ = np.linalg.lstsq(X, fl[oot], rcond=None)
            res = fl[oot] - X @ b
            d = b[0] - fl[tin].mean(); e = np.std(res) * np.sqrt(1 / tin.sum() + 1 / oot.sum())
            D.append(d); W.append(1 / e ** 2)
        if D:
            D, W = np.array(D), np.array(W)
            out[name] = dict(depth=float(np.sum(D * W) / W.sum()), err=float(1 / np.sqrt(W.sum())), npix=int(a.sum()), n=len(D))
    return out

results = {pl: [] for pl in PLS}
apdepth = {pl: {} for pl in PLS}
for s in SECTORS:
    tp = load_tpf(s)
    t = tp['t']
    ny, nx = tp['F'].shape[1:]
    # target pixel position from WCS
    tx, ty = tp['wcs'].all_world2pix([[324.3894166666666 + 120.399e-3 / 3600 / np.cos(np.radians(-77.979)) * 0, -77.97914722222222]], 0)[0]
    for pl in PLS:
        e = E[pl]; P, t0, dur = e['P'], e['t0'], e['dur']
        n0, n1 = int(np.ceil((t.min() - t0) / P)), int(np.floor((t.max() - t0) / P))
        diffs, directs, eps = [], [], []
        for n in range(n0, n1 + 1):
            r = diff_one(tp, t0 + n * P, dur)
            if r is None: continue
            diffs.append(r[0]); directs.append(r[1]); eps.append(t0 + n * P)
        if not diffs: continue
        D = np.mean(diffs, 0); Dir = np.mean(directs, 0)
        # fake epochs for empirical noise (same number of events averaged)
        msk = all_transit_mask(t, 2.0)
        fakes = []
        tries = 0
        while len(fakes) < 40 and tries < 2000:
            tries += 1
            tc = rng.uniform(t.min() + 0.5, t.max() - 0.5)
            if msk[np.abs(t - tc) < 0.75 * dur + 0.3].any(): continue
            r = diff_one(tp, tc, dur)
            if r is not None: fakes.append(r[0])
        fakes = np.array(fakes)
        err1 = np.std(fakes, 0)  # single-event per-pixel noise
        errD = err1 / np.sqrt(len(diffs))
        errD = np.where(errD > 0, errD, np.nanmedian(errD))
        # direct image fit
        pd_, pde = gauss_fit(Dir, np.sqrt(np.abs(Dir)) + 1, tx, ty)
        sig = abs(pd_[3])
        pdf, pdfe = gauss_fit(D, errD, pd_[1], pd_[2], sig=sig)
        # Monte-Carlo the offset using fake-epoch noise realisations (resample fake means)
        mc = []
        for _ in range(60):
            idx = rng.choice(len(fakes), len(diffs), replace=True)
            noise = fakes[idx].mean(0)
            # rescale: subtract mean fake to keep zero-mean
            noise = noise - fakes.mean(0)
            pp, _ = gauss_fit(D + noise, errD * np.sqrt(2), pdf[1], pdf[2], sig=sig)
            mc.append(pp[1:3])
        mc = np.array(mc)
        # pixel -> sky offset (arcsec) of difference centroid relative to direct centroid
        sky_dir = tp['wcs'].all_pix2world([[pd_[1], pd_[2]]], 0)[0]
        sky_dif = tp['wcs'].all_pix2world([[pdf[1], pdf[2]]], 0)[0]
        dra = (sky_dif[0] - sky_dir[0]) * np.cos(np.radians(sky_dir[1])) * 3600
        ddec = (sky_dif[1] - sky_dir[1]) * 3600
        # error ellipse from MC in pixel space -> sky via local Jacobian
        J = np.zeros((2, 2)); h = 0.01
        for k in range(2):
            dp = np.array([pd_[1], pd_[2]]); dp2 = dp.copy(); dp2[k] += h
            a = tp['wcs'].all_pix2world([dp], 0)[0]; b2 = tp['wcs'].all_pix2world([dp2], 0)[0]
            J[0, k] = (b2[0] - a[0]) * np.cos(np.radians(a[1])) * 3600 / h; J[1, k] = (b2[1] - a[1]) * 3600 / h
        Cpix = np.cov(mc.T) if len(mc) > 5 else np.diag(pdfe[1:3] ** 2)
        Csky = J @ Cpix @ J.T
        snr_diff = pdf[0] / pdfe[0]
        rec = dict(sector=s, n_events=len(diffs), dra=float(dra), ddec=float(ddec), cov=Csky.tolist(),
                   sig_ra=float(np.sqrt(Csky[0, 0])), sig_dec=float(np.sqrt(Csky[1, 1])), amp=float(pdf[0]), amp_snr=float(snr_diff),
                   direct_xy=[float(pd_[1]), float(pd_[2])], diff_xy=[float(pdf[1]), float(pdf[2])], psf_sigma_px=float(sig))
        results[pl].append(rec)
        np.savez(OUT / 'data' / f'diffimg_{pl}_s{s}.npz', D=D, Dir=Dir, errD=errD, ap=tp['ap'])
        apdepth[pl][s] = aperture_phot_depths(tp, (pd_[1], pd_[2]), eps, dur)
        print(f"s{s} {pl:4s} n={len(diffs):2d} offset RA {dra:+6.1f}+-{rec['sig_ra']:.1f}  Dec {ddec:+6.1f}+-{rec['sig_dec']:.1f} arcsec  diffSNR {snr_diff:.1f}", flush=True)
# combine per planet (inverse-covariance weighted)
summary = {}
for pl in PLS:
    rows = [r for r in results[pl] if np.isfinite(r['dra'])]
    Wsum = np.zeros((2, 2)); ws = np.zeros(2)
    for r in rows:
        C = np.array(r['cov']); Wi = np.linalg.inv(C)
        Wsum += Wi; ws += Wi @ np.array([r['dra'], r['ddec']])
    Cm = np.linalg.inv(Wsum); m = Cm @ ws
    chi2_0 = float(m @ Wsum @ m)
    # scatter chi2 about mean
    chi2_s = float(sum((np.array([r['dra'], r['ddec']]) - m) @ np.linalg.inv(np.array(r['cov'])) @ (np.array([r['dra'], r['ddec']]) - m) for r in rows))
    summary[pl] = dict(dra=float(m[0]), ddec=float(m[1]), sig_ra=float(np.sqrt(Cm[0, 0])), sig_dec=float(np.sqrt(Cm[1, 1])), cov=Cm.tolist(),
                       offset=float(np.hypot(*m)), chi2_offset_vs_zero=chi2_0, scatter_chi2=chi2_s, dof=2 * len(rows) - 2, n_sectors=len(rows))
    print(pl, json.dumps({k: round(v, 2) if isinstance(v, float) else v for k, v in summary[pl].items() if k != 'cov'}))
json.dump(dict(per_sector=results, summary=summary, aperture_depths=apdepth), open(OUT / 'data' / 'tpf_vet.json', 'w'), indent=1)
