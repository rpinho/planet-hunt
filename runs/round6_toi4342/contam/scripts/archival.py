"""Archival-imaging background-star test for TOI-4342.
usage: python -I archival.py <downloads_dir> <out_json> <out_png>"""
import sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import curve_fit
from scipy.ndimage import map_coordinates
dl, out_json, out_png = sys.argv[1], sys.argv[2], sys.argv[3]

SID = 6355915466181029376
EPOCH_NOW = Time('2026-10-13').decimalyear
TMAG_TGT = 11.0318  # TIC v8.2
gaia = Table.read(f'{dl}/gaia_260arcsec.ecsv'); gaia.rename_columns(gaia.colnames, [c.lower() for c in gaia.colnames])
tm = Table.read(f'{dl}/2mass_psc_260arcsec.ecsv')
tg = gaia[gaia['source_id'] == SID][0]
RA0, DE0, PMRA, PMDE = float(tg['ra']), float(tg['dec']), float(tg['pmra']), float(tg['pmdec'])

def prop(ra, dec, pmra, pmdec, t):
    """linear propagation from 2016.0 (pmra includes cos dec), mas/yr; parallax ignored (<=16 mas)"""
    dt = t - 2016.0
    pmra = np.nan_to_num(np.asarray(pmra, float)); pmdec = np.nan_to_num(np.asarray(pmdec, float))
    dec2 = dec + pmdec * dt / 3.6e6
    ra2 = ra + pmra * dt / 3.6e6 / np.cos(np.deg2rad(dec))
    return ra2, dec2

def tmag(G, c):
    c = np.asarray(c, float); G = np.asarray(G, float)
    t = G - 0.00522555*c**3 + 0.0891337*c**2 - 0.633923*c + 0.0324473
    return np.where(np.isfinite(c), t, G - 0.43)
gaia['T'] = tmag(gaia['phot_g_mean_mag'], gaia['bp_rp'])
# 2MASS <-> Gaia match (at 2000.55)
ra_g, de_g = prop(np.array(gaia['ra']), np.array(gaia['dec']), gaia['pmra'].filled(0) if hasattr(gaia['pmra'],'filled') else gaia['pmra'], gaia['pmdec'].filled(0) if hasattr(gaia['pmdec'],'filled') else gaia['pmdec'], 2000.55)
cg = SkyCoord(ra_g, de_g, unit='deg'); cm = SkyCoord(tm['RAJ2000'], tm['DEJ2000'], unit='deg')
idx, d2, _ = cg.match_to_catalog_sky(cm)
for b in 'JHK':
    v = np.array(tm[b+'mag'], float)[idx]; v[d2.arcsec > 1.5] = np.nan; gaia[b] = v
gaia['bp'] = np.array(gaia['phot_bp_mean_mag'], float); gaia['g'] = np.array(gaia['phot_g_mean_mag'], float); gaia['rp'] = np.array(gaia['phot_rp_mean_mag'], float)

IMAGES = [
 dict(key='SERC-J_1976', path='stsci_poss2ukstu_blue.fits', survey='DSS1 / UKST SERC-J (IIIaJ+GG395), plate 00Z6', band='B_J (proxy: Gaia BP)', cal='bp'),
 dict(key='AAO-SES-R_1992', path='stsci_poss2ukstu_red.fits', survey='DSS2 / AAO-SES (IIIaF+RG610), plate A0UD', band='R_F (proxy: Gaia G)', cal='g'),
 dict(key='SERC-I_1995', path='stsci_poss2ukstu_ir.fits', survey='DSS2 / UKST SERC-I (IVN+RG715), plate A2YG', band='I_N (proxy: Gaia RP)', cal='rp'),
 dict(key='2MASS-J_2000', path='2mass_J_000716s072_0267.fits.gz', survey='2MASS All-Sky Atlas 000716s072_0267', band='J', cal='J', date='2000-07-16'),
 dict(key='2MASS-H_2000', path='2mass_H_000716s072_0267.fits.gz', survey='2MASS All-Sky Atlas 000716s072_0267', band='H', cal='H', date='2000-07-16'),
 dict(key='2MASS-K_2000', path='2mass_K_000716s072_0267.fits.gz', survey='2MASS All-Sky Atlas 000716s072_0267', band='Ks', cal='K', date='2000-07-16'),
 dict(key='SkyMapper-r_2014', path='skymapper_r_20140928103735-29.fits', survey='SkyMapper DR4 short (5 s)', band='r (proxy: Gaia G)', cal='g'),
 dict(key='SkyMapper-i_2014', path='skymapper_i_20141003101808-29.fits', survey='SkyMapper DR4 main (100 s)', band='i (proxy: Gaia RP)', cal='rp'),
]

def g2d(xy, a, x0, y0, s, b):
    x, y = xy
    return a*np.exp(-((x-x0)**2+(y-y0)**2)/(2*s**2)) + b

def fit_star(img, x, y, r):
    ny, nx = img.shape
    x0, x1, y0, y1 = int(round(x-r)), int(round(x+r))+1, int(round(y-r)), int(round(y+r))+1
    if x0 < 0 or y0 < 0 or x1 > nx or y1 > ny: return None
    sub = img[y0:y1, x0:x1].astype(float)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    b0 = np.median(np.concatenate([sub[0], sub[-1], sub[:, 0], sub[:, -1]]))
    try:
        p, cov = curve_fit(g2d, (xx.ravel(), yy.ravel()), sub.ravel(), p0=[sub.max()-b0, x, y, r/3, b0], maxfev=4000)
    except Exception: return None
    if not (abs(p[1]-x) < r and abs(p[2]-y) < r and 0.3 < abs(p[3]) < r and p[0] > 0): return None
    return p[1], p[2], abs(p[3])*2.3548, p[0], p[4]

def aper(img, x, y, r, rin=None, rout=None):
    yy, xx = np.indices(img.shape)
    rr = np.hypot(xx-x, yy-y)
    m = rr <= r
    s = img[m].sum(); n = m.sum()
    if rin is not None:
        bg = np.median(img[(rr >= rin) & (rr <= rout)]); s -= bg*n
    return s

def epoch_of(h, im):
    if 'date' in im: return Time(im['date']).decimalyear
    if h.get('MJD-OBS'): return Time(h['MJD-OBS'], format='mjd').decimalyear
    return Time(h['DATE-OBS']).decimalyear

results = []; panels = {}
for im in IMAGES:
    try:
        hdu = fits.open(f"{dl}/{im['path']}")[0]
        img = np.array(hdu.data, float); h = hdu.header
        img[~np.isfinite(img)] = np.nanmedian(img)
        w = WCS(h)
        ep = epoch_of(h, im)
        scale = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600)
        # Gaia field stars at epoch
        rg, dg = prop(np.array(gaia['ra']), np.array(gaia['dec']), np.array(gaia['pmra'].filled(0) if hasattr(gaia['pmra'],'filled') else gaia['pmra']), np.array(gaia['pmdec'].filled(0) if hasattr(gaia['pmdec'],'filled') else gaia['pmdec']), ep)
        px, py = w.all_world2pix(rg, dg, 0)
        istgt = np.array(gaia['source_id']) == SID
        # first-pass FWHM guess
        fw_guess = 3.0/scale
        cs = SkyCoord(rg, dg, unit='deg'); sep_t = cs.separation(cs[istgt][0]).arcsec
        # isolation: no brighter-or-similar Gaia star within 4*fwhm (approx 12")
        Gm = np.array(gaia['phot_g_mean_mag'], float)
        iso = np.ones(len(gaia), bool)
        for i in range(len(gaia)):
            dd = cs.separation(cs[i]).arcsec
            iso[i] = not np.any((dd > 0.01) & (dd < 12) & (Gm < Gm[i]+3))
        meas = []
        for i in range(len(gaia)):
            if istgt[i] or sep_t[i] < 25 or not iso[i] or not (Gm[i] < 19.0): continue
            f = fit_star(img, px[i], py[i], max(3*fw_guess, 4))
            if f is not None: meas.append((i, px[i], py[i]) + f)
        meas = np.array(meas)
        # affine refinement predicted-pixel -> measured-pixel
        P = np.c_[meas[:, 1], meas[:, 2], np.ones(len(meas))]
        good = np.ones(len(meas), bool)
        for it in range(4):
            Ax, *_ = np.linalg.lstsq(P[good], meas[good, 3], rcond=None); Ay, *_ = np.linalg.lstsq(P[good], meas[good, 4], rcond=None)
            res = np.hypot(P@Ax - meas[:, 3], P@Ay - meas[:, 4])
            sig = 1.4826*np.median(res[good]); good = res < max(3*sig, 0.3)
        astrom_rms = float(np.sqrt(np.mean(res[good]**2))*scale)
        def topix(ra, dec):
            x, y = w.all_world2pix(ra, dec, 0); v = np.array([x, y, 1.0]); return float(v@Ax), float(v@Ay)
        # FWHM from moderately faint stars
        fwhms = meas[good, 5]
        fwhm_pix = float(np.median(fwhms)); fwhm = fwhm_pix*scale
        # target
        rt, dt_ = prop(RA0, DE0, PMRA, PMDE, ep); xt_pred, yt_pred = topix(rt, dt_)
        ft = fit_star(img, xt_pred, yt_pred, max(10/scale, 5))
        xt, yt = (ft[0], ft[1]) if ft else (xt_pred, yt_pred)
        off_meas_pred = float(np.hypot(xt-xt_pred, yt-yt_pred)*scale)
        rn, dn = prop(RA0, DE0, PMRA, PMDE, EPOCH_NOW); xn, yn = topix(rn, dn)
        disp = float(np.hypot(xn-xt, yn-yt)*scale)
        disp_pm = float(np.hypot(PMRA, PMDE)*(EPOCH_NOW-ep)/1000)
        # self-subtraction: residual = img(x) - img(2c-x)
        yy, xx = np.indices(img.shape).astype(float)
        def selfsub(a):
            rot = map_coordinates(a, [2*yt-yy, 2*xt-xx], order=3, mode='nearest')
            return a - rot
        R = selfsub(img)
        rap = max(fwhm_pix, 1.2)
        # noise: ring of apertures at same radius as current position
        rad = np.hypot(xn-xt, yn-yt); ang0 = np.arctan2(yn-yt, xn-xt)
        others = np.where((~istgt) & (sep_t < rad*scale+40))[0]
        pp = np.array([topix(rg[j], dg[j]) for j in others]) if len(others) else np.zeros((0, 2))
        vals = []
        for rf in (0.85, 1.0, 1.15):
            rr_ = rad*rf; nst = max(int(2*np.pi*rr_/rap), 8)
            for k in range(nst):
                a = ang0 + 2*np.pi*k/nst
                dang = np.angle(np.exp(1j*(a-ang0)))
                if abs(dang) < np.deg2rad(30) or abs(abs(dang)-np.pi) < np.deg2rad(30): continue
                xr, yr = xt+rr_*np.cos(a), yt+rr_*np.sin(a); xa, ya = 2*xt-xr, 2*yt-yr
                if len(pp) and (np.min(np.hypot(pp[:, 0]-xr, pp[:, 1]-yr)) < 2.5*fwhm_pix or np.min(np.hypot(pp[:, 0]-xa, pp[:, 1]-ya)) < 2.5*fwhm_pix): continue
                vals.append(aper(R, xr, yr, rap))
        vals = np.array(vals)
        sigma = float(1.4826*np.median(np.abs(vals-np.median(vals)))) if len(vals) >= 6 else np.nan
        resid_at_now = float(aper(R, xn, yn, rap))
        # local brightness at current position (saturation check)
        sky = float(np.median(img)); peak = float(np.percentile(img[int(yt)-2:int(yt)+3, int(xt)-2:int(xt)+3], 90))
        val_now = float(np.median(img[int(round(yn))-1:int(round(yn))+2, int(round(xn))-1:int(round(xn))+2]))
        frac_to_sat = (val_now-sky)/(peak-sky)
        # photometric calibration: aperture sum vs catalogue mag
        cal = np.array(gaia[im['cal']], float)
        cm_, cf_ = [], []
        for j in range(len(meas)):
            i = int(meas[j, 0])
            if not np.isfinite(cal[i]): continue
            F = aper(img, meas[j, 3], meas[j, 4], rap, 3*fwhm_pix, 5*fwhm_pix)
            if F > 0: cm_.append(cal[i]); cf_.append(np.log10(F))
        cm_, cf_ = np.array(cm_), np.array(cf_)
        m_tgt = {'bp': float(tg['phot_bp_mean_mag']), 'g': float(tg['phot_g_mean_mag']), 'rp': float(tg['phot_rp_mean_mag']), 'J': 9.832, 'H': 9.179, 'K': 9.018}[im['cal']]
        sel = cm_ > m_tgt + 2.5
        if sel.sum() < 5: sel = np.ones(len(cm_), bool)
        sl = np.ones(sel.sum(), bool)
        for it in range(4):
            coef = np.polyfit(cm_[sel][sl], cf_[sel][sl], 1)
            r_ = cf_[sel] - np.polyval(coef, cm_[sel]); s_ = 1.4826*np.median(np.abs(r_[sl])); sl = np.abs(r_) < 3*max(s_, 0.03)
        def F_of_m(m): return 10**np.polyval(coef, m)
        # injection of REAL field-star stamps (known catalogue mag) at the current position
        stamps = []
        hs = int(np.ceil(3*fwhm_pix))
        for j in range(len(meas)):
            i = int(meas[j, 0])
            if not np.isfinite(cal[i]) or cal[i] < m_tgt + 2.5: continue
            xc, yc = meas[j, 3], meas[j, 4]
            ix, iy = int(round(xc)), int(round(yc))
            if ix-hs-6 < 0 or iy-hs-6 < 0 or ix+hs+6 >= img.shape[1] or iy+hs+6 >= img.shape[0]: continue
            yy_, xx_ = np.indices(img.shape)
            rr2 = np.hypot(xx_-xc, yy_-yc); bg = np.median(img[(rr2 > 3*fwhm_pix) & (rr2 < 5*fwhm_pix)])
            # stamp resampled onto grid centred at current position
            gy, gx = np.mgrid[-hs:hs+1, -hs:hs+1]
            X0, Y0 = int(round(xn)), int(round(yn)); fx, fy = xn-X0, yn-Y0
            st = map_coordinates(img, [yc+gy-fy, xc+gx-fx], order=1, mode='nearest') - bg
            st[np.hypot(gx-fx, gy-fy) > 3*fwhm_pix] = 0
            inj = img.copy(); inj[Y0-hs:Y0+hs+1, X0-hs:X0+hs+1] += st
            exc = aper(selfsub(inj), xn, yn, rap) - resid_at_now
            stamps.append((cal[i], (resid_at_now+exc)/sigma))
        stamps = np.array(stamps); print(im['key'], 'stamps(mag,snr):', np.round(stamps[np.argsort(stamps[:,0])],2).tolist() if len(stamps) else [], 'sigma', sigma)
        m5 = m3 = None; lim_note = ''
        if len(stamps) >= 4:
            ok = stamps[:, 1] > 0.3
            cf2 = np.polyfit(stamps[ok, 0], np.log10(stamps[ok, 1]), 1)
            m5 = float((np.log10(5)-cf2[1])/cf2[0]); m3 = float((np.log10(3)-cf2[1])/cf2[0]); m5_raw = m5; print('uncapped m5', m5, 'm3', m3)
            faint = float(stamps[:, 0].max())
            det_faint = stamps[stamps[:, 1] > 5]
            lim_note = f'{len(stamps)} real-star stamps injected (mag {stamps[:,0].min():.1f}-{faint:.1f}); faintest recovered at 5 sigma: {det_faint[:,0].max() if len(det_faint) else None}'
            if m5 > faint + 0.3:
                lim_note += f'; 5-sigma fit extrapolates beyond faintest stamp ({faint:.2f}) - capped'
                m5 = min(m5, faint + 0.3); m3 = min(m3, faint + 0.3)
        resolved = disp > 1.5*fwhm
        sat_flag = frac_to_sat > 0.3
        if not resolved: lim_note += '; NOT a background test: shift to now < 1.5 FWHM (target still overlaps the current position)'
        if sat_flag: lim_note += '; current position lies in the target halo at >30% of peak (photographic/detector non-linearity): limit unreliable'
        if abs(resid_at_now/sigma) > 3: lim_note += f'; residual at current position = {resid_at_now/sigma:.1f} sigma (PSF asymmetry/defect, inspect)'
        cal_range = [float(cm_[sel][sl].min()), float(cm_[sel][sl].max())] if sel.sum() else None
        # background-colour correction: (band - T) median for field stars vs target
        bmT = np.array(gaia[im['cal']], float) - np.array(gaia['T'], float)
        fld = (~istgt) & np.isfinite(bmT) & (Gm < 19)
        bmT_bg = float(np.median(bmT[fld])); bmT_tgt = m_tgt - float(gaia['T'][istgt][0])
        d5 = None if m5 is None else m5 - m_tgt
        r = dict(survey=im['survey'], key=im['key'], band=im['band'], epoch=round(ep, 3), pixel_scale_arcsec=round(scale, 3),
                 fwhm_arcsec=round(fwhm, 2), n_astrom_stars=int(good.sum()), astrom_rms_arcsec=round(astrom_rms, 3),
                 target_measured_minus_predicted_arcsec=round(off_meas_pred, 2),
                 offset_of_target_from_now_arcsec=round(disp, 2), offset_from_pm_only_arcsec=round(disp_pm, 2),
                 offset_over_fwhm=round(disp/fwhm, 2),
                 current_position_brightness_frac_of_target_peak=round(frac_to_sat, 3),
                 resid_flux_at_current_pos_sigma=round(resid_at_now/sigma, 2) if sigma > 0 else None, n_ring_apertures=int(len(vals)),
                 target_mag_in_proxy=m_tgt, calibration_mag_range=cal_range, calib_slope_dlogF_dm=round(float(coef[0]), 3),
                 limit_mag_at_current_position=None if m5 is None else round(m5, 2),
                 limit_mag_at_current_position_3sigma=None if m3 is None else round(m3, 2),
                 delta_mag_limit=None if d5 is None else round(d5, 2), usable_for_background_test=bool(resolved and not sat_flag), cal_label={'bp':'BP','g':'G','rp':'RP','J':'J','H':'H','K':'Ks'}[im['cal']], injected_stamps_mag_snr=np.round(stamps,2).tolist() if len(stamps) else [], m5_uncapped=None if len(stamps)<4 else round(m5_raw,2), sigma_aperture=sigma, note=lim_note.strip('; '),
                 delta_T_limit_same_colour=None if d5 is None else round(d5, 2),
                 delta_T_limit_typical_field_colour=None if d5 is None else round(d5 - (bmT_bg - bmT_tgt), 2),
                 band_minus_T_target=round(bmT_tgt, 2), band_minus_T_field_median=round(bmT_bg, 2))
        results.append(r); print(json.dumps(r))
        panels[im['key']] = dict(img=img, w=w, xt=xt, yt=yt, xn=xn, yn=yn, xp=xt_pred, yp=yt_pred, R=R, scale=scale, ep=ep, fwhm=fwhm, Ax=Ax, Ay=Ay, rn=rn, dn=dn)
    except Exception as e:
        import traceback; traceback.print_exc()
        results.append(dict(key=im['key'], survey=im['survey'], error=repr(e)))
json.dump(results, open(out_json, 'w'), indent=1)

# ---------- figure ----------
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
keys = [('SERC-J_1976', 'img'), ('AAO-SES-R_1992', 'img'), ('2MASS-J_2000', 'img'), ('SkyMapper-i_2014', 'img'), ('SERC-J_1976', 'R'), ('AAO-SES-R_1992', 'R')]
fig, axs = plt.subplots(2, 3, figsize=(13.5, 9.4))
half = 40.0
for ax, (k, kind) in zip(axs.ravel(), keys):
    if k not in panels: ax.set_axis_off(); ax.set_title(f'{k}: failed'); continue
    p = panels[k]; s = p['scale']; hp = half/s
    a = p[kind]
    x0, x1, y0, y1 = int(p['xt']-hp), int(p['xt']+hp), int(p['yt']-hp), int(p['yt']+hp)
    cut = a[max(y0,0):y1, max(x0,0):x1]
    if kind == 'img':
        lo, hi = np.percentile(cut, [5, 99.7]); norm = lambda v: np.arcsinh((v-lo)/(0.05*(hi-lo)+1e-9)); show = norm(cut); cmap = 'gray_r'
        vmin, vmax = np.nanmin(show), np.nanmax(show)
    else:
        sd = 1.4826*np.median(np.abs(cut-np.median(cut))); show = cut; cmap = 'RdBu_r'; vmin, vmax = -8*sd, 8*sd
    ext = [(max(x0,0)-p['xt'])*s, (x1-p['xt'])*s, (max(y0,0)-p['yt'])*s, (y1-p['yt'])*s]
    ax.imshow(show, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, extent=ext, interpolation='nearest')
    # markers in arcsec offsets (pixel frame, axes scaled by pixel scale)
    ax.add_patch(Circle(((p['xp']-p['xt'])*s, (p['yp']-p['yt'])*s), 3.0, fill=False, ec='tab:orange', lw=1.5, label=f'Gaia-predicted at {p["ep"]:.1f}'))
    ax.plot((p['xn']-p['xt'])*s, (p['yn']-p['yt'])*s, marker='+', ms=16, mew=2, color='tab:red', label='position on 2026-10-13')
    # TESS 21" pixel aligned to RA/Dec, centred on current position
    corners = []
    for dx, dy in [(-10.5, -10.5), (10.5, -10.5), (10.5, 10.5), (-10.5, 10.5)]:
        rr = p['rn'] + dx/3600/np.cos(np.deg2rad(p['dn'])); dd = p['dn'] + dy/3600
        xx_, yy_ = p['w'].all_world2pix(rr, dd, 0); v = np.array([xx_, yy_, 1.0])
        corners.append(((v@p['Ax']-p['xt'])*s, (v@p['Ay']-p['yt'])*s))
    ax.add_patch(Polygon(corners, closed=True, fill=False, ec='tab:cyan', lw=1.3, ls='--', label='21" TESS pixel'))
    # north / east arrows
    x_, y_ = p['w'].all_world2pix([p['rn'], p['rn'], p['rn'] + 8/3600/np.cos(np.deg2rad(p['dn']))], [p['dn'], p['dn']+8/3600, p['dn']], 0)
    X = (np.c_[x_, y_, np.ones(3)]@p['Ax']-p['xt'])*s; Y = (np.c_[x_, y_, np.ones(3)]@p['Ay']-p['yt'])*s
    bx, by = -half*0.6, -half*0.72
    for (lab, j) in [('N', 1), ('E', 2)]:
        ax.annotate('', xy=(bx+X[j]-X[0], by+Y[j]-Y[0]), xytext=(bx, by), arrowprops=dict(arrowstyle='->', color='gold', lw=1.5))
        ax.text(bx+1.35*(X[j]-X[0]), by+1.35*(Y[j]-Y[0]), lab, color='gold', ha='center', va='center', fontsize=9)
    ax.set_xlim(-half, half); ax.set_ylim(-half, half)
    r = next(x for x in results if x.get('key') == k)
    ttl = f"{k.replace('_', ' ')}  (epoch {r['epoch']:.2f})" + ('' if kind == 'img' else '\n180-deg self-subtracted residual')
    ax.set_title(ttl, fontsize=10)
    lab2 = (f"5-sigma limit at +: d{r['cal_label']} >= {r['delta_mag_limit']}" if r['usable_for_background_test'] else
            ("not a background test: shift < 1.5 FWHM" if r['offset_over_fwhm'] < 1.5 else "+ inside target halo: no limit"))
    ax.text(0.02, 0.98, f"FWHM {r['fwhm_arcsec']:.1f}\"  shift to now {r['offset_of_target_from_now_arcsec']:.1f}\"\n" + lab2, transform=ax.transAxes, va='top', fontsize=8, color='k', bbox=dict(fc='white', alpha=0.75, ec='none'))
    ax.set_xlabel('offset (arcsec, image x)'); ax.set_ylabel('offset (arcsec, image y)')
axs[0, 0].legend(loc='lower right', fontsize=7, framealpha=0.85)
fig.suptitle('TOI-4342 (TIC 354944123): archival background-star check; Gaia DR3 pm = (120.3, -91.5) mas/yr', fontsize=12)
fig.tight_layout(); fig.savefig(out_png, dpi=130)
print('saved', out_png)
