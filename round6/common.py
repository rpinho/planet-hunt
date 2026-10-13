"""Shared loaders / models for round 6 (TOI-4342)."""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planet_hunt import config  # noqa
import numpy as np
from pathlib import Path
from astropy.io import fits
warnings.filterwarnings('ignore')
RAW = config.CACHE / 'r6' / 'raw'
OUT = Path(__file__).resolve().parent.parent / 'runs' / 'round6_toi4342'
TIC = 354944123
RSTAR, MSTAR, TEFF = 0.598175, 0.586597, 3883.0
RHO_TIC = 3.864
U = (0.32075, 0.2995)
# working ephemerides (refined later by fit_known.py -> data/ephem.json)
EPH = {
    'b': dict(P=5.538329, t0=1654.532426, dur=2.23/24),
    'c': dict(P=10.688608, t0=1659.348023, dur=2.67/24),
    'd18': dict(P=18.90829164266999, t0=1675.5347772070088, dur=1.887648681061937/24),
    'd2': dict(P=1.9576890755606433, t0=1658.9524010321534, dur=1.256898305349431/24),
}
SECTORS = [13, 27, 66, 67, 93, 94, 101, 102, 103]

def load_ephem():
    import json
    p = OUT / 'data' / 'ephem.json'
    e = {k: dict(v) for k, v in EPH.items()}
    if p.exists():
        for k, v in json.load(open(p)).items():
            e.setdefault(k, {}).update(v)
    return e

def sector_of(fn):
    import re
    m = re.search(r'-s0*(\d+)-', fn) or re.search(r's0*(\d+)-0000', fn) or re.search(r'-s0*(\d+)_', fn)
    return int(m.group(1))

def _read(path, fluxcol, qualcol='QUALITY', errcol=None, strict=True):
    with fits.open(path, memmap=False) as h:
        d = h[1].data
        t = np.array(d['TIME'], float)
        f = np.array(d[fluxcol], float)
        e = np.array(d[errcol], float) if errcol and errcol in d.columns.names else np.full_like(f, np.nan)
        q = np.array(d[qualcol]) if qualcol in d.columns.names else np.zeros(len(t), int)
        hdr = dict(h[0].header)
    ok = np.isfinite(t) & np.isfinite(f)
    if strict:
        ok &= (q == 0)
    t, f, e, q = t[ok], f[ok], e[ok], q[ok]
    med = np.nanmedian(f)
    return dict(t=t, f=f / med, e=e / med, q=q, hdr=hdr)

def load(kind, sector, strict=True, flux=None):
    """kind: '2min' (PDCSAP), '2min_sap', '20s', 'qlp', 'tspoc'."""
    files = sorted(os.listdir(RAW))
    if kind in ('2min', '2min_sap'):
        fn = [x for x in files if x.startswith('tess20') and x.endswith('-s_lc.fits') and f'-s{sector:04d}-' in x]
        col = 'SAP_FLUX' if kind == '2min_sap' else 'PDCSAP_FLUX'
        ec = 'SAP_FLUX_ERR' if kind == '2min_sap' else 'PDCSAP_FLUX_ERR'
    elif kind == '20s':
        fn = [x for x in files if x.endswith('fast-lc.fits') and f'-s{sector:04d}-' in x]
        col, ec = flux or 'PDCSAP_FLUX', 'PDCSAP_FLUX_ERR'
    elif kind == 'qlp':
        fn = [x for x in files if 'hlsp_qlp' in x and f'_s{sector:04d}-' in x]
        col, ec = flux or 'KSPSAP_FLUX', 'KSPSAP_FLUX_ERR'
    elif kind == 'tspoc':
        fn = [x for x in files if 'hlsp_tess-spoc' in x and f'-s{sector:04d}_' in x]
        col, ec = 'PDCSAP_FLUX', 'PDCSAP_FLUX_ERR'
    if not fn:
        return None
    if kind == 'qlp' and flux is None:
        names = fits.getheader(RAW / fn[0], 1)
        cn = [names[k] for k in names if k.startswith('TTYPE')]
        if 'KSPSAP_FLUX' not in cn:
            col, ec = 'DET_FLUX', 'DET_FLUX_ERR'
    d = _read(RAW / fn[0], col, errcol=ec, strict=strict)
    if kind == 'qlp' and not strict:
        pass
    if not np.isfinite(d['e']).any():
        d['e'][:] = 1.4826 * np.nanmedian(np.abs(np.diff(d['f']))) / np.sqrt(2)
    d['kind'], d['sector'], d['file'] = kind, sector, fn[0]
    return d

def in_tr(t, P, t0, dur, widen=1.0):
    ph = (t - t0 + 0.5 * P) % P - 0.5 * P
    return np.abs(ph) < 0.5 * dur * widen

def mask_all(t, eph, keys=('b', 'c', 'd18', 'd2'), widen=1.6, ttv=None):
    m = np.zeros(len(t), bool)
    for k in keys:
        e = eph[k]
        m |= in_tr(t, e['P'], e['t0'], e['dur'], widen)
    return m

def detrend(d, mask, window=0.75):
    from wotan import flatten
    f = d['f'].copy()
    ff = f.copy(); ff[mask] = np.nan
    _, trend = flatten(d['t'], ff, method='biweight', window_length=window, return_trend=True, break_tolerance=0.4, edge_cutoff=0.0)
    # interpolate trend across masked points
    good = np.isfinite(trend)
    if good.sum() > 10:
        trend = np.interp(d['t'], d['t'][good], trend[good])
    out = dict(d); out['fd'] = f / trend; out['trend'] = trend
    return out

def a_over_r(P, rho):
    G = 6.674e-8
    return ((G * rho * (P * 86400) ** 2) / (3 * np.pi)) ** (1 / 3)

def transit_model(t, tc, P, k, b, rho, u=U, exp_time=None, supersample=1):
    import batman
    p = batman.TransitParams()
    p.t0, p.per, p.rp = tc, P, k
    ar = a_over_r(P, rho)
    p.a = ar
    cosi = min(b / ar, 0.999)
    p.inc = np.degrees(np.arccos(cosi))
    p.ecc, p.w, p.u, p.limb_dark = 0., 90., list(u), 'quadratic'
    if exp_time and supersample > 1:
        m = batman.TransitModel(p, t, supersample_factor=supersample, exp_time=exp_time)
    else:
        m = batman.TransitModel(p, t)
    return m.light_curve(p)

def btjd_to_utc(btjd):
    from astropy.time import Time
    return Time(btjd + 2457000.0, format='jd', scale='tdb').utc.iso
