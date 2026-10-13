"""Future transit windows (ground-based), through 2027-12-31. tess-point: no TESS sectors 104-134."""
from common import *
import json
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_sun
import astropy.units as u
mc = json.load(open(OUT / 'data' / 'mcmc.json'))
SITES = {'CTIO/LCO-Chile': (-30.169, -70.806, 2200), 'SAAO/Sutherland': (-32.38, 20.81, 1798), 'SSO/Siding Spring': (-31.273, 149.071, 1165)}
star = SkyCoord(324.3894166666666 * u.deg, -77.97914722222222 * u.deg)
out = {}
for pl in ['d18', 'd2']:
    t0 = mc[pl]['t0']; P = mc[pl]['P']
    st0 = 0.5 * (t0[2] - t0[0]); sP = 0.5 * (P[2] - P[0])
    T14 = mc[pl]['T14_h'][1] / 24
    now = Time('2026-10-13').jd - 2457000; end = Time('2028-01-01').jd - 2457000
    n0 = int(np.ceil((now - t0[1]) / P[1])); n1 = int(np.floor((end - t0[1]) / P[1]))
    rows = []
    for n in range(n0, n1 + 1):
        tc = t0[1] + n * P[1]
        sig = np.sqrt(st0 ** 2 + (n * sP) ** 2) * 1440
        sig_tot = np.hypot(sig, 10.0)  # add ~10 min allowance for TTVs (observed O-C rms)
        row = dict(epoch=n, tc_btjd=round(tc, 4), tc_utc=Time(tc + 2457000, format='jd', scale='tdb').utc.iso[:16], sigma_ephem_min=round(sig, 1), sigma_with_ttv_min=round(sig_tot, 1), sites={})
        # observable if whole transit +- (30 min + 1 sigma) in astro night with alt > 25 deg
        pad = (30 + sig_tot) / 1440
        ts = Time(np.linspace(tc - T14 / 2 - pad, tc + T14 / 2 + pad, 40) + 2457000, format='jd', scale='tdb')
        for name, (lat, lon, h) in SITES.items():
            loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=h * u.m)
            fr = AltAz(obstime=ts, location=loc)
            sun = get_sun(ts).transform_to(fr).alt.deg
            alt = star.transform_to(fr).alt.deg
            ok = bool(np.all(sun < -12) and np.all(alt > 25))
            if ok:
                row['sites'][name] = dict(min_alt=round(float(alt.min()), 1), max_airmass=round(float(1 / np.sin(np.radians(alt.min()))), 2))
        rows.append(row)
    out[pl] = dict(t0=t0[1], P=P[1], T14_h=T14 * 24, rows=rows)
    good = [r for r in rows if r['sites']]
    print(pl, 'windows', len(rows), 'fully observable from >=1 site:', len(good))
    if pl == 'd18':
        for r in rows:
            print(f"  E{r['epoch']:4d} {r['tc_utc']} UTC  (BTJD {r['tc_btjd']})  +-{r['sigma_with_ttv_min']:.0f} min  ", ', '.join(f"{k} (X<={v['max_airmass']})" for k, v in r['sites'].items()) or '-')
json.dump(out, open(OUT / 'data' / 'future_windows.json', 'w'), indent=1)
