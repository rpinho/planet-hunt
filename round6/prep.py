"""Detrend every product per sector with all known/candidate transits masked; save npz."""
from common import *
import pickle
eph = load_ephem()
store = {}
for s in SECTORS:
    for kind in ['2min', '2min_sap', '20s', 'qlp', 'tspoc', '2min_loose']:
        if kind == '2min_loose':
            d = load('2min', s, strict=False)
            if d is None: continue
            # keep scattered-light / planet-search-exclude flagged cadences, drop hard flags
            hard = 1 | 2 | 4 | 8 | 32 | 128 | 1024
            ok = (d['q'] & hard) == 0
            for k in ['t', 'f', 'e', 'q']: d[k] = d[k][ok]
            d['kind'] = kind
        else:
            d = load(kind, s)
        if d is None: continue
        m = mask_all(d['t'], eph)
        dd = detrend(d, m, window=0.75)
        # 5-sigma clip of out-of-transit outliers (upward only + extreme down)
        r = dd['fd'] - 1
        mad = 1.4826 * np.nanmedian(np.abs(r - np.nanmedian(r)))
        keep = np.isfinite(r) & (r < 4 * mad) & ((r > -8 * mad) | m)
        for k in ['t', 'f', 'e', 'fd', 'trend', 'q']: dd[k] = dd[k][keep]
        dd['mad'] = mad
        dd.pop('hdr', None)
        store[(kind, s)] = dd
        print(kind, s, len(dd['t']), f'{mad*1e6:.0f} ppm', flush=True)
pickle.dump(store, open(OUT / 'data' / 'detrended.pkl', 'wb'))
