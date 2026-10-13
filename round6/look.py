"""Coverage of every predicted 18.9-d (and 1.96-d) window, every product."""
from common import *
import json
eph = load_ephem()
rows = []
for name in ['d18']:
    e = eph[name]
    for s in SECTORS:
        prods = {k: load(k, s) for k in ['2min', '20s', 'qlp', 'tspoc']}
        loose = load('2min', s, strict=False)
        tmin = min(p['t'].min() for p in prods.values() if p is not None)
        tmax = max(p['t'].max() for p in prods.values() if p is not None)
        n0, n1 = int(np.ceil((tmin - e['t0']) / e['P'] - 0.1)), int(np.floor((tmax - e['t0']) / e['P'] + 0.1))
        for n in range(n0, n1 + 1):
            tc = e['t0'] + n * e['P']
            row = dict(epoch=n, sector=s, tc=round(tc, 4), utc=btjd_to_utc(tc)[:16])
            for k, p in list(prods.items()) + [('2min_allq', loose)]:
                if p is None:
                    row[k] = None; continue
                w = np.abs(p['t'] - tc) < 0.5 * e['dur']
                cad = np.median(np.diff(p['t']))
                row[k] = round(w.sum() * cad / e['dur'], 2)
            # quality flags present in window
            w = np.abs(loose['t'] - tc) < 1.0 * e['dur']
            row['qflags'] = sorted(set(int(x) for x in loose['q'][w] if x))
            rows.append(row)
for r in rows: print(r)
json.dump(rows, open(OUT / 'data' / 'd18_windows_coverage.json', 'w'), indent=1)
