"""Per-transit mid-times for b and c (and every window of d18, d2-stacked per sector)."""
from fitlib import *
import json
store = load_store()
E = json.load(open(OUT / 'data' / 'ephem.json'))
rho = E['_rho']['rho']
res = {}
for pl in ['b', 'c', 'd18']:
    e = E[pl]; P, t0, dur = e['P'], e['t0'], e['dur']
    rows = []
    for s in SECTORS:
        kinds = ['2min'] if pl != 'd18' else ['2min', '2min_loose', '20s', 'qlp', 'tspoc', '2min_sap']
        for kind in kinds:
            d = store.get((kind, s))
            if d is None: continue
            t, f = d['t'], d['fd']
            n0, n1 = int(np.ceil((t.min() - t0) / P)), int(np.floor((t.max() - t0) / P))
            for n in range(n0, n1 + 1):
                tc0 = t0 + n * P
                half = max(3 * dur, 0.35)
                w = np.abs(t - tc0) < half
                if w.sum() < 20: continue
                cad = np.median(np.diff(t[w]))
                # require coverage of >= 50% of the transit
                cov = (np.abs(t[w] - tc0) < 0.5 * dur).sum() * cad / dur
                if cov < 0.3: continue
                tt, ff = t[w], f[w]
                ee = np.full(w.sum(), d['mad'])
                # other planets' transits in window -> remove model of others
                for o in ['b', 'c', 'd18', 'd2']:
                    if o == pl: continue
                    oe = E[o]
                    ff = ff - (transit_model(tt, oe['t0'], oe['P'], oe['k'], oe['b'], rho) - 1)
                span = 0.25 if pl == 'd18' else 0.1
                grid = tc0 + np.arange(-span, span + 1e-9, 0.5 / 1440)
                chi, amp, aerr = local_fit(tt, ff, ee, grid, (e['k'], e['b']), P, rho)
                # fixed at linear ephemeris
                chi0, amp0, aerr0 = local_fit(tt, ff, ee, np.array([tc0]), (e['k'], e['b']), P, rho)
                # null (no transit) chi2
                X = np.vstack([np.ones_like(tt), tt - tt.mean()]).T / ee[:, None]
                bnull, *_ = np.linalg.lstsq(X, ff / ee, rcond=None)
                chinull = np.sum((ff / ee - X @ bnull) ** 2)
                # only physical (positive) depths for timing
                cpos = np.where(amp > 0, chi, chinull)
                i = int(np.argmin(cpos))
                # red-noise beta from residuals at best
                m = transit_model(tt, grid[i], P, e['k'], e['b'], rho) - 1
                Xf = np.vstack([np.ones_like(tt), tt - tt.mean(), m]).T
                bb, *_ = np.linalg.lstsq(Xf / ee[:, None], ff / ee, rcond=None)
                r = ff - Xf @ bb
                beta = red_beta(tt, r, ee, dur)
                rchi = max(1.0, np.sum((r / ee) ** 2) / (len(tt) - 4))
                scale = rchi * beta ** 2
                dchi = (cpos - cpos[i]) / scale
                wpost = np.exp(-0.5 * dchi); wpost /= wpost.sum()
                cdf = np.cumsum(wpost)
                p16, p50, p84 = [float(np.interp(q, cdf, grid)) for q in (0.16, 0.5, 0.84)]
                inside = grid[dchi <= 1]
                lo, hi = (inside.min(), inside.max()) if len(inside) else (np.nan, np.nan)
                rows.append(dict(planet=pl, kind=kind, sector=s, epoch=n, tc_lin=tc0, tc=float(grid[i]),
                                 tc_lo=float(lo), tc_hi=float(hi), oc_min=float((grid[i] - tc0) * 1440),
                                 oc_err_min=float(0.5 * (hi - lo) * 1440) if np.isfinite(lo) else None,
                                 amp_best=float(amp[i]), amp_best_err=float(aerr[i] * np.sqrt(scale)),
                                 amp_lin=float(amp0[0]), amp_lin_err=float(aerr0[0] * np.sqrt(scale)),
                                 depth_lin_ppm=float(amp0[0] * e['k']**2 * 1e6 * 1.0),
                                 dchi2_det=float((chinull - cpos[i]) / scale), dchi2_lin=float((chinull - chi0[0]) / scale) * (1 if amp0[0] > 0 else -1),
                                 coverage=float(cov), post_oc_min=float((p50 - tc0) * 1440), post_err_min=float(0.5 * (p84 - p16) * 1440), post_tc=p50, beta=float(beta), rchi=float(rchi), at_edge=bool(i in (0, len(grid) - 1))))
    res[pl] = rows
    for r in rows:
        if r['kind'] in ('2min',) or pl == 'd18':
            print(pl, r['kind'], r['sector'], r['epoch'], f"O-C {r['oc_min']:+.1f}+-{r['oc_err_min'] or -1:.1f} min  amp_lin {r['amp_lin']:.2f}+-{r['amp_lin_err']:.2f} amp_best {r['amp_best']:.2f}  cov {r['coverage']:.2f} beta {r['beta']:.2f} dchi2det {r['dchi2_det']:.1f}" + (' EDGE' if r['at_edge'] else ''))
json.dump(res, open(OUT / 'data' / 'per_transit.json', 'w'), indent=1)
