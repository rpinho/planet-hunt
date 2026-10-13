"""Light-curve vetting for d18 and d2: odd/even, secondary, per-season, per-product, half-split,
background and centroid behaviour in transit, momentum dumps."""
from fitlib import *
import json
store = load_store()
E = json.load(open(OUT / 'data' / 'ephem.json'))
rho = E['_rho']['rho']
res = {}
def residual_lc(kind):
    lc = concat(store, kind)
    f = lc['f'].copy()
    return lc, f
def fit_amp(t, f, e, pl, tc_shift=0.0, epochs=None, phase_off=0.0):
    """Global linear amplitude of the fixed-shape model (others subtracted), local baselines already flat."""
    ep = E[pl]
    m = transit_model(t, ep['t0'] + phase_off * ep['P'] + tc_shift, ep['P'], ep['k'], ep['b'], rho) - 1
    sel = in_tr(t, ep['P'], ep['t0'] + phase_off * ep['P'], ep['dur'], 4)
    if epochs is not None:
        n = np.round((t - ep['t0'] - phase_off * ep['P']) / ep['P'])
        sel &= np.isin(n, epochs)
    if (np.abs(m[sel]) > 0).sum() < 5: return np.nan, np.nan
    w = 1 / e[sel] ** 2
    X = np.vstack([np.ones(sel.sum()), m[sel]]).T
    A = X.T @ (X * w[:, None]); bvec = X.T @ (f[sel] * w)
    c = np.linalg.solve(A, bvec); C = np.linalg.inv(A)
    r = f[sel] - X @ c
    s = max(1, np.sum(r ** 2 * w) / (sel.sum() - 2))
    return c[1], np.sqrt(C[1, 1] * s)
for pl in ['d18', 'd2']:
    R = {}
    ep = E[pl]
    for kind in ['2min', '2min_loose', '2min_sap', 'qlp', 'tspoc', '20s']:
        lc = concat(store, kind)
        if len(lc['t']) == 0: continue
        t, f, e = lc['t'], lc['f'].copy(), lc['e']
        for o in ['b', 'c', 'd18', 'd2']:
            if o == pl: continue
            oe = E[o]; f -= transit_model(t, oe['t0'], oe['P'], oe['k'], oe['b'], rho) - 1
        n = np.round((t - ep['t0']) / ep['P']).astype(int)
        a, s = fit_amp(t, f, e, pl)
        row = dict(amp=a, err=s, snr=a / s, depth_ppm=a * ep['k'] ** 2 * 1e6)
        if kind == '2min':
            ue = np.unique(n[in_tr(t, ep['P'], ep['t0'], ep['dur'])])
            odd = ue[ue % 2 != 0]; even = ue[ue % 2 == 0]
            ao, so = fit_amp(t, f, e, pl, epochs=odd); ae, se = fit_amp(t, f, e, pl, epochs=even)
            row['odd'] = (ao, so); row['even'] = (ae, se); row['oe_sigma'] = abs(ao - ae) / np.hypot(so, se)
            a5, s5 = fit_amp(t, f, e, pl, phase_off=0.5)
            row['secondary_amp'] = (a5, s5)
            # any-phase secondary: scan phase 0.1-0.9 with same shape
            scan = []
            for ph in np.arange(0.1, 0.9, 0.005):
                aa, ss = fit_amp(t, f, e, pl, phase_off=ph)
                scan.append((ph, aa / ss if np.isfinite(aa) else 0))
            scan = np.array(scan); i = np.argmax(np.abs(scan[:, 1]))
            row['max_any_phase'] = (float(scan[i, 0]), float(scan[i, 1]))
            row['any_phase_scan_sd'] = float(np.std(scan[:, 1]))
            # per-season and halves
            seas = {'Y1-2 (S13,S27)': [13, 27], 'Y5 (S66,S67)': [66, 67], 'Y7 (S93,S94)': [93, 94], 'Y8 (S101-103)': [101, 102, 103]}
            row['seasons'] = {}
            for lab, ss_ in seas.items():
                mm = np.isin(lc['s'], ss_)
                aa, sss = fit_amp(t[mm], f[mm], e[mm], pl)
                row['seasons'][lab] = (aa, sss)
            for lab, ss_ in {'first half (S13-S67)': [13, 27, 66, 67], 'second half (S93-S103)': [93, 94, 101, 102, 103]}.items():
                mm = np.isin(lc['s'], ss_); aa, sss = fit_amp(t[mm], f[mm], e[mm], pl); row['seasons'][lab] = (aa, sss)
            ch = [(v[0] - a) / v[1] for k_, v in row['seasons'].items() if k_.startswith('Y') and np.isfinite(v[0])]
            row['season_chi2'] = float(np.sum(np.square(ch))); row['season_dof'] = len(ch) - 1
            # duration test: amplitude with model of 0.5x, 0.75x, 1.25x, 1.5x duration (via density)
        R[kind] = row
        print(pl, kind, json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items() if k not in ('seasons',)}, default=lambda x: [round(float(y), 3) for y in x]))
        if 'seasons' in row:
            print('   seasons:', {k: f'{v[0]:.2f}+-{v[1]:.2f}' for k, v in row['seasons'].items()}, 'chi2', round(row['season_chi2'], 2), 'dof', row['season_dof'])
    res[pl] = R
json.dump(res, open(OUT / 'data' / 'vetlc.json', 'w'), indent=1, default=lambda x: [float(y) for y in x])
