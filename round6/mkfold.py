from fitlib import *
import json
store = load_store(); E = json.load(open(OUT / 'data' / 'ephem.json')); rho = E['_rho']['rho']
mc = json.load(open(OUT / 'data' / 'mcmc.json'))
lc = concat(store, '2min')
for pl in ['d18', 'd2']:
    f = lc['f'].copy()
    for o in ['b', 'c', 'd18', 'd2']:
        if o == pl: continue
        oe = E[o]; f -= transit_model(lc['t'], oe['t0'], oe['P'], oe['k'], oe['b'], rho) - 1
    P, t0 = mc[pl]['P'][1], mc[pl]['t0'][1]
    ph = (lc['t'] - t0 + 0.5 * P) % P - 0.5 * P
    s = np.abs(ph) < (0.3 if pl == 'd18' else 0.15)
    np.savez(OUT / 'tri' / f'fold_{pl}.npz', phase_d=ph[s], flux=f[s], err=lc['e'][s], period=P, depth=mc[pl]['depth_ppm'][1] * 1e-6, sectors=np.array(SECTORS))
    print(pl, s.sum())
