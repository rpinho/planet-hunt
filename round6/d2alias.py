from fitlib import *
import json
src = open('vetlc.py').read(); exec('def fit_amp' + src.split('def fit_amp')[1].split("for pl in ['d18', 'd2']:")[0])
store = load_store(); global E
E = json.load(open(OUT / 'data' / 'ephem.json')); rho = E['_rho']['rho']
lc = concat(store, '2min'); t, f, e = lc['t'], lc['f'].copy(), lc['e']
for o in ['b', 'c', 'd18']:
    oe = E[o]; f -= transit_model(t, oe['t0'], oe['P'], oe['k'], oe['b'], rho) - 1
out = {}
for off in [0, 1/3, 2/3, 0.25, 0.5, 0.75]:
    a, s = fit_amp(t, f, e, 'd2', phase_off=off); out[f'{off:.3f}'] = (a, s); print(f'phase {off:.3f}: amp {a:.2f}+-{s:.2f}')
ep = E['d2']; n = np.round((t - ep['t0']) / ep['P']).astype(int)
ue = np.unique(n[in_tr(t, ep['P'], ep['t0'], ep['dur'])])
for r in range(3):
    a, s = fit_amp(t, f, e, 'd2', epochs=ue[ue % 3 == r]); out[f'mod3={r}'] = (a, s); print(f'epochs mod 3 = {r}: {a:.2f}+-{s:.2f}')
json.dump(out, open(OUT / 'data' / 'd2_alias.json', 'w'), indent=1)
