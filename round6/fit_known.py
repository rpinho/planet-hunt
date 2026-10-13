"""Joint ML fit of b, c, d18, d2 (shared density) on 2-min PDCSAP; per-transit times for b and c."""
from fitlib import *
import json
store = load_store()
eph = load_ephem()
# reference epochs near middle of baseline to decorrelate t0,P
for pl in PL:
    e = eph[pl]; n = round((2900 - e['t0']) / e['P']); e['t0'] += n * e['P']
eph['b']['k'], eph['c']['k'], eph['d18']['k'], eph['d2']['k'] = 0.036, 0.037, 0.033, 0.017
lc = concat(store, '2min')
res, (t, f, e) = fit_joint(lc, eph, PL)
q = res.x
J = res.jac; cov = np.linalg.pinv(J.T @ J) * max(1, np.sum(res.fun ** 2) / (len(f) - len(q)))
err = np.sqrt(np.diag(cov))
out = {}
for i, pl in enumerate(PL):
    t0, P, k, b = q[4 * i:4 * i + 4]
    out[pl] = dict(t0=float(t0), t0_err=float(err[4*i]), P=float(P), P_err=float(err[4*i+1]), k=float(k), k_err=float(err[4*i+2]), b=float(b), b_err=float(err[4*i+3]),
                   depth=float(k**2), dur=eph[pl]['dur'])
    print(pl, f'P={P:.6f}+-{err[4*i+1]:.6f} t0={t0:.5f}+-{err[4*i]:.5f} k={k:.4f}+-{err[4*i+2]:.4f} b={b:.2f}+-{err[4*i+3]:.2f}')
print('rho', q[-1], '+-', err[-1], 'chi2/n', np.sum(res.fun**2), len(f))
out['_rho'] = dict(rho=float(q[-1]), err=float(err[-1]))
# durations from model
for pl in PL:
    o = out[pl]; ar = a_over_r(o['P'], q[-1]); inc = np.arccos(o['b'] / ar)
    T14 = o['P'] / np.pi * np.arcsin(np.sqrt((1 + o['k'])**2 - o['b']**2) / ar / np.sin(inc))
    o['dur'] = float(T14); print(pl, 'T14 h', T14 * 24)
json.dump(out, open(OUT / 'data' / 'ephem.json', 'w'), indent=1)
