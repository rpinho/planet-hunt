"""Do b and c show the near-2:1 TTV signal (super-period ~152.6 d)?"""
from common import *
import json
pt = json.load(open(OUT / 'data' / 'per_transit.json'))
E = json.load(open(OUT / 'data' / 'ephem.json'))
Pb, Pc = E['b']['P'], E['c']['P']
Psup = 1 / abs(2 / Pc - 1 / Pb)
print('super-period b-c (2:1)', Psup)
data = {}
for pl in ['b', 'c']:
    sel = [r for r in pt[pl] if r['kind'] == '2min' and r['coverage'] > 0.5]
    data[pl] = (np.array([r['epoch'] for r in sel]), np.array([r['post_tc'] for r in sel]), np.array([r['post_err_min'] for r in sel]) / 1440)
def design(pl, ep, tc, with_sin, Ps=Psup):
    cols = [np.ones_like(tc), ep.astype(float)]
    if with_sin:
        cols += [np.sin(2 * np.pi * tc / Ps), np.cos(2 * np.pi * tc / Ps)]
    return np.vstack(cols).T
res = {}
for pl in ['b', 'c']:
    ep, tc, e = data[pl]
    for ws in [False, True]:
        A = design(pl, ep, tc, ws) / e[:, None]
        c, *_ = np.linalg.lstsq(A, tc / e, rcond=None)
        chi = np.sum((tc / e - A @ c) ** 2)
        amp = np.hypot(c[2], c[3]) * 1440 if ws else 0
        res[(pl, ws)] = (chi, amp)
        print(pl, 'sinusoid' if ws else 'linear', f'chi2 {chi:.1f} n {len(tc)} amp {amp:.1f} min')
# periodogram of b residual TTVs over 40-1000 d
ep, tc, e = data['b']
Ps = np.exp(np.linspace(np.log(30), np.log(1500), 3000))
dchi = []
A0 = design('b', ep, tc, False) / e[:, None]; c0, *_ = np.linalg.lstsq(A0, tc / e, rcond=None); chi0 = np.sum((tc / e - A0 @ c0) ** 2)
for P in Ps:
    A = design('b', ep, tc, True, P) / e[:, None]; c, *_ = np.linalg.lstsq(A, tc / e, rcond=None)
    dchi.append(chi0 - np.sum((tc / e - A @ c) ** 2))
dchi = np.array(dchi)
print('b TTV periodogram best P', Ps[np.argmax(dchi)], 'dchi2', dchi.max(), 'at Psup', np.interp(Psup, Ps, dchi))
json.dump(dict(Psup=Psup, b_lin=res[('b', False)], b_sin=res[('b', True)], c_lin=res[('c', False)], c_sin=res[('c', True)],
               b_pgram_best=[float(Ps[np.argmax(dchi)]), float(dchi.max())]), open(OUT / 'data' / 'ttv_bc.json', 'w'), indent=1)
