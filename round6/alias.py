"""Which period does the data support? Matched-filter period scan + explicit alias tests."""
from fitlib import *
import json
store = load_store()
E = json.load(open(OUT / 'data' / 'ephem.json'))
rho = E['_rho']['rho']
lc = concat(store, '2min')
# include QLP-only stretches? keep 2-min only for scan (homogeneous); flagged S13 cadences via loose
lcl = concat(store, '2min_loose')
t, f, e = lcl['t'], lcl['f'].copy(), lcl['e']
# subtract b, c, d2 models
for o in ['b', 'c', 'd2']:
    oe = E[o]; f -= transit_model(t, oe['t0'], oe['P'], oe['k'], oe['b'], rho) - 1
d = E['d18']
dur = d['dur']
# --- 1. single-event matched filter: depth time series at 10-min steps using box of the transit duration
# compute local depth estimate (box, inverse variance) for each trial centre
def box_depth(tc, dur):
    w_in = np.abs(t - tc) < 0.35 * dur
    w_out = (np.abs(t - tc) > 0.75 * dur) & (np.abs(t - tc) < 0.75 * dur + 0.3)
    if w_in.sum() < 10 or w_out.sum() < 30: return np.nan, np.nan
    dep = np.mean(f[w_out]) - np.mean(f[w_in])
    err = np.std(f[w_out]) * np.sqrt(1 / w_in.sum() + 1 / w_out.sum())
    return dep, err
step = 10 / 1440
grid = np.arange(t.min(), t.max(), step)
grid = grid[np.min(np.abs(grid[:, None] - t[None, ::50]), axis=1) < 0.05] if False else grid
D = np.full(len(grid), np.nan); S = np.full(len(grid), np.nan)
ts = np.sort(t)
for i, tc in enumerate(grid):
    j0, j1 = np.searchsorted(ts, tc - 0.5), np.searchsorted(ts, tc + 0.5)
    if j1 - j0 < 60: continue
    D[i], S[i] = box_depth(tc, dur)
np.savez(OUT / 'data' / 'single_event.npz', grid=grid, D=D, S=S)
# --- 2. period scan: combine single-event depths phase-coherently
Ps = np.exp(np.arange(np.log(4.0), np.log(80.0), 2e-5))
good = np.isfinite(D) & (S > 0)
g, Dg, Sg = grid[good], D[good], S[good]
wv = 1 / Sg ** 2
best = []
t_ref = d['t0']
snr = np.zeros(len(Ps)); nwin = np.zeros(len(Ps), int)
for i, P in enumerate(Ps):
    # for each P, scan phase on coarse grid and keep best: too slow -> instead fix the epoch near the
    # observed strongest event family: evaluate phases on a step of 10 min
    ph = ((g - t_ref) / P) % 1.0
    nb = int(np.ceil(P / step))
    idx = np.floor(ph * nb).astype(int)
    sw = np.bincount(idx, wv, minlength=nb); sd = np.bincount(idx, wv * Dg, minlength=nb)
    with np.errstate(invalid='ignore', divide='ignore'):
        dep = sd / sw; s = dep * np.sqrt(sw)
    # each phase bin collects ~1 sample per transit; require >=3 distinct epochs
    k = np.nanargmax(s); snr[i] = s[k]
np.savez(OUT / 'data' / 'period_scan.npz', P=Ps, snr=snr)
order = np.argsort(-snr)
peaks = []
for j in order:
    if all(abs(Ps[j] / p - 1) > 0.002 for p, _ in peaks):
        peaks.append((Ps[j], snr[j]))
    if len(peaks) >= 15: break
print('top peaks (P, SNR):'); [print(f'  {p:.4f}  {s:.2f}') for p, s in peaks]
# --- 3. explicit aliases at the observed epoch
P0, T0 = d['P'], d['t0']
def windows(P, t0, offset=0.0):
    out = []
    n0, n1 = int(np.floor((t.min() - t0) / P)) - 1, int(np.ceil((t.max() - t0) / P)) + 1
    for n in range(n0, n1 + 1):
        tc = t0 + (n + offset) * P
        dep, err = box_depth(tc, dur)
        if np.isfinite(dep):
            cov = (np.abs(t - tc) < 0.5 * dur).sum() * (2 / 1440) / dur
            if cov > 0.4: out.append((n, tc, dep, err, cov))
    return out
alias = {}
for name, P, off in [('P', P0, 0), ('P/2 (new half-epochs)', P0, 0.5), ('P/3 (new third-epochs, +1/3)', P0, 1 / 3), ('P/3 (+2/3)', P0, 2 / 3),
                     ('P/4 (+1/4)', P0, .25), ('P/4 (+3/4)', P0, .75)]:
    w = windows(P, T0, off)
    if not w: alias[name] = dict(n=0); continue
    dd = np.array([x[2] for x in w]); ee = np.array([x[3] for x in w])
    mean = np.sum(dd / ee ** 2) / np.sum(1 / ee ** 2); err = 1 / np.sqrt(np.sum(1 / ee ** 2))
    alias[name] = dict(n=len(w), mean_ppm=mean * 1e6, err_ppm=err * 1e6, per_window=[(int(x[0]), round(x[1], 3), round(x[2] * 1e6), round(x[3] * 1e6)) for x in w])
    print(f'{name:32s} n={len(w):2d}  depth {mean*1e6:7.0f} +- {err*1e6:4.0f} ppm')
# 2P / 3P/2: split observed epochs by parity / mod 3
w = windows(P0, T0, 0)
ep = np.array([x[0] for x in w]); dd = np.array([x[2] for x in w]); ee = np.array([x[3] for x in w])
for m in [2, 3]:
    for r in range(m):
        s = (ep % m) == r
        if s.sum():
            mean = np.sum(dd[s] / ee[s] ** 2) / np.sum(1 / ee[s] ** 2); err = 1 / np.sqrt(np.sum(1 / ee[s] ** 2))
            alias[f'epochs mod {m} = {r}'] = dict(n=int(s.sum()), mean_ppm=mean * 1e6, err_ppm=err * 1e6, epochs=ep[s].tolist())
            print(f'epochs mod {m} = {r}: n={s.sum()} depth {mean*1e6:.0f} +- {err*1e6:.0f}  epochs {ep[s].tolist()}')
json.dump(dict(peaks=[(float(p), float(s)) for p, s in peaks], alias=alias), open(OUT / 'data' / 'alias.json', 'w'), indent=1)
