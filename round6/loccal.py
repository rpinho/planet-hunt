"""Calibrate per-sector template localisation with b+c (same pixels, same sector)."""
from common import *
import json
grid = np.arange(-110, 110.1, 2.0); n = len(grid); c0 = n // 2
nb = json.load(open(OUT / 'data' / 'neighbours_local.json'))
res = {}
def roll(m, di, dj):
    out = np.full_like(m, np.nan)
    src = m[max(0, di):n + min(0, di), max(0, dj):n + min(0, dj)]
    out[max(0, -di):n - max(0, di), max(0, -dj):n - max(0, dj)] = src
    return out
for pl, refs in [('d18', ['b', 'c']), ('d2', ['b', 'c']), ('b', ['c']), ('c', ['b'])]:
    tot = np.zeros((n, n)); offs = []
    for s in SECTORS:
        p = OUT / 'data' / f'locmap_{pl}_s{s}.npy'
        if not p.exists(): continue
        m = np.load(p)
        ref = sum(np.load(OUT / 'data' / f'locmap_{r}_s{s}.npy') for r in refs)
        i, j = np.unravel_index(np.argmin(ref), ref.shape)
        di, dj = i - c0, j - c0
        offs.append((s, float(grid[j]), float(grid[i])))
        mm = roll(m, di, dj)  # move reference best position to origin
        tot += np.nan_to_num(mm, nan=np.nanmax(m))
    tot -= tot.min()
    i0 = np.unravel_index(np.argmin(tot), tot.shape)
    XI, ETA = np.meshgrid(grid, grid)
    d0 = tot[c0, c0]
    r3 = float(np.max(np.hypot(XI[tot <= 11.8], ETA[tot <= 11.8])))
    nbr = []
    for x in nb:
        xi, eta = x['sep'] * np.sin(np.radians(x['pa'])), x['sep'] * np.cos(np.radians(x['pa']))
        ii = np.argmin(np.abs(grid - eta)); jj = np.argmin(np.abs(grid - xi))
        nbr.append(dict(source_id=x['source_id'], sep=x['sep'], tmag=x['tmag'], need1150=x['need1150'], need284=x['need284'], dchi2=float(tot[ii, jj])))
    res[pl] = dict(ref=refs, ref_offsets=offs, best_xi=float(XI[i0]), best_eta=float(ETA[i0]), dchi2_target=float(d0), max_sep_3sig=r3, neighbours=nbr)
    np.save(OUT / 'data' / f'locmap_cal_{pl}.npy', tot)
    mim = [x for x in nbr if (x['need1150'] if pl != 'd2' else x['need284']) < 1]
    print(pl, 'calibrated by', refs, f"best ({XI[i0]:+.0f},{ETA[i0]:+.0f}) dchi2(target) {d0:.1f}  3sig max sep {r3:.0f}\"",
          ' min dchi2 over mimic-capable neighbours:', round(min(x['dchi2'] for x in mim), 1) if mim else None, 'n_mimic', len(mim))
json.dump(res, open(OUT / 'data' / 'tpf_template_loc_cal.json', 'w'), indent=1)
