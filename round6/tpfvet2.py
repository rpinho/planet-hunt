"""Template-shift localisation: fit each sector's difference image with the sector's own
out-of-transit image (the empirical PRF, which needs no PSF model) shifted by a sky offset.
Sum chi^2 over sectors on a common sky grid; calibrate with planets b and c on the same pixels."""
from common import *
import json, glob
from astropy.wcs import WCS
from scipy.ndimage import shift as ndshift
PLS = ['b', 'c', 'd18', 'd2']
grid = np.arange(-110, 110.1, 2.0)
XI, ETA = np.meshgrid(grid, grid)  # arcsec, xi = +RA (east), eta = +Dec
nb = json.load(open(OUT / 'data' / 'neighbours_local.json')) if (OUT / 'data' / 'neighbours_local.json').exists() else []
def jac(wcs, x, y):
    J = np.zeros((2, 2)); h = 0.05
    a = wcs.all_pix2world([[x, y]], 0)[0]
    for k in range(2):
        p = [x, y]; p[k] += h
        b = wcs.all_pix2world([p], 0)[0]
        J[0, k] = (b[0] - a[0]) * np.cos(np.radians(a[1])) * 3600 / h; J[1, k] = (b[1] - a[1]) * 3600 / h
    return J
out = {}
for pl in PLS:
    tot = np.zeros_like(XI); per = []
    for s in SECTORS:
        p = OUT / 'data' / f'diffimg_{pl}_s{s}.npz'
        if not p.exists(): continue
        z = np.load(p); D, Dir, errD = z['D'], z['Dir'], z['errD']
        fn = glob.glob(str(RAW / f'*-s{s:04d}-*-s_tp.fits'))[0]
        w = WCS(fits.getheader(fn, 2))
        ny, nx = Dir.shape
        yy, xx = np.mgrid[0:ny, 0:nx]
        cy, cx = np.unravel_index(np.argmax(Dir), Dir.shape)
        Jinv = np.linalg.inv(jac(w, cx, cy))
        # use central 7x7 pixels around the star (edges carry little signal)
        sel = np.ones_like(Dir, bool)
        chi = np.zeros_like(XI)
        e = errD[sel]
        for i in range(XI.shape[0]):
            for j in range(XI.shape[1]):
                dx, dy = Jinv @ np.array([XI[i, j], ETA[i, j]])
                T = ndshift(Dir, (dy, dx), order=3, mode='constant', cval=0.0)[sel]
                X = np.vstack([T, np.ones_like(T)]).T / e[:, None]
                y = D[sel] / e
                bta, *_ = np.linalg.lstsq(X, y, rcond=None)
                r = y - X @ bta
                chi[i, j] = r @ r
        i0 = np.unravel_index(np.argmin(chi), chi.shape)
        dof = sel.sum() - 4
        scale = max(1.0, chi[i0] / dof)
        chi_s = chi / scale
        c0 = chi_s[XI.shape[0] // 2, XI.shape[1] // 2]
        per.append(dict(sector=s, best_xi=float(XI[i0]), best_eta=float(ETA[i0]), dchi2_target=float(c0 - chi_s.min()), rchi2=float(chi[i0] / dof)))
        tot += chi_s - chi_s.min()
        np.save(OUT / 'data' / f'locmap_{pl}_s{s}.npy', chi_s - chi_s.min())
    i0 = np.unravel_index(np.argmin(tot), tot.shape)
    d0 = tot[XI.shape[0] // 2, XI.shape[1] // 2] - tot.min()
    inside1 = tot - tot.min() <= 2.3
    r1 = float(np.max(np.hypot(XI[inside1], ETA[inside1])))
    inside3 = tot - tot.min() <= 11.8
    rec = dict(best_xi=float(XI[i0]), best_eta=float(ETA[i0]), dchi2_target=float(d0), radius_1sig=r1,
               max_sep_3sig=float(np.max(np.hypot(XI[inside3], ETA[inside3]))), per_sector=per)
    # neighbour positions
    rec['neighbours'] = []
    for n in nb:
        xi, eta = n['sep'] * np.sin(np.radians(n['pa'])), n['sep'] * np.cos(np.radians(n['pa']))
        ii = np.argmin(np.abs(grid - eta)); jj = np.argmin(np.abs(grid - xi))
        rec['neighbours'].append(dict(source_id=n['source_id'], sep=n['sep'], pa=n['pa'], tmag=n['tmag'], need1150=n['need1150'], need284=n['need284'], dchi2=float(tot[ii, jj] - tot.min())))
    out[pl] = rec
    np.save(OUT / 'data' / f'locmap_{pl}.npy', tot)
    print(pl, f"best ({rec['best_xi']:+.0f},{rec['best_eta']:+.0f})\"  dchi2(target)={d0:.1f}  1sig radius {r1:.0f}\"  3sig max sep {rec['max_sep_3sig']:.0f}\"",
          ' | per-sector dchi2_target:', [round(x['dchi2_target'], 1) for x in per], ' | min neighbour dchi2:', min(x['dchi2'] for x in rec['neighbours']))
json.dump(out, open(OUT / 'data' / 'tpf_template_loc.json', 'w'), indent=1)
