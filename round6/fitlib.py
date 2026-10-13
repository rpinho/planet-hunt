from common import *
import pickle
from scipy.optimize import least_squares, minimize_scalar

def load_store():
    return pickle.load(open(OUT / 'data' / 'detrended.pkl', 'rb'))

def concat(store, kind='2min', sectors=SECTORS):
    T, F, E, S = [], [], [], []
    for s in sectors:
        d = store.get((kind, s))
        if d is None: continue
        T.append(d['t']); F.append(d['fd']); E.append(np.full(len(d['t']), d['mad'])); S.append(np.full(len(d['t']), s))
    return dict(t=np.concatenate(T), f=np.concatenate(F), e=np.concatenate(E), s=np.concatenate(S))

PL = ['b', 'c', 'd18', 'd2']

def window_sel(t, e, half=2.5):
    return in_tr(t, e['P'], e['t0'], e['dur'], widen=2 * half)

def joint_model(t, theta, planets, rho, u=U):
    m = np.ones_like(t)
    for i, pl in enumerate(planets):
        t0, P, k, b = theta[4 * i: 4 * i + 4]
        m += transit_model(t, t0, P, k, b, rho, u) - 1
    return m

def fit_joint(lc, eph, planets, half=2.5):
    sel = np.zeros(len(lc['t']), bool)
    for pl in planets: sel |= window_sel(lc['t'], eph[pl], half)
    t, f, e = lc['t'][sel], lc['f'][sel], lc['e'][sel]
    x0 = []
    for pl in planets:
        ep = eph[pl]
        k0 = np.sqrt(ep.get('depth', 1e-3))
        x0 += [ep['t0'], ep['P'], ep.get('k', 0.035), ep.get('b', 0.3)]
    x0 += [RHO_TIC]
    lo, hi = [], []
    for pl in planets:
        ep = eph[pl]
        lo += [ep['t0'] - 0.1, ep['P'] - 0.01, 0.003, 0.0]; hi += [ep['t0'] + 0.1, ep['P'] + 0.01, 0.2, 1.0]
    lo += [0.5]; hi += [30.]
    x0 = np.clip(x0, np.array(lo) + 1e-9, np.array(hi) - 1e-9)
    res = least_squares(lambda q: (f - joint_model(t, q[:-1], planets, q[-1])) / e, x0, bounds=(lo, hi), x_scale='jac')
    return res, (t, f, e)

def local_fit(t, f, e, tc_grid, shape, P, rho, depth_free=True):
    """For each trial mid-time: linear baseline x (scale * transit). Returns chi2 array, best params.
    Model: f = c0 + c1*(t-tm) + A*(m-1), linear in (c0,c1,A); m from fixed shape."""
    k, b = shape
    tm = np.mean(t)
    chi = np.empty(len(tc_grid)); amps = np.empty(len(tc_grid)); aerr = np.empty(len(tc_grid))
    w = 1 / e
    for i, tc in enumerate(tc_grid):
        m = transit_model(t, tc, P, k, b, rho) - 1
        X = np.vstack([np.ones_like(t), t - tm, m]).T if depth_free else np.vstack([np.ones_like(t), t - tm]).T
        y = f - (0 if depth_free else m)
        Xw, yw = X * w[:, None], y * w
        beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        r = yw - Xw @ beta
        chi[i] = r @ r
        if depth_free:
            C = np.linalg.pinv(Xw.T @ Xw)
            amps[i] = beta[2]; aerr[i] = np.sqrt(C[2, 2])
        else:
            amps[i] = 1; aerr[i] = 0
    return chi, amps, aerr

def red_beta(t, r, e, dur, bins=(10, 20, 30, 45, 60)):
    """Time-averaging beta on residuals (minutes bins up to ~duration)."""
    betas = []
    for bm in bins:
        bw = bm / 1440
        if bw > dur * 1.2: break
        idx = np.floor((t - t.min()) / bw).astype(int)
        nb = np.bincount(idx); sb = np.bincount(idx, r)
        good = nb >= max(3, 0.6 * np.median(nb[nb > 0]))
        if good.sum() < 8: continue
        mb = sb[good] / nb[good]
        n = np.median(nb[good])
        exp = np.std(r) / np.sqrt(n) * np.sqrt(good.sum() / (good.sum() - 1))
        betas.append(np.std(mb) / exp)
    return max(1.0, np.median(betas)) if betas else 1.0
