"""Per-planet MCMC with free circular-orbit density (photo-eccentric / shape test), 2-min PDCSAP."""
from fitlib import *
import json, emcee, sys
store = load_store()
E = json.load(open(OUT / 'data' / 'ephem.json'))
rho0 = E['_rho']['rho']
lc = concat(store, '2min')
pls = sys.argv[1:] or ['d18', 'b', 'c', 'd2']
out = json.load(open(OUT / 'data' / 'mcmc_rhoprior.json')) if (OUT / 'data' / 'mcmc_rhoprior.json').exists() else {}
for pl in pls:
    e = E[pl]
    sel = in_tr(lc['t'], e['P'], e['t0'], e['dur'], widen=6)
    t, f, err = lc['t'][sel], lc['f'][sel].copy(), lc['e'][sel]
    for o in ['b', 'c', 'd18', 'd2']:
        if o == pl: continue
        oe = E[o]; f -= transit_model(t, oe['t0'], oe['P'], oe['k'], oe['b'], rho0) - 1
    def lnp(q):
        t0, P, k, b, lr, f0 = q
        if not (0.003 < k < 0.2 and 0 <= b < 1 + k and -1 < lr < 4 and abs(t0 - e['t0']) < 0.05 and abs(P - e['P']) < 0.002):
            return -np.inf
        m = transit_model(t, t0, P, k, b, np.exp(lr)) * f0
        return -0.5 * np.sum(((f - m) / err) ** 2) - 0.5 * ((np.exp(lr) - 3.86) / 0.36) ** 2
    p0 = np.array([e['t0'], e['P'], e['k'], min(e['b'], 0.8), np.log(rho0), 1.0])
    nw = 32
    P0 = p0 + np.array([1e-4, 1e-6, 1e-3, 0.05, 0.05, 1e-5]) * np.random.randn(nw, 6)
    P0[:, 3] = np.abs(P0[:, 3])
    S = emcee.EnsembleSampler(nw, 6, lnp)
    S.run_mcmc(P0, 7000, progress=False)
    ch = S.get_chain(discard=2000, thin=5, flat=True)
    rho = np.exp(ch[:, 4]); k = ch[:, 2]; b = ch[:, 3]
    ar = a_over_r(ch[:, 1], rho); inc = np.arccos(np.clip(b / ar, 0, 1))
    T14 = ch[:, 1] / np.pi * np.arcsin(np.clip(np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / ar / np.sin(inc), 0, 1)) * 24
    q = lambda x: [float(v) for v in np.percentile(x, [16, 50, 84])]
    res = dict(t0=q(ch[:, 0]), P=q(ch[:, 1]), k=q(k), b=q(b), rho_circ=q(rho), T14_h=q(T14), depth_ppm=q(k ** 2 * 1e6),
               rp_earth=q(k * RSTAR * 109.2), inc_deg=q(np.degrees(inc)), a_rs=q(ar),
               acc=float(np.mean(S.acceptance_fraction)), n=int(len(t)))
    # photo-eccentric: rho_circ/rho_star = ((1+e sin w)/sqrt(1-e^2))^3 ; g = (rho_c/rho_*)^(1/3)
    g = (rho / RHO_TIC) ** (1 / 3)
    res['g_ratio'] = q(g)
    np.save(OUT / 'data' / f'chain_rhoprior_{pl}.npy', ch)
    out[pl] = res
    print(pl, json.dumps({k_: [round(x, 5) for x in v] if isinstance(v, list) else v for k_, v in res.items()}), flush=True)
    json.dump(out, open(OUT / 'data' / 'mcmc_rhoprior.json', 'w'), indent=1)
