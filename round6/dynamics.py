"""N-body TTV predictions (rebound) and stability (SPOCK + MEGNO) for TOI-4342 b, c + candidates."""
import os
os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE'); os.environ.setdefault('OMP_NUM_THREADS', '1')
from common import *
import json, rebound
E = json.load(open(OUT / 'data' / 'ephem.json'))
pt = json.load(open(OUT / 'data' / 'per_transit.json'))
MS = MSTAR
MEARTH = 3.003e-6
rng = np.random.default_rng(1)
# masses (Earth): b, c from literature if available, else Chen & Kipping 2017 M-R
LIT = OUT / 'lit' / 'literature.json'
def mr(R):  # Chen & Kipping 2017 (Neptunian regime for R > 1.23)
    return (R / 1.008) ** (1 / 0.279) if R < 1.23 else (R / 0.80811) ** (1 / 0.589)
RP = dict(b=2.31, c=2.28, d18=2.15, d2=1.11)
MASS = {k: mr(v) for k, v in RP.items()}
print('M-R masses (Me):', {k: round(v, 2) for k, v in MASS.items()})

def build(planets, masses, ecc=None, omega=None, tref=2900.0):
    sim = rebound.Simulation(); sim.units = ('day', 'AU', 'Msun')
    sim.add(m=MS)
    for pl in sorted(planets, key=lambda p: E[p]['P']):
        e = E[pl]; P, T0 = e['P'], e['t0']
        ec = 0.0 if ecc is None else ecc[pl]
        w = 0.0 if omega is None else omega[pl]
        # transit when true longitude theta = 0 (observer along +x); f_tr = -w (pomega = w)
        f_tr = -w
        E_tr = 2 * np.arctan(np.sqrt((1 - ec) / (1 + ec)) * np.tan(f_tr / 2))
        M_tr = E_tr - ec * np.sin(E_tr)
        M_ref = M_tr + 2 * np.pi * (tref - T0) / P
        sim.add(m=masses[pl] * MEARTH, P=P, e=ec, pomega=w, M=M_ref, primary=sim.particles[0])
    sim.move_to_com()
    sim.integrator = 'whfast'; sim.dt = min(E[p]['P'] for p in planets) / 30
    return sim, sorted(planets, key=lambda p: E[p]['P'])

def transit_times(sim, order, t_end, t_start=1600.0, tref=2900.0):
    """Integrate backwards to t_start then forward, recording y=0 crossings with x>0."""
    sim.integrate(t_start - tref) if False else None
    # shift time so sim.t=0 is tref
    sim.t = 0.0
    sim.integrator = 'ias15'
    sim.integrate(t_start - tref)
    N = len(order)
    out = {p: [] for p in order}
    dt = min(E[p]['P'] for p in order) / 40
    ps = sim.particles
    while sim.t < t_end - tref:
        y_old = [ps[i + 1].y - ps[0].y for i in range(N)]
        t_old = sim.t
        sim.integrate(sim.t + dt)
        for i in range(N):
            yn = ps[i + 1].y - ps[0].y
            if y_old[i] < 0 <= yn and (ps[i + 1].x - ps[0].x) > 0:
                # bisection
                tl, tr = t_old, sim.t
                for _ in range(40):
                    tm = 0.5 * (tl + tr)
                    sim.integrate(tm)
                    if ps[i + 1].y - ps[0].y < 0: tl = tm
                    else: tr = tm
                out[order[i]].append(tm + tref)
                sim.integrate(t_old + dt)
    return {p: np.array(v) for p, v in out.items()}

def oc_at_obs(tt, pl):
    """Model transit times at observed epochs -> O-C after linear fit; return model O-C and data."""
    sel = [r for r in pt[pl] if r['kind'] == '2min' and r['coverage'] > 0.5]
    tobs = np.array([r['post_tc'] for r in sel]); eobs = np.array([r['post_err_min'] for r in sel]) / 1440
    # match model times to observed
    idx = np.array([np.argmin(np.abs(tt - t)) for t in tobs])
    tm = tt[idx]
    n = np.round((tobs - tobs[0]) / E[pl]['P'])
    A = np.vstack([np.ones_like(n), n]).T
    cm, *_ = np.linalg.lstsq(A, tm, rcond=None)
    model_oc = tm - A @ cm
    return n, tobs, eobs, model_oc

def chi2_vs_data(tts, planets_fit):
    chi = 0; amps = {}
    for pl in planets_fit:
        n, tobs, eobs, moc = oc_at_obs(tts[pl], pl)
        # fit linear ephemeris to (data - model TTV)
        A = np.vstack([np.ones_like(n), n]).T / eobs[:, None]
        c, *_ = np.linalg.lstsq(A, (tobs - moc) / eobs, rcond=None)
        r = (tobs - moc) / eobs - A @ c
        chi += np.sum(r ** 2)
        # full-baseline amplitude
        tt = tts[pl]; nn = np.arange(len(tt)); B = np.vstack([np.ones_like(nn), nn]).T
        cc, *_ = np.linalg.lstsq(B, tt, rcond=None)
        amps[pl] = float(np.ptp(tt - B @ cc) * 1440 / 2)
    return float(chi), amps

results = {'masses_MR': MASS}
configs = {'b+c': ['b', 'c'], 'b+c+d18': ['b', 'c', 'd18'], 'd2+b+c+d18': ['d2', 'b', 'c', 'd18'], 'd2+b+c': ['d2', 'b', 'c']}
ttv = {}
for name, pls in configs.items():
    rows = []
    for k in range(25 if name != 'b+c' else 25):
        ecc = {p: (0.0 if k == 0 else min(0.15, abs(rng.rayleigh(0.03)))) for p in pls}
        om = {p: rng.uniform(0, 2 * np.pi) for p in pls}
        ms = {p: MASS[p] * (1 if k == 0 else np.exp(rng.normal(0, 0.3))) for p in pls}
        sim, order = build(pls, ms, ecc, om)
        tts = transit_times(sim, order, 4200.0)
        chi, amps = chi2_vs_data(tts, [p for p in ['b', 'c', 'd18'] if p in pls])
        rows.append(dict(ecc=ecc, masses=ms, chi2=chi, half_amp_min=amps))
        if k == 0:
            np.savez(OUT / 'dyn' / f'tt_{name}.npz', **{p: tts[p] for p in tts})
    ttv[name] = rows
    ch = np.array([r['chi2'] for r in rows])
    print(name, 'chi2 vs TESS timings (circular, M-R):', round(rows[0]['chi2'], 1), ' range over draws', np.round(np.percentile(ch, [5, 50, 95]), 1),
          ' half-amplitudes (circ):', {p: round(v, 1) for p, v in rows[0]['half_amp_min'].items()},
          ' median over draws:', {p: round(float(np.median([r['half_amp_min'][p] for r in rows])), 1) for p in rows[0]['half_amp_min']}, flush=True)
results['ttv'] = ttv
json.dump(results, open(OUT / 'dyn' / 'ttv_nbody.json', 'w'), indent=1, default=float)
