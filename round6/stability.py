"""Stability: SPOCK (FeatureClassifier, 1e9 inner orbits) over random draws + MEGNO N-body checks."""
import os
os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE'); os.environ.setdefault('OMP_NUM_THREADS', '1')
from common import *
import json, rebound, time
from spock import FeatureClassifier
E = json.load(open(OUT / 'data' / 'ephem.json'))
MEARTH = 3.003e-6
def mr(R): return (R / 1.008) ** (1 / 0.279) if R < 1.23 else (R / 0.80811) ** (1 / 0.589)
RP = dict(b=2.31, c=2.28, d18=2.15, d2=1.11)
rng = np.random.default_rng(7)
model = FeatureClassifier()
def sim_for(pls, draw):
    sim = rebound.Simulation(); sim.units = ('day', 'AU', 'Msun'); sim.add(m=MSTAR)
    for pl in sorted(pls, key=lambda p: E[p]['P']):
        m = mr(RP[pl]) * (np.exp(rng.normal(0, 0.35)) if draw else 1)
        e = min(0.3, abs(rng.rayleigh(0.035))) if draw else 0.01
        inc = np.radians(rng.normal(0, 1.0)) if draw else 0.0
        sim.add(m=m * MEARTH, P=E[pl]['P'], e=e, pomega=rng.uniform(0, 2 * np.pi), M=rng.uniform(0, 2 * np.pi),
                inc=abs(inc), Omega=rng.uniform(0, 2 * np.pi))
    sim.move_to_com(); return sim
out = {}
for name, pls in {'b+c+d18': ['b', 'c', 'd18'], 'd2+b+c': ['d2', 'b', 'c'], 'd2+b+c+d18': ['d2', 'b', 'c', 'd18']}.items():
    ps = []
    for k in range(150):
        sim = sim_for(pls, draw=k > 0)
        try:
            ps.append(float(model.predict_stable(sim)))
        except Exception as ex:
            ps.append(float('nan'))
    ps = np.array(ps)
    out[name] = dict(p_nominal=float(ps[0]), p_median=float(np.nanmedian(ps)), p_mean=float(np.nanmean(ps)), p16_84=[float(x) for x in np.nanpercentile(ps, [16, 84])], frac_gt_0p5=float(np.nanmean(ps > 0.5)), n=int(np.isfinite(ps).sum()))
    print(name, out[name], flush=True)
# Hill spacing
def hill(p1, p2):
    a = lambda P: ((P / 365.25) ** 2 * MSTAR) ** (1 / 3)
    m1, m2 = mr(RP[p1]) * MEARTH, mr(RP[p2]) * MEARTH
    a1, a2 = a(E[p1]['P']), a(E[p2]['P'])
    RH = ((m1 + m2) / (3 * MSTAR)) ** (1 / 3) * (a1 + a2) / 2
    return (a2 - a1) / RH
out['mutual_hill'] = {f'{a}-{b}': float(hill(a, b)) for a, b in [('d2', 'b'), ('b', 'c'), ('c', 'd18')]}
print(out['mutual_hill'])
# MEGNO: 4-planet nominal + 8 draws, 2e5 orbits of the inner planet (~1e3 yr)
megno = []
for k in range(9):
    sim = sim_for(['d2', 'b', 'c', 'd18'], draw=k > 0)
    sim.integrator = 'whfast'; sim.dt = E['d2']['P'] / 25
    sim.init_megno()
    sim.exit_max_distance = 1.0
    t0 = time.time()
    try:
        sim.integrate(2e5 * E['d2']['P'], exact_finish_time=0)
        Y = sim.megno(); stat = 'ok'
    except rebound.Escape:
        Y = float('nan'); stat = 'escape'
    megno.append(dict(draw=k, megno=float(Y), status=stat, secs=time.time() - t0))
    print('MEGNO', megno[-1], flush=True)
out['megno_4pl'] = megno
json.dump(out, open(OUT / 'dyn' / 'stability.json', 'w'), indent=1)
