"""Round 8A: SPOCK (Tamayo et al. 2020) FeatureClassifier stability of TOI-4405 with the 3.60-d
candidate added to TOI-4405.01 (5.92 d) and .02 (1.09 d). Masses from Chen & Kipping (2017);
150 draws (mass scatter 0.35 dex-e, Rayleigh e ~ 0.035, inclination scatter 1 deg) + nominal.
    .r6venv/bin/python round8/stability.py
"""
import os
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import json
import numpy as np
import rebound
from spock import FeatureClassifier

MSTAR, RSTAR = 0.68, 0.687
MEARTH = 3.003e-6
PL = {".02": dict(P=1.09292034, depth=377e-6), "cand": dict(P=3.60082, depth=326e-6), ".01": dict(P=5.92142305, depth=1222e-6)}
for v in PL.values():
    v["R"] = np.sqrt(v["depth"]) * RSTAR * 109.076


def mr(R):
    return (R / 1.008) ** (1 / 0.279) if R < 1.23 else (R / 0.80811) ** (1 / 0.589)


rng = np.random.default_rng(4405)
model = FeatureClassifier()


def sim_for(names, draw):
    sim = rebound.Simulation(); sim.units = ("day", "AU", "Msun"); sim.add(m=MSTAR)
    for n in sorted(names, key=lambda k: PL[k]["P"]):
        m = mr(PL[n]["R"]) * (np.exp(rng.normal(0, 0.35)) if draw else 1)
        e = min(0.3, abs(rng.rayleigh(0.035))) if draw else 0.01
        inc = np.radians(rng.normal(0, 1.0)) if draw else 0.0
        sim.add(m=m * MEARTH, P=PL[n]["P"], e=e, pomega=rng.uniform(0, 2 * np.pi), M=rng.uniform(0, 2 * np.pi),
                inc=abs(inc), Omega=rng.uniform(0, 2 * np.pi))
    sim.move_to_com()
    return sim


out = {"radii": {k: v["R"] for k, v in PL.items()}, "masses": {k: mr(v["R"]) for k, v in PL.items()}}
names = [".02", "cand", ".01"]
ps = []
for k in range(151):
    try:
        ps.append(float(model.predict_stable(sim_for(names, k > 0))))
    except Exception:
        ps.append(float("nan"))
ps = np.array(ps)
out["spock"] = dict(p_nominal=float(ps[0]), p_median=float(np.nanmedian(ps)), p16_84=[float(x) for x in np.nanpercentile(ps, [16, 84])],
                    frac_gt_0p5=float(np.nanmean(ps > 0.5)), n=int(np.isfinite(ps).sum()))


def hill(a_, b_):
    a = lambda P: ((P / 365.25) ** 2 * MSTAR) ** (1 / 3)
    m1, m2 = mr(PL[a_]["R"]) * MEARTH, mr(PL[b_]["R"]) * MEARTH
    a1, a2 = a(PL[a_]["P"]), a(PL[b_]["P"])
    return (a2 - a1) / (((m1 + m2) / (3 * MSTAR)) ** (1 / 3) * (a1 + a2) / 2)


out["mutual_hill"] = {".02-cand": hill(".02", "cand"), "cand-.01": hill("cand", ".01")}
out["period_ratios"] = {"cand/.02": PL["cand"]["P"] / PL[".02"]["P"], ".01/cand": PL[".01"]["P"] / PL["cand"]["P"]}
print(json.dumps(out, indent=1))
json.dump(out, open("runs/round8_p1/toi4405_stability.json", "w"), indent=1)
