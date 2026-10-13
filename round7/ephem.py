"""Round 7: per-sector transit times -> linear ephemeris with errors (box shape fixed).
    .venv/bin/python round7/ephem.py <tic>   -> runs/round7_eumemic_weak/diag/<tic>_ephem.json"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent.parent)); sys.path.insert(0, str(Path(__file__).resolve().parent))
from planet_hunt import hostsearch as hs, multisector as ms
import revet as rv
import round5_hosts as r5

tic = int(sys.argv[1])
tg = pd.read_csv(rv.OUT / "targets.csv").drop_duplicates("tic").set_index("tic").loc[tic]
res = pd.read_csv(rv.OUT / "revet.csv").set_index("tic").loc[tic]
plans = {int(k): v for k, v in json.load(open(rv.OUT / "plans.json")).items()}
P, t0, dur = float(res.P_full), float(res.t0_full), float(res.dur_h) / 24
secs = rv.load_sectors(tic, plans[tic], int(tg.his_cut))
lc = hs.detrend(secs, hs.W_SHORT, lambda t: ms.in_tr(t, P, t0, dur, 1.6), hs.rotation(secs, lambda t: np.zeros(len(t), bool)))
beta = rv.sector_beta(lc, dur)
rows = []
for s in np.unique(lc["sector"]):
    sel = lc["sector"] == s
    t, y, w = lc["t"][sel], lc["flat"][sel], 1 / (lc["dy"][sel] * beta[sel]) ** 2
    n0 = np.round((np.median(t) - t0) / P)
    tc0 = t0 + n0 * P
    shifts = np.linspace(-0.08, 0.08, 321)
    chi = []
    for dt in shifts:
        ph = np.abs(ms.phase(t, P, tc0 + dt))
        intr = ph < 0.5 * dur
        oot = (ph > dur) & (ph < 3 * dur)
        lev = np.sum(w[oot] * y[oot]) / w[oot].sum()
        d = lev - np.sum(w[intr] * y[intr]) / w[intr].sum()
        chi.append(-(d ** 2) * w[intr].sum())          # -delta chi2 of a box of free depth
    chi = np.array(chi) - np.min(chi)
    i = int(np.argmin(chi))
    ok = chi <= 1.0
    lo, hi = shifts[ok].min(), shifts[ok].max()
    rows.append({"sector": int(s), "epoch": int(n0), "tc": float(tc0 + shifts[i]), "err": float(max((hi - lo) / 2, 0.002)),
                 "edge": bool(i in (0, len(shifts) - 1))})
d = pd.DataFrame(rows)
d = d[~d.edge]
A = np.column_stack([np.ones(len(d)), d.epoch]); W = np.diag(1 / d.err ** 2)
cov = np.linalg.inv(A.T @ W @ A); p = cov @ A.T @ W @ d.tc.values
resid = d.tc.values - A @ p; chi2 = float(np.sum(resid ** 2 / d.err ** 2))
scale = max(1.0, chi2 / max(len(d) - 2, 1))   # timing errors floored at 3 min: inflate to chi2/dof = 1
cov = cov * scale
# reference epoch at the information-weighted centre
nref = int(np.round(np.sum(d.epoch / d.err ** 2) / np.sum(1 / d.err ** 2)))
T_ref = p[0] + nref * p[1]; eT = float(np.sqrt(cov[0, 0] + nref ** 2 * cov[1, 1] + 2 * nref * cov[0, 1]))
out = {"tic": tic, "P": float(p[1]), "P_err": float(np.sqrt(cov[1, 1])), "T0_btjd": float(T_ref), "T0_err": eT,
       "T0_bjd": float(T_ref + 2457000), "chi2": chi2, "cov_scale": scale, "dof": int(len(d) - 2), "per_sector": rows}
# propagate to 2026-12-01 (BTJD ~ 4071) and 2027-06-01
for lab, tb in (("2026-12-01", 2461375.5 - 2457000), ("2027-06-01", 2461557.5 - 2457000)):
    n = np.round((tb - p[0]) / p[1]); var = cov[0, 0] + n ** 2 * cov[1, 1] + 2 * n * cov[0, 1]
    out[f"sigma_tc_{lab}_min"] = float(np.sqrt(var) * 1440)
Path(rv.OUT / "diag" / f"{tic}_ephem.json").write_text(json.dumps(out, indent=1))
print({k: v for k, v in out.items() if k != "per_sector"})
print(d.round(4).to_string())
