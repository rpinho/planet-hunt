"""Round 8A: fold of discovery vs hold-out data at the fixed discovery ephemeris.
    .r6venv/bin/python round8/plot_fold.py toi4405"""
import sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
key = sys.argv[1]
z = np.load(f"runs/round8_p1/{key}_fold.npz", allow_pickle=True)
r = json.load(open(f"runs/round8_p1/{key}_lcgate.json"))
ph, fl, er, role = z["phase_d"] * 24, z["flux"], z["err"], z["role"]
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
for a, sel, lab, d in ((ax[0], role == 1, "discovery data (SPOC 2-min)", r["disc"]), (ax[1], role == 0, "hold-out data (never searched)", r["hold"])):
    edges = np.arange(-4 * r["cand"]["dur_h"], 4 * r["cand"]["dur_h"] + 1e-9, r["cand"]["dur_h"] / 4)
    idx = np.digitize(ph[sel], edges)
    w = 1 / er[sel] ** 2
    xb = [np.average(ph[sel][idx == k], weights=w[idx == k]) for k in np.unique(idx) if (idx == k).sum() > 5]
    yb = [np.average(fl[sel][idx == k], weights=w[idx == k]) for k in np.unique(idx) if (idx == k).sum() > 5]
    eb = [1 / np.sqrt(w[idx == k].sum()) for k in np.unique(idx) if (idx == k).sum() > 5]
    a.errorbar(xb, (np.array(yb) - 1) * 1e6, np.array(eb) * 1e6, fmt="o", ms=4, color="k")
    a.axvspan(-r["cand"]["dur_h"] / 2, r["cand"]["dur_h"] / 2, color="0.9")
    a.axhline(0, color="0.6", lw=0.8)
    a.set_title(f"{lab}\n{d['depth_ppm']:.0f} ± {d['err_ppm']:.0f} ppm (z = {d['z']:.1f})", fontsize=9)
    a.set_xlabel("hours from predicted mid-transit")
ax[0].set_ylabel("relative flux (ppm)")
fig.suptitle(f"{r['cand']['name']}  P = {r['cand']['P']:.5f} d (ephemeris fixed before the hold-out data were loaded)", fontsize=10)
fig.tight_layout()
fig.savefig(f"runs/round8_p1/{key}_holdout_fold.png", dpi=120)
print("wrote", f"runs/round8_p1/{key}_holdout_fold.png")
