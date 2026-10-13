"""Round 7: diagnostic folds and per-sector depths for the deep-vet targets.

    .venv/bin/python round7/diag.py 119876352 144439305 ...
Writes runs/round7_eumemic_weak/diag/<tic>.png and <tic>.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from planet_hunt import config  # noqa: E402
from planet_hunt import hostsearch as hs  # noqa: E402
from planet_hunt import multisector as ms  # noqa: E402
import revet as rv  # noqa: E402

OUT = config.RUNS / "round7_eumemic_weak" / "diag"


def fold_bin(t, y, P, t0, dur, nb=60, span=3.0):
    ph = ms.phase(t, P, t0)
    m = np.abs(ph) < span * dur
    edges = np.linspace(-span * dur, span * dur, nb + 1)
    idx = np.digitize(ph[m], edges)
    x, yy, ee = [], [], []
    for k in range(1, nb + 1):
        s = idx == k
        if s.sum() < 5:
            continue
        x.append(ph[m][s].mean() * 24); yy.append(np.mean(y[m][s])); ee.append(np.std(y[m][s]) / np.sqrt(s.sum()))
    return np.array(x), (np.array(yy) - 1) * 1e6, np.array(ee) * 1e6


def main(tics):
    import matplotlib.pyplot as plt
    OUT.mkdir(parents=True, exist_ok=True)
    tg = pd.read_csv(rv.OUT / "targets.csv").drop_duplicates("tic").set_index("tic")
    rvt = pd.read_csv(rv.OUT / "revet.csv").set_index("tic")
    plans = {int(k): v for k, v in json.load(open(rv.OUT / "plans.json")).items()}
    import round5_hosts as r5
    for tic in tics:
        row, res = tg.loc[tic], rvt.loc[tic]
        P, t0, dur = float(res.P_full), float(res.t0_full), float(res.dur_h) / 24
        rstar, mstar = float(row.rstar), float(row.mstar)
        secs = rv.load_sectors(tic, plans[tic], int(row.his_cut))
        toi, ps = r5.load_catalogs()
        kl = hs.known_list(r5.known_ephemerides(tic, toi, ps), rstar, mstar)
        kmask = lambda t: hs.mask_known(t, kl)
        cmask = lambda t: kmask(t) | ms.in_tr(t, P, t0, dur, 1.5)
        rot = hs.rotation(secs, kmask)
        lc = hs.detrend(secs, hs.W_SHORT if P < hs.SPLIT_P else hs.W_LONG, cmask, rot=rot)  # candidate masked
        prov = {int(d["meta"]["sector"]): d["meta"].get("provenance") for d in secs}
        his = {int(d["meta"]["sector"]): d["meta"]["his"] for d in secs}
        beta = rv.sector_beta(lc, dur)
        rows = []
        for s in np.unique(lc["sector"]):
            sel = lc["sector"] == s
            d, e, n = rv.depth_at(lc, beta, P, t0, dur, sel)
            rows.append({"sector": int(s), "prov": prov.get(int(s)), "his": his.get(int(s)), "depth_ppm": d * 1e6,
                         "err_ppm": e * 1e6, "n_tr": n, "sigma_ppm": float(np.median(lc["dy"][sel]) * 1e6)})
        per = pd.DataFrame(rows)
        hsel = np.array([his.get(int(s), False) for s in lc["sector"]])
        st_all = rv.full_stats(lc, P, t0, dur, rstar, mstar, float(row.teff))
        # alias checks: depth at 2P (odd vs even) and P/2 (the would-be extra transits)
        d_half, e_half, _ = rv.depth_at(lc, beta, P, t0 + P / 2, dur)
        by_prov = {}
        for pv in per.prov.unique():
            g = per[(per.prov == pv) & per.depth_ppm.notna()]
            w = 1 / g.err_ppm ** 2
            by_prov[pv] = {"depth_ppm": float((g.depth_ppm * w).sum() / w.sum()), "err_ppm": float(1 / np.sqrt(w.sum())),
                           "n_sectors": int(len(g))}
        out = {"tic": int(tic), "P": P, "t0": t0, "dur_h": dur * 24, "stats_masked_detrend": st_all,
               "per_sector": rows, "depth_by_product": by_prov,
               "depth_at_phase_half_ppm": [d_half * 1e6, e_half * 1e6], "rotation": {k: rot.get(k) for k in ("active", "P_rot", "amp_ppm")}}
        (OUT / f"{tic}.json").write_text(json.dumps(out, default=rv._js, indent=1))
        fig, ax = plt.subplots(2, 3, figsize=(15, 8))
        for a, sel, lab in ((ax[0, 0], hsel, "his-era data"), (ax[0, 1], ~hsel, "data he did not have"),
                            (ax[0, 2], np.ones_like(hsel), "all data")):
            x, y, e = fold_bin(lc["t"][sel], lc["flat"][sel], P, t0, dur)
            a.errorbar(x, y, e, fmt="o", ms=3, color="k")
            a.axvspan(-dur * 12, dur * 12, color="C1", alpha=0.15)
            a.axhline(0, color="0.6"); a.axhline(-res.depth_full_ppm, color="C3", ls="--", lw=0.8)
            a.set_title(f"{lab}"); a.set_xlabel("hours from mid-transit"); a.set_ylabel("ppm")
        epoch = np.round((lc["t"] - t0) / P).astype(int)
        for k, c in ((0, "C0"), (1, "C3")):
            m = epoch % 2 == k
            x, y, e = fold_bin(lc["t"][m], lc["flat"][m], P, t0, dur)
            ax[1, 0].errorbar(x, y, e, fmt="o", ms=3, color=c, label="even" if k == 0 else "odd")
        ax[1, 0].legend(); ax[1, 0].set_title(f"odd/even {st_all['oe_sigma']:.1f} sigma"); ax[1, 0].axhline(0, color="0.6")
        x, y, e = fold_bin(lc["t"], lc["flat"], P, t0 + P / 2, dur)
        ax[1, 1].errorbar(x, y, e, fmt="o", ms=3, color="k"); ax[1, 1].axhline(0, color="0.6")
        ax[1, 1].set_title(f"phase 0.5: {d_half * 1e6:.0f} +- {e_half * 1e6:.0f} ppm")
        cols = {"SPOC-2min": "C2", "TESS-SPOC": "C0", "QLP": "C1"}
        for r in rows:
            ax[1, 2].errorbar(r["sector"], r["depth_ppm"], r["err_ppm"], fmt="o" if r["his"] else "s",
                              color=cols.get(r["prov"], "k"))
        ax[1, 2].axhline(0, color="0.6"); ax[1, 2].axhline(res.depth_full_ppm, color="C3", ls="--")
        ax[1, 2].set_title("per sector (o his, s new; green 2-min, blue TESS-SPOC, orange QLP)")
        ax[1, 2].set_xlabel("sector")
        fig.suptitle(f"TIC {tic}  P={P:.5f} d  depth {res.depth_full_ppm:.0f} ppm  SNR_full {res.snr_full:.1f}  "
                     f"z_hold {res.z_hold:.1f}  Teff {row.teff:.0f} R* {rstar:.2f} T={row.tmag:.2f}")
        fig.tight_layout(); fig.savefig(OUT / f"{tic}.png", dpi=90); plt.close(fig)
        print(tic, "by product", {k: (round(v["depth_ppm"]), round(v["err_ppm"]), v["n_sectors"]) for k, v in by_prov.items()},
              "| masked-detrend snr", round(st_all["snr"], 1), "depth", round(st_all["depth_ppm"]), "oe", round(st_all["oe_sigma"], 1),
              "sec", round(st_all["sec_sigma"], 1), "sec_any", round(st_all["sec_any_sigma"], 1), "season p", st_all["season_chi2_p"],
              "rho dur_ratio", round(st_all["dur_ratio"], 2), "rp", round(st_all["rp_earth"], 2), flush=True)


if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]])
