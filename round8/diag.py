"""Round 8, stage 6a: light-curve diagnostics for hold-out survivors (round-7 diag.py, keyed by signal id).

    .venv/bin/python round8/diag.py <signal id> ...      # e.g. E1_447558713_5.5577
Per signal, with the candidate masked in the detrend: folds (their data / new data / all), odd/even, phase 0.5,
per-sector and per-product depths, per-transit depths (chi2), alias depths (P/2, 2P, 3P, 1/2 and 3/2 phase),
a narrow BLS comparison of P vs 2P vs P/2, duration vs the TIC star, the ephemeris from per-sector mid-times.
Writes runs/round8_crowd/diag/<id>.json and .png.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from common import OUT, js
import holdout as H
from holdout import hs, ms, r5, r7

D = OUT / "diag"


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


def per_sector_t0(lc, beta, P, t0, dur, sel):
    """Best mid-time shift (grid +-0.6 dur) for one subset -> (tc, err) via box depth maximisation."""
    shifts = np.linspace(-0.6 * dur, 0.6 * dur, 61)
    dd = []
    for s in shifts:
        d, e, n = r7.depth_at(lc, beta, P, t0 + s, dur, sel)
        dd.append(d / e if e else np.nan)
    dd = np.array(dd)
    if not np.isfinite(dd).any() or np.nanmax(dd) < 3:
        return None
    i = int(np.nanargmax(dd))
    ok = dd > np.nanmax(dd) - 1          # delta(snr^2)... crude 1-sigma width in snr units
    width = (shifts[ok].max() - shifts[ok].min()) / 2 if ok.sum() > 1 else dur / 10
    epoch = np.round((np.median(lc["t"][sel]) - t0) / P)
    return float(t0 + epoch * P + shifts[i]), float(max(width, dur / 30)), int(epoch)


def run(sid: str):
    import matplotlib.pyplot as plt
    D.mkdir(parents=True, exist_ok=True)
    ho = pd.read_csv(OUT / "holdout.csv").set_index("id").loc[sid]
    tg = pd.read_csv(OUT / "targets.csv")
    tg["id"] = [H.sig_id(r) for r in tg.to_dict("records")]
    row = tg.set_index("id").loc[sid].to_dict()
    row["stream"], row["tic"] = sid.split("_")[0], int(sid.split("_")[1])
    tic = row["tic"]
    plans = H.load_plans(tg)
    sig = {"rule": row["rule"], "src_sectors": json.loads(row["src_sectors"]) if str(row["src_sectors"]).startswith("[") else []}
    secs = H.load_sectors(tic, plans[tic], sig)
    P, t0, dur = float(ho.P_full), float(ho.t0_full), float(ho.dur_h) / 24
    rstar, mstar, teff = float(row["rstar"]), float(row["mstar"]), float(row["teff"])
    toi, ps = r5.load_catalogs()
    kl = hs.known_list(r5.known_ephemerides(tic, toi, ps), rstar, mstar)
    kmask = lambda t: hs.mask_known(t, kl)
    cmask = lambda t: kmask(t) | ms.in_tr(t, P, t0, dur, 1.5)
    rot = hs.rotation(secs, kmask)
    lc = hs.detrend(secs, hs.W_SHORT if P < hs.SPLIT_P else hs.W_LONG, cmask, rot=rot)
    prov = {int(d["meta"]["sector"]): d["meta"].get("provenance") for d in secs}
    his = {int(d["meta"]["sector"]): d["meta"]["his"] for d in secs}
    beta = r7.sector_beta(lc, dur)
    rows = []
    for s in np.unique(lc["sector"]):
        sel = lc["sector"] == s
        d, e, n = r7.depth_at(lc, beta, P, t0, dur, sel)
        rows.append({"sector": int(s), "prov": prov.get(int(s)), "his": his.get(int(s)), "depth_ppm": d * 1e6,
                     "err_ppm": e * 1e6, "n_tr": n})
    per = pd.DataFrame(rows)
    st = r7.full_stats(lc, P, t0, dur, rstar, mstar, teff)
    hsel = np.array([his.get(int(s), False) for s in lc["sector"]])
    alias = {}
    for name, (PP, tt) in {"P": (P, t0), "phase0.5": (P, t0 + P / 2), "P/2_extra": (P, t0 + P / 2),
                           "2P_even": (2 * P, t0), "2P_odd": (2 * P, t0 + P), "3P_a": (3 * P, t0),
                           "3P_b": (3 * P, t0 + P), "3P_c": (3 * P, t0 + 2 * P),
                           "phase1/3": (P, t0 + P / 3), "phase2/3": (P, t0 + 2 * P / 3)}.items():
        d, e, n = r7.depth_at(lc, beta, PP, tt, dur)
        alias[name] = [d * 1e6, e * 1e6, n]
    # narrow-BLS power at P, 2P, P/2 on all data
    blsc = {}
    for name, f in (("P", 1), ("2P", 2), ("P/2", 0.5)):
        e = r7.narrow_bls(lc, np.ones(len(lc["t"]), bool), P * f, rstar, mstar, frac=0.002)
        blsc[name] = None if e is None else {"P": e["P"], "snr": e["bls_snr"], "dur_h": e["dur"] * 24}
    by_prov = {}
    for pv in per.prov.dropna().unique():
        g = per[(per.prov == pv) & per.depth_ppm.notna() & (per.err_ppm > 0)]
        if len(g):
            w = 1 / g.err_ppm ** 2
            by_prov[pv] = {"depth_ppm": float((g.depth_ppm * w).sum() / w.sum()), "err_ppm": float(1 / np.sqrt(w.sum())),
                           "n_sectors": int(len(g))}
    # per-sector mid-times -> linear ephemeris
    tcs = []
    for s in np.unique(lc["sector"]):
        r = per_sector_t0(lc, beta, P, t0, dur, lc["sector"] == s)
        if r:
            tcs.append((int(s),) + r)
    eph = None
    if len(tcs) >= 3:
        E = np.array([x[3] for x in tcs], float); T = np.array([x[1] for x in tcs]); W = 1 / np.array([x[2] for x in tcs]) ** 2
        A = np.vstack([np.ones_like(E), E]).T
        cov = np.linalg.inv(A.T @ (A * W[:, None]))
        sol = cov @ (A.T @ (W * T))
        res = T - A @ sol
        chi2 = float(np.sum(res ** 2 * W)); dof = max(len(T) - 2, 1)
        sc = max(1.0, chi2 / dof)
        eph = {"T0": float(sol[0]), "P": float(sol[1]), "eT0": float(np.sqrt(cov[0, 0] * sc)),
               "eP": float(np.sqrt(cov[1, 1] * sc)), "chi2_dof": chi2 / dof, "n": len(T),
               "oc_min": (res * 1440).round(1).tolist()}
    exp_dur = ms.central_duration(P, rstar, mstar) * 24 if hasattr(ms, "central_duration") else np.nan
    out = {"id": sid, "tic": tic, "P": P, "t0": t0, "dur_h": dur * 24, "exp_central_dur_h": float(np.atleast_1d(exp_dur)[0]),
           "stats_masked_detrend": st, "per_sector": rows, "depth_by_product": by_prov, "alias_depths_ppm": alias,
           "bls_compare": blsc, "ephemeris": eph, "rotation": {k: rot.get(k) for k in ("active", "P_rot", "amp_ppm")},
           "star": {"teff": teff, "rstar": rstar, "mstar": mstar, "tmag": row.get("tmag"), "ra": row.get("ra"),
                    "dec": row.get("dec"), "gaia": row.get("GAIA"), "contratio": row.get("contratio")}}
    (D / f"{sid}.json").write_text(json.dumps(out, default=js, indent=1))
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for a, sel, lab in ((ax[0, 0], hsel, "their data"), (ax[0, 1], ~hsel, "new data (hold-out)"),
                        (ax[0, 2], np.ones_like(hsel), "all data")):
        x, y, e = fold_bin(lc["t"][sel], lc["flat"][sel], P, t0, dur)
        a.errorbar(x, y, e, fmt="o", ms=3, color="k")
        a.axvspan(-dur * 12, dur * 12, color="C1", alpha=0.15)
        a.axhline(0, color="0.6"); a.axhline(-ho.depth_full_ppm, color="C3", ls="--", lw=0.8)
        a.set_title(lab); a.set_xlabel("hours from mid-transit"); a.set_ylabel("ppm")
    epoch = np.round((lc["t"] - t0) / P).astype(int)
    for k, c in ((0, "C0"), (1, "C3")):
        m = epoch % 2 == k
        x, y, e = fold_bin(lc["t"][m], lc["flat"][m], P, t0, dur)
        ax[1, 0].errorbar(x, y, e, fmt="o", ms=3, color=c, label="even" if k == 0 else "odd")
    ax[1, 0].legend(); ax[1, 0].set_title(f"odd/even {st['oe_sigma']:.1f} sigma"); ax[1, 0].axhline(0, color="0.6")
    x, y, e = fold_bin(lc["t"], lc["flat"], P, t0 + P / 2, dur)
    ax[1, 1].errorbar(x, y, e, fmt="o", ms=3, color="k"); ax[1, 1].axhline(0, color="0.6")
    ax[1, 1].set_title(f"phase 0.5: {alias['phase0.5'][0]:.0f} +- {alias['phase0.5'][1]:.0f} ppm")
    cols = {"SPOC-2min": "C2", "TESS-SPOC": "C0", "QLP": "C1"}
    for r in rows:
        ax[1, 2].errorbar(r["sector"], r["depth_ppm"], r["err_ppm"], fmt="o" if r["his"] else "s", color=cols.get(r["prov"], "k"))
    ax[1, 2].axhline(0, color="0.6"); ax[1, 2].axhline(ho.depth_full_ppm, color="C3", ls="--")
    ax[1, 2].set_title("per sector (o theirs, s new; green 2-min, blue TESS-SPOC, orange QLP)"); ax[1, 2].set_xlabel("sector")
    fig.suptitle(f"{sid}  TIC {tic}  P={P:.5f} d  depth {ho.depth_full_ppm:.0f} ppm  SNR_full {ho.snr_full:.1f}  "
                 f"z_hold {ho.z_hold:.1f}  Teff {teff:.0f} R* {rstar:.2f} T={row.get('tmag'):.2f}")
    fig.tight_layout(); fig.savefig(D / f"{sid}.png", dpi=80); plt.close(fig)
    print(sid, "| prod", {k: (round(v["depth_ppm"]), round(v["err_ppm"]), v["n_sectors"]) for k, v in by_prov.items()},
          "| snr", round(st["snr"], 1), "depth", round(st["depth_ppm"]), "oe", round(st["oe_sigma"], 1), "sec",
          round(st["sec_sigma"], 1), "sec_any", round(st["sec_any_sigma"], 1), "season p", round(st["season_chi2_p"], 3),
          "tr p", round(st["tr_chi2_p"], 3), "dur_ratio", round(st["dur_ratio"], 2), "rp", round(st["rp_earth"], 2),
          "| bls", {k: (round(v["snr"], 1) if v else None) for k, v in blsc.items()},
          "| ph0.5", [round(x) for x in alias["phase0.5"][:2]], "| eph", None if eph is None else
          (round(eph["P"], 6), f"{eph['eP']:.1e}", round(eph["chi2_dof"], 2)), flush=True)


if __name__ == "__main__":
    for s in sys.argv[1:]:
        try:
            run(s)
        except Exception as exc:
            import traceback
            traceback.print_exc()
