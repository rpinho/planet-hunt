"""Round 5 track A: deep vetting of triage survivors.

    .venv/bin/python round5_deep.py specs [--include TIC:P ...]   # vetting/round5_targets.json (+ on-target controls)
    .venv/bin/python deep_vet.py vetting/round5_targets.json --out runs/round5_hosts/deep
    .venv/bin/python round5_deep.py folds                          # fold + aperture npz per candidate
    .trivenv/bin/python run_triceratops.py ...                     # (command lines printed by `folds`)
    .venv/bin/python round5_deep.py plots                          # runs/round5_hosts/<tic>_<P>.png

For every survivor: a pixel-localisation spec (TESScut, Gaia DR3 scene, joint
PSF fit; known-planet transits excluded from the difference images), plus the
star's best known planet as an on-target control on the same pixels; the
folded 2-min light curve and the SPOC optimal apertures for TRICERATOPS.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from planet_hunt import config
from planet_hunt import hostsearch as hs
from planet_hunt import multisector as ms
import round5_hosts as r5

OUT = r5.OUT
DEEP = OUT / "deep"
SPEC = config.ROOT / "vetting" / "round5_targets.json"


def candidates(include=()):
    tri = pd.read_csv(OUT / "triage.csv")
    sel = tri[tri.status == "survivor"].copy()
    for item in include:
        tic, P = item.split(":")
        m = (tri.tic == int(tic)) & (np.abs(tri.period / float(P) - 1) < 0.002)
        sel = pd.concat([sel, tri[m]])
    return sel.drop_duplicates(["tic", "period"])


def _name(tic, P):
    return f"tic{int(tic)}_{P:.3f}".replace(".", "p")


def cmd_specs(a):
    sel = candidates(a.include)
    tg = pd.read_csv(OUT / "targets.csv")
    specs = []
    for s in sel.itertuples():
        row = tg[tg.tic == s.tic].iloc[0]
        kl = hs.known_list(json.loads(row["known"]), row["rstar"], row["mstar"])
        excl = [[k["period"], k["t0"], k["dur"] * 1.5] for k in kl if k["masked"]]
        base = {"tic": int(s.tic), "ra": float(row.ra), "dec": float(row.dec), "gaia_id": None,
                "tmag": float(row.tmag), "teff": float(row.teff) if np.isfinite(row.teff) else 3500.0,
                "rstar": float(row.rstar), "mstar": float(row.mstar), "max_sectors": a.max_sectors}
        specs.append(base | {"name": _name(s.tic, s.period), "role": "candidate", "period": float(s.period),
                             "t0": float(s.t0), "duration_h": float(s.duration_h), "depth_ppm": float(s.depth_ppm),
                             "exclude_windows": excl})
        # on-target control: the star's deepest masked known planet
        ctl = [k for k in json.loads(row["known"]) if k.get("t0") and np.isfinite(k.get("t0") or np.nan)
               and np.isfinite(k.get("duration_h") or np.nan) and k.get("disp") not in ("FP", "FA")]
        if ctl:
            k = max(ctl, key=lambda k: k.get("depth_ppm") or 0)
            P, dur = float(k["period"]), float(k["duration_h"]) / 24
            excl_c = [[P_, t_, d_] for (P_, t_, d_) in excl if abs(P_ / P - 1) > 1e-3]
            excl_c.append([float(s.period), float(s.t0), float(s.duration_h) / 24 * 1.5])
            specs.append(base | {"name": _name(s.tic, P) + "_ctl", "role": "control_on", "period": P,
                                 "t0": float(k["t0"]), "duration_h": float(k["duration_h"]),
                                 "depth_ppm": float(k.get("depth_ppm") or 1000), "exclude_windows": excl_c,
                                 "control_of": _name(s.tic, s.period), "known_name": k["name"]})
    SPEC.parent.mkdir(exist_ok=True)
    SPEC.write_text(json.dumps({"targets": specs}, indent=1))
    print(f"{len(specs)} specs ({len(sel)} candidates) -> {SPEC}")


def fold_lc(tic, P, t0, dur, known_json, rstar, mstar, teff):
    secs = r5.load_lc(int(tic))
    kl = hs.known_list(json.loads(known_json), rstar, mstar)
    km = lambda t: hs.mask_known(t, kl, 1.5)
    rot = hs.rotation(secs, km)
    W = hs.W_SHORT if P < hs.SPLIT_P else hs.W_LONG
    lc = hs.detrend(secs, W, lambda t: km(t) | ms.in_tr(t, P, t0, dur, 1.6), rot)
    lc = hs.subset(lc, ~km(lc["t"]))
    return secs, lc, rot


def cmd_folds(a):
    sel = candidates(a.include)
    tg = pd.read_csv(OUT / "targets.csv")
    DEEP.mkdir(parents=True, exist_ok=True)
    for s in sel.itertuples():
        row = tg[tg.tic == s.tic].iloc[0]
        dur = s.duration_h / 24
        secs, lc, rot = fold_lc(s.tic, s.period, s.t0, dur, row["known"], row.rstar, row.mstar, row.teff)
        ph = ms.phase(lc["t"], s.period, s.t0)
        near = np.abs(ph) < max(3 * dur, 0.25)
        name = _name(s.tic, s.period)
        np.savez(DEEP / f"{name}_fold.npz", phase_d=ph[near], flux=lc["flat"][near], err=lc["dy"][near],
                 period=s.period, depth=s.depth_ppm * 1e-6, sectors=np.array(sorted(set(lc["sector"].astype(int)))))
        np.savez(DEEP / f"{name}_apertures.npz", **{f"s{d['meta']['sector']}": d["aperture"] for d in secs})
        print(f".trivenv/bin/python run_triceratops.py --tic {int(s.tic)} --fold {DEEP / (name + '_fold.npz')} "
              f"--lc-dir {DEEP / (name + '_apertures.npz')} --runs 3 --out {DEEP} --name triceratops_{name} "
              f"--local-cutouts {config.CACHE / 'tesscut' / f'{row.ra:.5f}_{row.dec:+.5f}'}")


def cmd_plots(a):
    import matplotlib.pyplot as plt
    tri = pd.read_csv(OUT / "triage.csv")
    tg = pd.read_csv(OUT / "targets.csv")
    picks = [p.split(":") for p in a.pick]
    for tic, P in picks:
        s = tri[(tri.tic == int(tic)) & (np.abs(tri.period / float(P) - 1) < 0.002)].iloc[0]
        row = tg[tg.tic == s.tic].iloc[0]
        dur = s.duration_h / 24
        secs, lc, rot = fold_lc(s.tic, s.period, s.t0, dur, row["known"], row.rstar, row.mstar, row.teff)
        res = json.load(open(OUT / "search" / f"{int(s.tic)}.json"))
        fig = plt.figure(figsize=(12, 8))
        gs = fig.add_gridspec(3, 3)
        ph = ms.phase(lc["t"], s.period, s.t0) * 24
        near = np.abs(ph) < 6 * s.duration_h
        ax = fig.add_subplot(gs[0, :2])
        ax.plot(ph[near], (lc["flat"][near] - 1) * 1e6, ",", c="0.75")
        b = np.arange(-6 * s.duration_h, 6 * s.duration_h + 1e-9, max(s.duration_h / 6, 0.1))
        idx = np.digitize(ph[near], b)
        w = 1 / lc["dy"][near] ** 2
        xb = [np.average(ph[near][idx == k]) for k in np.unique(idx)]
        yb = [np.average((lc["flat"][near][idx == k] - 1) * 1e6, weights=w[idx == k]) for k in np.unique(idx)]
        eb = [1e6 / np.sqrt(w[idx == k].sum()) for k in np.unique(idx)]
        ax.errorbar(xb, yb, eb, fmt="o", c="#1f6fb4", ms=4)
        ax.set_ylim(-4 * s.depth_ppm, 3 * s.depth_ppm)
        ax.axhline(0, c="k", lw=0.5); ax.axhline(-s.depth_ppm, c="r", lw=0.5, ls=":")
        ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("ppm")
        ax.set_title(f"TIC {int(s.tic)} (TOI-{int(s.host_toi) if s.host_toi == s.host_toi else '-'} host): "
                     f"P {s.period:.5f} d, {s.depth_ppm:.0f} ppm, SNR {s.snr:.1f}, {int(s.n_transits)} transits, "
                     f"Rp {s.rp_earth:.2f} Re", fontsize=10)
        # odd / even
        ax2 = fig.add_subplot(gs[0, 2])
        ep = np.round((lc["t"] - s.t0) / s.period).astype(int)
        for par, col in ((0, "#1f6fb4"), (1, "#d95f02")):
            m = near & (ep % 2 == par)
            idx2 = np.digitize(ph[m], b)
            w2 = 1 / lc["dy"][m] ** 2
            ax2.plot([np.average(ph[m][idx2 == k]) for k in np.unique(idx2)],
                     [np.average((lc["flat"][m][idx2 == k] - 1) * 1e6, weights=w2[idx2 == k]) for k in np.unique(idx2)],
                     "o-", ms=3, c=col, label="even" if par == 0 else "odd")
        ax2.legend(fontsize=8, frameon=False); ax2.set_title(f"odd/even {s.oe_sigma:.1f} sigma", fontsize=9)
        ax2.set_ylim(-4 * s.depth_ppm, 3 * s.depth_ppm)
        # per-sector depths
        ax3 = fig.add_subplot(gs[1, :2])
        psd = json.loads(s.per_sector)
        ax3.errorbar([p["sector"] for p in psd], [p["depth_ppm"] for p in psd], [p["err_ppm"] for p in psd],
                     fmt="o", c="#1f6fb4", ms=4)
        ax3.axhline(s.depth_ppm, c="r", lw=0.7); ax3.axhline(0, c="k", lw=0.5)
        ax3.set_xlabel("sector"); ax3.set_ylabel("depth (ppm)")
        ax3.set_title(f"per-sector depth (season p={s.season_chi2_p:.2g}, crowding slope "
                      f"{s.crowd_slope_sigma:.1f} sigma)", fontsize=9)
        # phase 0.5
        ax4 = fig.add_subplot(gs[1, 2])
        ph2 = ms.phase(lc["t"], s.period, s.t0 + 0.5 * s.period) * 24
        n2 = np.abs(ph2) < 6 * s.duration_h
        idx4 = np.digitize(ph2[n2], b)
        w4 = 1 / lc["dy"][n2] ** 2
        ax4.plot([np.average(ph2[n2][idx4 == k]) for k in np.unique(idx4)],
                 [np.average((lc["flat"][n2][idx4 == k] - 1) * 1e6, weights=w4[idx4 == k]) for k in np.unique(idx4)],
                 "o", ms=3, c="0.3")
        ax4.set_ylim(-4 * s.depth_ppm, 3 * s.depth_ppm); ax4.set_title(f"phase 0.5 (sec {s.sec_sigma:.1f} sigma)", fontsize=9)
        # periodogram peaks: per-transit depths vs time
        ax5 = fig.add_subplot(gs[2, :])
        intr = ms.in_tr(lc["t"], s.period, s.t0, dur)
        oot = ~ms.in_tr(lc["t"], s.period, s.t0, dur, 3)
        lev = np.average(lc["flat"][oot], weights=1 / lc["dy"][oot] ** 2)
        es, ds, de = [], [], []
        for e in np.unique(ep[intr]):
            m = intr & (ep == e)
            ww = 1 / lc["dy"][m] ** 2
            if m.sum() >= 5:
                es.append(np.median(lc["t"][m])); ds.append((lev - np.average(lc["flat"][m], weights=ww)) * 1e6)
                de.append(1e6 / np.sqrt(ww.sum()))
        ax5.errorbar(es, ds, de, fmt="o", ms=2, c="0.4", elinewidth=0.5)
        ax5.axhline(s.depth_ppm, c="r", lw=0.7); ax5.axhline(0, c="k", lw=0.5)
        ax5.set_xlabel("BTJD"); ax5.set_ylabel("per-transit depth (ppm)")
        ax5.set_ylim(np.percentile(ds, 1) - 2 * s.depth_ppm, np.percentile(ds, 99) + 2 * s.depth_ppm)
        fig.tight_layout()
        path = OUT / f"{_name(s.tic, s.period)}.png"
        fig.savefig(path, dpi=110)
        plt.close(fig)
        print(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("specs"); p.add_argument("--include", nargs="*", default=[])
    p.add_argument("--max-sectors", type=int, default=24)
    p = sub.add_parser("folds"); p.add_argument("--include", nargs="*", default=[])
    p = sub.add_parser("plots"); p.add_argument("--pick", nargs="+", required=True, help="TIC:period")
    a = ap.parse_args()
    {"specs": cmd_specs, "folds": cmd_folds, "plots": cmd_plots}[a.cmd](a)


if __name__ == "__main__":
    main()
