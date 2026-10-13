"""Figure 1 of the TIC 144439305 research note, built from the independent re-reduction
(round7/rnaas_check.py: raw HLSP FITS + own TESScut photometry), not from the pipeline cache.

    python round7/rnaas_fig.py --raw <HLSP fits dir> --cut <tesscut dir> --fit <recompute.json> --out runs/round7_rnaas
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rnaas_check as rc  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({"font.size": 8.5})


def binned(t, r, width_min=10, span_h=2.5):
    ph = rc.phase(t) * rc.P0 * 24
    s = np.abs(ph) < span_h
    edges = np.arange(-span_h, span_h + 1e-9, width_min / 60)
    idx = np.digitize(ph[s], edges)
    xs, ys, es = [], [], []
    for i in np.unique(idx):
        m = idx == i
        if m.sum() < 5:
            continue
        xs.append(ph[s][m].mean()); ys.append(r[s][m].mean() * 1e6); es.append(r[s][m].std() / np.sqrt(m.sum()) * 1e6)
    return np.array(xs), np.array(ys), np.array(es)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--cut", required=True)
    ap.add_argument("--fit", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    fit = json.load(open(a.fit))

    lcs = {}
    for p in sorted(glob.glob(f"{a.raw}/*_lc.fits") + glob.glob(f"{a.raw}/*_llc.fits")):
        s, prov, t, f = rc.load_hlsp(p)
        lcs[s] = rc.detrend(t, f)[:2]
    cuts = {}
    for p in sorted(glob.glob(f"{a.cut}/*_astrocut.fits")):
        s, prov, t, f = rc.load_cut(p)
        cuts[s] = rc.detrend(t, f)[:2]

    sets = [("TESS-SPOC S2, 29, 69 (used by the public search)", [lcs[s] for s in (2, 29, 69)], "k"),
            ("QLP S96, 103, 104 (not used by it)", [lcs[s] for s in (96, 103, 104)], "tab:orange"),
            ("TESScut S105, 106 (no light-curve products)", [cuts[s] for s in (105, 106)], "tab:blue")]

    # limb-darkened model at the posterior median (rnaas_check.py fit)
    import batman
    L = fit["ld_fit"]
    k, b, rho = L["k"][1], L["b"][1], L["rho_cgs"][1]
    aRs = (6.674e-8 * rho * (rc.P0 * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    prm = batman.TransitParams()
    prm.t0, prm.per, prm.rp, prm.a = 0.0, rc.P0, k, aRs
    prm.inc, prm.ecc, prm.w = np.degrees(np.arccos(b / aRs)), 0.0, 90.0
    prm.limb_dark, prm.u = "quadratic", [0.25, 0.35]
    tm = np.linspace(-2.6, 2.6, 600) / 24
    model = (batman.TransitModel(prm, tm, supersample_factor=5, exp_time=10 / 1440).light_curve(prm) - 1) * 1e6

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw=dict(width_ratios=[1.25, 1]))
    off = 0
    for lab, parts, col in sets:
        t = np.concatenate([p[0] for p in parts]); r = np.concatenate([p[1] for p in parts])
        x, y, e = binned(t, r)
        ax[0].errorbar(x, y - off, e, fmt="o", ms=2.5, color=col, elinewidth=0.7)
        ax[0].plot(tm * 24, model - off, "-", color="r", lw=1)
        ax[0].text(-2.45, -off + 450, lab, fontsize=7.5, color=col)
        off += 1600
    ax[0].set_xlim(-2.5, 2.5)
    ax[0].set_xlabel("Hours from mid-transit")
    ax[0].set_ylabel("Relative flux (ppm, offset)")
    ax[0].text(0.98, 0.02, "(a)", transform=ax[0].transAxes, ha="right")

    order = sorted(int(s) for s in fit["tesscut_sector"])
    pos = {s: i for i, s in enumerate(order)}
    for s, v in fit["depth_sector"].items():
        s = int(s)
        mk, col = ("o", "k") if v["prov"] == "TESS-SPOC" else ("s", "tab:orange")
        ax[1].errorbar(pos[s] - 0.15, v["depth_ppm"], v["err_ppm"], fmt=mk, color=col, ms=4)
    for s, v in fit["tesscut_sector"].items():
        s = int(s)
        col = "tab:blue" if s >= 105 else "0.45"
        ax[1].errorbar(pos[s] + 0.15, v["depth_ppm"], v["err_ppm"], fmt="D", mfc="none", color=col, ms=4)
    ax[1].set_xticks(range(len(order)), [str(s) for s in order])
    ax[1].axhline(fit["depth_all6"]["depth_ppm"], color="r", lw=0.8)
    ax[1].axhline(0, color="0.6", lw=0.6)
    ax[1].set_xlabel("TESS Sector")
    ax[1].set_ylabel("Depth at the ephemeris (ppm)")
    ax[1].set_ylim(-150, 1350)
    ax[1].text(0.98, 0.02, "(b)", transform=ax[1].transAxes, ha="right")
    fig.tight_layout()
    out = Path(a.out)
    fig.savefig(out / "fig1.pdf")
    fig.savefig(out / "fig1.png", dpi=200)
    print("ok")


if __name__ == "__main__":
    main()
