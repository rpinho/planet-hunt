"""Is an FFI signal also in the star's SPOC 2-minute light curves?

    .venv/bin/python twomin_check.py --tic 259233660 --period 4.213463 --t0 1686.714262 \
        --duration-h 1.68 --depth-ppm 820 --mask 0.448956:1684.47826:0.73392

SPOC 2-minute apertures are small and PDCSAP is crowding-corrected, so a dip
on the target keeps its depth there; a dip from a neighbour that only the
larger FFI apertures include shrinks or vanishes. Known signals on the star
can be masked (--mask P:t0:duration_h). Reports the folded depth, SNR vs the
SNR expected at the FFI depth, per-year depths, odd/even and secondary, and
a local period/phase scan to make sure a slightly wrong ephemeris is not the
reason for a non-detection. Writes runs/round4_vet/twomin_<tic>.json/.png.
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from planet_hunt import config, lcvet
from planet_hunt import multisector as ms
from crossyear import spoc_2min, box_fit

OUT = config.RUNS / "round4_vet"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tic", type=int, required=True)
    ap.add_argument("--period", type=float, required=True)
    ap.add_argument("--t0", type=float, required=True)
    ap.add_argument("--duration-h", type=float, required=True)
    ap.add_argument("--depth-ppm", type=float, required=True)
    ap.add_argument("--mask", nargs="*", default=[])
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    P, t0, dur = args.period, args.t0, args.duration_h / 24
    masks = [tuple(float(v) for v in m.split(":")) for m in args.mask]
    masks = [(p, t, d / 24) for p, t, d in masks]

    def mask_fn(t):
        m = ms.in_tr(t, P, t0, dur, 1.6)
        for p, tt, d in masks:
            m |= ms.in_tr(t, p, tt, d, 1.6)
        return m
    lcs = spoc_2min(args.tic)
    lc = ms.stitch(list(lcs.values()), transit_mask=mask_fn)
    keep = np.ones(len(lc["t"]), bool)
    for p, tt, d in masks:
        keep &= ~ms.in_tr(lc["t"], p, tt, d, 1.6)
    lc = {k: (v[keep] if isinstance(v, np.ndarray) and len(v) == len(keep) else v) for k, v in lc.items()}
    bs = lcvet.box_stats(lc, P, t0, dur)
    w = 1 / lc["dy"] ** 2
    exp_snr = args.depth_ppm / bs["err_ppm"]
    # local scan: periods within the FFI period's uncertainty, full phase
    scan = []
    for dp in np.linspace(-3e-4, 3e-4, 13) * P:
        for sh in np.linspace(-0.5 * P, 0.5 * P, int(P / (dur / 3))):
            d, e = box_fit(lc["t"], lc["flat"], w, P + dp, t0 + sh, dur)
            if np.isfinite(d):
                scan.append((P + dp, t0 + sh, d / e, d))
    scan = np.array(scan)
    k = int(np.argmax(scan[:, 2]))
    res = {"tic": args.tic, "n_sectors": len(lcs), "sectors": sorted(lcs),
           "crowdsap": {int(s): v["meta"]["crowdsap"] for s, v in lcs.items()},
           "fold_at_ffi_ephemeris": {k_: bs[k_] for k_ in bs if k_ != "per_sector"},
           "expected_snr_at_ffi_depth": float(exp_snr),
           "scan_best": {"period": float(scan[k, 0]), "t0": float(scan[k, 1]), "snr": float(scan[k, 2]),
                         "depth_ppm": float(scan[k, 3] * 1e6), "n_trials": int(len(scan))},
           "scan_snr_at_ffi_ephemeris": float(bs["snr"])}
    # robustness: per-epoch depths, drop the deepest 5% of epochs, leave one season out
    intr = ms.in_tr(lc["t"], P, t0, dur)
    oot = ~ms.in_tr(lc["t"], P, t0, dur, 3.0)
    lev = np.average(lc["flat"][oot], weights=w[oot])
    ep = np.round((lc["t"] - t0) / P).astype(int)
    rows = []
    for e in np.unique(ep[intr]):
        m = intr & (ep == e)
        if m.sum() >= 0.5 * dur / (2 / 1440):
            rows.append((int(e), float(lev - np.average(lc["flat"][m], weights=w[m])), float(1 / np.sqrt(w[m].sum()))))
    rows = np.array(rows)
    def snr_of(keep_epochs):
        m = intr & np.isin(ep, keep_epochs)
        d = lev - np.average(lc["flat"][m], weights=w[m]); e_ = 1 / np.sqrt(w[m].sum())
        return float(d * 1e6), float(e_ * 1e6), float(d / e_)
    order = np.argsort(rows[:, 1] / rows[:, 2])[::-1]
    n_drop = max(1, int(round(0.05 * len(rows))))
    clip = snr_of(rows[order[n_drop:], 0].astype(int))
    loso = []
    for i, idx in enumerate(ms.seasons(lc["t"])):
        eps = np.unique(ep[idx])
        keep = rows[~np.isin(rows[:, 0], eps), 0].astype(int)
        loso.append({"season_dropped": i + 1, "sectors": sorted(set(int(x) for x in lc["sector"][idx])),
                     "depth_ppm": snr_of(keep)[0], "snr": snr_of(keep)[2]})
    from scipy import stats
    z = (rows[:, 1] - np.average(rows[:, 1], weights=1 / rows[:, 2] ** 2)) / rows[:, 2]
    res["robustness"] = {"n_epochs": int(len(rows)), "drop_top5pct": {"n_dropped": n_drop, "depth_ppm": clip[0],
                                                                       "err_ppm": clip[1], "snr": clip[2]},
                         "leave_one_season_out": loso,
                         "per_epoch_chi2_p": float(stats.chi2.sf(np.sum(z ** 2), len(z) - 1)),
                         "max_epoch_sigma": float(np.max(np.abs(z)))}
    per_sec = {}
    for s_ in np.unique(lc["sector"]):
        m = intr & (lc["sector"] == s_)
        if m.sum() > 10:
            per_sec[int(s_)] = (float((lev - np.average(lc["flat"][m], weights=w[m])) * 1e6), float(1e6 / np.sqrt(w[m].sum())))
    res["per_sector"] = per_sec
    # folded light curve near transit, for statistical validation (TRICERATOPS)
    near = np.abs(ms.phase(lc["t"], P, t0)) < 3 * dur
    np.savez(OUT / f"twomin_{args.tic}_fold.npz", phase_d=ms.phase(lc["t"], P, t0)[near],
             flux=lc["flat"][near], err=lc["dy"][near], period=P, depth=bs["depth_ppm"] * 1e-6,
             sectors=np.array(sorted(lcs)))
    (OUT / f"twomin_{args.tic}.json").write_text(json.dumps(res, indent=1, default=float))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ph = ms.phase(lc["t"], P, t0) * 24
    m = np.abs(ph) < 6
    so = np.argsort(ph[m])
    tb, yb, eb = ms.bin_weighted(ph[m][so] / 24, lc["flat"][m][so], lc["dy"][m][so], minutes=15)
    o = np.argsort(tb)
    ax[0].errorbar(tb[o] * 24, (yb[o] - 1) * 1e6, eb[o] * 1e6, fmt="k.", ms=3)
    ax[0].axhline(-args.depth_ppm, color="r", ls=":", label=f"FFI depth {args.depth_ppm:.0f} ppm")
    ax[0].axvspan(-args.duration_h / 2, args.duration_h / 2, color="r", alpha=0.08)
    ax[0].set_xlabel("hours from FFI mid-transit"); ax[0].set_ylabel("ppm"); ax[0].legend()
    ax[0].set_title(f"SPOC 2-min, {len(lcs)} sectors: {bs['depth_ppm']:.0f} ± {bs['err_ppm']:.0f} ppm "
                    f"(expected SNR {exp_snr:.0f})", fontsize=10)
    ax[1].hist(scan[:, 2], bins=60, color="0.6")
    ax[1].axvline(bs["snr"], color="r", label="FFI ephemeris")
    ax[1].set_xlabel("box SNR over local period/phase scan"); ax[1].legend()
    fig.tight_layout(); fig.savefig(OUT / f"twomin_{args.tic}.png", dpi=110); plt.close(fig)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
