"""Does a few-transit signal from one sector persist in the star's other years?

    .venv/bin/python crossyear.py --tic 11528365 --period 8.171535 --t0 4267.422972 \
        --duration-h 1.944 --depth-ppm 814 --discovery-sector 107

For a signal found in one sector (here: 3 transits in sector 107), the
period is known to ~0.1-0.5%, too loose to fold years of data at the
catalogue ephemeris. This script:
  1. collects every sector of the star: SPOC 2-minute PDCSAP light curves where
     they exist, otherwise a fixed-aperture light curve built from a TESScut
     cutout (pixelvet scene fit + lcvet.cutout_lightcurve);
  2. measures the signal in each sector ALONE at the discovery ephemeris
     propagated with its uncertainty (a 1-D scan of mid-time inside the
     propagated window), and checks the new, independent sectors;
  3. scans a fine period grid around the candidate with all sectors together
     (coherent box fit, fixed duration), and compares the best stitched depth
     and SNR with what a constant transit would give.
Writes runs/round4_vet/crossyear_<tic>.json and .png.
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from planet_hunt import config
from planet_hunt import lcvet, pixelvet as pv
from planet_hunt import multisector as ms

OUT = config.RUNS / "round4_vet"


ARCHIVE = "https://archive.stsci.edu/missions/tess/tid"


def _lc_path(tic: int, sector: int):
    """Cached SPOC 2-min light curve, or download it from the static TESS
    archive tree (tid/sNNNN/xxxx/xxxx/xxxx/xxxx/). None if the star had no
    2-minute data in that sector."""
    import re
    import urllib.error
    import urllib.request
    hit = sorted(config.MAST_DIR.glob(f"tess*-s{sector:04d}-{tic:016d}-*-s_lc.fits"))
    if hit:
        return hit[0]
    t = f"{tic:016d}"
    url = f"{ARCHIVE}/s{sector:04d}/{t[0:4]}/{t[4:8]}/{t[8:12]}/{t[12:16]}/"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                names = sorted(set(re.findall(r"tess\d+-s\d{4}-\d{16}-\d{4}-s_lc\.fits", r.read().decode())))
            break
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
        except Exception as exc:
            print(f"[2min] S{sector} listing attempt {attempt + 1}: {exc!r}"[:120], flush=True)
    else:
        return None
    if not names:
        return None
    path = config.MAST_DIR / names[0]
    with urllib.request.urlopen(url + names[0], timeout=180) as r:
        path.write_bytes(r.read())
    return path


def spoc_2min(tic: int, max_sector: int = 115) -> dict[int, dict]:
    """SPOC 2-minute PDCSAP light curves of a star, every sector that has one."""
    import concurrent.futures as cf
    from astropy.io import fits
    with cf.ThreadPoolExecutor(8) as ex:
        paths = [p for p in ex.map(lambda s: _lc_path(tic, s), range(1, max_sector + 1)) if p is not None]
    print(f"[2min] TIC {tic}: {len(paths)} SPOC 2-min light curves", flush=True)
    out = {}
    for path in paths:
        with fits.open(path) as hd:
            d = hd[1].data
            s = int(hd[0].header["SECTOR"])
            crowd = float(hd[1].header.get("CROWDSAP", np.nan) or np.nan)
            t = np.asarray(d["TIME"], float)
            f = np.asarray(d["PDCSAP_FLUX"], float)
            q = np.asarray(d["QUALITY"])
        ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
        if ok.sum() < 100:
            continue
        out[s] = {"time": t[ok], "flux": f[ok] / np.nanmedian(f[ok]), "cx": np.zeros(ok.sum()),
                  "cy": np.zeros(ok.sum()), "meta": {"sector": s, "provenance": "SPOC-2min", "crowdsap": crowd}}
    return out


def tesscut_lcs(ra, dec, gaia_id, tmag, P, t0, dur, skip) -> dict[int, dict]:
    scene = pv.gaia_scene(ra, dec, gaia_id=gaia_id, tmag=tmag)
    out = {}
    for s, path in pv.fetch_cutouts(ra, dec):
        if s in skip:
            continue
        cut = pv.load_cutout(path)
        # scene fit needs an out-of-transit image; ephemeris here only sets which
        # frames are averaged, so a loose ephemeris is harmless
        sd = pv.prepare_sector(cut, scene, P, t0, dur)
        if sd is None:
            sd = pv.prepare_sector(cut, scene, P, t0 + P / 2, dur)
        if sd is None:
            continue
        out[s] = lcvet.cutout_lightcurve(cut, sd)
    return out


def box_fit(t, y, w, P, t0, dur):
    intr = ms.in_tr(t, P, t0, dur)
    oot = ~ms.in_tr(t, P, t0, dur, 3.0)
    if intr.sum() < 3:
        return np.nan, np.nan
    lev = np.average(y[oot], weights=w[oot])
    return float(lev - np.average(y[intr], weights=w[intr])), float(1 / np.sqrt(w[intr].sum()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tic", type=int, required=True)
    ap.add_argument("--period", type=float, required=True)
    ap.add_argument("--t0", type=float, required=True)
    ap.add_argument("--duration-h", type=float, required=True)
    ap.add_argument("--depth-ppm", type=float, required=True)
    ap.add_argument("--discovery-sector", type=int, required=True)
    ap.add_argument("--dp-frac", type=float, default=0.01, help="half-width of the period scan (fraction)")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    from astroquery.mast import Catalogs
    row = Catalogs.query_criteria(catalog="TIC", ID=args.tic).to_pandas().iloc[0]
    P, t0, dur = args.period, args.t0, args.duration_h / 24
    lcs = spoc_2min(args.tic)
    lcs.update(tesscut_lcs(float(row.ra), float(row.dec), row.GAIA, float(row.Tmag), P, t0, dur, set(lcs)))
    clean = {}
    for s, d in sorted(lcs.items()):
        c = ms.clean_sector(d, window=0.75, transit_mask=lambda t: ms.in_tr(t, P, t0, dur, 1.5))
        if c:
            # red noise: inflate by the scatter of transit-length bins vs white noise
            c["cad_arr"] = np.full(len(c["t"]), c["cadence"])
            tmp = {"t": c["t"], "flat": c["flat"], "dy": np.full(len(c["t"]), c["sigma"]),
                   "cad": c["cad_arr"], "sector": np.full(len(c["t"]), s)}
            oot = ~ms.in_tr(c["t"], P, t0, dur, 2.0)
            beta = float(np.median(ms.red_noise_beta(tmp, dur, oot)))
            c["beta"] = beta
            c["dy"] = np.full(len(c["t"]), c["sigma"] * beta)
            clean[s] = c
    res = {"tic": args.tic, "sectors": {}, "inputs": vars(args)}
    # 1. discovery sector alone: refine P, t0 and their errors with a local grid
    d = clean[args.discovery_sector]
    w = 1 / d["dy"] ** 2
    Ps = np.linspace(P * 0.985, P * 1.015, 301)
    t0s = t0 + np.linspace(-0.2, 0.2, 121)
    best = (-np.inf, P, t0)
    snr_grid = np.zeros((len(Ps), len(t0s)))
    for i, p in enumerate(Ps):
        for j, tt in enumerate(t0s):
            dd, ee = box_fit(d["t"], d["flat"], w, p, tt, dur)
            snr_grid[i, j] = dd / ee if np.isfinite(dd) else -np.inf
    i, j = np.unravel_index(np.argmax(snr_grid), snr_grid.shape)
    snr_max = snr_grid[i, j]
    good = snr_grid ** 2 > snr_max ** 2 - 1.0     # delta chi2 < 1: the 1-sigma region
    pi_, tj_ = np.where(good)
    P_d, t0_d = Ps[i], t0s[j]
    sP = float(max(np.ptp(Ps[pi_]) / 2, Ps[1] - Ps[0]))
    st0 = float(max(np.ptp(t0s[tj_]) / 2, t0s[1] - t0s[0]))
    res["discovery"] = {"period": float(P_d), "period_err": sP, "t0": float(t0_d), "t0_err": st0,
                        "snr": float(snr_max)}
    # 2. each other sector alone: scan mid-time inside the propagated window
    for s, c in clean.items():
        w = 1 / c["dy"] ** 2
        n = np.round((np.median(c["t"]) - t0_d) / P_d)
        win = float(np.hypot(st0, abs(n) * sP))
        tc = t0_d + n * P_d
        exp_err = box_fit(c["t"], c["flat"], w, P_d, tc, dur)[1]
        at_pred = box_fit(c["t"], c["flat"], w, P_d, tc, dur)
        scan = []
        if win < P_d / 2:
            for sh in np.linspace(-min(3 * win, P_d / 2), min(3 * win, P_d / 2), 241):
                dd, ee = box_fit(c["t"], c["flat"], w, P_d, tc + sh, dur)
                scan.append((float(sh), dd, ee))
        best_s = max(scan, key=lambda x: x[1] / x[2] if np.isfinite(x[1]) else -np.inf) if scan else None
        res["sectors"][int(s)] = {
            "provenance": c["provenance"], "cadence_min": c["cadence"] * 1440, "noise_ppm": c["sigma"] * 1e6,
            "red_noise_beta": c["beta"],
            "cycles_from_discovery": int(n), "ephemeris_window_d": win,
            "depth_at_prediction_ppm": at_pred[0] * 1e6 if np.isfinite(at_pred[0]) else None,
            "err_ppm": at_pred[1] * 1e6 if np.isfinite(at_pred[1]) else None,
            "expected_snr": float(args.depth_ppm * 1e-6 / exp_err) if np.isfinite(exp_err) else None,
            "best_in_window": None if best_s is None else {"shift_d": best_s[0], "depth_ppm": best_s[1] * 1e6,
                                                           "snr": best_s[1] / best_s[2]},
            "window_trials": int(len(scan) and (6 * win) / (dur / 3))}
    # 3. coherent period scan over all sectors
    t = np.concatenate([c["t"] for c in clean.values()])
    y = np.concatenate([c["flat"] for c in clean.values()])
    dy = np.concatenate([c["dy"] for c in clean.values()])
    w = 1 / dy ** 2
    base = t.max() - t.min()
    dP = dur / 4 * P_d / base
    grid = np.arange(P_d * (1 - args.dp_frac), P_d * (1 + args.dp_frac), dP)
    out = []
    for p in grid:
        # best t0 near the discovery transit for this period
        best_p = (-np.inf, np.nan, np.nan, np.nan)
        for sh in np.linspace(-1.5 * st0 - dur, 1.5 * st0 + dur, 13):
            dd, ee = box_fit(t, y, w, p, t0_d + sh, dur)
            if np.isfinite(dd) and dd / ee > best_p[0]:
                best_p = (dd / ee, dd, ee, t0_d + sh)
        out.append((p,) + best_p)
    out = np.array(out)
    k = int(np.nanargmax(out[:, 1]))
    # SNR expected for a constant transit at the discovery depth with all data
    exp = box_fit(t, y, w, out[k, 0], out[k, 4], dur)[1]
    res["stitched"] = {"best_period": float(out[k, 0]), "best_t0": float(out[k, 4]), "snr": float(out[k, 1]),
                       "depth_ppm": float(out[k, 2] * 1e6), "err_ppm": float(out[k, 3] * 1e6),
                       "expected_snr_if_constant": float(args.depth_ppm * 1e-6 / exp), "n_periods": int(len(grid)),
                       "baseline_d": float(base),
                       "snr_without_discovery_sector": None}
    # does the best stitched ephemeris survive without the discovery sector?
    m = ~np.isin(np.concatenate([np.full(len(c["t"]), s) for s, c in clean.items()]), [args.discovery_sector])
    dd, ee = box_fit(t[m], y[m], w[m], out[k, 0], out[k, 4], dur)
    res["stitched"]["snr_without_discovery_sector"] = float(dd / ee) if np.isfinite(dd) else None
    # the top 5 periods
    order = np.argsort(out[:, 1])[::-1]
    picks = []
    for i in order:
        if all(abs(out[i, 0] - q) > 5 * dP for q in picks):
            picks.append(out[i, 0])
        if len(picks) == 5:
            break
    res["stitched"]["top_periods"] = [{"period": float(p), "snr": float(out[np.argmin(abs(out[:, 0] - p)), 1])}
                                      for p in picks]
    # noise level of the scan: SNR distribution over the grid
    res["stitched"]["scan_snr_p99"] = float(np.nanpercentile(out[:, 1], 99))
    path = OUT / f"crossyear_{args.tic}.json"
    path.write_text(json.dumps(res, indent=1, default=float))
    plot(res, out, clean, args, path.with_suffix(".png"))
    print(json.dumps({k: v for k, v in res.items() if k != "sectors"}, indent=1, default=float))
    for s, v in res["sectors"].items():
        print(s, v)


def plot(res, scan, clean, args, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    ax[0].plot(scan[:, 0], scan[:, 1], "k-", lw=0.6)
    ax[0].axvline(res["discovery"]["period"], color="r", ls=":")
    ax[0].set_xlabel("period (d)"); ax[0].set_ylabel("box SNR (all sectors)")
    ax[0].set_title("coherent scan, all years", fontsize=10)
    s = res["stitched"]
    for a, (P, t0, title) in zip(ax[1:], ((res["discovery"]["period"], res["discovery"]["t0"],
                                            f"S{args.discovery_sector} alone"),
                                           (s["best_period"], s["best_t0"], "other sectors at best stitched P"))):
        for sec, c in clean.items():
            if ("alone" in title) != (sec == args.discovery_sector):
                continue
            ph = ms.phase(c["t"], P, t0) * 24
            m = np.abs(ph) < 8
            a.plot(ph[m], (c["flat"][m] - 1) * 1e6, ".", ms=1.5, alpha=0.4, label=f"S{sec}")
        a.set_xlabel("hours"); a.set_ylabel("ppm"); a.set_ylim(-4000, 3000); a.legend(fontsize=7, markerscale=5)
        a.set_title(title, fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()
