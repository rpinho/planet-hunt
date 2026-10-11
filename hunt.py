"""Step 2+3: search one TESS sector for transits and cross-match the hits.

    .venv/bin/python hunt.py --sector 90 --n 150
    .venv/bin/python hunt.py --sector 92 --n 300 --tmag 12 --controls 0

What it does
  1. picks ~N bright FGK dwarfs with SPOC 2-min light curves in the sector
     (plus up to --controls known-TOI hosts as positive controls, labelled)
  2. downloads their PDCSAP light curves from MAST (~2 MB each, cached)
  3. detrends, runs BLS (0.5 d .. min(15 d, baseline/2)), measures SNR and
     the vetting statistics, flags likely eclipsing binaries
  4. labels every signal against the ExoFOP TOI list + NASA confirmed planets
Outputs runs/sector_NNNN/{targets.csv, results.csv, summary.json, lc/*.npz}
Then run:  .venv/bin/python report.py --sector N
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import time

import numpy as np
import pandas as pd

from planet_hunt import config
from planet_hunt.catalogs import (crossmatch, load_confirmed, load_ctoi, load_toi,
                                  period_match)
from planet_hunt.data import (detrend, download_many, in_transit,
                              load_lightcurve, select_targets)
from planet_hunt.search import analyze


def search_one(args):
    """Worker: full analysis of one star (runs in a separate process).

    Iterative search: find the strongest BLS signal, then cut its transits
    out of the light curve and search again, up to `max_signals` times.
    This finds second planets, and recovers a real planet when a one-off
    systematic dip won the first round. The search stops after the first
    round that is not a threshold-crossing event (TCE); that row is kept
    because "the same period came back" is itself a vetting test.
    """
    tic, path, outdir, max_signals = args
    rows = []
    try:
        raw, meta = load_lightcurve(path)
        for k in range(max_signals):
            if len(raw) < 1000:
                break
            res = analyze(detrend(raw, meta), raw, meta)
            st, b, lc = res["stats"], res["bls"], res["lc"]
            # Save arrays so report.py can redraw figures without re-running BLS.
            np.savez_compressed(
                os.path.join(outdir, f"{tic}_{k + 1}.npz"), time=lc.time,
                flux=lc.flux.astype(np.float32), trend=lc.trend.astype(np.float32),
                flat=lc.flat.astype(np.float32), periods=b["periods"].astype(np.float32),
                power=b["power"].astype(np.float32))
            row = {"tic": tic, "signal": k + 1, **st,
                   **{f"star_{key}": v for key, v in lc.star.items()},
                   "baseline_d": float(lc.time.max() - lc.time.min()),
                   "n_points": len(lc.time), "t_start": float(lc.time.min()),
                   "t_end": float(lc.time.max())}
            row["eb_flags"] = "; ".join(st["eb_flags"])
            row["caution_flags"] = "; ".join(st["caution_flags"])
            rows.append(row)
            if not st["tce"]:
                break
            cut = in_transit(raw.time.value, st["period"], st["t0"],
                             st["duration_h"] / 24, widen=2.0)
            raw = raw[~cut]
    except Exception as exc:  # one bad file must not kill the hunt
        rows.append({"tic": tic, "signal": len(rows) + 1, "error": repr(exc)[:200]})
    return rows or [{"tic": tic, "signal": 1, "error": "too few good cadences"}]


def _add_flag(res, i, msg):
    old = res.at[i, "eb_flags"] if isinstance(res.at[i, "eb_flags"], str) else ""
    res.at[i, "eb_flags"] = "; ".join(filter(None, [old, msg]))


def flag_repeats(res: pd.DataFrame) -> pd.DataFrame:
    """If the same period (or 2x / 0.5x) comes back AFTER its transits were
    cut out, the "dips" were never transits: it is continuous variability
    (spots, pulsation, ellipsoidal/contact binary) or wide eclipses.
    Flag the first signal and demote the repeat."""
    for tic, g in res.groupby("tic"):
        if len(g) < 2:
            continue
        first = g[g["signal"] == 1]
        if first.empty or pd.isna(first["period"].iloc[0]):
            continue
        p1, snr1, d1 = first[["period", "snr", "depth_ppm"]].iloc[0]
        for i, r in g[g["signal"] > 1].iterrows():
            # comparable strength required: a weak echo at the same period is
            # just the imperfectly-masked edges of a real transit
            if (period_match(r["period"], p1, tol=0.02) is not None
                    and r["snr"] >= 0.5 * snr1 and r["depth_ppm"] >= 0.5 * d1):
                _add_flag(res, i, "repeat of signal 1")
                res.at[i, "tce"] = False
                _add_flag(res, first.index[0], "persists after masking (variable star / EB?)")
    return res


def flag_shared_ephemeris(res: pd.DataFrame, min_others: int = 2,
                          tol: float = 0.1) -> pd.DataFrame:
    """Flag signals whose transits line up with signals on OTHER stars.

    Real planets on different stars do not share transit times. If the
    predicted transits of this signal and of >= `min_others` other stars'
    signals coincide one-to-one (every in-sector transit of each has a partner
    within `tol` days), the dip is almost certainly a spacecraft systematic
    (momentum dump, scattered light, downlink ramp), not a planet.
    """
    ok = res.dropna(subset=["period", "t0"])
    t_lo, t_hi = ok["t_start"].min(), ok["t_end"].max()
    times, owner = {}, {}
    for r in ok.itertuples():
        n = np.arange(np.floor((t_lo - r.t0) / r.period), np.ceil((t_hi - r.t0) / r.period) + 1)
        tt = r.t0 + n * r.period
        times[r.Index] = tt[(tt >= t_lo) & (tt <= t_hi)]
        owner[r.Index] = r.tic

    def aligned(a, b):
        if len(a) == 0 or len(b) == 0:
            return False
        d = np.abs(a[:, None] - b[None, :])
        return d.min(axis=1).max() < tol and d.min(axis=0).max() < tol

    # count OTHER stars (not other signals of the same star) with aligned transits
    shared = {i: len({owner[j] for j, other in times.items()
                      if owner[j] != owner[i] and aligned(tt, other)})
              for i, tt in times.items()}
    res["shared_ephemeris"] = pd.Series(shared).reindex(res.index).fillna(0).astype(int)
    for i in res.index[res["shared_ephemeris"] >= min_others]:
        _add_flag(res, i, f"systematic? same transit times as "
                          f"{res.at[i, 'shared_ephemeris']} other stars")
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sector", type=int, default=90)
    ap.add_argument("--n", type=int, default=150, help="number of stars to search")
    ap.add_argument("--tmag", type=float, default=11.0, help="faintest TESS magnitude")
    ap.add_argument("--controls", type=int, default=15,
                    help="how many known-TOI hosts to include as positive controls")
    ap.add_argument("--max-signals", type=int, default=2,
                    help="signals to search per star (mask & re-search)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    args = ap.parse_args()
    t_start = time.time()

    run = config.RUNS / f"sector_{args.sector:04d}"
    (run / "lc").mkdir(parents=True, exist_ok=True)

    toi, conf, ctoi = load_toi(), load_confirmed(), load_ctoi()
    targets = select_targets(args.sector, args.n, tmag_max=args.tmag,
                             n_controls=args.controls, toi=toi, seed=args.seed)
    targets.to_csv(run / "targets.csv", index=False)
    print(f"[hunt] {len(targets)} targets ({targets['control'].sum()} TOI controls)", flush=True)

    paths = download_many(list(targets["dataURL"]))
    jobs = [(int(r.tic), str(paths[r.dataURL]), str(run / "lc"), args.max_signals)
            for r in targets.itertuples() if r.dataURL in paths]
    mb = sum(os.path.getsize(j[1]) for j in jobs) / 1e6
    print(f"[hunt] {len(jobs)} light curves on disk ({mb:.0f} MB); searching "
          f"with {args.workers} workers...", flush=True)

    rows = []
    with cf.ProcessPoolExecutor(args.workers) as ex:
        for k, r in enumerate(ex.map(search_one, jobs, chunksize=2), 1):
            rows.extend(r)
            if k % 20 == 0 or k == len(jobs):
                print(f"[bls] {k}/{len(jobs)}", flush=True)

    res = pd.DataFrame(rows).reset_index(drop=True)
    if "error" not in res:
        res["error"] = np.nan
    res = res.merge(targets[["tic", "Tmag", "Teff", "logg", "rad", "control"]],
                    on="tic", how="left")
    # Catalog labels for every star (searched or failed), so we can also ask
    # "did we recover the known TOIs among our targets?"
    labels = [crossmatch(int(r.tic), float(r.period) if pd.notna(getattr(r, "period", np.nan)) else np.nan,
                         toi, conf, ctoi) for r in res.itertuples()]
    res["catalog_status"] = [lab["status"] for lab in labels]
    res["catalog_label"] = [lab["label"] for lab in labels]
    res["tce"] = res["tce"].fillna(False).astype(bool)
    res = flag_repeats(res)
    res = flag_shared_ephemeris(res)
    res = res.sort_values(["tce", "snr"], ascending=False).reset_index(drop=True)
    res.to_csv(run / "results.csv", index=False)

    tce = res[res["tce"]]
    clean = tce[tce["eb_flags"].fillna("") == ""]
    summary = {
        "sector": args.sector, "n_targets": int(len(targets)),
        "n_searched": int(res.loc[res["error"].isna(), "tic"].nunique()),
        "n_failed": int(res.loc[res["error"].notna(), "tic"].nunique()),
        "n_controls": int(targets["control"].sum()),
        "n_tce": int(len(tce)), "n_tce_no_eb_flags": int(len(clean)),
        "tce_by_status": tce["catalog_status"].value_counts().to_dict(),
        "controls_recovered": int(tce[(tce["control"] == True)  # noqa: E712
                                      & tce["catalog_status"].isin(
                                          ["known_planet", "known_toi", "known_fp", "known_ctoi"])]["tic"].nunique()),
        "download_mb": round(mb, 1), "runtime_min": round((time.time() - t_start) / 60, 1),
        "settings": vars(args), "date": time.strftime("%Y-%m-%d"),
    }
    (run / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))
    cols = ["tic", "signal", "Tmag", "period", "depth_ppm", "duration_h", "snr", "n_transits",
            "eb_flags", "caution_flags", "catalog_label"]
    with pd.option_context("display.width", 250, "display.max_colwidth", 60):
        print(tce[cols].head(30).to_string(index=False))


if __name__ == "__main__":
    main()
