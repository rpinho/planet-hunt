"""Step 3b: the cheapest real vetting test -- does the signal come back in
the star's OTHER TESS sectors?

A real planet transits every orbit, so it must reappear at the same period
and depth whenever TESS looks again. Noise and most systematics do not.

    .venv/bin/python followup.py --sector 90               # all clean non-catalogued TCEs
    .venv/bin/python followup.py --sector 90 --tic 278896737 --signal 1

For every other SPOC 2-min sector (up to --max-sectors) it
  * runs a narrow BLS within +/-1.5% of the candidate period and measures
    the depth/SNR there, and
  * computes the SNR the signal SHOULD have there given that sector's noise;
then it stitches all sectors and folds them near the candidate period
(the decisive test: a real planet's SNR grows ~sqrt(N sectors)), so
"faded where it should have been obvious" is distinguishable from
"too noisy to tell". Writes runs/sector_NNNN/followup/<tic>_<signal>.json
(report.py shows it on the candidate's card).
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from planet_hunt import config
from planet_hunt.data import CleanLC, detrend, download, load_lightcurve
from planet_hunt.search import bls_search, followup_targets, transit_stats

import lightkurve as lk  # noqa: E402


def check_signal(tic: int, period: float, depth_ppm: float, duration_h: float,
                 this_sector: int, this_snr: float, max_sectors: int = 6) -> dict:
    sr = lk.search_lightcurve(f"TIC {tic}", author="SPOC", exptime=120)
    tab = sr.table.to_pandas() if len(sr) else pd.DataFrame()
    sectors = []
    if len(tab):
        tab["sector"] = tab["mission"].str.extract(r"Sector (\d+)").astype(int)
        tab = tab.drop_duplicates("sector")
        # this sector + the nearest others in time (best period continuity)
        tab = tab.iloc[np.argsort(np.abs(tab["sector"] - this_sector))].head(max_sectors + 1)
        sectors = sorted(zip(tab["sector"], tab["dataURI"]))
    per_sector, stitched_parts = [], []
    for sector, uri in sectors:
        try:
            raw, meta = load_lightcurve(download(uri))
            lc0 = detrend(raw, meta)
            # stitch only sectors within ~1 year: across multi-year gaps the
            # trial-period grid needed to bridge the gap explodes
            if abs(sector - this_sector) <= 14:
                stitched_parts.append((sector, raw, meta, lc0))
            if sector == this_sector:
                continue
            b = bls_search(lc0.time, lc0.flat, pmin=period * 0.985, pmax=period * 1.015,
                           n_chunks=1)
            lc = detrend(raw, meta, transit=(b["period"], b["t0"], b["duration"]))
            st = transit_stats(lc, b["period"], b["t0"], b["duration"])
            # SNR the candidate SHOULD have here: same depth, this sector's noise
            # Same recipe as the detection SNR: the lower of a white-noise and a
            # red-noise (transit-timescale scatter) estimate. White noise alone
            # over-promises on spotted/granulating stars and would call real
            # signals "faded".
            n_in_expected = st["n_in"] * (duration_h / st["duration_h"])
            exp_white = depth_ppm / (st["sigma_ppm"] / np.sqrt(max(n_in_expected, 1)))
            exp_red = depth_ppm / (st["sigma_dur_ppm"] / np.sqrt(max(st["n_transits"], 1)))
            expected = min(exp_white, exp_red)
            per_sector.append({"sector": int(sector), "period": b["period"],
                               "depth_ppm": st["depth_ppm"], "snr": st["snr"],
                               "n_transits": st["n_transits"], "expected_snr": float(expected)})
            print(f"  S{sector}: best P in +/-1.5% = {b['period']:.4f} depth={st['depth_ppm']:.0f} ppm "
                  f"SNR={st['snr']:.1f} (expected ~{expected:.0f})", flush=True)
        except Exception as exc:
            per_sector.append({"sector": int(sector), "error": repr(exc)[:150]})

    # Stitch every sector and fold: a real planet keeps one exact period for
    # years, so its SNR grows ~sqrt(number of sectors). Noise does not add up.
    stitched = None
    if len(stitched_parts) >= 2:
        def stitch(lcs):
            t = np.concatenate([p.time for p in lcs])
            f = np.concatenate([p.flat for p in lcs])
            o = np.argsort(t)
            return CleanLC(tic=tic, sector=0, time=t[o], flux=f[o], flux_err=f[o] * 0,
                           trend=f[o] * 0 + 1, flat=f[o], star=lcs[0].star)
        # 1) find the exact period on the first-pass (unmasked) flattening
        lc_all = stitch([p[3] for p in stitched_parts])
        b = bls_search(lc_all.time, lc_all.flat, pmin=period * 0.99, pmax=period * 1.01,
                       n_chunks=1, frequency_factor=0.5)
        # 2) re-flatten every sector with those transits masked (otherwise the
        #    1-day filter eats ~1/3 of a 4-hour transit) and measure
        eph = (b["period"], b["t0"], b["duration"])
        lc_all = stitch([detrend(raw, meta, transit=eph) for _, raw, meta, _ in stitched_parts])
        st = transit_stats(lc_all, *eph)
        used = {p[0] for p in stitched_parts}
        exp = np.sqrt(this_snr ** 2 + sum(r.get("expected_snr", 0) ** 2 for r in per_sector
                                          if r["sector"] in used))
        stitched = {"n_sectors": len(stitched_parts), "sectors": sorted(int(x) for x in used),
                    "period": b["period"], "depth_ppm": st["depth_ppm"], "snr": st["snr"],
                    "n_transits": st["n_transits"], "expected_snr": float(exp),
                    "oe_sigma": st["oe_sigma"]}
        print(f"  stitched sectors {stitched['sectors']}: P={b['period']:.5f} "
              f"depth={st['depth_ppm']:.0f} ppm SNR={st['snr']:.1f} (expected ~{exp:.0f} if real)",
              flush=True)

    seen = [r for r in per_sector if r.get("snr", 0) >= 7
            and abs(r["period"] / period - 1) < 0.003 and 0.5 < r["depth_ppm"] / depth_ppm < 2]
    should = [r for r in per_sector if r.get("expected_snr", 0) >= 10]
    if not per_sector:
        verdict = "only one SPOC 2-min sector: cannot check yet"
    elif stitched is None:   # other sectors exist but are years away: per-sector test
        if seen:
            verdict = f"REAPPEARS in {len(seen)}/{len(per_sector)} other sector(s): worth a closer look"
        elif should:
            verdict = (f"NOT SEEN in {len(should)} other sector(s) where it should have had "
                       "SNR >= 10: most likely not a real periodic transit")
        else:
            verdict = "inconclusive: other sectors too noisy to tell"
    elif (stitched["snr"] >= 0.6 * stitched["expected_snr"] and stitched["snr"] >= 10
          and abs(stitched["period"] / period - 1) < 0.003
          and 0.5 < stitched["depth_ppm"] / depth_ppm < 2):
        verdict = (f"PERSISTS across {stitched['n_sectors']} stitched sectors "
                   f"(SNR {stitched['snr']:.1f}): survives this test, worth a closer look")
    elif stitched["expected_snr"] >= 12 and stitched["snr"] < 0.6 * stitched["expected_snr"]:
        verdict = (f"FADES when {stitched['n_sectors']} sectors are stitched (SNR "
                   f"{stitched['snr']:.1f} vs ~{stitched['expected_snr']:.0f} expected for a real "
                   "transit): most likely noise, not a planet")
    else:
        verdict = "inconclusive: other sectors too noisy to tell"
    return {"tic": tic, "period": period, "depth_ppm": depth_ppm,
            "other_sectors": per_sector, "stitched": stitched, "verdict": verdict}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sector", type=int, default=90)
    ap.add_argument("--tic", type=int)
    ap.add_argument("--signal", type=int, default=1)
    ap.add_argument("--max-sectors", type=int, default=6,
                    help="other sectors to download per star (~2 MB each)")
    ap.add_argument("--borderline", action="store_true",
                    help="also re-check signals whose only red flag is one transit at a gap edge")
    args = ap.parse_args()

    run = config.RUNS / f"sector_{args.sector:04d}"
    (run / "followup").mkdir(exist_ok=True)
    res = pd.read_csv(run / "results.csv")
    if args.tic:
        todo = res[(res["tic"] == args.tic) & (res["signal"] == args.signal)]
    else:
        todo = followup_targets(res, borderline=args.borderline)
    for r in todo.itertuples():
        print(f"TIC {r.tic} signal {r.signal}: P={r.period:.4f} d, {r.depth_ppm:.0f} ppm", flush=True)
        out = check_signal(int(r.tic), r.period, r.depth_ppm, r.duration_h,
                           args.sector, r.snr, args.max_sectors)
        print("  ->", out["verdict"], flush=True)
        (run / "followup" / f"{r.tic}_{r.signal}.json").write_text(
            json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    main()
