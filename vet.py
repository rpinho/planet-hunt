"""Step 3c: is the dip on the target star? (centroid shift + neighbours)

    .venv/bin/python vet.py --sector 106
    .venv/bin/python vet.py --sector 106 --tic 123456789 --signal 1

For the same signals followup.py checks (clean, not catalogued), runs
  * the two light-curve-level "on target?" tests from planet_hunt/centroid.py:
    flux-weighted centroid shift in transit (this sector's SPOC file) and TIC
    neighbours within 2.5 pixels bright enough to produce the depth, and
  * a per-transit consistency test (do the individual dips agree?), and
  * the SPOC TCE check from planet_hunt/spoc.py: did NASA's own pipeline
    raise a threshold-crossing event at this period?
Writes runs/sector_NNNN/vet/<tic>_<signal>.json (shown by the reports).
"""
from __future__ import annotations

import argparse
import json

import pandas as pd

from planet_hunt import config
from planet_hunt.centroid import centroid_test, centroid_verdict, neighbour_check
import numpy as np

from planet_hunt.search import followup_targets, per_transit_check
from planet_hunt.spoc import spoc_tces


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--tic", type=int)
    ap.add_argument("--signal", type=int, default=1)
    ap.add_argument("--borderline", action="store_true",
                    help="also re-check signals whose only red flag is one transit at a gap edge")
    args = ap.parse_args()

    run = config.RUNS / f"sector_{args.sector:04d}"
    (run / "vet").mkdir(exist_ok=True)
    res = pd.read_csv(run / "results.csv")
    tg = pd.read_csv(run / "targets.csv").set_index("tic")
    if args.tic:
        todo = res[(res["tic"] == args.tic) & (res["signal"] == args.signal)]
    else:
        todo = followup_targets(res, borderline=args.borderline)
    for r in todo.itertuples():
        t = tg.loc[r.tic]
        path = config.MAST_DIR / str(t["dataURL"]).rsplit("/", 1)[-1]
        out = {"tic": int(r.tic), "signal": int(r.signal)}
        try:
            c = centroid_test(path, r.period, r.t0, r.duration_h, r.depth_ppm)
        except Exception as exc:
            c = {"ok": False, "why": repr(exc)[:150]}
        c["verdict"] = centroid_verdict(c)
        out["centroid"] = c
        try:
            out["neighbours"] = neighbour_check(int(r.tic), float(t["ra"]), float(t["dec"]),
                                                float(t["Tmag"]), float(r.depth_ppm))
        except Exception as exc:
            out["neighbours"] = {"error": repr(exc)[:150]}
        try:
            z = np.load(run / "lc" / f"{r.tic}_{r.signal}.npz")
            out["transits"] = per_transit_check(z["time"], z["flat"].astype(float), r.period,
                                                r.t0, r.duration_h / 24)
        except Exception as exc:
            out["transits"] = {"ok": False, "verdict": f"failed: {exc!r}"[:150]}
        try:
            out["spoc"] = spoc_tces(int(r.tic), args.sector, float(r.period))
        except Exception as exc:
            out["spoc"] = {"verdict": f"SPOC DV lookup failed: {exc!r}"[:150]}
        nb = out["neighbours"]
        print(f"TIC {r.tic} #{r.signal}: centroid -> {c['verdict']}; neighbours able to fake "
              f"{r.depth_ppm:.0f} ppm: {nb.get('n_able', '?')} of {nb.get('n_within', '?')}; "
              f"{out['spoc']['verdict']}; transits: {out['transits']['verdict']}", flush=True)
        (run / "vet" / f"{r.tic}_{r.signal}.json").write_text(json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    main()
