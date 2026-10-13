"""Do TESS-SPOC and QLP light curves give the same transit depth?

    .venv/bin/python product_calibration.py runs/mdwarf/results.csv

For every known planet / TOI re-found by mdwarf_hunt.py whose star has both
TESS-SPOC and QLP sectors cached, measure the depth separately in each
product (same ephemeris, same detrending) and report the QLP/TESS-SPOC ratio.
On-target planets should give ~1 if both products correct dilution the same
way; a candidate whose ratio is far outside the known planets' spread is
depending on the aperture, i.e. on a neighbour. Extra (tic, signal) pairs to
compare can be given with --also.
Writes runs/round4_vet/product_calibration.json.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from planet_hunt import config, ffi, lcvet
from planet_hunt import multisector as ms


def ratio_for(tic: int, P: float, t0: float, dur: float) -> dict | None:
    secs = ffi.load_star(tic)
    by = {}
    for prov in ("TESS-SPOC", "QLP"):
        sub = [d for d in secs if d["meta"].get("provenance") == prov]
        if len(sub) < 2:
            return None
        lc = ms.stitch(sub, transit_mask=lambda t: ms.in_tr(t, P, t0, dur, 1.6))
        if not lc:
            return None
        bs = lcvet.box_stats(lc, P, t0, dur)
        by[prov] = (bs["depth_ppm"], bs["err_ppm"], len(sub))
    (d1, e1, n1), (d2, e2, n2) = by["TESS-SPOC"], by["QLP"]
    if not (d1 > 3 * e1):
        return None
    r = d2 / d1
    er = abs(r) * np.hypot(e1 / d1, e2 / d2)
    return {"tic": tic, "spoc_ppm": d1, "spoc_err": e1, "n_spoc": n1, "qlp_ppm": d2, "qlp_err": e2,
            "n_qlp": n2, "ratio": r, "ratio_err": er}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results")
    ap.add_argument("--also", nargs="*", default=[], help="tic:signal pairs (candidates)")
    args = ap.parse_args()
    res = pd.read_csv(args.results)
    known = res[res["tce"] & res["catalog_status"].isin(["known_planet", "known_toi"])]
    rows = []
    for r in known.itertuples():
        x = ratio_for(int(r.tic), r.period, r.t0, r.duration_h / 24)
        if x:
            x.update({"label": r.catalog_label, "signal": int(r.signal), "contratio": float(r.contratio)})
            rows.append(x)
            print(f"{r.catalog_label:45s} QLP/SPOC = {x['ratio']:.2f} ± {x['ratio_err']:.2f}", flush=True)
    cands = []
    for pair in args.also:
        tic, sig = (int(v) for v in pair.split(":"))
        r = res[(res["tic"] == tic) & (res["signal"] == sig)].iloc[0]
        x = ratio_for(tic, r.period, r.t0, r.duration_h / 24)
        if x:
            x["signal"] = sig
            cands.append(x)
            print(f"candidate TIC {tic} #{sig}: QLP/SPOC = {x['ratio']:.2f} ± {x['ratio_err']:.2f}", flush=True)
    rr = np.array([x["ratio"] for x in rows]); ee = np.array([x["ratio_err"] for x in rows])
    w = 1 / ee ** 2
    summ = {"n": len(rows), "weighted_mean": float(np.sum(w * rr) / w.sum()),
            "median": float(np.median(rr)), "scatter": float(1.4826 * np.median(np.abs(rr - np.median(rr)))),
            "chi2_about_mean": float(np.sum(w * (rr - np.sum(w * rr) / w.sum()) ** 2))}
    for c in cands:
        c["sigma_from_known"] = float((c["ratio"] - summ["weighted_mean"]) /
                                      np.hypot(c["ratio_err"], summ["scatter"]))
    out = config.RUNS / "round4_vet" / "product_calibration.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"known": rows, "summary": summ, "candidates": cands}, indent=1, default=float))
    print(json.dumps(summ, indent=1), json.dumps(cands, indent=1, default=float))


if __name__ == "__main__":
    main()
