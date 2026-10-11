"""Step 1: prove the pipeline works by recovering a known planet.

Default target: pi Mensae c (TOI-144.01, TIC 261136679), a 2 Earth-radius
planet on a 6.27-day orbit around a naked-eye star -- the first planet TESS
announced (Huang et al. 2018). Its ~270 ppm transit is shallow enough to be a
real test of the detrending + BLS, but the star is bright enough that a
single sector should recover it at high SNR.

    .venv/bin/python validate.py                       # pi Men c, sector 1
    .venv/bin/python validate.py --tic 150428135 --planet "TOI-700 c" --sector 13

Outputs runs/validation/<slug>/{result.json, figure.jpg, arrays.npz}
"""
from __future__ import annotations

import argparse
import json
import urllib.parse
import urllib.request

import numpy as np

from planet_hunt import config
from planet_hunt.data import download, load_lightcurve, detrend
from planet_hunt.search import analyze
from planet_hunt.plots import candidate_figure, fig_to_jpeg

import lightkurve as lk  # noqa: E402  (after config: caches redirected)

# Fallback published values (Huang et al. 2018 / NASA Exoplanet Archive)
FALLBACK = {"pi Men c": {"period": 6.26784, "depth_ppm": 267.6, "duration_h": 2.952,
                         "rp_earth": 2.02}}


def published(planet: str) -> dict:
    """Published period/depth/duration from the NASA Exoplanet Archive."""
    q = ("select pl_name,pl_orbper,pl_trandep,pl_trandur,pl_rade from pscomppars "
         f"where pl_name = '{planet}'")
    url = f"{config.EXOARCHIVE_TAP}?query={urllib.parse.quote(q)}&format=json"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            row = json.loads(r.read())[0]
        return {"period": row["pl_orbper"],
                "depth_ppm": (row["pl_trandep"] or np.nan) * 1e4,  # % -> ppm
                "duration_h": row["pl_trandur"], "rp_earth": row["pl_rade"],
                "source": "NASA Exoplanet Archive (pscomppars)"}
    except Exception as exc:
        print(f"[validate] archive query failed ({exc}); using fallback values")
        return {**FALLBACK.get(planet, {}), "source": "fallback (Huang+2018)"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tic", type=int, default=261136679)
    ap.add_argument("--planet", default="pi Men c")
    ap.add_argument("--sector", type=int, default=1)
    args = ap.parse_args()

    # Find the SPOC 2-minute light curve with lightkurve's MAST search.
    sr = lk.search_lightcurve(f"TIC {args.tic}", author="SPOC", exptime=120,
                              sector=args.sector)
    if len(sr) == 0:
        raise SystemExit(f"No SPOC 2-min light curve for TIC {args.tic} in sector {args.sector}")
    print(sr)
    path = download(sr.table["dataURI"][0])

    raw, meta = load_lightcurve(path)
    first = detrend(raw, meta)
    out = analyze(first, raw, meta)
    st, bls, lc = out["stats"], out["bls"], out["lc"]
    pub = published(args.planet)

    def rel(a, b):
        return abs(a / b - 1) if b and np.isfinite(b) else np.nan

    checks = {
        "period_rel_err": rel(st["period"], pub.get("period")),
        "depth_rel_err": rel(st["depth_ppm"], pub.get("depth_ppm")),
        "duration_rel_err": rel(st["duration_h"], pub.get("duration_h")),
    }
    recovered = bool(checks["period_rel_err"] < 0.01 and st["snr"] >= 7.1)
    verdict = "RECOVERED" if recovered else "NOT RECOVERED"

    slug = args.planet.replace(" ", "_")
    run = config.RUNS / "validation" / slug
    run.mkdir(parents=True, exist_ok=True)
    title = (f"Validation: {args.planet} (TIC {args.tic}), sector {args.sector}: {verdict}")
    (run / "figure.jpg").write_bytes(fig_to_jpeg(candidate_figure(lc, bls, st, title)))
    clean = {k: (float(v) if isinstance(v, (float, np.floating)) else v)
             for k, v in st.items()}
    result = {"planet": args.planet, "tic": args.tic, "sector": args.sector,
              "found": clean, "published": pub, "checks": checks,
              "recovered": recovered, "star": lc.star}
    (run / "result.json").write_text(json.dumps(result, indent=2, default=float))

    print(f"\n=== {title} ===")
    print(f"{'':12s}{'found':>12s}{'published':>12s}{'rel. err':>10s}")
    for key, lab in (("period", "period d"), ("depth_ppm", "depth ppm"),
                     ("duration_h", "duration h")):
        k2 = {"period": "period_rel_err", "depth_ppm": "depth_rel_err",
              "duration_h": "duration_rel_err"}[key]
        print(f"{lab:12s}{st[key]:12.4f}{pub.get(key, np.nan):12.4f}{checks[k2]:10.1%}")
    print(f"SNR {st['snr']:.1f} (white {st['snr_white']:.1f}, red {st['snr_red']:.1f}), "
          f"SDE {st['sde']:.1f}, transits {st['n_transits']}, "
          f"Rp ~ {st['rp_earth']:.2f} Re (published {pub.get('rp_earth')})")
    print(f"flags: EB={st['eb_flags']} caution={st['caution_flags']}")
    print(f"saved -> {run}")


if __name__ == "__main__":
    main()
