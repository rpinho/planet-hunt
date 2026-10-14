"""Round 8A: list every TESS product for each P1 candidate (MAST, live) and the TESScut sectors.

    .venv/bin/python round8/inventory.py           -> runs/round8_p1/inventory.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from planet_hunt import config  # noqa: E402
from planet_hunt import pixelvet as pv  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from targets import CANDS  # noqa: E402

OUT = config.RUNS / "round8_p1"
OUT.mkdir(parents=True, exist_ok=True)


def retry(fn, n=6):
    for k in range(n):
        try:
            return fn()
        except Exception as exc:
            print(f"  retry {k + 1}: {exc!r}"[:160], flush=True)
            time.sleep(15 * (k + 1))
    raise RuntimeError("MAST unavailable")


def tic_coords(tic):
    from astroquery.mast import Catalogs
    r = retry(lambda: Catalogs.query_criteria(catalog="TIC", ID=tic))
    r = r[0]
    return {k: (float(r[k]) if r[k] is not None and str(r[k]) != "--" else None)
            for k in ("ra", "dec", "Tmag", "Teff", "rad", "mass", "logg", "GAIAmag", "Vmag", "Jmag", "Hmag", "Kmag",
                      "plx", "d", "contratio", "e_rad", "e_mass", "e_Teff", "gaiabp", "gaiarp")} | {
        "GAIA": str(r["GAIA"]), "lumclass": str(r["lumclass"])}


def products(tic):
    from astroquery.mast import Observations
    obs = retry(lambda: Observations.query_criteria(target_name=str(tic), obs_collection=["TESS", "HLSP"]))
    if len(obs) == 0:
        return []
    prods = retry(lambda: Observations.get_product_list(obs))
    rows = []
    for r in prods:
        fn = str(r["productFilename"])
        if not fn.endswith(".fits"):
            continue
        kind = None
        if fn.startswith("tess") and fn.endswith("-s_lc.fits"):
            kind = "SPOC-2min"
        elif fn.endswith("fast-lc.fits"):
            kind = "SPOC-20s"
        elif fn.startswith("tess") and fn.endswith("-s_tp.fits"):
            kind = "SPOC-2min-TPF"
        elif "hlsp_qlp" in fn and fn.endswith("_llc.fits"):
            kind = "QLP"
        elif "hlsp_tess-spoc" in fn and fn.endswith("_lc.fits"):
            kind = "TESS-SPOC"
        elif "hlsp_tglc" in fn and fn.endswith("_llc.fits"):
            kind = "TGLC"
        elif "hlsp_gsfc-eleanor-lite" in fn:
            kind = "eleanor-lite"
        if kind is None:
            continue
        import re
        m = (re.search(r"-s0*(\d+)-", fn) or re.search(r"_s0*(\d+)[-_]", fn) or re.search(r"-s0*(\d+)_", fn))
        s = int(m.group(1)) if m else None
        rows.append({"kind": kind, "sector": s, "file": fn, "uri": str(r["dataURI"]),
                     "size": int(r["size"]) if str(r["size"]) != "--" else None})
    seen, out = set(), []
    for r in rows:
        if r["file"] in seen:
            continue
        seen.add(r["file"]); out.append(r)
    return out


def main():
    inv = {}
    for key, c in CANDS.items():
        print(f"== {key} TIC {c['tic']}", flush=True)
        tic = tic_coords(c["tic"])
        ra, dec = c["ra"] or tic["ra"], c["dec"] or tic["dec"]
        pr = products(c["tic"])
        cut = retry(lambda: pv.available_sectors(ra, dec))
        by = {}
        for r in pr:
            by.setdefault(r["kind"], set()).add(r["sector"])
        summ = {k: sorted(v) for k, v in by.items()}
        print("  ", {k: v for k, v in summ.items()}, "TESScut:", cut, flush=True)
        inv[key] = {"tic": c["tic"], "ra": ra, "dec": dec, "tic_row": tic, "products": pr, "summary": summ,
                    "tesscut": cut, "queried_utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}
    json.dump(inv, open(OUT / "inventory.json", "w"), indent=1)


if __name__ == "__main__":
    main()
