"""Did NASA's own pipeline see it? SPOC threshold-crossing events (TCEs).

SPOC runs its Transiting Planet Search on every 2-minute star, and for every
signal above its 7.1-sigma threshold it writes a Data Validation (DV) report
to MAST: `*_dvr.xml` (plus PDF summaries), per sector ("s0107-s0107") and for
multi-sector runs ("s0092-s0107"). DV files exist only for stars with at
least one TCE. So, once a sector's DV products are public:

  * our period matches a SPOC TCE  -> independent confirmation of the signal
    (not of a planet: TCEs include eclipsing binaries and systematics);
  * the star has no TCE in any run covering the sector -> SPOC's search,
    which is more sensitive than this kit's BLS, did not find it either:
    weak evidence against.

The DV XML is small (~100 kB) and cached in cache/dv/.
"""
from __future__ import annotations

import re
import urllib.request

from . import config
from .catalogs import period_match

DV_DIR = config.CACHE / "dv"


def _dv_products(tic: int):
    from astroquery.mast import Observations, conf
    conf.timeout = 120   # a stuck MAST request must not hang the run
    obs = Observations.query_criteria(obs_collection="TESS", target_name=str(tic),
                                      dataproduct_type="timeseries", provenance_name="SPOC")
    if len(obs) == 0:
        return []
    prods = Observations.get_product_list(obs)
    out = []
    for uri, fn in zip(prods["dataURI"], prods["productFilename"]):
        m = re.search(r"-s(\d{4})-s(\d{4})-\d{16}-\d+_dvr\.xml$", str(fn))
        if m:
            out.append((int(m.group(1)), int(m.group(2)), str(uri), str(fn)))
    return sorted(set(out))


def _parse_dvr(text: str) -> list[dict]:
    tces = []
    for block in re.findall(r"<dv:planetResults(.*?)</dv:planetResults>", text, flags=re.S):
        head = block[:block.find(">")]
        num = re.search(r'planetNumber="(\d+)"', head)
        toi = re.search(r'toiId="([^"]+)"', head)
        per = re.search(r'orbitalPeriodDays"[^>]*value="([0-9.eE+-]+)"', block)
        dep = re.search(r'transitDepthPpm"[^>]*value="([0-9.eE+-]+)"', block)
        mes = re.search(r'maxMultipleEventSigma="([0-9.eE+-]+)"', block)
        tces.append({"n": int(num.group(1)) if num else None,
                     "period": float(per.group(1)) if per else float("nan"),
                     "depth_ppm": float(dep.group(1)) if dep else float("nan"),
                     "mes": float(mes.group(1)) if mes else float("nan"),
                     "toi": toi.group(1) if toi else ""})
    return tces


def spoc_tces(tic: int, sector: int, period: float) -> dict:
    """SPOC TCEs for this star from DV runs that include `sector`."""
    DV_DIR.mkdir(parents=True, exist_ok=True)
    runs = [p for p in _dv_products(tic) if p[0] <= sector <= p[1]]
    found = []
    for s0, s1, uri, fn in runs:
        path = DV_DIR / fn
        if not path.exists():
            with urllib.request.urlopen(config.MAST_DOWNLOAD + uri, timeout=120) as r:
                path.write_bytes(r.read())
        for t in _parse_dvr(path.read_text(errors="ignore")):
            t["run"] = f"S{s0}" if s0 == s1 else f"S{s0}-{s1}"
            t["match"] = period_match(period, t["period"], tol=0.01)
            found.append(t)
    hits = [t for t in found if t["match"] is not None]
    if hits:
        h = hits[0]
        verdict = (f"SPOC TCE matches: P={h['period']:.4f} d, {h['depth_ppm']:.0f} ppm in DV run "
                   f"{h['run']}{h['match']}: NASA's pipeline sees it too")
    elif found:
        others = "; ".join(f"P={t['period']:.3f} d ({t['run']})" for t in found[:4])
        verdict = f"SPOC raised TCE(s) on this star, but at other periods: {others}"
    else:
        verdict = (f"no SPOC TCE for this star in any DV run covering sector {sector}: "
                   "NASA's (more sensitive) search did not flag it")
    return {"n_runs": len(runs), "tces": found, "match": bool(hits), "verdict": verdict}
