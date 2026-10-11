"""Known-object catalogs and the cross-match that labels each candidate.

Sources (both public):
  * ExoFOP-TESS TOI list (every TESS Object of Interest + its TFOPWG
    disposition: PC/APC = planet candidate, CP/KP = confirmed/known planet,
    FP = false positive, FA = false alarm).
  * NASA Exoplanet Archive, Planetary Systems Composite table (confirmed
    planets), queried over TAP for every planet with a TIC ID.
"""
from __future__ import annotations

import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

from . import config


def _fetch(url: str, path, max_age_days: float = 3.0) -> pd.DataFrame:
    """Download a CSV once and reuse it for `max_age_days`."""
    if not (path.exists() and time.time() - path.stat().st_mtime < max_age_days * 86400):
        print(f"[catalog] downloading {path.name} ...", flush=True)
        with urllib.request.urlopen(url, timeout=180) as r:
            path.write_bytes(r.read())
    return pd.read_csv(path, comment="#", low_memory=False)


def load_toi() -> pd.DataFrame:
    df = _fetch(config.TOI_URL, config.CATALOG_DIR / "toi.csv")
    df = df.rename(columns={"TIC ID": "tic", "Period (days)": "period",
                            "TFOPWG Disposition": "tfop", "TESS Disposition": "tess_disp",
                            "Depth (ppm)": "depth_ppm", "Duration (hours)": "duration_h",
                            "Planet Name": "planet_name", "Comments": "comments"})
    df["tic"] = df["tic"].astype("int64")
    return df


def toi_coverage() -> dict[int, int]:
    """Number of TOIs whose ExoFOP `Sectors` list includes each sector."""
    toi = _fetch(config.TOI_URL, config.CATALOG_DIR / "toi.csv")
    cov: dict[int, int] = {}
    for s in toi["Sectors"].dropna().astype(str):
        for x in s.split(","):
            x = x.strip()
            if x.isdigit():
                cov[int(x)] = cov.get(int(x), 0) + 1
    return cov


CTOI_URL = "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv"


def load_ctoi() -> pd.DataFrame:
    """Community TOIs (candidates submitted by outside teams/citizen scientists)."""
    df = _fetch(CTOI_URL, config.CATALOG_DIR / "ctoi.csv")
    df = df.rename(columns={"TIC ID": "tic", "Period (days)": "period",
                            "User Disposition": "user_disp",
                            "TFOPWG Disposition": "tfop", "Promoted to TOI": "promoted"})
    df["tic"] = df["tic"].astype("int64")
    return df


def load_confirmed() -> pd.DataFrame:
    q = ("select pl_name,hostname,tic_id,pl_orbper,pl_trandep,pl_trandur,pl_rade "
         "from pscomppars where tic_id is not null")
    url = f"{config.EXOARCHIVE_TAP}?query={urllib.parse.quote(q)}&format=csv"
    df = _fetch(url, config.CATALOG_DIR / "confirmed_planets.csv")
    df["tic"] = pd.to_numeric(df["tic_id"].str.replace("TIC", "").str.strip(),
                              errors="coerce").astype("Int64")
    return df.dropna(subset=["tic"])


# Ratios BLS commonly lands on instead of the true period.
ALIASES = {1.0: "", 2.0: " (BLS at 2x)", 0.5: " (BLS at 1/2)",
           3.0: " (BLS at 3x)", 1 / 3: " (BLS at 1/3)", 0.25: " (BLS at 1/4)"}


def period_match(p_found: float, p_known: float, tol: float = 0.01) -> str | None:
    """Return an alias note if p_found ~ r * p_known for a common ratio r."""
    if not (np.isfinite(p_known) and p_known > 0):
        return None
    for r, note in ALIASES.items():
        if abs(p_found / (r * p_known) - 1) < tol:
            return note
    return None


def crossmatch(tic: int, period: float, toi: pd.DataFrame,
               conf: pd.DataFrame, ctoi: pd.DataFrame | None = None) -> dict:
    """Label one BLS signal against the catalogs.

    status is one of
      known_planet   period matches a confirmed planet / TOI dispositioned CP or KP
      known_toi      period matches a TOI that is still a planet candidate
      known_fp       period matches a TOI dispositioned FP/FA (EB, systematics...)
      known_ctoi     period matches a Community TOI (not yet promoted to TOI)
      host_other     star hosts known TOI(s)/planet(s) but at a different period
                     -> the BLS peak is usually a residual of the known signal,
                        a systematic, or (rarely) an additional planet
      new            nothing in either catalog for this star
    """
    rows = toi[toi["tic"] == tic]
    planets = conf[conf["tic"] == tic]
    for _, p in planets.iterrows():
        note = period_match(period, p["pl_orbper"])
        if note is not None:
            return {"status": "known_planet", "label": f"Known planet {p['pl_name']}{note}",
                    "match": p["pl_name"]}
    for _, r in rows.iterrows():
        note = period_match(period, r["period"])
        if note is None:
            continue
        name = f"TOI-{r['TOI']}"
        disp = str(r["tfop"]) if pd.notna(r["tfop"]) else str(r["tess_disp"])
        if disp in ("CP", "KP"):
            pn = r["planet_name"] if pd.notna(r["planet_name"]) else ""
            return {"status": "known_planet",
                    "label": f"Known planet {name} {pn} ({disp}){note}".replace("  ", " "),
                    "match": name}
        if disp in ("FP", "FA"):
            return {"status": "known_fp",
                    "label": f"Known false positive {name} ({disp}){note}", "match": name}
        return {"status": "known_toi", "label": f"Known TOI {name} ({disp}){note}",
                "match": name}
    crow = ctoi[ctoi["tic"] == tic] if ctoi is not None else ctoi
    if crow is not None:
        for _, c in crow.iterrows():
            note = period_match(period, c["period"])
            if note is not None:
                disp = c["tfop"] if pd.notna(c["tfop"]) else c["user_disp"]
                return {"status": "known_ctoi",
                        "label": f"Known community TOI CTOI-{c['CTOI']} ({disp}){note}",
                        "match": f"CTOI-{c['CTOI']}"}
    if len(rows) or len(planets) or (crow is not None and len(crow)):
        known = [f"TOI-{t} P={p:.3f}d" for t, p in zip(rows["TOI"], rows["period"])]
        known += [f"{n} P={p:.3f}d" if np.isfinite(p) else f"{n} (no transit period)"
                  for n, p in zip(planets["pl_name"], planets["pl_orbper"].astype(float))]
        if crow is not None:
            known += [f"CTOI-{c} P={p:.3f}d" for c, p in zip(crow["CTOI"], crow["period"])]
        return {"status": "host_other",
                "label": "Known host, DIFFERENT period: " + "; ".join(known[:4]),
                "match": ""}
    return {"status": "new", "label": "NEW: not in TOI, CTOI or confirmed-planet catalogs",
            "match": ""}
