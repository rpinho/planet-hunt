"""Round 8: live catalogue snapshots for exclusion (stage 1) and the novelty gate (stage 7). Read-only.

    .venv/bin/python round8/catalogs.py fetch     # ExoFOP TOI + CTOI, NEA pscomppars, SPOC TCE bulk (S1-107)
    .venv/bin/python round8/catalogs.py build     # -> known.csv: (tic, P, name, kind) for alias matching

Files go to <data>/r8crowd/novelty/ with a fetched_utc.txt stamp.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime as dt
import glob
import json
import re
import urllib.request

import numpy as np
import pandas as pd

from common import DATA, EUM, SRC, download

NOV = DATA / "novelty"
NOV.mkdir(parents=True, exist_ok=True)
TCE_DIR = NOV / "tce"
# The link was ~1 MB/s shared with another agent, which had already pulled the same public SPOC TCE bulk
# tables (S1-107, ~500 MB); we read its copies READ-ONLY instead of downloading them twice.
TCE_SHARED = DATA.parent / "r8" / "novelty" / "tce" / "csv"
URLS = {
    "toi_live.csv": "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv",
    "ctoi_live.csv": "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv",
    "pscomppars.csv": ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
                       "select+pl_name,tic_id,pl_orbper,ra,dec+from+pscomppars&format=csv"),
    "nea_toi.csv": ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
                    "select+toi,tid,pl_orbper,tfopwg_disp+from+toi&format=csv"),
}
TCE_INDEX = "https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html"


def cmd_fetch(a):
    for name, url in URLS.items():
        ok = download(url, NOV / name, min_bytes=1000)
        print(name, ok, flush=True)
    TCE_DIR.mkdir(exist_ok=True)
    idx = NOV / "tce_index.html"
    download(TCE_INDEX, idx, min_bytes=1000)
    links = sorted(set(re.findall(r'(?:https://archive\.stsci\.edu)?(/missions/tess/catalogs/tce/[^"\'\s>]+\.csv)',
                                  idx.read_text())))
    print(len(links), "TCE files", flush=True)
    with cf.ThreadPoolExecutor(4) as ex:
        res = list(ex.map(lambda l: download("https://archive.stsci.edu" + l, TCE_DIR / l.split("/")[-1],
                                             min_bytes=200), links))
    print("TCE fetched", sum(res), "of", len(res))
    (NOV / "fetched_utc.txt").write_text(dt.datetime.utcnow().isoformat() + "Z\n")


def _tces() -> pd.DataFrame:
    out = []
    files = {}
    for dd in (TCE_SHARED, TCE_DIR):
        for f in glob.glob(str(dd / "*.csv")):
            files[f.split("/")[-1]] = f
    print(len(files), "TCE tables")
    for f in sorted(files.values()):
        try:
            d = pd.read_csv(f, comment="#", low_memory=False)
        except Exception:
            continue
        cols = {c.lower(): c for c in d.columns}
        tc = cols.get("ticid") or cols.get("tic_id")
        pc = cols.get("tce_period")
        if tc is None or pc is None:
            continue
        sec = re.search(r"-(s\d{4}-s\d{4})_", f).group(1)
        out.append(pd.DataFrame({"tic": d[tc], "P": d[pc], "name": "TCE " + sec + " " + d[tc].astype(str), "kind": "TCE"}))
    return pd.concat(out, ignore_index=True)


def eumemic_lists() -> pd.DataFrame:
    import importlib.util
    from common import config
    spec = importlib.util.spec_from_file_location("r7select", config.ROOT / "round7" / "select.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    h = m.his_candidates()
    return pd.DataFrame({"tic": h.tic, "P": h.P, "name": "Eumemic " + h.his_list.astype(str), "kind": "eumemic_list"})


def community_claims() -> pd.DataFrame:
    rows = []
    d = SRC / "Th0rgal_open-exoplanet-hunt"
    for f in glob.glob(str(d / "candidates/*/certificate.json")):
        c = json.load(open(f))
        tic = int(re.search(r"TIC-(\d+)", f).group(1))
        P = next((c[k] for k in ("period_d", "best_period_d", "period") if k in c), np.nan)
        if isinstance(P, dict):
            P = P.get("value", np.nan)
        rows.append((tic, P, "Th0rgal candidate/watchlist", "community_claim"))
    try:
        ext = json.load(open(d / "catalog/external_discoveries.json"))
        ext = ext if isinstance(ext, list) else ext.get("discoveries", ext.get("entries", []))
        for e in ext:
            tic = e.get("tic_id") or e.get("tic")
            if tic:
                rows.append((int(str(tic).replace("TIC", "").strip()), e.get("period_d", np.nan),
                             "Th0rgal external: " + str(e.get("source", e.get("discoverer", "")))[:60], "community_claim"))
    except Exception as exc:
        print("th0rgal external:", exc)
    a = pd.read_csv(SRC / "Aidenn8_tess-planet-search/results/all_signals.csv")
    for r in a[a.verdict.isin(["candidate", "weak candidate"])].itertuples():
        rows.append((r.tic, r.period, f"Aiden Huang {r.verdict}", "community_claim"))
    for t in (370612437, 252721481, 350883383):
        rows.append((t, np.nan, "Paras Chopra named/preregistered", "community_claim"))
    for t in (156733696, 345494583, 180791301):
        rows.append((t, np.nan, "Avijit Thawani preregistered (retracted)", "community_claim"))
    for f in glob.glob(str(SRC / "moosesteve21_PlanetCandidates/TIC*")):
        rows.append((int(re.search(r"TIC(\d+)", f).group(1)), np.nan, "moosesteve21 named candidate", "community_claim"))
    v = pd.read_csv(SRC / "alexvaldex_nasa-tess-planet-search/results/candidates.csv")
    for r in v[v.status == "candidate"].itertuples():
        rows.append((r.tic, r.period, "Alejandro Valdez candidate", "community_claim"))
    for t in (240488148, 4206066, 417732194, 264151965, 69912980, 451641911, 26432942, 150064768, 404160133,
              374098817, 233729835, 230202008, 38396135):
        rows.append((t, np.nan, "named on GitHub/Zenodo/PHT (community map)", "community_claim"))
    return pd.DataFrame(rows, columns=["tic", "P", "name", "kind"])


def cmd_build(a):
    toi = pd.read_csv(NOV / "toi_live.csv", comment="#", low_memory=False)
    toi = pd.DataFrame({"tic": toi["TIC ID"], "P": pd.to_numeric(toi["Period (days)"], errors="coerce"),
                        "name": "TOI-" + toi["TOI"].astype(str) + " " + toi["TFOPWG Disposition"].astype(str), "kind": "TOI"})
    ctoi = pd.read_csv(NOV / "ctoi_live.csv", comment="#", low_memory=False)
    pc = [c for c in ctoi.columns if c.lower().startswith("period")][0]
    tc = [c for c in ctoi.columns if c.upper().startswith("TIC")][0]
    nc = [c for c in ctoi.columns if c.upper().startswith("CTOI")][0]
    ctoi = pd.DataFrame({"tic": ctoi[tc], "P": pd.to_numeric(ctoi[pc], errors="coerce"),
                         "name": "CTOI-" + ctoi[nc].astype(str), "kind": "CTOI"})
    ps = pd.read_csv(NOV / "pscomppars.csv")
    ps = pd.DataFrame({"tic": pd.to_numeric(ps["tic_id"].astype(str).str.replace("TIC", "").str.strip(), errors="coerce"),
                       "P": ps["pl_orbper"], "name": ps["pl_name"], "kind": "confirmed"})
    parts = [toi, ctoi, ps, _tces(), eumemic_lists(), community_claims()]
    k = pd.concat(parts, ignore_index=True)
    k["tic"] = pd.to_numeric(k["tic"], errors="coerce")
    k = k.dropna(subset=["tic"])
    k["tic"] = k["tic"].astype("int64")
    k.to_csv(NOV / "known.csv", index=False)
    print(k.kind.value_counts().to_string())


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["fetch", "build"])
    a = p.parse_args()
    {"fetch": cmd_fetch, "build": cmd_build}[a.cmd](a)
