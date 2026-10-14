"""Round 8: MAST light-curve product inventory (SPOC 2-min, TESS-SPOC, QLP) for any TIC list.

    .venv/bin/python round8/inventory.py <csv with a tic column> [--workers 4]

Reads the shared cache/mast_tables/lcprod.csv read-only as a seed, writes ONLY our own
cache/r8/lcprod_r8.csv (another agent may append to the shared one). Network errors are
retried with back-off for up to ~15 min per chunk (the Mini may switch networks).
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import threading
import time

import pandas as pd

from common import R8C, config

PATH = R8C / "lcprod_r8.csv"
COLS = ["tic", "provenance_name", "sequence_number", "t_exptime", "t_min", "t_max", "dataURL"]
_lock = threading.Lock()


def load_all() -> pd.DataFrame:
    parts = []
    for p in (config.TABLE_DIR / "lcprod.csv", PATH):
        if p.exists():
            parts.append(pd.read_csv(p, usecols=lambda c: c in COLS, low_memory=False))
    d = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=COLS)
    d["tic"] = d["tic"].astype("int64")
    return d.drop_duplicates(["tic", "provenance_name", "sequence_number", "dataURL"])


def _query(ids: list[str]) -> pd.DataFrame | None:
    from astroquery.mast import Observations, conf
    conf.timeout = 300
    t0, k = time.time(), 0
    while True:
        try:
            obs = Observations.query_criteria(obs_collection=["TESS", "HLSP"], dataproduct_type="timeseries",
                                              target_name=ids, provenance_name=["SPOC", "TESS-SPOC", "QLP"])
            break
        except Exception as exc:
            k += 1
            if time.time() - t0 > 900:
                print(f"[mast] giving up on a chunk: {exc!r}"[:150], flush=True)
                return None
            time.sleep(min(60, 5 * k))
    o = obs["target_name", "provenance_name", "sequence_number", "t_exptime", "t_min", "t_max",
            "dataURL"].to_pandas() if len(obs) else pd.DataFrame(columns=["target_name"] + COLS[1:])
    o["tic"] = pd.to_numeric(o["target_name"], errors="coerce")
    o = o.dropna(subset=["tic"])
    url = o["dataURL"].astype(str)
    o = o[url.str.endswith("_lc.fits") | url.str.endswith("_llc.fits")]
    o = o[~((o["provenance_name"] == "SPOC") & (o["t_exptime"] < 100))]
    missing = set(int(t) for t in ids) - set(o["tic"].astype("int64"))
    if missing:
        o = pd.concat([o, pd.DataFrame({"tic": sorted(missing)})], ignore_index=True)
    o["tic"] = o["tic"].astype("int64")
    return o[COLS]


def inventory(tics, workers=4, chunk=150) -> pd.DataFrame:
    have = load_all()
    done = set(have.tic)
    todo = sorted({int(t) for t in tics} - done)
    print(f"[mast] {len(todo)} stars to inventory ({len(done)} known)", flush=True)
    chunks = [[str(t) for t in todo[i:i + chunk]] for i in range(0, len(todo), chunk)]
    n = 0
    with cf.ThreadPoolExecutor(workers) as ex:
        for o in ex.map(_query, chunks):
            n += 1
            if o is None:
                continue
            with _lock:
                o.to_csv(PATH, mode="a", header=not PATH.exists(), index=False)
            if n % 10 == 0:
                print(f"[mast] {n}/{len(chunks)} chunks", flush=True)
    d = load_all()
    return d[d.tic.isin({int(t) for t in tics})]


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("csv")
    p.add_argument("--workers", type=int, default=4)
    a = p.parse_args()
    t = pd.read_csv(a.csv, usecols=["tic"]).tic.dropna().astype("int64").unique()
    d = inventory(t, a.workers)
    print(len(d), "rows for", d.tic.nunique(), "stars")
