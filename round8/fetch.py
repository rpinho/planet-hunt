"""Round 8A: download every light-curve product (QLP, TESS-SPOC, SPOC 2-min LC, TGLC/eleanor-lite
for sectors with nothing else) and every TESScut 11x11 cutout for the P1 candidates.

    .venv/bin/python round8/fetch.py [--workers 4] [--only toi1846 ...]

Files: cache/r8/raw/<tic>/<file> (light curves), cache/tesscut/<ra>_<dec>/ (cutouts, shared with
earlier rounds). Downloads retry with backoff on connection errors and resume partial files.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from planet_hunt import config  # noqa: E402
from planet_hunt import pixelvet as pv  # noqa: E402

OUT = config.RUNS / "round8_p1"
RAW = config.CACHE / "r8" / "raw"
MIN_FREE_GB = 27.0


def free_gb():
    return shutil.disk_usage("/").free / 1e9


def curl(url, path, tries=8):
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(path.suffix + ".part")
    for k in range(tries):
        r = subprocess.run(["curl", "-s", "-f", "-L", "-C", "-", "--retry", "3", "--retry-all-errors",
                            "--connect-timeout", "30", "-m", "900", "-o", str(part), url],
                           capture_output=True, text=True)
        if r.returncode == 0 and part.exists() and part.stat().st_size > 2000:
            part.rename(path)
            return True
        if r.returncode == 22 and "416" not in r.stderr:   # HTTP error (not 'range not satisfiable')
            part.unlink(missing_ok=True)
        time.sleep(min(120, 10 * 2 ** k))
    return False


def want(p, summary):
    k = p["kind"]
    if k in ("QLP", "TESS-SPOC", "SPOC-2min"):
        return True
    if k in ("TGLC", "eleanor-lite"):
        # only for sectors with no QLP / TESS-SPOC / 2-min light curve
        have = set(summary.get("QLP", [])) | set(summary.get("TESS-SPOC", [])) | set(summary.get("SPOC-2min", []))
        return p["sector"] not in have
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--no-cutouts", action="store_true")
    ap.add_argument("--tpf", action="store_true", help="only the SPOC 2-min target pixel files")
    ap.add_argument("--sectors", type=int, nargs="*", help="restrict to these sectors")
    a = ap.parse_args()
    inv = json.load(open(OUT / "inventory.json"))
    jobs = []
    for key, v in inv.items():
        if a.only and key not in a.only:
            continue
        for p in v["products"]:
            if a.tpf and p["kind"] != "SPOC-2min-TPF":
                continue
            if a.sectors and p["sector"] not in a.sectors:
                continue
            if (a.tpf and p["kind"] == "SPOC-2min-TPF") or (not a.tpf and want(p, v["summary"])):
                dest = RAW / str(v["tic"]) / p["file"]
                if not dest.exists():
                    jobs.append((config.MAST_DOWNLOAD + p["uri"], dest))
    print(f"{len(jobs)} light-curve files to fetch; free {free_gb():.1f} GB", flush=True)
    with cf.ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(curl, u, d): d for u, d in jobs}
        for i, f in enumerate(cf.as_completed(futs)):
            d = futs[f]
            print(f"[{i + 1}/{len(jobs)}] {'ok ' if f.result() else 'FAIL'} {d.name}", flush=True)
            if free_gb() < MIN_FREE_GB:
                print("DISK GUARD", flush=True)
                break
    if a.no_cutouts or a.tpf:
        return
    for key, v in inv.items():
        if a.only and key not in a.only:
            continue
        if free_gb() < MIN_FREE_GB:
            print("DISK GUARD", flush=True)
            break
        for attempt in range(4):
            got = pv.fetch_cutouts(v["ra"], v["dec"], v["tesscut"], workers=a.workers)
            missing = sorted(set(v["tesscut"]) - {s for s, _ in got})
            print(f"{key}: {len(got)} cutouts, missing {missing}; free {free_gb():.1f} GB", flush=True)
            if not missing:
                break
            time.sleep(60 * (attempt + 1))


if __name__ == "__main__":
    main()
