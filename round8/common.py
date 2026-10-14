"""Round 8 shared helpers: paths, a network-switch-tolerant downloader, cache accounting.

The Mini may change networks mid-run (Starlink <-> Spectrum), so every download
retries on connection errors with a long back-off (up to ~15 min in total) and
resumes partial files (curl -C -). A 404 is final.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "round7"))
from planet_hunt import config  # noqa: E402

OUT = config.RUNS / "round8_crowd"
OUT.mkdir(parents=True, exist_ok=True)
DATA = Path.home() / "dev/planet-hunt-data/r8crowd"
SRC = DATA / "src"                         # read-only clones of the community repos
EUM = Path.home() / "dev/planet-hunt-data/r7/eumemic_src/repo"   # eumemic/exoplanets @ da722613 (still HEAD 10-13)
R8C = config.CACHE / "r8crowd"                # our own MAST inventory etc. (never the shared lcprod.csv)
R8C.mkdir(parents=True, exist_ok=True)
MIN_FREE_GB = 25.0
CACHE_GROWTH_CAP_GB = 12.0
BASELINE = OUT / "cache_baseline_kb.txt"   # du -sk of the shared dirs we add to, at start


def free_gb() -> float:
    return shutil.disk_usage("/").free / 1e9


def _du_kb(p: Path) -> int:
    if not p.exists():
        return 0
    r = subprocess.run(["du", "-sk", str(p)], capture_output=True, text=True)
    try:
        return int(r.stdout.split()[0])
    except Exception:
        return 0


GROW_DIRS = ["ffi", "r5/lc", "tesscut", "r8crowd", "mast"]


def cache_growth_gb() -> float:
    if not BASELINE.exists():
        BASELINE.write_text("\n".join(f"{d} {_du_kb(config.CACHE / d)}" for d in GROW_DIRS))
    base = dict(line.split() for line in BASELINE.read_text().splitlines())
    return sum(_du_kb(config.CACHE / d) - int(base.get(d, 0)) for d in GROW_DIRS) / 1e6


def guard() -> str | None:
    """None if we may keep downloading, else the reason to stop."""
    if free_gb() < MIN_FREE_GB:
        return f"free disk {free_gb():.1f} GB < {MIN_FREE_GB}"
    return None


def download(url: str, path: Path, min_bytes: int = 5000, max_wait_s: float = 900) -> bool:
    """curl with resume; retries on network errors for up to max_wait_s; False on 404."""
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    t0, k = time.time(), 0
    while True:
        r = subprocess.run(["curl", "-s", "-L", "-f", "-C", "-", "--connect-timeout", "20", "-m", "600",
                            "-o", str(part), "-w", "%{http_code}", url], capture_output=True, text=True)
        code = r.stdout.strip()[-3:]
        if r.returncode == 0 and part.exists() and part.stat().st_size >= min_bytes:
            part.rename(path)
            return True
        if code in ("404", "403", "400") or r.returncode == 22 and code.startswith("4") and code != "416":
            part.unlink(missing_ok=True)
            return False
        if code == "416":            # resume past end: restart clean
            part.unlink(missing_ok=True)
        k += 1
        if time.time() - t0 > max_wait_s:
            part.unlink(missing_ok=True)
            return False
        time.sleep(min(60, 5 * k))


def js(o):
    import numpy as np
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, float) and o != o:
        return None
    return str(o)
