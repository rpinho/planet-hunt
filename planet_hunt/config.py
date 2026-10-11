"""Project paths and cache setup.

Import this module BEFORE astropy / astroquery / lightkurve. It points every
cache those libraries use at directories inside the project (./cache), so the
starter kit never writes to ~/.astropy or ~/.lightkurve and you can wipe all
downloaded data by deleting ./cache.
"""
from __future__ import annotations

import os
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "cache"
MAST_DIR = CACHE / "mast"            # downloaded light-curve FITS files
TABLE_DIR = CACHE / "mast_tables"    # cached MAST observation queries (slow)
CATALOG_DIR = CACHE / "catalogs"     # TOI list, confirmed planets, TIC rows
RUNS = ROOT / "runs"                 # one sub-folder per hunt / validation
REPORTS = ROOT / "reports"           # final HTML reports

# astropy (and therefore astroquery) uses $XDG_*_HOME/astropy when it exists.
_XDG_CACHE = CACHE / "xdg_cache"
_XDG_CONFIG = CACHE / "xdg_config"
for d in (MAST_DIR, TABLE_DIR, CATALOG_DIR, RUNS, REPORTS,
          _XDG_CACHE / "astropy", _XDG_CONFIG / "astropy"):
    d.mkdir(parents=True, exist_ok=True)
os.environ["XDG_CACHE_HOME"] = str(_XDG_CACHE)
os.environ["XDG_CONFIG_HOME"] = str(_XDG_CONFIG)
os.environ.setdefault("MPLBACKEND", "Agg")   # headless plotting

# lightkurve prints a noisy warning about an optional sub-module; silence it.
warnings.filterwarnings("ignore", message=".*tpfmodel submodule.*")
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*invalid value.*")

# MAST public download endpoint for a "mast:TESS/product/..." data URI.
MAST_DOWNLOAD = "https://mast.stsci.edu/api/v0.1/Download/file?uri="
TOI_URL = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
EXOARCHIVE_TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
