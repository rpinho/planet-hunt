"""planet_hunt: a small, readable TESS transit-search starter kit.

Modules
-------
config    paths + cache redirection (import first)
data      target selection, downloads, light-curve cleaning/detrending
search    Box Least Squares search, transit statistics, vetting flags
catalogs  ExoFOP TOI list + NASA Exoplanet Archive cross-match
plots     the diagnostic figures used by validate.py and report.py
"""
from . import config  # noqa: F401  (side effect: cache redirection)
