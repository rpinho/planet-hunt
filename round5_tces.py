"""SPOC TCEs (every DV run, single and multi-sector) for each round-5 survivor."""
import json, urllib.request
import pandas as pd
from planet_hunt import config, spoc
from planet_hunt.catalogs import period_match
import round5_hosts as r5

tri = pd.read_csv(r5.OUT / "triage.csv")
sel = tri[tri.status.isin(["survivor", "flagged (rotation-saturated star)"])]
out = {}
for tic, g in sel.groupby("tic"):
    try:
        prods = spoc._dv_products(int(tic))
    except Exception as exc:
        out[int(tic)] = {"error": repr(exc)[:200]}; continue
    tces = []
    for s0, s1, uri, fn in prods:
        path = spoc.DV_DIR / fn
        if not path.exists():
            try:
                with urllib.request.urlopen(config.MAST_DOWNLOAD + uri, timeout=120) as r:
                    path.write_bytes(r.read())
            except Exception:
                continue
        for t in spoc._parse_dvr(path.read_text(errors="ignore")):
            t["run"] = f"S{s0}-{s1}"
            tces.append(t)
    res = {"n_runs": len(prods), "longest_run": (lambda r: f"S{r[0]}-{r[1]}" if r else None)(max(prods, key=lambda p: p[1] - p[0], default=None)), "candidates": {}}
    for P in g.period:
        hits = [t | {"match": period_match(P, t["period"], tol=0.01)} for t in tces
                if period_match(P, t["period"], tol=0.01) is not None]
        res["candidates"][f"{P:.5f}"] = hits
    res["all_tce_periods"] = sorted({round(t["period"], 4) for t in tces})
    out[int(tic)] = res
    print(tic, res["n_runs"], res["longest_run"], {k: [(h["run"], round(h["period"], 4), h["mes"]) for h in v] for k, v in res["candidates"].items()}, res["all_tce_periods"], flush=True)
(r5.OUT / "spoc_tces.json").write_text(json.dumps(out, indent=1, default=str))
