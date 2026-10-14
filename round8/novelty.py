"""Round 8, stage 7: live, read-only novelty gate for survivors.

    .venv/bin/python round8/novelty.py <tic> [<tic> ...] [--period P ...]   # -> runs/round8_crowd/novelty/<tic>.json

Per TIC (and Gaia DR3 id, and a 3' cone at the period or an alias):
  ExoFOP TOI + CTOI (fresh download), ExoFOP target page (imaging/spectra/notes), NASA Exoplanet Archive
  (pscomppars, toi), SPOC TCE bulk tables S1-107 + live MAST DV-product query, arXiv API (metadata),
  Zenodo search, Planet Hunters TESS Talk search, GitHub code search, every community clone (grep, incl. .gz),
  Eumemic's full results_public (raw/vetted/candidate tables), T16 (ApJS MRT tables 1-2) and RAVEN tables.
Nothing is written anywhere but this repo's runs/.
"""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import time
import urllib.parse

import numpy as np
import pandas as pd

from common import DATA, EUM, OUT, SRC, download

NOV = DATA / "novelty"
OTHER = DATA.parent / "r8" / "novelty"         # another agent's public-table copies, read-only (T16, RAVEN)
GH = str(__import__("pathlib").Path.home() / "dev/.bin/gh")
ALIASES = [1, 2, 3, 0.5, 1 / 3, 1.5, 2 / 3]


def _get(url, timeout=60):
    r = subprocess.run(["curl", "-s", "-L", "-m", str(timeout), "--retry", "4", "--retry-all-errors", url],
                       capture_output=True, text=True)
    return r.stdout


def alias_hit(P, Ps):
    for p in Ps:
        if not (np.isfinite(p) and p > 0 and np.isfinite(P)):
            continue
        for a in ALIASES:
            if abs(P / (p * a) - 1) < 0.003:
                return float(p), a
    return None


def sep_arcmin(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    return np.degrees(np.arccos(np.clip(c, -1, 1))) * 60


def check(tic: int, P: float, ra: float, dec: float, gaia: str | None, toi, ctoi) -> dict:
    out = {"tic": tic, "P": P, "checked_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    # ExoFOP TOI/CTOI: star, and 3' cone at the period/alias
    t = toi[toi["TIC ID"] == tic]
    out["toi_same_star"] = [f"TOI-{r['TOI']} P={r['Period (days)']} {r['TFOPWG Disposition']}" for _, r in t.iterrows()]
    tcol = next(cc for cc in ctoi.columns if cc.upper().startswith("TIC"))
    c = ctoi[pd.to_numeric(ctoi[tcol], errors="coerce") == tic]
    out["ctoi_same_star"] = [str(r.get("CTOI", "")) + f" P={r.get('Period (days)')}" for _, r in c.iterrows()]
    near = []
    for name, tab, rc, dc, pc in (("TOI", toi, "RA", "Dec", "Period (days)"), ("CTOI", ctoi, "RA", "Dec", "Period (days)")):
        if rc not in tab:
            continue
        rr = pd.to_numeric(tab[rc], errors="coerce")
        if rr.isna().all():
            from astropy.coordinates import Angle
            rr = tab[rc].apply(lambda s: Angle(str(s), unit="hourangle").deg if ":" in str(s) else np.nan)
            dd = tab[dc].apply(lambda s: Angle(str(s), unit="deg").deg if ":" in str(s) else np.nan)
        else:
            dd = pd.to_numeric(tab[dc], errors="coerce")
        s = sep_arcmin(ra, dec, rr.to_numpy(float), dd.to_numpy(float))
        for i in np.where(s < 3.0)[0]:
            p = pd.to_numeric(tab.iloc[i][pc], errors="coerce")
            if alias_hit(P, [p]):
                near.append(f"{name} {tab.iloc[i].iloc[0]} at {s[i] * 60:.0f}\" P={p}")
    out["toi_ctoi_cone3_alias"] = near
    # ExoFOP target page
    js = _get(f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json")
    try:
        j = json.loads(js)
        out["exofop_page"] = {k: len(j.get(k, []) or []) for k in ("imaging", "spectroscopy", "tseries", "notes",
                                                                    "planet_parameters", "tois", "ctois", "files")}
    except Exception:
        out["exofop_page"] = "unreadable"
    # NEA
    q = urllib.parse.quote(f"select pl_name,pl_orbper from pscomppars where tic_id='TIC {tic}'")
    out["nea_pscomppars"] = _get(f"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query={q}&format=csv").strip().splitlines()[1:]
    q = urllib.parse.quote(f"select toi,pl_orbper,tfopwg_disp from toi where tid={tic}")
    out["nea_toi"] = _get(f"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query={q}&format=csv").strip().splitlines()[1:]
    # TCEs (bulk tables) + live DV products
    k = pd.read_csv(NOV / "known.csv")
    kk = k[(k.tic == tic)]
    out["known_table_same_star"] = [f"{r.kind}: {r['name']} P={r.P}" for _, r in kk.iterrows()][:20]
    try:
        from astroquery.mast import Observations
        o = Observations.query_criteria(target_name=str(tic), obs_collection="TESS", dataproduct_type="timeseries")
        prods = Observations.get_product_list(o) if len(o) else []
        out["mast_dv_products"] = int(sum("dvr" in str(x) or "dvs" in str(x) or "dvm" in str(x)
                                         for x in (prods["productFilename"] if len(prods) else [])))
    except Exception as exc:
        out["mast_dv_products"] = f"error {exc!r}"[:80]
    # arXiv (metadata)
    xs = _get(f"http://export.arxiv.org/api/query?search_query=all:%22TIC%20{tic}%22&max_results=20")
    out["arxiv_tic"] = xs.count("<entry>")
    if gaia:
        xs = _get(f"http://export.arxiv.org/api/query?search_query=all:%22{gaia}%22&max_results=20")
        out["arxiv_gaia"] = xs.count("<entry>")
    # Zenodo
    try:
        z = json.loads(_get(f"https://zenodo.org/api/records?q=%22TIC%20{tic}%22&size=20") or "{}")
        out["zenodo"] = [h["metadata"]["title"] for h in z.get("hits", {}).get("hits", [])]
        z = json.loads(_get(f"https://zenodo.org/api/records?q=%22{tic}%22&size=20") or "{}")
        out["zenodo_bare"] = [h["metadata"]["title"] for h in z.get("hits", {}).get("hits", [])]
    except Exception as exc:
        out["zenodo"] = f"error {exc!r}"[:80]
    # PHT Talk (project 7929)
    try:
        tk = json.loads(_get(f"https://talk.zooniverse.org/searches?http_cache=true&section=project-7929&query=TIC%20{tic}&page_size=20") or "{}")
        out["pht_talk"] = [(c.get("discussion_title") or c.get("title") or "")[:80] for c in tk.get("searches", [])]
        tk = json.loads(_get(f"https://talk.zooniverse.org/searches?http_cache=true&section=project-7929&query={tic}&page_size=20") or "{}")
        out["pht_talk_bare"] = [(c.get("discussion_title") or c.get("title") or "")[:80] for c in tk.get("searches", [])]
    except Exception as exc:
        out["pht_talk"] = f"error {exc!r}"[:80]
    # GitHub code search (read-only API)
    env = dict(__import__("os").environ, GH_CONFIG_DIR=str(__import__("pathlib").Path.home() / ".config/gh"))
    r = subprocess.run([GH, "api", "-X", "GET", "search/code", "-f", f"q={tic}", "--jq",
                        ".total_count, (.items[]? | .repository.full_name + \":\" + .path)"], capture_output=True, text=True, env=env)
    out["github_code"] = r.stdout.strip().splitlines()[:30] if r.returncode == 0 else f"error {r.stderr[:80]}"
    time.sleep(7)    # code search: 10 requests / minute
    # community clones (grep, including gzipped tables)
    hits = []
    for d in sorted(SRC.iterdir()):
        r = subprocess.run(["grep", "-rIl", "--exclude-dir=.git", str(tic), str(d)], capture_output=True, text=True)
        hits += [p.replace(str(SRC) + "/", "") for p in r.stdout.splitlines()]
        for gz in d.rglob("*.gz"):
            if ".git" in gz.parts:
                continue
            try:
                if str(tic).encode() in gzip.open(gz).read():
                    hits.append(str(gz).replace(str(SRC) + "/", ""))
            except Exception:
                pass
    out["community_repos"] = hits
    eh = []
    for f in sorted((EUM / "results_public").glob("*")):
        try:
            d = pd.read_csv(f, low_memory=False)
        except Exception:
            continue
        col = "tic" if "tic" in d else ("TIC" if "TIC" in d else None)
        if col is None:
            continue
        g = d[pd.to_numeric(d[col], errors="coerce") == tic]
        if len(g):
            pcol = next((c for c in ("period", "P_d") if c in g), None)
            snr = next((c for c in ("snr", "MES") if c in g), None)
            eh.append(f"{f.name}: " + ", ".join(f"P={r[pcol]:.5f}" + (f" {snr}={r[snr]:.2f}" if snr else "")
                                                for _, r in g.iterrows()))
    out["eumemic_tables"] = eh
    # T16 + RAVEN (public tables, other agent's copies read only)
    t16 = []
    for f in ("t16/apjsae5b6ct1_mrt.txt", "t16/apjsae5b6ct2_mrt.txt"):
        p = OTHER / f
        if p.exists():
            for line in open(p, errors="ignore"):
                if f" {tic} " in line[19:32] or line[20:30].strip() == str(tic):
                    t16.append(f.split("/")[-1] + ": " + line[:60].strip())
    out["t16"] = t16 if (OTHER / "t16").exists() else "tables not available"
    rv = []
    for f in ("raven/vet_table.csv", "raven/val_table.csv", "raven/nsfp09_table.csv", "raven/val_rp8_table.csv"):
        p = OTHER / f
        if p.exists():
            d = pd.read_csv(p, usecols=lambda c: c in ("ticid", "bls_per", "toi", "ctoi"), low_memory=False)
            g = d[d.ticid == tic]
            rv += [f"{f}: P={r.bls_per}" for _, r in g.iterrows()]
    out["raven"] = rv
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--targets", default=str(OUT / "novelty_targets.csv"),
                   help="csv with tic, period, ra, dec, GAIA")
    a = p.parse_args()
    nd = OUT / "novelty"
    nd.mkdir(exist_ok=True)
    download("https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv", nd / "toi_gate.csv", 1000)
    download("https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv", nd / "ctoi_gate.csv", 1000)
    toi = pd.read_csv(nd / "toi_gate.csv", comment="#", low_memory=False)
    ctoi = pd.read_csv(nd / "ctoi_gate.csv", comment="#", low_memory=False)
    tg = pd.read_csv(a.targets)
    for r in tg.itertuples():
        res = check(int(r.tic), float(r.period), float(r.ra), float(r.dec),
                    str(int(r.GAIA)) if pd.notna(getattr(r, "GAIA", np.nan)) else None, toi, ctoi)
        (nd / f"{int(r.tic)}_{float(r.period):.4f}.json").write_text(json.dumps(res, indent=1, default=str))
        print(int(r.tic), {k: (len(v) if isinstance(v, list) else v) for k, v in res.items() if k not in ("tic", "P")},
              flush=True)


if __name__ == "__main__":
    main()
