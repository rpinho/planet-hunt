"""Round 8, stage 1: one table of every unvetted / below-threshold public signal we can test.

    .venv/bin/python round8/sources.py eumemic     # -> runs/round8_crowd/src_eumemic.csv
    .venv/bin/python round8/sources.py community   # -> runs/round8_crowd/src_community.csv
    .venv/bin/python round8/sources.py merge       # -> signals.csv (deduplicated, excluded rows kept apart)

Common columns: stream, src_file, src_row, author, tic, period, t0_btjd (NaN if unknown),
dur_h, depth_ppm, snr_src, n_sec_src, rule (how to label the sectors the source used),
max_sec_src, note.

Eumemic streams (github.com/eumemic/exoplanets @ da722613, the HEAD on 2026-10-13):
  E1  his raw BLS tables, SNR 6-9 band, the round-7 table cuts (stage1.csv, 35,407 signals),
      minus the 120 round 7 already tested;
  E2  signals that entered his LEO-Vetter stage (SNR >= 9) and failed at most one SOFT test
      (noise-like or stellar-density tests, never secondary / odd-even / V-shape / radius),
      plus those that failed none but were not promoted (MES < 7.1, chi2 or centroid);
  E3  his full-sky catalogue rows he lists but does not count: tier C (MES < 10, grazing or RUWE),
      tier L (classifier < 0.8, MES >= 9) and "likely false alarm".  Novelty here is weaker: he
      has published the signal (uncounted); a survivor confirms his listing.
"""
from __future__ import annotations

import argparse
import glob
import json
import re

import numpy as np
import pandas as pd

from common import EUM, OUT, SRC, config

R7 = config.ROOT.parent / "planet-hunt-r7-eumemic-weak" / "runs" / "round7_eumemic_weak"   # read-only
SOFT = {"unphysical transit orbit", "inconsistent transit SNRs", "events not unique in phased light curve",
        "events not unique in local light curve", "not enough valid transits", "too many transits near gaps",
        "bad transit model fit", "dominated by single event"}
TBJD = 2457000.0


def r7select():
    """round7/select.py loaded under another name (it would shadow the stdlib 'select')."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("r7select", config.ROOT / "round7" / "select.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _std(df, **kw):
    out = pd.DataFrame({k: (df[v] if isinstance(v, str) and v in df else v) for k, v in kw.items()})
    return out


def eumemic() -> pd.DataFrame:
    rows = []
    # ---- E1
    s1 = pd.read_csv(R7 / "stage1.csv")
    tested = pd.read_csv(R7 / "targets.csv")[["tic", "period"]]
    key = lambda d: d.tic.astype(str) + "_" + np.round(d.period, 4).astype(str)
    s1 = s1[~key(s1).isin(set(key(tested)))]
    e1 = _std(s1, tic="tic", period="period", t0_btjd="t0", dur_h="duration_h",
              snr_src="snr", n_sec_src="n_sectors", sde="sde", rho="rho")
    e1["depth_ppm"] = s1["depth"].values * 1e6
    e1["stream"] = "E1"
    e1["src_file"] = "results_public/all_signals_" + s1["src"].map(
        {"km": "", "toih": "toi_hosts", "f1": "faint_t1", "f2": "faint_t2", "tce": "tce", "long": "long",
         "longtoi": "long_toi_hosts"}).values + ".csv.gz"
    e1["src_file"] = e1.src_file.str.replace("_.csv", ".csv")
    e1["src_row"] = "rank " + s1["rank"].astype(str).values
    e1["rule"] = "eum_early"
    e1["max_sec_src"] = 105
    rows.append(e1)
    # ---- E2
    v = []
    for f in sorted(glob.glob(str(EUM / "results_public/vetted_signals*.csv.gz"))):
        d = pd.read_csv(f)
        d["src_file"] = "results_public/" + f.split("/")[-1]
        v.append(d)
    v = pd.concat(v, ignore_index=True)
    fails = v.fails.fillna("").str.split("; ").apply(lambda xs: [x for x in xs if x])
    soft = fails.apply(lambda xs: len(xs) == 0 or (len(xs) == 1 and xs[0] in SOFT))
    v = v[soft & (v.known_status == "new") & (v.period >= 0.5) & (v.n_good_events >= 3)]
    e2 = _std(v, tic="tic", period="period", t0_btjd="t0", dur_h="duration_h", snr_src="snr",
              n_sec_src="n_sectors", sde="sde", rho="rho", src_file="src_file")
    e2["depth_ppm"] = v["depth"].values * 1e6
    e2["MES"] = v["MES"].values
    e2["stream"] = "E2"
    e2["src_row"] = ("rank " + v["rank"].astype(str) + " fails=" + v.fails.fillna("none")).values
    e2["rule"] = "eum_early"
    e2["max_sec_src"] = 105
    rows.append(e2)
    # ---- E3
    c = pd.read_csv(EUM / "results_public/candidates_allsky.csv")
    c = c[c.tier.isin(["C", "L", "likely false alarm"])]
    e3 = pd.DataFrame({"tic": c.TIC, "period": c.P_d, "t0_btjd": c.T0_BJD - TBJD, "dur_h": c.T14_h,
                       "depth_ppm": c.depth_ppm, "snr_src": c.MES, "MES": c.MES, "n_sec_src": np.nan,
                       "src_file": "results_public/candidates_allsky.csv",
                       "src_row": "tier " + c.tier.astype(str) + " reason=" + c.c_reason.fillna("").astype(str)})
    e3["stream"] = "E3"
    e3["rule"] = "eum_allsky"
    e3["max_sec_src"] = 104
    rows.append(e3)
    out = pd.concat(rows, ignore_index=True)
    out["author"] = "Eumemic (eumemic/exoplanets)"
    return out


# ----------------------------------------------------------------------------
# community tables (clones under <data>/r8crowd/src, read-only)
# ----------------------------------------------------------------------------
CLAIMED = {   # named / preregistered / retracted by their authors: never re-tested as leads
    "paras": {370612437, 252721481, 350883383},
    "avijit": {156733696, 345494583, 180791301},
}


def _git_head(d):
    import subprocess
    r = subprocess.run(["git", "-C", str(d), "log", "-1", "--format=%h"], capture_output=True, text=True)
    return r.stdout.strip()


def community() -> pd.DataFrame:
    rows = []
    # ---- C1 Aiden Huang: 1,279 long-observed M dwarfs, SPOC 2-min only, explicit sectors per star
    d = SRC / "Aidenn8_tess-planet-search"
    a = pd.read_csv(d / "results/all_signals.csv")
    a = a[a.verdict.isin(["not vetted", "below threshold"])]
    secs = {}
    for t in a.tic.unique():
        f = d / f"results/search/{t}.json"
        if f.exists():
            secs[t] = json.load(open(f)).get("sectors", [])
    c1 = pd.DataFrame({"tic": a.tic, "period": a.period, "t0_btjd": a.t0_btjd, "dur_h": a.duration_h,
                       "depth_ppm": a.depth_ppm, "snr_src": a.snr_bls, "sde": a.sde, "n_sec_src": a.n_sectors,
                       "src_row": "signal " + a.signal.astype(str) + " verdict=" + a.verdict})
    c1["src_sectors"] = [json.dumps(secs.get(t, [])) for t in a.tic]
    c1["rule"] = "explicit"
    c1["stream"] = "C1"
    c1["author"] = f"Aiden Huang (Aidenn8/tess-planet-search @{_git_head(d)})"
    c1["src_file"] = "results/all_signals.csv"
    rows.append(c1)
    # ---- C2/C3 Paras Chopra and Avijit Thawani (same pipeline): phase-2 signals with no hard flag
    for key, repo, st in (("paras", "paraschopra_tess-planet-candidates", "C2"),
                          ("avijit", "avijit-thawani_tess-niyati-avijit", "C3")):
        d = SRC / repo
        t = pd.read_csv(d / "results/phase2/triage.csv")
        t = t[t.hard_flags.isna() & ~t.tic.isin(CLAIMED[key])]
        c = pd.DataFrame({"tic": t.tic, "period": t.period, "t0_btjd": t.t0, "dur_h": t.duration_hr,
                          "depth_ppm": t.depth_ppm, "snr_src": t.snr, "sde": t.sde, "n_sec_src": t.n_sectors_used,
                          "src_row": "signal " + t.signal.astype(str) + " soft=" + t.soft_flags.fillna("none").astype(str)
                          + " survivor=" + t.survivor.astype(str)})
        c["src_sectors"] = ""
        c["rule"] = "lk_ffi"          # TESS-SPOC > QLP via lightkurve, rows with NaN errors dropped
        c["max_sec_src"] = 104
        c["stream"] = st
        who = "Paras Chopra" if key == "paras" else "Avijit Thawani"
        c["author"] = f"{who} ({repo.replace('_', '/', 1)} @{_git_head(d)})"
        c["src_file"] = "results/phase2/triage.csv"
        rows.append(c)
    # ---- C4 Thomas Marchand registry: best signal of every scanned star; REJECTED rows only
    d = SRC / "Th0rgal_open-exoplanet-hunt"
    reg = [json.loads(l) for l in open(d / "registry/scanned_targets.jsonl")]
    soft_g = {"gate2_multi_sector_baseline", "gate3_combined_bls_snr", "gate4_leave_one_epoch_out",
              "gate9_per_transit_coherence"}       # noise-like gates; odd/even, secondary, plausibility,
    reg = [r for r in reg if r.get("verdict") == "REJECTED_BY_FALSIFICATION" and r.get("best_period_d")   # centroid,
           and set(r.get("failed_gates", [])) <= soft_g]                                                     # rotation = hard
    c4 = pd.DataFrame({"tic": [r["tic_id"] for r in reg], "period": [r["best_period_d"] for r in reg],
                       "t0_btjd": [r.get("t0_btjd", np.nan) for r in reg], "dur_h": np.nan,
                       "depth_ppm": [r.get("weighted_dip_ppm", np.nan) for r in reg],
                       "snr_src": [r.get("combined_bls_snr", np.nan) for r in reg],
                       "n_sec_src": [len(r.get("sectors", [])) for r in reg],
                       "src_sectors": [json.dumps(r.get("sectors", [])) for r in reg],
                       "src_row": ["failed=" + ",".join(g.replace("gate", "g").split("_")[0] for g in r.get("failed_gates", []))
                                   + f" loeo_mean_z={r.get('loeo_mean_z')}" for r in reg]})
    c4["rule"] = "explicit"
    c4["stream"] = "C4"
    c4["author"] = f"Thomas Marchand (Th0rgal/open-exoplanet-hunt @{_git_head(d)})"
    c4["src_file"] = "registry/scanned_targets.jsonl"
    rows.append(c4)
    # ---- C5 Alejandro Valdez bot: status low_score (not vetted further), explicit 1-3 sectors
    d = SRC / "alexvaldex_nasa-tess-planet-search"
    v = pd.read_csv(d / "results/candidates.csv")
    v = v[v.status == "low_score"]
    hard = v["flags"].fillna("").str.contains("too_large|secondary|odd_even|v_shaped")
    v = v[~hard]
    c5 = pd.DataFrame({"tic": v.tic, "period": v.period, "t0_btjd": v.t0_btjd, "dur_h": v.duration_hr,
                       "depth_ppm": v.depth_ppm, "snr_src": v.snr, "sde": v.sde,
                       "n_sec_src": v.sectors.astype(str).str.split().str.len(),
                       "src_sectors": [json.dumps([int(x) for x in str(s).split()]) for s in v.sectors],
                       "src_row": v["id"] + " score=" + v.score.astype(str) + " flags=" + v["flags"].fillna("").astype(str)})
    c5["rule"] = "explicit"
    c5["stream"] = "C5"
    c5["author"] = f"Alejandro Valdez (alexvaldex/nasa-tess-planet-search @{_git_head(d)})"
    c5["src_file"] = "results/candidates.csv"
    rows.append(c5)
    # ---- C6 Marc Balcells ExoNit: single-sector strong candidates (L0, 'do not submit'); P to 0.01 d
    d = SRC / "marc252_exonit-results"
    x = pd.read_csv(d / "candidates/strong_candidates.csv")
    x = x.drop_duplicates(["tic_id", "period_days"])
    c6 = pd.DataFrame({"tic": x.tic_id, "period": x.period_days, "t0_btjd": np.nan, "dur_h": np.nan,
                       "depth_ppm": x.depth_ppm, "snr_src": x.snr, "n_sec_src": 1,
                       "src_sectors": [json.dumps([int(s[1:])]) for s in x.sector],
                       "src_row": x.run_id + " verdict=" + x.verdict.astype(str) + " eb=" + x.eb_flags.astype(str)})
    c6["rule"] = "explicit"
    c6["mode"] = "window"
    c6["stream"] = "C6"
    c6["author"] = f"Marc Balcells (marc252/exonit-results @{_git_head(d)})"
    c6["src_file"] = "candidates/strong_candidates.csv"
    rows.append(c6)
    # ---- C7 Bruno Rahle: 186 unclassified Sector-1 dossiers; DISCOVERY hypothesis (S1 only)
    d = SRC / "brahle_exoplanet-candidates"
    rr = []
    for f in sorted((d / "candidates").glob("*/evidence/ephemeris.json")):
        e = json.load(open(f))
        disc = next((h for h in e.get("hypotheses", []) if h.get("name") == "discovery"), None)
        sig = (disc or {}).get("signal") or e.get("ephemeris", {})
        cj = json.load(open(f.parent.parent / "candidate.json"))
        sec = int(re.search(r"-s(\d+)$", cj["id"]).group(1))
        rr.append({"tic": int(cj["target"]["tic"]), "period": sig.get("period_d"),
                   "t0_btjd": sig.get("epoch_bjd_tdb", np.nan) - TBJD, "dur_h": sig.get("duration_h"),
                   "depth_ppm": sig.get("depth_ppm"), "snr_src": np.nan, "n_sec_src": 1,
                   "src_sectors": json.dumps([sec]), "src_row": cj["id"] + " adopted=" + str(e.get("chosen_hypothesis"))})
    c7 = pd.DataFrame(rr)
    c7["rule"] = "explicit"
    c7["mode"] = "window"
    c7["stream"] = "C7"
    c7["author"] = f"Bruno Rahle (brahle/exoplanet-candidates @{_git_head(d)})"
    c7["src_file"] = "candidates/<id>/evidence/ephemeris.json"
    rows.append(c7)
    out = pd.concat(rows, ignore_index=True)
    out["mode"] = out["mode"].fillna("fixed") if "mode" in out else "fixed"
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["eumemic", "community"])
    a = p.parse_args()
    if a.cmd == "eumemic":
        d = eumemic()
        d.to_csv(OUT / "src_eumemic.csv", index=False)
        print(d.groupby("stream").agg(n=("tic", "size"), stars=("tic", "nunique")))
    if a.cmd == "community":
        d = community()
        d.to_csv(OUT / "src_community.csv", index=False)
        print(d.groupby("stream").agg(n=("tic", "size"), stars=("tic", "nunique"), author=("author", "first")).to_string())


if __name__ == "__main__":
    main()
