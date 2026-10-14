"""Round 8, stages 3-5: hold-out test of crowd signals on data their searches did not use.

    .venv/bin/python round8/holdout.py fetch   [--workers 4]   # light curves for targets.csv
    .venv/bin/python round8/holdout.py run     [--workers 8]   # per-signal JSON in runs/round8_crowd/holdout/
    .venv/bin/python round8/holdout.py inject  [--workers 8]   # completeness (runs/round8_crowd/inject/)
    .venv/bin/python round8/holdout.py collect                 # holdout.csv, noise expectation

PRE-REGISTERED (written and committed before any round-8 hold-out result was computed, 2026-10-13):

Data per star: SPOC 2-min > TESS-SPOC > QLP per sector (QLP files with empty error columns kept), through
the latest products on MAST (2-min S107, QLP S104; no TESS-SPOC/QLP light curves exist yet after S104).
Each sector is labelled THEIRS (data the source's search could have used) or NEW, by a rule fixed per source:
  eum_early  (Eumemic's raw tables, E1/E2): theirs = sector <= 105, except QLP files with empty errors in
             S94-97/S101-104 (his reader dropped them).  Identical to round 7.
  eum_allsky (Eumemic full-sky catalogue, E3): theirs = sector <= 104, except any QLP file with empty errors
             and QLP-only S74-79 (missing from his full-sky index, his issue #9).
  lk_ffi     (Paras Chopra / Avijit Thawani pipeline, C2/C3: TESS-SPOC > QLP via lightkurve, rows with
             non-finite errors dropped): theirs = sector <= 104 with a TESS-SPOC file, or a QLP file whose
             errors are not empty; a sector seen only at 2-min is NEW; when unsure, THEIRS (conservative).
  explicit   (C1 Aiden Huang, C4 Marchand, C5 Valdez, C6 ExoNit, C7 Rahle): theirs = the sectors the source
             lists for that star/signal.
Mode, fixed by the span of THEIR data (not by any result):
  FIXED  (their data span >= 150 d): the round-7 procedure unchanged -
    1. narrow BLS (source period +-0.5%) on THEIR data alone -> ephemeris;
    2. depth at that ephemeris in NEW data alone -> z_hold (red-noise errors);
    3. narrow BLS on all data -> SNR_full, odd/even, secondary, season tests;
    4. survivor = z_hold >= 2.5 AND hold-out depth within 3 sigma of their-data depth AND SNR_full >= 8 within
       0.1% of the source period AND |odd/even| < 3 AND secondary < 3 sigma AND season chi2 p > 0.001.
    Noise: 300 random ephemerides (P within +-5% but > 0.5% away, random T0) in the same NEW data.
  WINDOW (their data span < 150 d, e.g. one sector: the period cannot be propagated over years):
    1. narrow BLS (source period +-max(0.5%, 0.006 d/P)) on THEIR data -> P_their, depth_their;
    2. narrow BLS on NEW data alone in P_their +-0.5%, T0 free -> z_hold = depth/err at the best peak;
    3. as FIXED step 3 (window +-0.5%);
    4. survivor = z_hold >= 4 AND z_hold above every one of 60 null windows (same width, P_their x (1 +- U[3%,25%]),
       same NEW data) AND depth consistent (< 3 sigma) AND SNR_full >= 8 within 0.5% AND the FIXED-mode
       odd/even, secondary and season cuts.
    Noise: the 60 null windows give an empirical look-elsewhere p for every signal.
Both modes also run the whole procedure on the inverted light curve (flux -> 2 - flux).
Completeness: each target gets a transit with its source depth injected at another period (x0.7-1.3) and the
identical procedure is run.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

from common import OUT, config, download, free_gb, guard, cache_growth_gb, CACHE_GROWTH_CAP_GB, js

sys.path.insert(0, str(config.ROOT))
import revet as r7  # noqa: E402  (round7/revet.py)
import round5_hosts as r5  # noqa: E402
from planet_hunt import ffi  # noqa: E402
from planet_hunt import hostsearch as hs  # noqa: E402
from planet_hunt import multisector as ms  # noqa: E402

HO = OUT / "holdout"
INJ = OUT / "inject"
R5 = config.CACHE / "r5"
LC2 = R5 / "lc"
DROPPED_QLP = set(range(94, 98)) | set(range(101, 105))
QLP_TIMING = set(range(99, 105))
SPAN_WINDOW_D = 150.0
NULL_FIXED = 300
NULL_WINDOW = 60


def _dl(url, path, retries=4):        # resilient replacement for round5's downloader
    return download(url, Path(path), min_bytes=10000)


r5._download = _dl


# --------------------------------------------------------------------------
# plans and fetch
# --------------------------------------------------------------------------
def plan_star(tic: int, two: pd.DataFrame, prod: pd.DataFrame) -> list[dict]:
    rows = {}
    for r in two[two.tic == tic].itertuples():
        rows[int(r.sector)] = {"sector": int(r.sector), "prov": "SPOC-2min", "fname": r.fname, "alts": ["SPOC-2min"]}
    p = prod[prod.tic == tic].dropna(subset=["dataURL"])
    for prov in ("TESS-SPOC", "QLP"):
        for r in p[p.provenance_name == prov].itertuples():
            s = int(r.sequence_number)
            if s <= 0 or s > 200:
                continue
            if s not in rows:
                rows[s] = {"sector": s, "prov": prov, "url": r.dataURL, "alts": [prov]}
            elif prov not in rows[s]["alts"]:
                rows[s]["alts"].append(prov)
    return [rows[s] for s in sorted(rows)]


def fetch_ffi(tic, s, url, prov):
    out = ffi._npz_path(tic, s)
    if out.exists():
        return True
    out.parent.mkdir(parents=True, exist_ok=True)
    ffi.TMP_DIR.mkdir(parents=True, exist_ok=True)
    tmp = ffi.TMP_DIR / f"r8c_{tic}_{s}_{prov}.fits"
    if not download(config.MAST_DOWNLOAD + url, tmp):
        return False
    try:
        z = ffi.extract_fits(tmp, prov)
        if len(z["time"]) < 50:
            return False
        meta = z.pop("meta")
        part = out.with_name(out.stem + ".r8c.part")
        with open(part, "wb") as fh:
            np.savez_compressed(fh, **z, **{f"meta_{k}": v for k, v in meta.items()})
        part.rename(out)
        return True
    except Exception:
        return False
    finally:
        tmp.unlink(missing_ok=True)


def fetch_target(tic: int, plan: list[dict]) -> dict:
    two = [p["fname"] for p in plan if p["prov"] == "SPOC-2min"]
    st = {"tic": tic}
    if two and not (LC2 / f"{tic}.npz").exists():
        st["2min"] = r5.fetch_star(tic, two)["status"]
    n_ok = sum(fetch_ffi(tic, p["sector"], p["url"], p["prov"]) for p in plan if p["prov"] != "SPOC-2min")
    st["ffi_ok"], st["ffi_planned"] = n_ok, sum(p["prov"] != "SPOC-2min" for p in plan)
    return st


def load_plans(tg: pd.DataFrame) -> dict:
    path = OUT / "plans.json"
    plans = {int(k): v for k, v in json.load(open(path)).items()} if path.exists() else {}
    todo = [int(t) for t in tg.tic.unique() if int(t) not in plans]
    if todo:
        from inventory import load_all
        prod = load_all()
        prod = prod[prod.tic.isin(set(todo))]
        two = pd.read_csv(R5 / "twomin_index.csv.gz")
        two = two[two.tic.isin(set(todo))].drop_duplicates(["tic", "sector"])
        for t in todo:
            plans[t] = plan_star(t, two, prod)
        json.dump(plans, open(path, "w"))
    return plans


def cmd_fetch(a):
    tg = pd.read_csv(OUT / "targets.csv")
    plans = load_plans(tg)
    order = list(dict.fromkeys(int(t) for t in tg.tic))       # rank order
    print(f"{len(order)} stars, {sum(len(plans[t]) for t in order)} sector files planned; free {free_gb():.1f} GB",
          flush=True)
    t0 = time.time()
    with cf.ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(fetch_target, t, plans[t]): t for t in order}
        for k, f in enumerate(cf.as_completed(futs)):
            try:
                r = f.result()
            except Exception as exc:
                r = {"tic": futs[f], "err": repr(exc)[:100]}
            msg = f"[{k + 1}/{len(order)}] {r} free {free_gb():.1f} GB {time.time() - t0:.0f}s"
            if (k + 1) % 25 == 0:
                g = cache_growth_gb()
                msg += f" cache+{g:.2f} GB"
                if g > CACHE_GROWTH_CAP_GB:
                    print(msg + " CACHE CAP: stopping", flush=True)
                    for x in futs:
                        x.cancel()
                    break
            print(msg, flush=True)
            if guard():
                print("GUARD:", guard(), flush=True)
                for x in futs:
                    x.cancel()
                break


# --------------------------------------------------------------------------
# load + label
# --------------------------------------------------------------------------
def theirs(meta: dict, alts: list[str], sig: dict) -> bool:
    s, prov, empty = int(meta["sector"]), meta.get("provenance"), bool(meta.get("empty_err"))
    rule = sig["rule"]
    if rule == "eum_early":
        return s <= 105 and not (prov == "QLP" and empty and s in DROPPED_QLP)
    if rule == "eum_allsky":
        if s > 104:
            return False
        if prov == "QLP" and (empty or 74 <= s <= 79):
            return False
        return True
    if rule == "lk_ffi":
        if s > 104:
            return False
        if prov == "TESS-SPOC" or "TESS-SPOC" in alts:
            return True
        if prov == "QLP":
            return not empty
        return "QLP" in alts          # 2-min chosen: theirs if a QLP file exists (unknown errors: conservative)
    if rule == "explicit":
        return s in set(sig["src_sectors"])
    raise ValueError(rule)


def load_sectors(tic: int, plan: list[dict], sig: dict) -> list[dict]:
    secs = []
    alts = {p["sector"]: p.get("alts", [p["prov"]]) for p in plan}
    if (LC2 / f"{tic}.npz").exists():
        for d in r5.load_lc(tic):
            d["meta"]["provenance"] = "SPOC-2min"
            d["meta"]["empty_err"] = False
            secs.append(d)
    have = {int(d["meta"]["sector"]) for d in secs}
    for p in plan:
        if p["sector"] in have or p["prov"] == "SPOC-2min":
            continue
        path = ffi._npz_path(tic, p["sector"])
        if not path.exists():
            continue
        z = np.load(path)
        d = {"time": z["time"], "flux": z["flux"], "cx": z["cx"], "cy": z["cy"]}
        meta = {k[5:]: z[k].item() for k in z.files if k.startswith("meta_")}
        meta["empty_err"] = bool(not np.isfinite(z["flux_err"]).any())
        d["meta"] = meta
        secs.append(d)
    for d in secs:
        m = d["meta"]
        m["his"] = bool(theirs(m, alts.get(int(m["sector"]), []), sig))
    return sorted(secs, key=lambda d: d["meta"]["sector"])


# --------------------------------------------------------------------------
# window mode
# --------------------------------------------------------------------------
def z_at(lc, beta, e, sel):
    d, err, n = r7.depth_at(lc, beta, e["P"], e["t0"], e["dur"], sel)
    return d, err, n, (d / err if err else np.nan)


def run_window(lc, his_sel, new_sel, P0, rstar, mstar, teff, rng=None, n_null=NULL_WINDOW):
    w0 = max(0.005, 0.006 / P0)
    e_his = r7.narrow_bls(lc, his_sel, P0, rstar, mstar, frac=w0)
    if e_his is None:
        return {"error": "no their data"}
    beta = r7.sector_beta(lc, e_his["dur"])
    d_h, e_h, n_h, z_h = z_at(lc, beta, e_his, his_sel)
    out = {"mode": "window",
           "his": {**e_his, "depth_ppm": d_h * 1e6, "err_ppm": e_h * 1e6, "snr": z_h, "n_tr": n_h}}
    if new_sel.sum() > 100:
        e_new = r7.narrow_bls(lc, new_sel, e_his["P"], rstar, mstar, frac=0.005)
    else:
        e_new = None
    if e_new is None:
        out["hold"] = {"depth_ppm": np.nan, "err_ppm": np.nan, "z": np.nan, "n_tr": 0, "consistency_sigma": np.nan}
    else:
        d_n, e_n, n_n, z_n = z_at(lc, beta, e_new, new_sel)
        out["hold"] = {**{f"new_{k}": v for k, v in e_new.items()}, "depth_ppm": d_n * 1e6, "err_ppm": e_n * 1e6,
                       "z": z_n, "n_tr": n_n, "dP_frac": e_new["P"] / e_his["P"] - 1,
                       "consistency_sigma": float((d_n - d_h) / np.hypot(e_n, e_h)) if np.isfinite(d_n) else np.nan}
    e_all = r7.narrow_bls(lc, np.ones(len(lc["t"]), bool), e_his["P"], rstar, mstar, frac=0.005)
    st = r7.full_stats(lc, e_all["P"], e_all["t0"], e_all["dur"], rstar, mstar, teff)
    out["full"] = {**e_all, **{k: st[k] for k in ("depth_ppm", "depth_err_ppm", "snr", "snr_white", "beta_med",
                                                   "n_transits", "oe_sigma", "sec_sigma", "sec_any_sigma",
                                                   "season_chi2_p", "tr_chi2_p", "rp_earth", "dur_ratio",
                                                   "sin_ratio", "shape", "tr_max_over_mean")},
                   "season_depths": st["season_depths"]}
    if rng is not None and e_new is not None:
        zs = []
        for _ in range(n_null):
            Pn = e_his["P"] * (1 + rng.choice([-1, 1]) * rng.uniform(0.03, 0.25))
            if Pn < 0.3:
                continue
            en = r7.narrow_bls(lc, new_sel, Pn, rstar, mstar, frac=0.005)
            if en is None:
                continue
            zs.append(z_at(lc, beta, en, new_sel)[3])
        zs = np.array(zs, float)
        out["null"] = {"n": int(np.isfinite(zs).sum()), "z_max": float(np.nanmax(zs)) if len(zs) else None,
                       "z_mean": float(np.nanmean(zs)) if len(zs) else None,
                       "z_std": float(np.nanstd(zs)) if len(zs) else None,
                       "p_z": float(np.mean(zs >= out["hold"]["z"])) if len(zs) and np.isfinite(out["hold"]["z"]) else None,
                       "rate_z4": float(np.mean(zs >= 4)) if len(zs) else None, "zs": zs.round(2).tolist()}
    out["pass"] = passes_window(out, P0)
    return out


def passes_window(o, P0):
    h, f, nl = o.get("hold", {}), o.get("full", {}), o.get("null", {})
    z = np.nan_to_num(h.get("z"), nan=-9)
    c = {"z_hold": bool(z >= 4),
         "above_nulls": bool(nl.get("z_max") is None or z > nl["z_max"]) if nl else True,
         "consistent": bool(abs(np.nan_to_num(h.get("consistency_sigma"), nan=9)) < 3),
         "snr_full": bool(np.nan_to_num(f.get("snr"), nan=0) >= 8 and abs(f.get("P", 0) / o["his"]["P"] - 1) < 0.005),
         "oe": bool(abs(np.nan_to_num(f.get("oe_sigma"), nan=0)) < 3),
         "sec": bool(np.nan_to_num(f.get("sec_sigma"), nan=0) < 3),
         "season": bool(np.nan_to_num(f.get("season_chi2_p"), nan=1) > 0.001)}
    c["all"] = all(c.values())
    return c


# --------------------------------------------------------------------------
# per signal
# --------------------------------------------------------------------------
def sig_id(row) -> str:
    return f"{row['stream']}_{int(row['tic'])}_{float(row['period']):.4f}"


def holdout_signal(row: dict, plan: list[dict], inject: dict | None = None) -> dict:
    tic = int(row["tic"])
    rstar, mstar, teff = float(row["rstar"]), float(row["mstar"]), float(row["teff"])
    P0 = float(row["period"])
    sig = {"rule": row["rule"], "src_sectors": json.loads(row["src_sectors"]) if isinstance(row.get("src_sectors"), str)
           and row["src_sectors"].startswith("[") else []}
    secs = load_sectors(tic, plan, sig)
    if not secs:
        return {"tic": tic, "error": "no data"}
    if inject is not None:
        secs = hs.inject_transit(secs, inject, rstar, mstar, teff)
        P0 = inject["period"]
    toi, ps = r5.load_catalogs()
    kl = hs.known_list(r5.known_ephemerides(tic, toi, ps), rstar, mstar)
    kmask = lambda t: hs.mask_known(t, kl)
    rot = hs.rotation(secs, kmask)
    lc = hs.detrend(secs, hs.W_SHORT if P0 < hs.SPLIT_P else hs.W_LONG, kmask, rot=rot)
    if not lc:
        return {"tic": tic, "error": "detrend failed"}
    lc = hs.subset(lc, ~kmask(lc["t"]))
    sec_his = {int(d["meta"]["sector"]): d["meta"]["his"] for d in secs}
    his_sel = np.array([sec_his.get(int(s), False) for s in lc["sector"]])
    new_sel = ~his_sel
    info = [{"sector": int(d["meta"]["sector"]), "prov": d["meta"].get("provenance"), "his": d["meta"]["his"],
             "empty_err": d["meta"]["empty_err"]} for d in secs]
    span = float(np.ptp(lc["t"][his_sel])) if his_sel.sum() > 1 else 0.0
    mode = "fixed" if span >= SPAN_WINDOW_D else "window"
    res = {"id": sig_id(row), "tic": tic, "row": {k: (v.item() if hasattr(v, "item") else v) for k, v in row.items()},
           "sectors": info, "n_his_sectors": int(sum(i["his"] for i in info)),
           "n_new_sectors": int(sum(not i["his"] for i in info)), "his_span_d": span, "mode": mode,
           "n_qlp_timing": int(sum(i["prov"] == "QLP" and i["sector"] in QLP_TIMING for i in info)),
           "rotation": {k: rot.get(k) for k in ("active", "P_rot", "saturated") if k in rot},
           "known_masked": [k["name"] for k in kl if k["masked"]]}
    if his_sel.sum() < 100:
        res["error"] = "no their-era data on disk"
        return res
    rng = np.random.default_rng(tic % (2 ** 32))
    if mode == "fixed":
        real = r7.run_procedure(lc, his_sel, new_sel, P0, rstar, mstar, teff)
        if "error" in real:
            res["error"] = real["error"]
            return res
        beta, e_his = real.pop("_beta"), real.pop("_e_his")
        real["mode"] = "fixed"
        if inject is None:
            zs, cons = r7.null_ephemerides(lc, beta, his_sel, new_sel, e_his["P"], e_his["dur"], rng, n=NULL_FIXED)
            if len(zs):
                z_obs = real["hold"]["z"]
                res["null"] = {"n": int(np.isfinite(zs).sum()),
                               "p_z": float(np.mean(zs >= z_obs)) if np.isfinite(z_obs) else None,
                               "rate_z25_cons": float(np.mean((zs >= 2.5) & (np.abs(cons) < 3))),
                               "z_std": float(np.nanstd(zs)), "z_mean": float(np.nanmean(zs))}
    else:
        real = run_window(lc, his_sel, new_sel, P0, rstar, mstar, teff, rng=rng)
        if "error" in real:
            res["error"] = real["error"]
            return res
        if "null" in real:
            res["null"] = real.pop("null")
            res["null"]["rate_pass"] = res["null"].pop("rate_z4")
    res["real"] = real
    if inject is not None:
        res["inject"] = inject
        return res
    inv = dict(lc)
    inv["flat"] = 2.0 - lc["flat"]
    if mode == "fixed":
        ir = r7.run_procedure(inv, his_sel, new_sel, P0, rstar, mstar, teff)
        ir.pop("_beta", None), ir.pop("_e_his", None)
    else:
        ir = run_window(inv, his_sel, new_sel, P0, rstar, mstar, teff, rng=np.random.default_rng(tic % 7919), n_null=20)
        ir.pop("null", None)
    res["inverted"] = {k: ir.get(k) for k in ("his", "hold", "pass")}
    res["inverted"]["full_snr"] = ir.get("full", {}).get("snr")
    return res


def _worker(args):
    row, plan = args
    path = HO / f"{sig_id(row)}.json"
    if path.exists():
        return path.stem, "cached"
    t0 = time.time()
    try:
        r = holdout_signal(row, plan)
        r["elapsed_s"] = time.time() - t0
        path.write_text(json.dumps(r, default=js))
        rr = r.get("real", {})
        if not rr:
            return path.stem, r.get("error")
        return path.stem, (f"{r['mode']} z_hold {rr['hold'].get('z', np.nan):.2f} snr_full "
                           f"{rr['full'].get('snr', np.nan):.1f} pass {rr['pass']['all']} "
                           f"inv {r.get('inverted', {}).get('pass', {}).get('all')} {time.time() - t0:.0f}s")
    except Exception as exc:
        (HO / f"{path.stem}.err").write_text(traceback.format_exc())
        return path.stem, f"error {exc!r}"[:150]


def cmd_run(a):
    HO.mkdir(parents=True, exist_ok=True)
    tg = pd.read_csv(OUT / "targets.csv")
    plans = load_plans(tg)
    if a.only:
        tg = tg[tg.tic.isin(a.only)]
    if a.fetched:     # only stars whose fetch finished (fetch.log completion lines)
        import re
        done = {int(m) for m in re.findall(r"\{'tic': (\d+)", "".join(p.read_text() for p in OUT.glob("fetch*.log")))}
        tg = tg[tg.tic.isin(done)]
    jobs = [(r, plans.get(int(r["tic"]), [])) for r in tg.to_dict("records")]
    jobs = [j for j in jobs if not (HO / f"{sig_id(j[0])}.json").exists()]
    print(f"{len(jobs)} signals to run, {a.workers} processes", flush=True)
    with cf.ProcessPoolExecutor(a.workers) as ex:
        for k, (name, msg) in enumerate(ex.map(_worker, jobs, chunksize=1)):
            print(f"[{k + 1}/{len(jobs)}] {name}: {msg}", flush=True)


def _inj_worker(args):
    row, plan, inj = args
    path = INJ / f"{inj['id']}.json"
    if path.exists():
        return inj["id"], "cached"
    try:
        r = holdout_signal(row, plan, inject=inj)
        path.write_text(json.dumps(r, default=js))
        rr = r.get("real", {})
        return inj["id"], f"{r.get('mode')} z_hold {rr.get('hold', {}).get('z', np.nan):.2f} pass {rr.get('pass', {}).get('all')}"
    except Exception as exc:
        (INJ / f"{inj['id']}.err").write_text(traceback.format_exc())
        return inj["id"], f"error {exc!r}"[:150]


def cmd_inject(a):
    INJ.mkdir(parents=True, exist_ok=True)
    tg = pd.read_csv(OUT / "targets.csv")
    plans = load_plans(tg)
    rng = np.random.default_rng(8)
    import re
    done = {int(m) for m in re.findall(r"\{'tic': (\d+)", "".join(p.read_text() for p in OUT.glob("fetch*.log")))}
    jobs = []
    for i, r in enumerate(tg.to_dict("records")):
        if a.n and i >= a.n:
            break
        P = float(r["period"]) * float(rng.choice([-1, 1]) * rng.uniform(0.08, 0.3) + 1)
        P = float(np.clip(P, 0.55, 29.5))
        depth = max(float(r["depth_ppm"]), 10.0) * 1e-6
        rp = float(np.sqrt(depth) * r["rstar"] * 109.1)
        inj = {"id": f"inj_{sig_id(r)}", "period": P, "rp_earth": rp, "b": 0.3,
               "phase": float(rng.uniform(0, 1)), "depth_src_ppm": depth * 1e6}
        if int(r["tic"]) in done:      # draws are made for every row, so ids/periods do not depend on timing
            jobs.append((r, plans.get(int(r["tic"]), []), inj))
    jobs = [j for j in jobs if not (INJ / f"{j[2]['id']}.json").exists()]
    print(f"{len(jobs)} injections", flush=True)
    with cf.ProcessPoolExecutor(a.workers) as ex:
        for k, (name, msg) in enumerate(ex.map(_inj_worker, jobs, chunksize=1)):
            print(f"[{k + 1}/{len(jobs)}] {name}: {msg}", flush=True)


def cmd_collect(a):
    from scipy import stats
    rows = []
    for f in sorted(HO.glob("*.json")):
        r = json.load(open(f))
        row = r.get("row", {})
        base = {"id": f.stem, "tic": r["tic"], "stream": row.get("stream"), "author": row.get("author"),
                "P_src": row.get("period"), "snr_src": row.get("snr_src"), "E_hold": row.get("E_hold"),
                "teff": row.get("teff"), "rstar": row.get("rstar"), "tmag": row.get("tmag"), "mode": r.get("mode"),
                "n_his": r.get("n_his_sectors"), "n_new": r.get("n_new_sectors")}
        if "real" not in r:
            rows.append({**base, "error": r.get("error"), "verdict": "NO DATA"})
            continue
        rr = r["real"]
        h, H, F = rr["his"], rr["hold"], rr["full"]
        nl = r.get("null", {})
        cons = H.get("consistency_sigma")
        if rr["pass"]["all"]:
            verdict = "SURVIVOR"
        elif cons is not None and np.isfinite(cons) and cons < -3:
            verdict = "NOT REPRODUCED"
        else:
            verdict = "INCONCLUSIVE"
        if r["mode"] == "fixed":
            p_noise = nl.get("rate_z25_cons")
        else:   # P(a pure-noise z is the max of n+1 exchangeable window draws AND >= 4)
            zs = np.asarray(nl.get("zs", []), float)
            zs = zs[np.isfinite(zs)]
            if len(zs):
                q = (float(np.sum(zs >= 4)) + 0.5) / (len(zs) + 1)     # never exactly 0
                p_noise = (1.0 / (len(zs) + 1)) * (1 - (1 - q) ** (len(zs) + 1))
            else:
                p_noise = None
        inv = r.get("inverted", {})
        rows.append({**base,
                     "snr_his_ours": h.get("snr"), "depth_his_ppm": h.get("depth_ppm"), "err_his_ppm": h.get("err_ppm"),
                     "depth_hold_ppm": H.get("depth_ppm"), "err_hold_ppm": H.get("err_ppm"), "z_hold": H.get("z"),
                     "z_hold_expected": (h["depth_ppm"] / H["err_ppm"]) if H.get("err_ppm") else np.nan,
                     "n_tr_hold": H.get("n_tr"), "consistency_sigma": cons,
                     "P_full": F.get("P"), "t0_full": F.get("t0"), "dur_h": F.get("dur", np.nan) * 24,
                     "snr_full": F.get("snr"), "depth_full_ppm": F.get("depth_ppm"), "rp_full": F.get("rp_earth"),
                     "n_tr_full": F.get("n_transits"), "oe": F.get("oe_sigma"), "sec": F.get("sec_sigma"),
                     "season_p": F.get("season_chi2_p"), "null_p": nl.get("p_z"), "null_zmax": nl.get("z_max"),
                     "p_noise": p_noise, "inv_pass": (inv.get("pass") or {}).get("all"),
                     "inv_z_hold": (inv.get("hold") or {}).get("z"), "n_qlp_timing": r.get("n_qlp_timing"),
                     "rot_active": (r.get("rotation") or {}).get("active"),
                     "pass_flags": ";".join(k for k, v in rr["pass"].items() if not v and k != "all"),
                     "verdict": verdict})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "holdout.csv", index=False)
    print(df.verdict.value_counts().to_string())
    ok = df[df.verdict != "NO DATA"]
    for mode, g in ok.groupby("mode"):
        z = g.z_hold.to_numpy(float)
        zi = g.inv_z_hold.to_numpy(float)
        print(f"[{mode}] n={len(g)} survivors={int((g.verdict == 'SURVIVOR').sum())} expected from noise "
              f"(sum p_noise)={g.p_noise.sum():.2f}  inverted passing={int(g.inv_pass.fillna(False).sum())}/{g.inv_pass.notna().sum()}")
        print(f"   z_hold mean {np.nanmean(z):.2f} std {np.nanstd(z):.2f}; inverted mean {np.nanmean(zi):.2f} std {np.nanstd(zi):.2f}"
              + (f"; KS vs N(0,1) p {stats.kstest(z[np.isfinite(z)], 'norm').pvalue:.3g}" if mode == "fixed" else ""))
    cols = ["id", "mode", "snr_src", "z_hold", "z_hold_expected", "null_zmax", "consistency_sigma", "snr_full",
            "rp_full", "n_tr_full", "verdict", "pass_flags"]
    print(ok.sort_values("z_hold", ascending=False)[cols].head(40).round(2).to_string())


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("fetch"); s.add_argument("--workers", type=int, default=4)
    s = sp.add_parser("run"); s.add_argument("--workers", type=int, default=8); s.add_argument("--only", type=int, nargs="*")
    s.add_argument("--fetched", action="store_true")
    s = sp.add_parser("inject"); s.add_argument("--workers", type=int, default=8); s.add_argument("--n", type=int, default=0)
    sp.add_parser("collect")
    a = p.parse_args()
    {"fetch": cmd_fetch, "run": cmd_run, "inject": cmd_inject, "collect": cmd_collect}[a.cmd](a)


if __name__ == "__main__":
    os.chdir(config.ROOT)
    main()
