"""Round 7, stage 2: re-measure Eumemic's weak signals with every sector, and a hold-out test.

    .venv/bin/python round7/revet.py fetch  [--workers 6]   # 2-min + FFI light curves for targets.csv
    .venv/bin/python round7/revet.py run    [--workers 6]   # per-target JSON in runs/round7_eumemic_weak/revet/
    .venv/bin/python round7/revet.py collect                # revet.csv + survivors + noise count

Data per star: SPOC 2-min PDCSAP where it exists, else TESS-SPOC FFI, else QLP
(the same preference as Eumemic). QLP files whose error column is empty are KEPT
(his reader dropped them: ~900k light curves in S94-97 and S101-104).

Each sector is labelled
  his  = data he could have searched (sector <= his cutoff, not an empty-error QLP file)
  new  = data he did not have (empty-error QLP in S94-97/S101-104, or sectors after his cutoff).

Pre-registered test (fixed before any result was looked at, 2026-10-13):
  1. narrow BLS (his period +-0.5%) on the HIS data alone -> ephemeris E_his;
  2. HOLD-OUT: depth at E_his in the NEW data alone -> z_hold (red-noise errors);
  3. narrow BLS on ALL data -> SNR_full, plus odd/even, secondary, season tests;
  4. survivor = z_hold >= 2.5, hold-out depth consistent with the his-data depth
     (< 3 sigma apart), SNR_full >= 8 within 0.1% of his period, |odd/even| < 3,
     secondary < 3 sigma, season chi2 p > 0.001.
Noise:
  a. 300 random ephemerides per star (P within +-5% of his but > 0.5% away, random
     T0) measured in the same NEW data -> empirical p-value of z_hold, and the
     expected number of survivors from noise = sum of per-star null pass rates;
  b. the WHOLE procedure (1-4) on the inverted light curve of every star.
QLP S99 (second half) to S104 carry a known timing error (minutes; MAST notice),
recorded per target as n_qlp_timing.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from planet_hunt import config  # noqa: E402
from planet_hunt import ffi  # noqa: E402
from planet_hunt import hostsearch as hs  # noqa: E402
from planet_hunt import multisector as ms  # noqa: E402

OUT = config.RUNS / "round7_eumemic_weak"
RV = OUT / "revet"
R5 = config.CACHE / "r5"
LC2 = R5 / "lc"
ARCHIVE = "https://archive.stsci.edu/missions/tess/tid"
DROPPED_QLP = set(range(94, 98)) | set(range(101, 105))
QLP_TIMING = set(range(99, 105))
MIN_FREE_GB = 20.0
NULL_N = 300


def free_gb():
    import shutil
    return shutil.disk_usage("/").free / 1e9


# --------------------------------------------------------------------------
# fetch
# --------------------------------------------------------------------------
def _curl(url, path, retries=3):
    for k in range(retries):
        r = subprocess.run(["curl", "-s", "-f", "-L", "-m", "300", "-o", str(path), url],
                           capture_output=True, text=True)
        if r.returncode == 0 and path.exists() and path.stat().st_size > 5000:
            return True
        path.unlink(missing_ok=True)
        if r.returncode == 22:
            return False
        time.sleep(4 * (k + 1))
    return False


def plan_star(tic: int, twomin: pd.DataFrame, prod: pd.DataFrame) -> list[dict]:
    """One product per sector: SPOC 2-min > TESS-SPOC > QLP."""
    rows = {}
    for r in twomin[twomin.tic == tic].itertuples():
        rows[int(r.sector)] = {"sector": int(r.sector), "prov": "SPOC-2min", "fname": r.fname}
    p = prod[(prod.tic == tic)].dropna(subset=["dataURL"])
    for prov in ("TESS-SPOC", "QLP"):
        for r in p[p.provenance_name == prov].itertuples():
            s = int(r.sequence_number)
            if s not in rows:
                rows[s] = {"sector": s, "prov": prov, "url": r.dataURL}
    return [rows[s] for s in sorted(rows)]


def fetch_ffi(tic, s, url, prov):
    out = ffi._npz_path(tic, s)
    if out.exists():
        return True
    out.parent.mkdir(parents=True, exist_ok=True)
    ffi.TMP_DIR.mkdir(parents=True, exist_ok=True)
    tmp = ffi.TMP_DIR / f"r7_{tic}_{s}_{prov}.fits"
    ok = _curl(config.MAST_DOWNLOAD + url, tmp)
    if not ok:
        return False
    try:
        z = ffi.extract_fits(tmp, prov)
        if len(z["time"]) < 50:
            return False
        meta = z.pop("meta")
        np.savez_compressed(out, **z, **{f"meta_{k}": v for k, v in meta.items()})
        return True
    except Exception:
        return False
    finally:
        tmp.unlink(missing_ok=True)


def fetch_target(tic: int, plan: list[dict]) -> dict:
    import round5_hosts as r5
    two = [p["fname"] for p in plan if p["prov"] == "SPOC-2min"]
    st = {"tic": tic}
    if two and not (LC2 / f"{tic}.npz").exists():
        st["2min"] = r5.fetch_star(tic, two)["status"]
    n_ok = 0
    for p in plan:
        if p["prov"] == "SPOC-2min":
            continue
        n_ok += fetch_ffi(tic, p["sector"], p["url"], p["prov"])
    st["ffi_ok"] = n_ok
    st["ffi_planned"] = sum(p["prov"] != "SPOC-2min" for p in plan)
    return st


def cmd_fetch(a):
    tg = pd.read_csv(OUT / "targets.csv")
    two = pd.read_csv(R5 / "twomin_index.csv.gz")
    two = two[two.tic.isin(set(tg.tic))].drop_duplicates(["tic", "sector"])
    prod = pd.read_csv(config.TABLE_DIR / "lcprod.csv")
    prod = prod[prod.tic.isin(set(tg.tic))]
    plans = {int(t): plan_star(int(t), two, prod) for t in tg.tic.unique()}
    json.dump(plans, open(OUT / "plans.json", "w"))
    print(f"{len(plans)} stars, {sum(len(v) for v in plans.values())} sector files planned; free {free_gb():.1f} GB",
          flush=True)
    t_start = time.time()
    with cf.ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(fetch_target, t, p): t for t, p in plans.items()}
        for k, f in enumerate(cf.as_completed(futs)):
            try:
                r = f.result()
            except Exception as exc:
                r = {"tic": futs[f], "err": repr(exc)[:100]}
            print(f"[{k + 1}/{len(plans)}] {r} free {free_gb():.1f} GB {time.time() - t_start:.0f}s", flush=True)
            if free_gb() < MIN_FREE_GB:
                print("DISK GUARD: < 20 GB free, stopping", flush=True)
                for g in futs:
                    g.cancel()
                break


# --------------------------------------------------------------------------
# load
# --------------------------------------------------------------------------
def load_sectors(tic: int, plan: list[dict], his_cut: int) -> list[dict]:
    import round5_hosts as r5
    secs = []
    if (LC2 / f"{tic}.npz").exists():
        for d in r5.load_lc(tic):
            d["meta"]["provenance"] = "SPOC-2min"
            d["meta"]["empty_err"] = False
            secs.append(d)
    have = {d["meta"]["sector"] for d in secs}
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
        s = int(m["sector"])
        dropped = m.get("provenance") == "QLP" and m["empty_err"] and s in DROPPED_QLP
        m["his"] = bool(s <= his_cut and not dropped)
        m["new_reason"] = "" if m["his"] else ("dropped_qlp" if dropped else "after_cut")
    return sorted(secs, key=lambda d: d["meta"]["sector"])


# --------------------------------------------------------------------------
# measurement helpers
# --------------------------------------------------------------------------
def sector_beta(lc, dur):
    oot = np.ones(len(lc["t"]), bool)
    return ms.red_noise_beta(lc, dur, oot)


def depth_at(lc, beta, P, t0, dur, sel=None):
    """Weighted depth and red-noise error at one ephemeris (subset sel)."""
    t, y = lc["t"], lc["flat"]
    if sel is None:
        sel = np.ones(len(t), bool)
    ph = np.abs(ms.phase(t, P, t0))
    intr = sel & (ph < 0.5 * dur)
    oot = sel & (ph > dur) & (ph < 3 * dur)
    if intr.sum() < 3 or oot.sum() < 10:
        return np.nan, np.nan, 0
    w = 1 / (lc["dy"] * beta) ** 2
    lev = np.sum(w[oot] * y[oot]) / w[oot].sum()
    d = lev - np.sum(w[intr] * y[intr]) / w[intr].sum()
    e = np.sqrt(1 / w[intr].sum() + 1 / w[oot].sum())
    ep = np.unique(np.round((t[intr] - t0) / P))
    return float(d), float(e), int(len(ep))


def narrow_bls(lc, sel, P0, rstar, mstar, frac=0.005):
    """Coherent BLS on 10-min bins of the subset, periods within +-frac of P0."""
    from astropy.timeseries import BoxLeastSquares
    t, y, dy = lc["t"][sel], lc["flat"][sel], lc["dy"][sel]
    if len(t) < 100:
        return None
    tb, yb, eb = ms.bin_weighted(t, y, dy, minutes=10.0)
    T = tb.max() - tb.min()
    periods = ms.period_grid(P0 * (1 - frac), P0 * (1 + frac), T, rstar, mstar, oversample=3.0)
    if len(periods) < 5:
        periods = np.linspace(P0 * (1 - frac), P0 * (1 + frac), 200)
    durs = ms.durations_for(P0, rstar, mstar)
    durs = durs[durs < 0.5 * P0 * (1 - frac)]
    bls = BoxLeastSquares(tb, yb, eb)
    r = bls.power(periods, durs, objective="snr")
    i = int(np.nanargmax(r.power))
    return {"P": float(r.period[i]), "t0": float(r.transit_time[i]), "dur": float(r.duration[i]),
            "bls_snr": float(r.depth_snr[i]), "n_periods": len(periods)}


def full_stats(lc, P, t0, dur, rstar, mstar, teff, sel=None):
    sub = hs.subset(lc, sel) if sel is not None else lc
    s = ms.transit_stats(sub, P, t0, dur, rstar, mstar, teff)
    s.pop("per_transit", None)
    return s


def run_procedure(lc, his_sel, new_sel, P0, rstar, mstar, teff):
    """Steps 1-3 of the pre-registered test on one light curve (real or inverted)."""
    out = {}
    e_his = narrow_bls(lc, his_sel, P0, rstar, mstar)
    if e_his is None:
        return {"error": "no his data"}
    beta = sector_beta(lc, e_his["dur"])
    d_h, e_h, n_h = depth_at(lc, beta, e_his["P"], e_his["t0"], e_his["dur"], his_sel)
    out["his"] = {**e_his, "depth_ppm": d_h * 1e6, "err_ppm": e_h * 1e6, "snr": d_h / e_h if e_h else np.nan,
                  "n_tr": n_h}
    if new_sel.sum() > 50:
        d_n, e_n, n_n = depth_at(lc, beta, e_his["P"], e_his["t0"], e_his["dur"], new_sel)
    else:
        d_n, e_n, n_n = np.nan, np.nan, 0
    out["hold"] = {"depth_ppm": d_n * 1e6, "err_ppm": e_n * 1e6, "z": d_n / e_n if e_n else np.nan, "n_tr": n_n,
                   "consistency_sigma": float((d_n - d_h) / np.hypot(e_n, e_h)) if np.isfinite(d_n) else np.nan}
    allsel = np.ones(len(lc["t"]), bool)
    e_all = narrow_bls(lc, allsel, P0, rstar, mstar)
    st = full_stats(lc, e_all["P"], e_all["t0"], e_all["dur"], rstar, mstar, teff)
    out["full"] = {**e_all, **{k: st[k] for k in ("depth_ppm", "depth_err_ppm", "snr", "snr_white", "beta_med",
                                                   "n_transits", "oe_sigma", "sec_sigma", "sec_any_sigma",
                                                   "season_chi2_p", "tr_chi2_p", "rp_earth", "dur_ratio",
                                                   "sin_ratio", "shape", "tr_max_over_mean")},
                   "season_depths": st["season_depths"]}
    out["pass"] = passes(out, P0)
    out["_beta"] = beta
    out["_e_his"] = e_his
    return out


def passes(o, P0) -> dict:
    h, f = o.get("hold", {}), o.get("full", {})
    c = {
        "z_hold": bool(np.nan_to_num(h.get("z"), nan=-9) >= 2.5),
        "consistent": bool(abs(np.nan_to_num(h.get("consistency_sigma"), nan=9)) < 3),
        "snr_full": bool(np.nan_to_num(f.get("snr"), nan=0) >= 8 and abs(f.get("P", 0) / P0 - 1) < 0.001),
        "oe": bool(abs(np.nan_to_num(f.get("oe_sigma"), nan=0)) < 3),
        "sec": bool(np.nan_to_num(f.get("sec_sigma"), nan=0) < 3),
        "season": bool(np.nan_to_num(f.get("season_chi2_p"), nan=1) > 0.001),
    }
    c["all"] = all(c.values())
    return c


def null_ephemerides(lc, beta, his_sel, new_sel, P_his, dur, rng, n=NULL_N):
    """z_hold at random ephemerides (same star, same new data, same duration)."""
    if new_sel.sum() < 50:
        return np.array([]), np.array([])
    t = lc["t"]
    zs, cons = [], []
    for _ in range(n):
        while True:
            dp = rng.uniform(-0.05, 0.05)
            if abs(dp) > 0.005:
                break
        P = P_his * (1 + dp)
        t0 = t.min() + rng.uniform(0, P)
        d_n, e_n, _ = depth_at(lc, beta, P, t0, dur, new_sel)
        d_h, e_h, _ = depth_at(lc, beta, P, t0, dur, his_sel)
        zs.append(d_n / e_n if e_n else np.nan)
        cons.append((d_n - d_h) / np.hypot(e_n, e_h) if np.isfinite(d_n + d_h) else np.nan)
    return np.array(zs), np.array(cons)


# --------------------------------------------------------------------------
# per target
# --------------------------------------------------------------------------
def revet_target(row: dict, plan: list[dict], inject: dict | None = None) -> dict:
    import round5_hosts as r5
    tic = int(row["tic"])
    rstar, mstar, teff = float(row["rstar"]), float(row["mstar"]), float(row["teff"])
    P0 = float(row["period"])
    secs = load_sectors(tic, plan, int(row["his_cut"]))
    if not secs:
        return {"tic": tic, "error": "no data"}
    if inject is not None:
        secs = hs.inject_transit(secs, inject, rstar, mstar, teff)
        P0 = inject["period"]
    toi, ps = r5.load_catalogs()
    known = r5.known_ephemerides(tic, toi, ps)
    kl = hs.known_list(known, rstar, mstar)
    kmask = lambda t: hs.mask_known(t, kl)
    rot = hs.rotation(secs, kmask)
    window = hs.W_SHORT if P0 < hs.SPLIT_P else hs.W_LONG
    lc = hs.detrend(secs, window, kmask, rot=rot)
    if not lc:
        return {"tic": tic, "error": "detrend failed"}
    keep = ~kmask(lc["t"])
    lc = hs.subset(lc, keep)
    sec_his = {int(d["meta"]["sector"]): d["meta"]["his"] for d in secs}
    his_sel = np.array([sec_his.get(int(s), False) for s in lc["sector"]])
    new_sel = ~his_sel
    info = [{"sector": int(d["meta"]["sector"]), "prov": d["meta"].get("provenance"), "his": d["meta"]["his"],
             "new_reason": d["meta"]["new_reason"], "empty_err": d["meta"]["empty_err"]} for d in secs]
    res = {"tic": tic, "row": {k: (v.item() if hasattr(v, "item") else v) for k, v in row.items()},
           "sectors": info, "n_his_sectors": int(sum(i["his"] for i in info)),
           "n_new_sectors": int(sum(not i["his"] for i in info)),
           "n_qlp_timing": int(sum(i["prov"] == "QLP" and i["sector"] in QLP_TIMING for i in info)),
           "rotation": {k: rot.get(k) for k in ("active", "P_rot", "saturated") if k in rot},
           "known_masked": [k["name"] for k in kl if k["masked"]]}
    real = run_procedure(lc, his_sel, new_sel, P0, rstar, mstar, teff)
    if "error" in real:
        res["error"] = real["error"]
        return res
    if inject is not None:
        real.pop("_beta"); real.pop("_e_his")
        res["inject"] = inject
        res["real"] = real
        return res
    beta, e_his = real.pop("_beta"), real.pop("_e_his")
    res["real"] = real
    rng = np.random.default_rng(tic % (2 ** 32))
    zs, cons = null_ephemerides(lc, beta, his_sel, new_sel, e_his["P"], e_his["dur"], rng)
    if len(zs):
        z_obs = real["hold"]["z"]
        res["null"] = {"n": int(np.isfinite(zs).sum()),
                       "p_z": float(np.mean(zs >= z_obs)) if np.isfinite(z_obs) else None,
                       "rate_z25": float(np.mean(zs >= 2.5)),
                       "rate_z25_cons": float(np.mean((zs >= 2.5) & (np.abs(cons) < 3))),
                       "z_std": float(np.nanstd(zs)), "z_mean": float(np.nanmean(zs))}
    inv = dict(lc)
    inv["flat"] = 2.0 - lc["flat"]
    ir = run_procedure(inv, his_sel, new_sel, P0, rstar, mstar, teff)
    ir.pop("_beta", None); ir.pop("_e_his", None)
    res["inverted"] = ir
    return res


def _js(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return str(o)


def _worker(args):
    row, plan = args
    path = RV / f"{int(row['tic'])}_{row['period']:.4f}.json"
    if path.exists():
        return path.stem, "cached"
    t0 = time.time()
    try:
        r = revet_target(row, plan)
        r["elapsed_s"] = time.time() - t0
        path.write_text(json.dumps(r, default=_js))
        rr = r.get("real", {})
        return path.stem, (f"z_hold {rr.get('hold', {}).get('z', np.nan):.2f} snr_full "
                           f"{rr.get('full', {}).get('snr', np.nan):.1f} pass {rr.get('pass', {}).get('all')} "
                           f"inv {r.get('inverted', {}).get('pass', {}).get('all')} {time.time() - t0:.0f}s"
                           if rr else r.get("error"))
    except Exception as exc:
        (RV / f"{path.stem}.err").write_text(traceback.format_exc())
        return path.stem, f"error {exc!r}"[:150]


def _inj_worker(args):
    row, plan, inj = args
    d = RV.parent / "inject"
    path = d / f"{inj['id']}.json"
    if path.exists():
        return inj["id"], "cached"
    try:
        r = revet_target(row, plan, inject=inj)
        path.write_text(json.dumps(r, default=_js))
        rr = r.get("real", {})
        return inj["id"], f"his snr {rr.get('his', {}).get('snr', np.nan):.1f} z_hold {rr.get('hold', {}).get('z', np.nan):.2f} pass {rr.get('pass', {}).get('all')}"
    except Exception as exc:
        (d / f"{inj['id']}.err").write_text(traceback.format_exc())
        return inj["id"], f"error {exc!r}"[:150]


def cmd_inject(a):
    """Inject a transit with HIS signal's depth and duration scale at another period
    into each target's raw flux; run the identical procedure (completeness)."""
    d = RV.parent / "inject"
    d.mkdir(parents=True, exist_ok=True)
    tg = pd.read_csv(OUT / "targets.csv")
    plans = {int(k): v for k, v in json.load(open(OUT / "plans.json")).items()}
    rng = np.random.default_rng(7)
    jobs = []
    for i, r in enumerate(tg.to_dict("records")):
        P = float(r["period"]) * float(rng.choice([-1, 1]) * rng.uniform(0.08, 0.3) + 1)
        P = float(np.clip(P, 0.55, 29.5))
        rp = float(np.sqrt(max(r["depth"], 1e-5)) * r["rstar"] * 109.1)
        inj = {"id": f"inj_{int(r['tic'])}_{i:03d}", "period": P, "rp_earth": rp, "b": 0.3,
               "phase": float(rng.uniform(0, 1)), "depth_eum_ppm": float(r["depth"]) * 1e6}
        jobs.append((r, plans.get(int(r["tic"]), []), inj))
    print(f"{len(jobs)} injections", flush=True)
    with cf.ProcessPoolExecutor(a.workers) as ex:
        for k, (name, msg) in enumerate(ex.map(_inj_worker, jobs, chunksize=1)):
            print(f"[{k + 1}/{len(jobs)}] {name}: {msg}", flush=True)


def cmd_run(a):
    RV.mkdir(parents=True, exist_ok=True)
    tg = pd.read_csv(OUT / "targets.csv")
    plans = {int(k): v for k, v in json.load(open(OUT / "plans.json")).items()}
    if a.only:
        tg = tg[tg.tic.isin(a.only)]
    jobs = [(r, plans.get(int(r["tic"]), [])) for r in tg.to_dict("records")]
    print(f"{len(jobs)} targets, {a.workers} processes", flush=True)
    with cf.ProcessPoolExecutor(a.workers) as ex:
        for k, (name, msg) in enumerate(ex.map(_worker, jobs, chunksize=1)):
            print(f"[{k + 1}/{len(jobs)}] {name}: {msg}", flush=True)


def cmd_collect(a):
    rows = []
    for f in sorted(RV.glob("*.json")):
        r = json.load(open(f))
        if "real" not in r:
            rows.append({"name": f.stem, "tic": r["tic"], "error": r.get("error")})
            continue
        rr, row = r["real"], r["row"]
        h, H, F = rr["his"], rr["hold"], rr["full"]
        exp_z = h["depth_ppm"] / H["err_ppm"] if H.get("err_ppm") else np.nan
        cons = H.get("consistency_sigma")
        if rr["pass"]["all"]:
            verdict = "SURVIVOR"
        elif cons is not None and np.isfinite(cons) and cons < -3:
            verdict = "REFUTED"
        else:
            verdict = "INCONCLUSIVE"
        inv = r.get("inverted", {})
        nl = r.get("null", {})
        rows.append({
            "name": f.stem, "tic": r["tic"], "src": row["src"], "P_his": row["period"], "snr_eum": row["snr"],
            "n_sec_eum": row["n_sectors"], "n_his_ours": r["n_his_sectors"], "n_new_ours": r["n_new_sectors"],
            "teff": row["teff"], "rstar": row["rstar"], "tmag": row["tmag"],
            "snr_his_ours": h["snr"], "depth_his_ppm": h["depth_ppm"], "err_his_ppm": h["err_ppm"],
            "depth_hold_ppm": H["depth_ppm"], "err_hold_ppm": H["err_ppm"], "z_hold": H["z"], "z_hold_expected": exp_z,
            "n_tr_hold": H["n_tr"], "consistency_sigma": cons,
            "P_full": F["P"], "t0_full": F["t0"], "dur_h": F["dur"] * 24, "snr_full": F["snr"],
            "depth_full_ppm": F["depth_ppm"], "rp_full": F["rp_earth"], "n_tr_full": F["n_transits"],
            "oe": F["oe_sigma"], "sec": F["sec_sigma"], "season_p": F["season_chi2_p"],
            "null_p": nl.get("p_z"), "null_rate": nl.get("rate_z25_cons"), "null_zstd": nl.get("z_std"),
            "inv_pass": inv.get("pass", {}).get("all"), "inv_z_hold": inv.get("hold", {}).get("z"),
            "inv_snr_full": inv.get("full", {}).get("snr"),
            "n_qlp_timing": r["n_qlp_timing"], "rot_active": r["rotation"].get("active"),
            "pass_flags": ";".join(k for k, v in rr["pass"].items() if not v and k != "all"),
            "verdict": verdict})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "revet.csv", index=False)
    ok = df.dropna(subset=["z_hold"]) if "z_hold" in df else df
    print(df.verdict.value_counts(dropna=False).to_string())
    print(f"targets with a hold-out measurement: {len(ok)}")
    print(f"expected survivors from noise (sum of per-star null pass rates of z_hold>=2.5 & consistent): "
          f"{ok.null_rate.sum():.2f}  (Gaussian 0.62% x N = {0.0062 * len(ok):.2f})")
    print(f"inverted light curves passing the whole procedure: {int(df.inv_pass.fillna(False).sum())} of "
          f"{df.inv_pass.notna().sum()}")
    from scipy import stats
    z = ok.z_hold.to_numpy(float)
    print(f"hold-out z: mean {np.nanmean(z):.2f}, std {np.nanstd(z):.2f}, n(z>=2.5) {(z >= 2.5).sum()}, "
          f"n(z>=2) {(z >= 2).sum()}; KS vs N(0,1) p {stats.kstest(z[np.isfinite(z)], 'norm').pvalue:.3g}")
    zi = df.inv_z_hold.to_numpy(float)
    print(f"inverted hold-out z: mean {np.nanmean(zi):.2f}, std {np.nanstd(zi):.2f}, n(z>=2.5) {(zi >= 2.5).sum()}")
    cols = ["name", "src", "snr_eum", "snr_his_ours", "z_hold", "z_hold_expected", "consistency_sigma",
            "snr_full", "rp_full", "n_tr_full", "null_p", "verdict", "pass_flags"]
    print(df.sort_values("z_hold", ascending=False)[cols].head(25).round(2).to_string())


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("fetch"); s.add_argument("--workers", type=int, default=6)
    s = sp.add_parser("run"); s.add_argument("--workers", type=int, default=6); s.add_argument("--only", type=int, nargs="*")
    sp.add_parser("collect")
    s = sp.add_parser("inject"); s.add_argument("--workers", type=int, default=6)
    a = p.parse_args()
    {"fetch": cmd_fetch, "run": cmd_run, "collect": cmd_collect, "inject": cmd_inject}[a.cmd](a)


if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parent.parent)
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    main()
