"""Method B: hidden planets in archival HARPS radial velocities (RVBank).

    .venv/bin/python rv_hunt.py --validate            # HD 32564 first (the benchmark)
    .venv/bin/python rv_hunt.py --min-nights 50       # every RVBank star with >= 50 nights
    .venv/bin/python rv_hunt.py --star HD10180        # one star

Idea credited to Mehran Jalali's public thread
(https://x.com/mehran__jalali/status/2108339123669111033): three planets
of 9-13 Earth masses (2.6, 11.1, 23.6 d, plus a 49.1 d candidate) around
HD 32564 in 191 nights of public HARPS data from 2009-2016.

For each star: nightly RVs (SERVAL, NZP-corrected; pre/post 2015 fibre
upgrade as separate instruments) -> GLS + iterative pre-whitening with joint
Keplerian fits, stopping when the next peak has FAP > --fap (Baluev) ->
vetting (activity indicators, aliases, coherence across the two halves of the
data) -> cross-match with the NASA Exoplanet Archive (pscomppars + ps).
Outputs runs/rv/{results.csv, stars.csv, star/<name>.npz|json}.
See planet_hunt/rv.py for the details of each step.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import time

import numpy as np
import pandas as pd

from planet_hunt import config
from planet_hunt import rv

RUN = config.RUNS / "rv"


def analyse_star(job):
    name, star_rows, info, planets, fap_max, max_signals = job
    out_rows, meta = [], {"name": name}
    try:
        n = rv.nightly(star_rows)
        if n.empty:
            meta["error"] = "no usable RVs"
            return out_rows, meta
        meta.update(n_nights=len(n), n_pre=int((n["inst"] == 0).sum()),
                    n_post=int((n["inst"] == 1).sum()))
        if len(n) < 20:
            return out_rows, meta
        T = float(n["bjd"].max() - n["bjd"].min())
        model, steps = rv.prewhiten(n, fap_max=fap_max, max_signals=max_signals)
        model = rv.try_eccentric(model)
        errs = model.param_errors()
        r, _ = model.residuals()
        meta.update(T_days=T, rms_raw=float(np.std(n["rv"])), rms_resid=float(np.std(r)),
                    jitter={str(k): v for k, v in model.jit.items()}, eccentric=model.eccentric,
                    slope_ms_per_yr=float(model.lin[-1]),
                    steps=[{k: v for k, v in s.items() if k not in ("freq", "pow")} for s in steps])
        teff = info.get("Teff", np.nan)
        mstar = rv.mass_from_teff(teff)
        absG = (info.get("Gmag", np.nan) - 5 * np.log10(info.get("Dist", np.nan) / 10)
                if info.get("Dist", 0) and np.isfinite(info.get("Dist", np.nan)) else np.nan)
        evolved = bool(np.isfinite(absG) and absG < 3.0 and np.isfinite(teff) and teff < 5600)
        pg = rv.indicator_periodograms(n, T)
        earlier = []
        for k, s in enumerate(model.signals):
            e = s["sec"] ** 2 + s["ses"] ** 2
            sig_curve = rv.keplerian(model.t, s["P"], s["K"], s["sec"], s["ses"], s["lam"], model.tref)
            act = rv.activity_check(pg, s["P"], sig_curve, T)
            coh = rv.coherence_check(model, k)
            alias = rv.alias_flags(s["P"], earlier, T)
            earlier.append(s["P"])
            step = steps[k] if k < len(steps) else {}
            perr = errs[k]["P_err"] if k < len(errs) else np.nan
            tol = float(np.clip(3 * perr / s["P"], 0.05, 0.25)) if np.isfinite(perr) else 0.05
            known = rv.match_known(s["P"], planets, tol=tol) if len(planets) else None
            status = ("known_planet" if known else
                      "host_other" if len(planets) else "new")
            red = [f for f in act["flags"] if "same period" in f or "correlates" in f]
            amber = [f for f in act["flags"] if f not in red]
            red += [f for f in alias if "harmonic" not in f]
            amber += [f for f in alias if "harmonic" in f]
            flags = red
            if not coh["coherent"]:
                flags.append(f"not coherent: K = {coh['K_first']:.1f} vs {coh['K_second']:.1f} m/s in "
                             f"the two halves ({coh['diff_sigma']:.1f} sigma)")
            kerr = errs[k]["K_err"] if k < len(errs) else np.nan
            if e > 0.5 and np.isfinite(kerr) and kerr > s["K"] / 3:
                flags.append(f"unconstrained eccentric fit (e = {e:.2f}, K = {s['K']:.1f} ± {kerr:.1f} m/s): "
                             "a few outlying nights, not an orbit")
            if evolved:
                flags.append("evolved star (oscillations/granulation/rotation can mimic signals)")
            dl = step.get("alias_dlnL") or {}
            ambiguous = [f"{nm} at {step['alias_periods'][nm]:.3f} d (dlnL {v:.1f})"
                         for nm, v in dl.items() if v < 5]
            if ambiguous:
                amber.append("alias ambiguity: " + "; ".join(ambiguous))
            prot_lo, prot_hi = (8, 150) if (np.isfinite(teff) and teff < 4000) else (8, 60)
            if prot_lo < s["P"] < prot_hi:
                amber.append(f"period within the typical rotation range of this star type ({prot_lo}-{prot_hi} d)")
            out_rows.append({
                "star": name, "signal": k + 1, "period": s["P"],
                "period_err": errs[k]["P_err"] if k < len(errs) else np.nan,
                "K": s["K"], "K_err": errs[k]["K_err"] if k < len(errs) else np.nan, "ecc": e,
                "msini_earth": rv.min_mass_earth(s["K"], s["P"], mstar, e), "mstar_est": mstar,
                "fap": step.get("fap", np.nan), "gls_power": step.get("power", np.nan),
                "alias_choice": step.get("alias_choice", ""), "alias_ambiguous": "; ".join(ambiguous),
                "n_nights": len(n), "T_days": T, "Teff": teff, "Vmag": info.get("Vmag", np.nan),
                "dist_pc": info.get("Dist", np.nan), "SpType": info.get("SpType", ""),
                "status": status, "known_match": known or "",
                "known_planets": "; ".join(f"{p.pl_name} P={p.pl_orbper:.2f}" for p in
                                           planets.drop_duplicates("pl_name").itertuples()
                                           if np.isfinite(p.pl_orbper))[:300],
                "coh_K1": coh["K_first"], "coh_K2": coh["K_second"], "coh_sigma": coh["diff_sigma"],
                "coherent": coh["coherent"], "flags": "; ".join(flags), "cautions": "; ".join(amber),
                "activity_best": "; ".join(f"{c}:{v['best_period']:.1f}d(FAP {v['best_fap']:.0e})"
                                           for c, v in act["indicators"].items()
                                           if v["best_fap"] < 1e-3),
            })
        # arrays for figures
        (RUN / "star").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(RUN / "star" / f"{name}.npz", t=model.t, rv=model.rv, err=model.err,
                            inst=model.inst, sigma=model.sigma(), lin=model.lin, tref=model.tref,
                            signals=json.dumps(model.signals),
                            pg_freq=steps[0]["freq"] if steps else np.zeros(1),
                            **{f"pg{s['k']}": s["pow"].astype(np.float32) for s in steps})
    except Exception as exc:
        meta["error"] = repr(exc)[:200]
    return out_rows, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-nights", type=int, default=50)
    ap.add_argument("--fap", type=float, default=1e-3)
    ap.add_argument("--max-signals", type=int, default=6)
    ap.add_argument("--star", nargs="*")
    ap.add_argument("--validate", action="store_true", help="only HD 32564")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    args = ap.parse_args()
    t0 = time.time()
    RUN.mkdir(parents=True, exist_ok=True)

    df = rv.add_rhk(rv.load_rvbank())
    print(f"[rv] R'HK attached to {int(df['RHK'].notna().sum())} of {len(df)} RVs", flush=True)
    lst = rv.load_star_list().set_index("Name")
    known = rv.load_known_planets()
    good = df[(df["Flag"] == 0) & np.isfinite(df["DRVmlcnzp"])]
    nights = good.assign(night=np.floor(good["BJD"] - rv.LA_SILLA_NOON)).groupby("Name")["night"].nunique()
    if args.validate:
        names = ["HD32564"]
    elif args.star:
        names = args.star
    else:
        names = sorted(nights[nights >= args.min_nights].index)
    print(f"[rv] {len(df)} RVs, {df['Name'].nunique()} stars; {int((nights >= args.min_nights).sum())} "
          f"with >= {args.min_nights} nights; analysing {len(names)}", flush=True)

    groups = dict(tuple(good[good["Name"].isin(names)].groupby("Name")))
    ids = rv.simbad_ids(list(names))
    jobs = []
    for nm in names:
        if nm not in groups:
            continue
        info = lst.loc[nm].to_dict() if nm in lst.index else {}
        info.setdefault("Name", nm)
        star = pd.Series({"Name": nm, "ra": info.get("ra", np.nan), "dec": info.get("dec", np.nan)})
        planets = rv.host_planets(star, known, ids.get(nm)) if np.isfinite(star["ra"]) else known.iloc[0:0]
        jobs.append((nm, groups[nm], info, planets, args.fap, args.max_signals))

    rows, metas = [], []
    with cf.ProcessPoolExecutor(args.workers) as ex:
        for k, (r, m) in enumerate(ex.map(analyse_star, jobs, chunksize=1), 1):
            rows.extend(r)
            metas.append(m)
            if k % 50 == 0 or k == len(jobs):
                print(f"[rv] {k}/{len(jobs)} stars ({time.time() - t0:.0f} s)", flush=True)
    tag = "validation" if args.validate else ("custom" if args.star else "all")
    res = pd.DataFrame(rows)
    res.to_csv(RUN / f"results_{tag}.csv", index=False)
    (RUN / f"stars_{tag}.json").write_text(json.dumps(metas, indent=1, default=float))
    if len(res):
        print(res[["star", "signal", "period", "K", "msini_earth", "fap", "status", "known_match",
                   "flags"]].head(40).to_string(), flush=True)
        print(res["status"].value_counts().to_string(), flush=True)
    print(f"[rv] done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
