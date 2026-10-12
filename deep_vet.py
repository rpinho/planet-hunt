"""Deep vetting of a few candidates from TESScut pixels (round 4).

    .venv/bin/python deep_vet.py vetting/round4_targets.json            # all targets
    .venv/bin/python deep_vet.py vetting/round4_targets.json --only tic167723381
    .venv/bin/python deep_vet.py vetting/round4_targets.json --assemble vetting/round4_verdicts.json

Each target in the JSON file: name, tic, ra, dec, gaia_id, tmag, teff, rstar,
mstar, period, t0 (BTJD), duration_h, depth_ppm, role (candidate /
control_on = known planet / control_off = known off-target false positive),
optional max_sectors (evenly spaced subset, to cap downloads) and sectors.

Per target (planet_hunt/pixelvet.py, planet_hunt/lcvet.py):
  1. TESScut 11x11 cutouts for every (or a capped subset of) sector;
  2. Gaia DR3 scene, PSF scene fit per sector, transit difference image with
     null-ephemeris noise;
  3. joint PSF localisation over all sectors, the chi^2 of "the dip is on
     star j" for every Gaia star within 90", and the same per sector;
  4. depth in rings of pixels vs the prediction for each source, and the
     eclipse depth each neighbour would need (flux budget);
  5. fixed-aperture (same pixel rule every sector) depths per sector and per
     year/camera, and the same for the pipeline light curves (TESS-SPOC and
     QLP separately) when they are cached;
  6. odd/even, secondary eclipse, and a limb-darkened transit fit with free
     stellar density, compared with the TIC density.
Writes runs/round4_vet/<name>.json and runs/round4_vet/<name>.png;
--assemble merges them with hand-written verdicts into runs/round4_vet.json.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np

from planet_hunt import config
from planet_hunt import lcvet, pixelvet as pv
from planet_hunt import multisector as ms

OUT = config.RUNS / "round4_vet"


def pick_sectors(avail: list[int], spec: dict) -> list[int]:
    s = spec.get("sectors") or avail
    s = [x for x in s if x in avail]
    n = spec.get("max_sectors")
    if n and len(s) > n:
        s = [s[i] for i in np.unique(np.linspace(0, len(s) - 1, n).round().astype(int))]
    return s


def per_sector_sources(sectors, scene, hyp_rows: list[int]) -> list[dict]:
    """Per sector: chi^2 of each candidate source relative to the best of them."""
    rows = []
    xi = scene["xi"].to_numpy()[hyp_rows]
    eta = scene["eta"].to_numpy()[hyp_rows]
    for sd in sectors:
        chi, amps = pv._chi2_at([sd], xi, eta, per_sector=True)
        A, eA = amps[0]
        it = hyp_rows.index(sd["it"]) if sd["it"] in hyp_rows else None
        best = int(np.argmin(chi))
        rows.append({"sector": sd["sector"], "best_source_row": int(hyp_rows[best]),
                     "dchi2_target": float(chi[it] - chi[best]) if it is not None else None,
                     "dchi2": {int(hyp_rows[j]): float(chi[j] - chi[best]) for j in range(len(hyp_rows))},
                     "snr_at_target": float(A[it] / eA[it]) if it is not None else None})
    return rows


def run_target(spec: dict) -> dict:
    t_start = time.time()
    name = spec["name"]
    P, t0, dur = float(spec["period"]), float(spec["t0"]), float(spec["duration_h"]) / 24
    depth = float(spec["depth_ppm"]) * 1e-6
    avail = pv.available_sectors(spec["ra"], spec["dec"])
    want = pick_sectors(avail, spec)
    print(f"[{name}] {len(want)} of {len(avail)} TESScut sectors: {want}", flush=True)
    cuts = pv.fetch_cutouts(spec["ra"], spec["dec"], want)
    scene = pv.gaia_scene(spec["ra"], spec["dec"], gaia_id=spec.get("gaia_id"), tmag=spec.get("tmag"))
    excl = [tuple(w) for w in spec.get("exclude_windows", [])]
    sectors, lcs, dropped = [], [], []
    for s, path in cuts:
        try:
            cut = pv.load_cutout(path)
            sd = pv.prepare_sector(cut, scene, P, t0, dur, exclude_windows=excl)
        except Exception as exc:
            print(f"[{name}] S{s}: failed {exc!r}"[:160], flush=True)
            continue
        if sd is None:
            continue
        if not sd["quality_ok"]:
            dropped.append({"sector": s, "target_frac": sd["target_frac"], "bad_dist_px": sd["bad_dist_px"],
                            "scene_resid": sd["fit"]["frac_resid"]})
            print(f"[{name}] S{s}: dropped by quality gate {dropped[-1]}", flush=True)
            continue
        sectors.append(sd)
        lcs.append(lcvet.cutout_lightcurve(cut, sd))
        del cut
    print(f"[{name}] {len(sectors)} sectors with transits ({time.time() - t_start:.0f} s)", flush=True)
    out = {"spec": spec, "n_sectors": len(sectors), "sectors_used": [sd["sector"] for sd in sectors],
           "sectors_dropped": dropped}
    tgt = scene[scene["is_target"]].iloc[0]
    out["gaia_target"] = {"source_id": int(tgt["source_id"]), "gmag": float(tgt["Gmag"]),
                          "tmag_gaia": float(tgt["tmag"]), "ruwe": float(tgt["RUWE"]) if tgt["RUWE"] == tgt["RUWE"] else None}
    near = scene[(scene["sep"] <= 63) & ~scene["is_target"]]
    out["gaia_neighbours_63as"] = [
        {"source_id": int(r.source_id), "sep_arcsec": round(float(r.sep), 1), "pa_deg": round(float(r.pa_deg)),
         "gmag": round(float(r.Gmag), 2), "tmag": round(float(r.tmag), 2),
         "bp_rp": None if r.bp_rp != r.bp_rp else round(float(r.bp_rp), 2),
         "plx_mas": None if r.Plx != r.Plx else round(float(r.Plx), 2),
         "delta_t": round(float(r.tmag - tgt["tmag"]), 2)} for r in near.itertuples()]

    # localisation
    L = pv.localize(sectors, scene)
    out["localization"] = {k: v for k, v in L.items() if k not in ("grid",)}
    hyp_rows = [int(np.flatnonzero(scene["source_id"].to_numpy() == h["source_id"])[0])
                for h in L["hypotheses"][:6]]
    it = int(np.flatnonzero(scene["is_target"].to_numpy())[0])
    if it not in hyp_rows:
        hyp_rows.append(it)
    out["per_sector_sources"] = per_sector_sources(sectors, scene, hyp_rows)
    ps = out["per_sector_sources"]
    out["per_sector_summary"] = {
        "n": len(ps), "target_best": int(sum(r["best_source_row"] == it for r in ps)),
        "target_disfavoured_dchi2_gt_4": int(sum((r["dchi2_target"] or 0) > 4 for r in ps)),
        "best_counts": {int(k): int(v) for k, v in zip(*np.unique([r["best_source_row"] for r in ps],
                                                                   return_counts=True))}}
    out["source_rows"] = {int(j): {"source_id": int(scene.iloc[j]["source_id"]),
                                   "sep": float(scene.iloc[j]["sep"]),
                                   "tmag": float(scene.iloc[j]["tmag"])} for j in hyp_rows}
    # ring test, flux budget
    out["rings"] = pv.ring_test(sectors, scene, hyp_rows)
    out["flux_budget"] = pv.flux_budget(sectors, scene, depth)[:12]
    # fixed-aperture depths
    psd = pv.per_sector_depths(sectors)
    out["fixed_aperture_sectors"] = psd
    out["fixed_aperture_const"] = lcvet.chi2_const(psd)
    by = {}
    for key in ("camera", "year"):
        groups = {}
        for r in psd:
            g = int(r[key]) if key == "camera" else int(np.floor(r["year"]))
            groups.setdefault(g, []).append(r)
        by[key] = {str(g): lcvet.chi2_const(v) for g, v in sorted(groups.items())}
    cad_groups = {}
    for r in psd:
        cad_groups.setdefault(f"{r['cadence_min']:.0f}min", []).append(r)
    by["cadence"] = {g: lcvet.chi2_const(v) for g, v in cad_groups.items()}
    out["fixed_aperture_by"] = by
    from scipy import stats
    yrs = np.array([r["year"] for r in psd]); dd = np.array([r["depth_ppm"] for r in psd])
    ee = np.array([r["err_ppm"] for r in psd])
    if len(yrs) > 3:
        w = 1 / ee ** 2
        A = np.vstack([np.ones_like(yrs), yrs - yrs.mean()]).T
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        coef = cov @ (A.T @ (w * dd))
        out["fixed_aperture_trend_ppm_per_yr"] = [float(coef[1]), float(np.sqrt(cov[1, 1]))]

    # fixed-aperture stitched light curve -> box stats + transit fit
    mask_fn = lambda t: ms.in_tr(t, P, t0, dur, 1.6)
    lc = ms.stitch(lcs, transit_mask=mask_fn)
    out["fixed_lc"] = lcvet.box_stats(lc, P, t0, dur)
    rho_tic = 1.41 * float(spec["mstar"]) / float(spec["rstar"]) ** 3
    try:
        out["transit_fit"] = lcvet.fit_transit(lc, P, t0, dur, depth, float(spec["teff"]), rho_tic,
                                               profile=spec["role"] == "candidate")
    except Exception as exc:
        out["transit_fit"] = {"error": repr(exc)[:200]}
    # pipeline light curves (TESS-SPOC / QLP), if cached by mdwarf_hunt fetch
    from planet_hunt import ffi
    pipe = ffi.load_star(int(spec["tic"]))
    if pipe:
        res = {}
        for prov in sorted({d["meta"].get("provenance", "") for d in pipe}):
            sub = [d for d in pipe if d["meta"].get("provenance") == prov]
            plc = ms.stitch(sub, transit_mask=mask_fn)
            if plc:
                bs = lcvet.box_stats(plc, P, t0, dur)
                res[prov] = {k: bs[k] for k in ("depth_ppm", "err_ppm", "snr", "per_sector", "oe_sigma",
                                                "sec_sigma")}
        out["pipeline_lc"] = res
    # the dip prefers a neighbour: look at that star's own light curve
    top = L["hypotheses"][0]
    if not top["is_target"] and L["target_dchi2"] > 25:
        try:
            out["source_followup"] = source_followup(spec, scene, sectors, cuts, top, P, t0, dur)
        except Exception as exc:
            out["source_followup"] = {"error": repr(exc)[:200]}
    out["elapsed_s"] = time.time() - t_start
    plot_target(name, spec, scene, sectors, L, out)
    return out


def source_followup(spec, scene, sectors, cuts, hyp, P, t0, dur) -> dict:
    """Fixed-aperture light curve centred on the preferred neighbour, its TIC
    stellar parameters, and the transit fit with ITS density."""
    from astroquery.mast import Catalogs
    j = int(np.flatnonzero(scene["source_id"].to_numpy() == hyp["source_id"])[0])
    tic = Catalogs.query_criteria(catalog="TIC", GAIA=str(hyp["source_id"])).to_pandas()
    star = tic.iloc[0] if len(tic) else None
    info = {"source_id": hyp["source_id"], "sep_arcsec": hyp["sep_arcsec"], "tmag": hyp["tmag"]}
    if star is not None:
        info.update({"tic": int(star["ID"]), "teff": float(star["Teff"]), "rstar": float(star["rad"]),
                     "mstar": float(star["mass"]), "dist_pc": float(star["d"])})
    by_sector = {s: p for s, p in cuts}
    lcs = []
    for sd in sectors:
        cut = pv.load_cutout(by_sector[sd["sector"]])
        lcs.append(lcvet.cutout_lightcurve(cut, sd, npix=2, row=j))
        del cut
    lc = ms.stitch(lcs, transit_mask=lambda t: ms.in_tr(t, P, t0, dur, 1.6))
    info["fixed_lc"] = lcvet.box_stats(lc, P, t0, dur)
    if star is not None and np.isfinite(info.get("rstar", np.nan)):
        rho = 1.41 * info["mstar"] / info["rstar"] ** 3
        tf = lcvet.fit_transit(lc, P, t0, dur, max(hyp["implied_depth"], 1e-3), info["teff"], rho, profile=True)
        tf.pop("binned", None)
        info["transit_fit"] = tf
        k = tf["k"][1]
        info["radius_earth"] = float(k * info["rstar"] * 109.1)
    return info


# --------------------------------------------------------------------------
def plot_target(name, spec, scene, sectors, L, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    fig, ax = plt.subplots(2, 3, figsize=(16, 10))
    best = max(sectors, key=lambda s: float(np.nansum(np.clip(s["diff"] / s["err"], 0, None)[
        np.abs(np.arange(11)[:, None] - s["tpix"][1]) + np.abs(np.arange(11)[None, :] - s["tpix"][0]) < 3])))
    for a, img, title in ((ax[0, 0], best["oot"], f"S{best['sector']} out-of-transit (sky removed)"),
                          (ax[0, 1], best["diff"] / best["err"], f"S{best['sector']} difference image / noise")):
        im = a.imshow(img if "noise" in title else np.log10(np.clip(img - best["fit"]["c0"], 1, None)),
                      origin="lower", cmap="viridis" if "noise" not in title else "RdBu_r",
                      vmin=None if "noise" not in title else -5, vmax=None if "noise" not in title else 5)
        fig.colorbar(im, ax=a, fraction=0.046)
        sel = scene["rel_flux"] > 10 ** (-0.4 * 6)
        px = best["px"][sel.to_numpy()] + best["fit"]["dx"]
        py = best["py"][sel.to_numpy()] + best["fit"]["dy"]
        a.scatter(px, py, s=8 + 60 * scene["rel_flux"][sel] ** 0.5, facecolors="none", edgecolors="w")
        a.plot(best["tpix"][0], best["tpix"][1], "r+", ms=14, mew=2)
        bx = best["tpix"] + best["J"] @ np.array([L["xi"], L["eta"]])
        a.plot(bx[0], bx[1], "mx", ms=12, mew=2)
        a.set_xlim(-0.5, 10.5); a.set_ylim(-0.5, 10.5); a.set_title(title, fontsize=10)
    # localisation map (sky offsets)
    a = ax[0, 2]
    g = np.array(L["grid"]["xi"]); chi = np.array(L["grid"]["chi"])
    a.contourf(g, g, np.clip(chi, 0, 60), levels=[0, 2.3, 6.2, 11.8, 25, 60], cmap="Greys_r")
    s2 = scene[scene["sep"] < 80]
    a.scatter(s2["xi"], s2["eta"], s=10 + 80 * s2["rel_flux"] ** 0.5, c="tab:orange", edgecolors="k", zorder=3)
    for r in s2.itertuples():
        if r.rel_flux > 0.01 or r.is_target:
            a.annotate(f"T{r.tmag:.1f}", (r.xi, r.eta), fontsize=7, xytext=(3, 3), textcoords="offset points")
    a.plot(0, 0, "r+", ms=14, mew=2)
    a.plot(L["xi"], L["eta"], "mx", ms=12, mew=2)
    cov = np.array(L["cov_jackknife"]) if np.all(np.isfinite(L["cov_jackknife"])) else np.array(L["cov"])
    if np.all(np.isfinite(cov)):
        w_, v_ = np.linalg.eigh(cov)
        ang = np.degrees(np.arctan2(v_[1, 1], v_[0, 1]))
        a.add_patch(Ellipse((L["xi"], L["eta"]), 2 * 2 * np.sqrt(max(w_[1], 0)), 2 * 2 * np.sqrt(max(w_[0], 0)),
                            angle=ang, fill=False, color="m"))
    if len(L["jackknife"]):
        jk = np.array(L["jackknife"]); a.plot(jk[:, 0], jk[:, 1], "m.", ms=3, alpha=0.5)
    a.set_xlim(60, -60); a.set_ylim(-60, 60); a.set_aspect("equal")
    a.set_xlabel("arcsec East"); a.set_ylabel("arcsec North")
    a.set_title(f"source localisation, {L['n_sectors']} sectors: offset {L['offset_arcsec']:.1f}\"; "
                f"target dchi2 {L['target_dchi2']:.0f}", fontsize=9)
    # rings
    a = ax[1, 0]
    rc = np.array(out["rings"]["ring_corrected_depth"])
    a.errorbar(np.arange(len(rc)) + 1, rc[:, 0] * 1e6, rc[:, 1] * 1e6, fmt="ko", label="corrected as if on target")
    a.axhline(spec["depth_ppm"], color="grey", ls=":", label="catalogue depth")
    a.set_xticks(np.arange(len(rc)) + 1, ["1", "2-4", "5-9", "10-16", "17-25"])
    a.set_xlabel("pixel ring (ranked by target PSF)"); a.set_ylabel("depth (ppm)")
    a.set_title("depth vs aperture ring (target-normalised)", fontsize=10); a.legend(fontsize=8)
    # per-sector depths
    a = ax[1, 1]
    psd = out["fixed_aperture_sectors"]
    for cam in sorted({r["camera"] for r in psd}):
        rr = [r for r in psd if r["camera"] == cam]
        a.errorbar([r["year"] for r in rr], [r["depth_ppm"] for r in rr], [r["err_ppm"] for r in rr],
                   fmt="o", ms=4, label=f"TESScut fixed, cam {cam}")
    for prov, mk in (("TESS-SPOC", "s"), ("QLP", "^")):
        pr = out.get("pipeline_lc", {}).get(prov)
        if pr:
            yrs = [2000 + (r["t_mid"] + 2457000 - 2451545) / 365.25 for r in pr["per_sector"]]
            a.errorbar(np.array(yrs) + 0.03, [r["depth_ppm"] for r in pr["per_sector"]],
                       [r["err_ppm"] for r in pr["per_sector"]], fmt=mk, ms=3, alpha=0.5, label=prov)
    a.axhline(spec["depth_ppm"], color="grey", ls=":")
    med_e = np.median([r["err_ppm"] for r in psd])
    a.set_ylim(np.percentile([r["depth_ppm"] for r in psd], 5) - 3 * med_e,
               np.percentile([r["depth_ppm"] for r in psd], 95) + 3 * med_e)
    a.set_xlabel("year"); a.set_ylabel("depth (ppm)"); a.legend(fontsize=7)
    a.set_title("per-sector depth: fixed TESScut aperture vs pipeline products", fontsize=10)
    # folded light curve + model
    a = ax[1, 2]
    tf = out.get("transit_fit", {})
    if "binned" in tf:
        x = np.array(tf["binned"]["x"]) * 24; y = np.array(tf["binned"]["y"]); e = np.array(tf["binned"]["e"])
        o = np.argsort(x)
        nb = max(1, len(x) // 80)
        xb = [x[o][i:i + nb].mean() for i in range(0, len(x), nb)]
        yb = [np.average(y[o][i:i + nb], weights=1 / e[o][i:i + nb] ** 2) for i in range(0, len(x), nb)]
        a.plot(x, (y - 1) * 1e6, ".", color="0.8", ms=2)
        a.plot(xb, (np.array(yb) - 1) * 1e6, "ko", ms=3)
        bf = tf["best"]
        xx = np.linspace(x.min(), x.max(), 400) / 24
        u1, u2 = tf["u"]
        m = lcvet.transit_model(xx - bf["dt0"], np.full(len(xx), 2 / 1440), spec["period"], bf["k"], bf["b"],
                                bf["rho"], u1, u2)
        a.plot(xx * 24, (m - 1) * 1e6, "r-", label=f"b={tf['b'][1]:.2f}, rho={tf['rho_cgs'][1]:.1f} "
                                                   f"(TIC {tf['rho_tic_cgs']:.1f}) g/cc")
        a.legend(fontsize=8)
    fl = out["fixed_lc"]
    a.set_title(f"fixed-aperture fold: {fl['depth_ppm']:.0f}±{fl['err_ppm']:.0f} ppm, odd/even "
                f"{fl['oe_sigma']:.1f}σ, sec {fl['sec_sigma']:.1f}σ", fontsize=9)
    a.set_xlabel("hours from mid-transit"); a.set_ylabel("ppm")
    fig.suptitle(f"{name} (TIC {spec['tic']}) P={spec['period']} d, {spec['role']}", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=110)
    plt.close(fig)


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("targets")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--assemble", help="verdicts JSON: merge per-target results into runs/round4_vet.json")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    specs = json.load(open(args.targets))["targets"]
    if args.assemble:
        verdicts = json.load(open(args.assemble))
        allres = {"generated": time.strftime("%Y-%m-%d"), "notes": verdicts.get("notes", ""), "targets": []}
        for s in specs:
            p = OUT / f"{s['name']}.json"
            if not p.exists():
                continue
            r = json.load(open(p))
            r.pop("transit_fit", {}).get("binned") if False else None
            if "transit_fit" in r:
                r["transit_fit"].pop("binned", None)
            r["localization"].pop("jackknife", None)
            r["verdict"] = verdicts.get("verdicts", {}).get(s["name"])
            r["plot"] = str((OUT / f"{s['name']}.png").relative_to(config.ROOT))
            allres["targets"].append(r)
        (config.RUNS / "round4_vet.json").write_text(json.dumps(_jsonable(allres), indent=1))
        print(f"wrote {config.RUNS / 'round4_vet.json'}")
        return
    for s in specs:
        if args.only and s["name"] not in args.only:
            continue
        r = run_target(s)
        (OUT / f"{s['name']}.json").write_text(json.dumps(_jsonable(r), indent=1))
        L = r["localization"]
        print(f"[{s['name']}] offset {L['offset_arcsec']:.1f}\" ± {L['offset_err_arcsec']:.1f}; target dchi2 "
              f"{L['target_dchi2']:.1f}; best {L['hypotheses'][0]['sep_arcsec']:.1f}\" T{L['hypotheses'][0]['tmag']:.1f} "
              f"(depth {L['hypotheses'][0]['implied_depth']:.3f}); per-sector target best "
              f"{r['per_sector_summary']['target_best']}/{r['per_sector_summary']['n']}; fixed-ap depth "
              f"{r['fixed_aperture_const'].get('mean', np.nan):.0f} (chi2 p {r['fixed_aperture_const'].get('p', np.nan):.2g}); "
              f"rho {r['transit_fit'].get('rho_cgs', ['?'] * 3)[1]} vs TIC {r['transit_fit'].get('rho_tic_cgs', '?')}; "
              f"{r['elapsed_s']:.0f} s", flush=True)


if __name__ == "__main__":
    main()
