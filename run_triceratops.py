"""Statistical validation of a transit candidate with TRICERATOPS.

    # TRICERATOPS has heavy, pinned dependencies: give it its own environment
    uv venv --python 3.12 .trivenv && uv pip install --python .trivenv/bin/python triceratops
    .trivenv/bin/python run_triceratops.py --tic 259233660 \
        --fold runs/round4_vet/twomin_259233660_fold.npz --lc-dir cache/mast --runs 3

Giacalone et al. (2021): computes the probability of each astrophysical
scenario (planet on target, eclipsing binary on target, hierarchical EB,
planet or EB on a nearby star, background planet/EB) from the folded light
curve, the TIC/Gaia neighbours inside the photometric apertures and stellar
priors. FPP = 1 - P(planet on target, incl. hierarchical); NFPP = probability
that the signal comes from a nearby resolved star. The published validation
bar is FPP < 0.015 and NFPP < 0.001; FPP < 0.5 is the usual bar for a
candidate. Without a high-resolution contrast curve, unresolved companions
are only constrained by priors, so FPP is an upper-end estimate.

Apertures come from the SPOC 2-minute light-curve files (APERTURE extension,
bit 2 = optimal aperture). The folded light curve (phase in days, flux,
error) is written by twomin_check.py. This script imports nothing from the
planet_hunt package so it can run in TRICERATOPS's own environment.
Writes runs/round4_vet/triceratops_<tic>.json.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np


def spoc_apertures(lc_dir: str, tic: int, sectors) -> tuple[list, list]:
    from astropy.io import fits
    aps, used = [], []
    for s in sectors:
        hit = sorted(glob.glob(f"{lc_dir}/tess*-s{int(s):04d}-{tic:016d}-*-s_lc.fits"))
        if not hit:
            continue
        with fits.open(hit[0]) as hd:
            img = hd["APERTURE"].data
            h = hd["APERTURE"].header
            c0, r0 = int(h["CRVAL1P"]), int(h["CRVAL2P"])
        yy, xx = np.nonzero(img & 2)
        aps.append(np.column_stack([c0 + xx, r0 + yy]))
        used.append(int(s))
    return aps, used


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tic", type=int, required=True)
    ap.add_argument("--fold", required=True)
    ap.add_argument("--lc-dir", required=True)
    ap.add_argument("--runs", type=int, default=3)
    ap.add_argument("--n", type=int, default=1_000_000)
    ap.add_argument("--max-sectors", type=int, default=12, help="apertures from this many sectors")
    ap.add_argument("--bin-min", type=float, default=4.0)
    ap.add_argument("--out", default="runs/round4_vet")
    args = ap.parse_args()
    import triceratops.triceratops as tr
    z = np.load(args.fold)
    ph, fl, er = z["phase_d"], z["flux"], z["err"]
    P, depth = float(z["period"]), float(z["depth"])
    sectors = [int(s) for s in z["sectors"]]
    pick = [sectors[i] for i in np.unique(np.linspace(0, len(sectors) - 1, args.max_sectors).round().astype(int))]
    aps, used = spoc_apertures(args.lc_dir, args.tic, pick)
    # bin the fold (TRICERATOPS evaluates every sample on every point)
    edges = np.arange(ph.min(), ph.max() + 1e-9, args.bin_min / 1440)
    idx = np.digitize(ph, edges)
    w = 1 / er ** 2
    xb, yb, eb = [], [], []
    for k in np.unique(idx):
        m = idx == k
        if m.sum() < 3:
            continue
        xb.append(np.average(ph[m], weights=w[m])); yb.append(np.average(fl[m], weights=w[m]))
        eb.append(1 / np.sqrt(w[m].sum()))
    xb, yb, eb = np.array(xb), np.array(yb), np.array(eb)
    yb = yb / np.median(yb[np.abs(xb) > 0.1])
    target = tr.target(ID=args.tic, sectors=np.array(used), background_population_source="gaia")
    target.calc_depths(tdepth=depth, all_ap_pixels=aps)
    stars = target.stars
    runs = []
    for i in range(args.runs):
        target.calc_probs(time=xb, flux_0=yb, flux_err_0=float(np.median(eb)), P_orb=P, N=args.n,
                          exptime=args.bin_min / 1440, verbose=0)
        probs = target.probs
        runs.append({"FPP": float(target.FPP), "NFPP": float(target.NFPP),
                     "scenarios": {f"{r.ID}:{r.scenario}": float(r.prob) for r in probs.itertuples()
                                   if r.prob > 1e-4}})
        print(f"run {i + 1}: FPP={target.FPP:.4f} NFPP={target.NFPP:.5f}", flush=True)
    fpp = np.array([r["FPP"] for r in runs]); nfpp = np.array([r["NFPP"] for r in runs])
    out = {"tic": args.tic, "period": P, "depth": depth, "aperture_sectors": used, "n_bins": int(len(xb)),
           "FPP_mean": float(fpp.mean()), "FPP_std": float(fpp.std()), "NFPP_mean": float(nfpp.mean()),
           "NFPP_std": float(nfpp.std()), "runs": runs, "n_stars_in_search": int(len(stars)),
           "stars_with_flux_in_aperture": [{"ID": str(r.ID), "Tmag": float(r.Tmag), "sep": float(r.sep),
                                            "tdepth": float(r.tdepth)} for r in stars.itertuples()
                                           if getattr(r, "fluxratio", 0) > 0][:15],
           "note": "no contrast curve: unresolved companions constrained by priors only"}
    Path(args.out).mkdir(parents=True, exist_ok=True)
    (Path(args.out) / f"triceratops_{args.tic}.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps({k: v for k, v in out.items() if k != "runs"}, indent=1, default=str))


if __name__ == "__main__":
    main()
