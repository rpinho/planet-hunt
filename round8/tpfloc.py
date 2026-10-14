"""Round 8A pixel localisation from SPOC 2-min target pixel files (generalised round-6 method;
difference imaging in the sense of Bryson et al. 2013, with an empirical PRF).

For the candidate AND the host's known planet (the on-target control, same pixels):
  * per sector: mean in-transit difference image (local linear per-pixel baseline per event),
    out-of-transit direct image, and a per-pixel noise image from fake epochs (same number of
    events averaged, placed away from every known/candidate transit);
  * template-shift localisation: each difference image is fitted with the sector's own direct
    image (the star's empirical PRF) shifted to every sky offset on a 2" grid (+-90"); chi^2
    surfaces are rescaled by their reduced chi^2 and summed over sectors;
  * depth vs aperture (core, r<=1.5, 2.5, 3.5 px, SPOC aperture) relative to the control.
Neighbours: Gaia DR3 within 110"; each gets the summed delta-chi2 at its position and the
eclipse depth it would need (all of its flux in the aperture).

    .r6venv/bin/python round8/tpfloc.py toi6008 [--grid 90]
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from planet_hunt import config  # noqa: E402
from planet_hunt import pixelvet as pv  # noqa: E402
from targets import CANDS  # noqa: E402

OUT = config.RUNS / "round8_p1" / "tpf"
RAWS = [config.CACHE / "r8" / "raw", config.CACHE / "r6" / "raw"]


def tpf_files(tic):
    fs = []
    for r in RAWS:
        fs += glob.glob(str(r / str(tic) / f"tess*-{tic:016d}-*-s_tp.fits")) + glob.glob(str(r / f"tess*-{tic:016d}-*-s_tp.fits"))
    by = {}
    for f in sorted(fs):
        s = int(Path(f).name.split("-s")[1][:4])
        by.setdefault(s, f)
    return by


def load_tpf(fn):
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(fn, memmap=False) as h:
        d = h[1].data
        t = np.array(d["TIME"], float)
        F = np.array(d["FLUX"], float)
        q = np.array(d["QUALITY"])
        ap = np.array(h[2].data)
        w = WCS(h[2].header)
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1)
    return dict(t=t[ok], F=F[ok], ap=ap, wcs=w)


def in_tr(t, P, t0, half):
    return np.abs(((t - t0) / P + 0.5) % 1 - 0.5) * P < half


def diff_one(t, F, tc, dur, excl, nmin=6):
    w = np.abs(t - tc) < 0.75 * dur + max(0.25, 1.5 * dur)
    tin = (np.abs(t - tc) < 0.35 * dur) & ~excl
    oot = w & (np.abs(t - tc) > 0.75 * dur) & ~excl
    if tin.sum() < nmin or oot.sum() < 4 * nmin or (t[oot] < tc).sum() < nmin or (t[oot] > tc).sum() < nmin:
        return None
    X = np.vstack([np.ones(oot.sum()), t[oot] - tc]).T
    Y = F[oot].reshape(oot.sum(), -1)
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    Xi = np.vstack([np.ones(tin.sum()), t[tin] - tc]).T
    diff = (Xi @ beta - F[tin].reshape(tin.sum(), -1)).mean(0).reshape(F.shape[1:])
    return diff, beta[0].reshape(F.shape[1:])


def jac(wcs, x, y):
    J = np.zeros((2, 2)); h = 0.05
    a = wcs.all_pix2world([[x, y]], 0)[0]
    for k in range(2):
        p = [x, y]; p[k] += h
        b = wcs.all_pix2world([p], 0)[0]
        J[0, k] = (b[0] - a[0]) * np.cos(np.radians(a[1])) * 3600 / h
        J[1, k] = (b[1] - a[1]) * 3600 / h
    return J


def aperture_depths(t, F, ap, centre, epochs, dur, excl):
    ny, nx = F.shape[1:]
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = np.hypot(xx - centre[0], yy - centre[1])
    aps = {"core": r <= 0.75, "r1.5": r <= 1.5, "r2.5": r <= 2.5, "r3.5": r <= 3.5, "SPOC": (ap & 2) > 0}
    out = {}
    for name, a in aps.items():
        if a.sum() == 0:
            continue
        fl = F[:, a].sum(1); fl = fl / np.median(fl)
        D, W = [], []
        for tc in epochs:
            w = np.abs(t - tc) < 0.75 * dur + max(0.25, 1.5 * dur)
            tin = (np.abs(t - tc) < 0.35 * dur) & ~excl
            oot = w & (np.abs(t - tc) > 0.75 * dur) & ~excl
            if tin.sum() < 6 or oot.sum() < 24:
                continue
            X = np.vstack([np.ones(oot.sum()), t[oot] - tc]).T
            b, *_ = np.linalg.lstsq(X, fl[oot], rcond=None)
            res = fl[oot] - X @ b
            D.append(b[0] - fl[tin].mean()); W.append(1 / (np.std(res) ** 2 * (1 / tin.sum() + 1 / oot.sum())))
        if D:
            D, W = np.array(D), np.array(W)
            out[name] = dict(depth=float(np.sum(D * W) / W.sum()), err=float(1 / np.sqrt(W.sum())), npix=int(a.sum()))
    return out


def main():
    from scipy.ndimage import shift as ndshift
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("key")
    ap_.add_argument("--grid", type=float, default=90.0)
    ap_.add_argument("--step", type=float, default=2.0)
    ap_.add_argument("--nfake", type=int, default=40)
    a = ap_.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    c = CANDS[a.key]
    lg = json.load(open(config.RUNS / "round8_p1" / f"{a.key}_lcgate.json"))
    inv = json.load(open(config.RUNS / "round8_p1" / "inventory.json"))[a.key]
    stars = json.load(open(config.RUNS / "round8_p1" / "stars.json"))[a.key]
    known = lg["known_masked"]
    sigs = {"cand": dict(P=c["P"], t0=c["t0"], dur=c["dur_h"] / 24, depth=lg["all"]["depth_ppm"] * 1e-6)}
    if known:
        k = known[0]
        sigs["ctrl"] = dict(P=k["P"], t0=k["t0"], dur=k["dur_h"] / 24, depth=None, name=k["name"])
    rng = np.random.default_rng(c["tic"] % 2 ** 32)
    grid = np.arange(-a.grid, a.grid + 0.1, a.step)
    XI, ETA = np.meshgrid(grid, grid)
    files = tpf_files(c["tic"])
    print(f"{a.key}: {len(files)} TPF sectors {sorted(files)}", flush=True)
    tot = {s: np.zeros_like(XI) for s in sigs}
    per = {s: [] for s in sigs}
    apd = {s: {} for s in sigs}
    for sec, fn in sorted(files.items()):
        tp = load_tpf(fn)
        t, F = tp["t"], tp["F"]
        allmask = np.zeros(len(t), bool)
        for kk in known:
            allmask |= in_tr(t, kk["P"], kk["t0"], 0.75 * kk["dur_h"] / 24 + 0.02)
        allmask |= in_tr(t, c["P"], c["t0"], 0.75 * c["dur_h"] / 24 + 0.02)
        for name, sg in sigs.items():
            P, t0, dur = sg["P"], sg["t0"], sg["dur"]
            excl = np.zeros(len(t), bool)                 # other signals' transits
            for kk in known:
                if abs(kk["P"] / P - 1) > 1e-4:
                    excl |= in_tr(t, kk["P"], kk["t0"], 0.75 * kk["dur_h"] / 24 + 0.02)
            if name == "ctrl":
                excl |= in_tr(t, c["P"], c["t0"], 0.75 * c["dur_h"] / 24 + 0.02)
            n0, n1 = int(np.ceil((t.min() - t0) / P)), int(np.floor((t.max() - t0) / P))
            diffs, directs, eps = [], [], []
            for n in range(n0, n1 + 1):
                r = diff_one(t, F, t0 + n * P, dur, excl)
                if r is not None:
                    diffs.append(r[0]); directs.append(r[1]); eps.append(t0 + n * P)
            if len(diffs) < 1:
                continue
            D = np.mean(diffs, 0); Dir = np.mean(directs, 0)
            fakes, tries = [], 0
            while len(fakes) < a.nfake and tries < 50 * a.nfake:
                tries += 1
                fk = []
                for _ in range(len(diffs)):
                    for _t in range(30):
                        tc = rng.uniform(t.min() + 0.3, t.max() - 0.3)
                        if not allmask[np.abs(t - tc) < 0.75 * dur + 0.1].any():
                            r = diff_one(t, F, tc, dur, np.zeros(len(t), bool))
                            if r is not None:
                                fk.append(r[0]); break
                if len(fk) == len(diffs):
                    fakes.append(np.mean(fk, 0))
            if len(fakes) < 10:
                continue
            fakes = np.array(fakes)
            errD = np.std(fakes, 0)
            errD = np.where(errD > 0, errD, np.nanmedian(errD))
            ny, nx = Dir.shape
            cy, cx = np.unravel_index(np.argmax(Dir), Dir.shape)
            tx, ty = tp["wcs"].all_world2pix([[inv["ra"], inv["dec"]]], 0)[0]
            Jinv = np.linalg.inv(jac(tp["wcs"], cx, cy))
            chi = np.zeros_like(XI)
            e = errD.ravel()
            y = D.ravel() / e
            for i in range(XI.shape[0]):
                for j in range(XI.shape[1]):
                    # template = direct image shifted so that the star moves by (xi, eta) from the TARGET position
                    dx, dy = Jinv @ np.array([XI[i, j], ETA[i, j]])
                    T = ndshift(Dir, (dy, dx), order=3, mode="constant", cval=0.0).ravel()
                    X = np.vstack([T, np.ones_like(T)]).T / e[:, None]
                    b, *_ = np.linalg.lstsq(X, y, rcond=None)
                    r = y - X @ b
                    chi[i, j] = r @ r
            i0 = np.unravel_index(np.argmin(chi), chi.shape)
            dof = D.size - 4
            scale = max(1.0, chi[i0] / dof)
            cs = chi / scale - (chi / scale).min()
            c0 = float(cs[XI.shape[0] // 2, XI.shape[1] // 2])
            amp_snr = float(np.sum(D * Dir / errD ** 2) / np.sqrt(np.sum(Dir ** 2 / errD ** 2)))
            per[name].append(dict(sector=sec, n_events=len(diffs), best_xi=float(XI[i0]), best_eta=float(ETA[i0]),
                                  dchi2_target=c0, rchi2=float(chi[i0] / dof), diff_snr=amp_snr))
            tot[name] += cs
            apd[name][sec] = aperture_depths(t, F, tp["ap"], (tx, ty), eps, dur, excl)
            print(f"  s{sec} {name}: n={len(diffs)} best ({XI[i0]:+.0f},{ETA[i0]:+.0f})\" dchi2_target {c0:.1f} "
                  f"rchi2 {chi[i0] / dof:.2f} diffSNR {amp_snr:.1f}", flush=True)
    scene = pv.gaia_scene(inv["ra"], inv["dec"], radius_arcsec=200, gaia_id=stars.get("gaia_dr3"), tmag=c["tmag"])
    res = {"key": a.key, "sectors": sorted(files), "grid_arcsec": a.grid, "step": a.step}
    for name in sigs:
        T = tot[name]
        if not per[name]:
            continue
        i0 = np.unravel_index(np.argmin(T), T.shape)
        d0 = float(T[XI.shape[0] // 2, XI.shape[1] // 2] - T.min())
        in3 = T - T.min() <= 11.8
        rec = {"best_xi": float(XI[i0]), "best_eta": float(ETA[i0]), "dchi2_target": d0,
               "max_sep_3sig": float(np.max(np.hypot(XI[in3], ETA[in3]))),
               "radius_1sig": float(np.max(np.hypot(XI[T - T.min() <= 2.3], ETA[T - T.min() <= 2.3]))),
               "per_sector": per[name], "n_sectors": len(per[name])}
        nbs = []
        depth = sigs["cand"]["depth"]
        for r in scene.itertuples():
            if r.is_target or r.sep > min(a.grid, 110):
                continue
            ii = int(np.argmin(np.abs(grid - r.eta))); jj = int(np.argmin(np.abs(grid - r.xi)))
            need = depth / r.rel_flux * (1 + 0) if r.rel_flux > 0 else np.inf
            nbs.append({"source_id": int(r.source_id), "sep": float(r.sep), "pa": float(r.pa_deg), "tmag": float(r.tmag),
                        "need_depth_cand": float(min(need, 99)), "can_mimic_cand": bool(need < 1),
                        "dchi2": float(T[ii, jj] - T.min())})
        rec["neighbours"] = sorted(nbs, key=lambda x: x["sep"])
        mim = [x for x in nbs if x["can_mimic_cand"]]
        rec["min_dchi2_mimic_capable"] = float(min([x["dchi2"] for x in mim], default=np.nan))
        rec["closest_mimic_capable"] = min(mim, key=lambda x: x["dchi2"]) if mim else None
        # depth vs aperture, combined over sectors
        comb = {}
        for apn in ("core", "r1.5", "r2.5", "r3.5", "SPOC"):
            d_ = [v[apn]["depth"] for v in apd[name].values() if apn in v]
            e_ = [v[apn]["err"] for v in apd[name].values() if apn in v]
            if d_:
                w = 1 / np.array(e_) ** 2
                comb[apn] = dict(depth=float(np.sum(np.array(d_) * w) / w.sum()), err=float(1 / np.sqrt(w.sum())))
        rec["aperture_depths"] = comb
        res[name] = rec
        np.save(OUT / f"{a.key}_{name}_locmap.npy", T)
        print(f"{name}: best ({rec['best_xi']:+.0f},{rec['best_eta']:+.0f})\" dchi2(target) {d0:.1f}; 3-sigma max sep "
              f"{rec['max_sep_3sig']:.0f}\"; min dchi2 of a mimic-capable neighbour {rec['min_dchi2_mimic_capable']:.1f}", flush=True)
    if "ctrl" in res and "cand" in res:
        # calibrated test: is the candidate consistent with the position where the known (on-target)
        # planet's own signal localises on the same pixels? (removes the PRF/WCS systematic offset)
        Tc = tot["cand"]
        ii = int(np.argmin(np.abs(grid - res["ctrl"]["best_eta"]))); jj = int(np.argmin(np.abs(grid - res["ctrl"]["best_xi"])))
        res["cand"]["dchi2_at_ctrl_position"] = float(Tc[ii, jj] - Tc.min())
        off = (res["ctrl"]["best_xi"], res["ctrl"]["best_eta"])
        for nb in res["cand"]["neighbours"]:
            x = nb["sep"] * np.sin(np.radians(nb["pa"])) + off[0]; y = nb["sep"] * np.cos(np.radians(nb["pa"])) + off[1]
            i2 = int(np.argmin(np.abs(grid - y))); j2 = int(np.argmin(np.abs(grid - x)))
            nb["dchi2_ctrl_frame"] = float(Tc[i2, j2] - Tc.min())
        mim = [x for x in res["cand"]["neighbours"] if x["can_mimic_cand"]]
        res["cand"]["min_dchi2_mimic_capable_ctrl_frame"] = float(min([x["dchi2_ctrl_frame"] for x in mim], default=np.nan))
        print(f"candidate dchi2 at the control's position {res['cand']['dchi2_at_ctrl_position']:.1f}; "
              f"min mimic-capable neighbour (control frame) {res['cand']['min_dchi2_mimic_capable_ctrl_frame']:.1f}")
        rr = {}
        for apn in ("core", "r1.5", "r2.5", "r3.5"):
            ca, cs_ = res["cand"]["aperture_depths"].get(apn), res["cand"]["aperture_depths"].get("SPOC")
            ka, ks = res["ctrl"]["aperture_depths"].get(apn), res["ctrl"]["aperture_depths"].get("SPOC")
            if ca and cs_ and ka and ks:
                r1 = ca["depth"] / cs_["depth"]; r2 = ka["depth"] / ks["depth"]
                e1 = abs(r1) * np.hypot(ca["err"] / ca["depth"], cs_["err"] / cs_["depth"])
                e2 = abs(r2) * np.hypot(ka["err"] / ka["depth"], ks["err"] / ks["depth"])
                rr[apn] = dict(cand_over_ctrl=r1 / r2, err=abs(r1 / r2) * np.hypot(e1 / r1, e2 / r2))
        res["aperture_ratio_vs_ctrl"] = rr
        print("aperture ratio (cand/SPOC) / (ctrl/SPOC):", {k: f"{v['cand_over_ctrl']:.2f}+-{v['err']:.2f}" for k, v in rr.items()})
    (OUT / f"{a.key}_tpfloc.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
