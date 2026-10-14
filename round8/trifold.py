"""Round 8: fold + 3x3 FFI apertures for TRICERATOPS (round-7 trifold.py, keyed by signal id).

    .venv/bin/python round8/trifold.py <signal id>
Writes runs/round8_crowd/tri/<id>_fold.npz and <id>_apertures.npz and prints the run_triceratops.py command.
"""
from __future__ import annotations

import glob
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from common import OUT, config
import holdout as H
from holdout import hs, ms, r5

T = OUT / "tri"


def main(sid: str):
    from astropy.io import fits
    from astropy.wcs import WCS
    T.mkdir(parents=True, exist_ok=True)
    dj = json.load(open(OUT / "diag" / f"{sid}.json"))
    spec = [s for s in json.load(open(OUT / "deepvet_targets.json"))["targets"] if s["name"] == sid][0]
    tg = pd.read_csv(OUT / "targets.csv")
    tg["id"] = [H.sig_id(r) for r in tg.to_dict("records")]
    row = tg.set_index("id").loc[sid].to_dict()
    tic = int(dj["tic"])
    plans = H.load_plans(tg)
    sig = {"rule": row["rule"], "src_sectors": json.loads(row["src_sectors"]) if str(row["src_sectors"]).startswith("[") else []}
    secs = H.load_sectors(tic, plans[tic], sig)
    eph = dj.get("ephemeris") or {}
    P, t0, dur = float(eph.get("P", dj["P"])), float(dj["t0"]), float(dj["dur_h"]) / 24
    kl = hs.known_list(r5.known_ephemerides(tic, *r5.load_catalogs()), row["rstar"], row["mstar"])
    km = lambda t: hs.mask_known(t, kl)
    rot = hs.rotation(secs, km)
    lc = hs.detrend(secs, hs.W_SHORT if P < hs.SPLIT_P else hs.W_LONG, lambda t: km(t) | ms.in_tr(t, P, t0, dur, 1.6), rot)
    ph = ms.phase(lc["t"], P, t0)
    near = np.abs(ph) < max(3 * dur, 0.25 * P if P < 1 else 0.25)
    depth = float(dj["stats_masked_detrend"]["depth_ppm"]) * 1e-6
    np.savez(T / f"{sid}_fold.npz", phase_d=ph[near], flux=lc["flat"][near], err=lc["dy"][near], period=P,
             depth=depth, sectors=np.array(sorted(set(lc["sector"].astype(int)))))
    cands = glob.glob(str(config.CACHE / "tesscut" / f"{spec['ra']:.5f}_*"))
    cdir = Path(cands[0])
    aps = {}
    for f in sorted(cdir.glob("*_astrocut.fits")):
        s = int(re.search(r"tess-s(\d{4})-", f.name).group(1))
        with fits.open(f) as hd:
            h1 = hd[1].header
            w = WCS(hd[2].header)
            x, y = w.world_to_pixel_values(spec["ra"], spec["dec"])
            c0, r0 = int(h1["1CRV4P"]), int(h1["2CRV4P"])
        xi, yi = int(round(float(x))), int(round(float(y)))
        aps[f"s{s}"] = np.array([[c0 + xi + dx, r0 + yi + dy] for dx in (-1, 0, 1) for dy in (-1, 0, 1)], int)
    np.savez(T / f"{sid}_apertures.npz", **aps)
    print(f"apertures for sectors {sorted(int(k[1:]) for k in aps)}; fold {near.sum()} points; depth {depth * 1e6:.0f} ppm")
    print(f".trivenv/bin/python run_triceratops.py --tic {tic} --fold {T / f'{sid}_fold.npz'} "
          f"--lc-dir {T / f'{sid}_apertures.npz'} --runs 5 --out {T} --name triceratops_{sid} --local-cutouts {cdir}")


if __name__ == "__main__":
    main(sys.argv[1])
