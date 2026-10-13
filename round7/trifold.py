"""Round 7: fold + FFI apertures for TRICERATOPS on an FFI-only star.

    .venv/bin/python round7/trifold.py <tic>
Writes runs/round7_eumemic_weak/tri/<tic>_fold.npz and <tic>_apertures.npz and prints
the run_triceratops.py command. Apertures: a 3x3 pixel box centred on the
target's pixel in each TESScut cutout (the TESS-SPOC / QLP apertures of a T~12
star are 2x2 to 3x3), so neighbours are, if anything, over-included.
"""
from __future__ import annotations

import glob
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from planet_hunt import config  # noqa: E402
from planet_hunt import hostsearch as hs  # noqa: E402
from planet_hunt import multisector as ms  # noqa: E402
import revet as rv  # noqa: E402

OUT = config.RUNS / "round7_eumemic_weak" / "tri"


def main(tic: int):
    from astropy.io import fits
    from astropy.wcs import WCS
    import round5_hosts as r5
    OUT.mkdir(parents=True, exist_ok=True)
    tg = pd.read_csv(rv.OUT / "targets.csv").drop_duplicates("tic").set_index("tic")
    res = pd.read_csv(rv.OUT / "revet.csv").set_index("tic").loc[tic]
    spec = [s for s in json.load(open(Path(__file__).parent / "deepvet_targets.json"))["targets"] if s["tic"] == tic][0]
    row = tg.loc[tic]
    plans = {int(k): v for k, v in json.load(open(rv.OUT / "plans.json")).items()}
    P, t0, dur = float(res.P_full), float(res.t0_full), float(res.dur_h) / 24
    secs = rv.load_sectors(tic, plans[tic], int(row.his_cut))
    toi, ps = r5.load_catalogs()
    kl = hs.known_list(r5.known_ephemerides(tic, toi, ps), row.rstar, row.mstar)
    km = lambda t: hs.mask_known(t, kl)
    rot = hs.rotation(secs, km)
    lc = hs.detrend(secs, hs.W_SHORT if P < hs.SPLIT_P else hs.W_LONG, lambda t: km(t) | ms.in_tr(t, P, t0, dur, 1.6), rot)
    ph = ms.phase(lc["t"], P, t0)
    near = np.abs(ph) < max(3 * dur, 0.25 * P if P < 1 else 0.25)
    np.savez(OUT / f"{tic}_fold.npz", phase_d=ph[near], flux=lc["flat"][near], err=lc["dy"][near], period=P,
             depth=float(res.depth_full_ppm) * 1e-6, sectors=np.array(sorted(set(lc["sector"].astype(int)))))
    cdir = config.CACHE / "tesscut" / f"{spec['ra']:.5f}_{spec['dec']:+.5f}".replace("_+", "_")
    if not cdir.exists():
        cdir = Path(glob.glob(str(config.CACHE / "tesscut" / f"{spec['ra']:.5f}_*"))[0])
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
    np.savez(OUT / f"{tic}_apertures.npz", **aps)
    print(f"apertures for sectors {sorted(int(k[1:]) for k in aps)}; fold {near.sum()} points")
    print(f".trivenv/bin/python run_triceratops.py --tic {tic} --fold {OUT / f'{tic}_fold.npz'} "
          f"--lc-dir {OUT / f'{tic}_apertures.npz'} --runs 5 --out {OUT} --name triceratops_{tic} --local-cutouts {cdir}")


if __name__ == "__main__":
    main(int(sys.argv[1]))
