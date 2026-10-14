"""Round 8A: TRICERATOPS (Giacalone et al. 2021) inputs and MOLUSC (Wood et al. 2021) companion
limits for each P1 candidate.

    .venv/bin/python round8/tri.py prep <key>        # fold + apertures + contrast curves (CSV arcsec,dmag)
    round8/tri.py molusc <key>                        # prints the MOLUSC command (Gaia contrast + RUWE + AO)
    round8/tri.py cmds <key>                          # prints the TRICERATOPS commands for every configuration

Configurations (each 3 runs x 1e6 draws unless stated):
  base     : no contrast curve, no clearing (round-5-like)
  cc       : every ExoFOP contrast curve
  cc_clear : + Gaia neighbours excluded by the pixel localisation set to zero depth
  cc_clear_molusc : + MOLUSC bound-companion population (Gaia contrast + RUWE + the deepest AO/speckle curve)
"""
from __future__ import annotations

import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from planet_hunt import config  # noqa: E402
from planet_hunt import pixelvet as pv  # noqa: E402
from targets import CANDS  # noqa: E402

OUT = config.RUNS / "round8_p1"
TRI = OUT / "tri"
ROOT = config.ROOT


def prep(key):
    TRI.mkdir(parents=True, exist_ok=True)
    c = CANDS[key]
    z = np.load(OUT / f"{key}_fold.npz", allow_pickle=True)
    cad = z["cad"] * 1440
    keep = cad <= (2.5 if (cad <= 2.5).sum() > 1000 else 11)
    ph, fl, er, sec = z["phase_d"][keep], z["flux"][keep], z["err"][keep], z["sector"][keep]
    lg = json.load(open(OUT / f"{key}_lcgate.json"))
    # depth for TRICERATOPS: the limb-darkened fit's k^2-free transit depth is not what it wants; it
    # wants the observed depth -> box depth of the data actually folded
    depth = max(lg["all"]["depth_ppm"] * 1e-6, 1e-5)
    np.savez(TRI / f"{key}_fold.npz", phase_d=ph, flux=fl, err=er, period=c["P"], depth=depth,
             sectors=np.array(sorted(set(sec.astype(int)))))
    print(f"{key}: fold {keep.sum()} points (cadence <= {cad[keep].max():.1f} min), sectors {sorted(set(sec.astype(int)))}, depth {depth*1e6:.0f} ppm")
    inv = json.load(open(OUT / "inventory.json"))[key]
    raw = config.CACHE / "r8" / "raw" / str(c["tic"])
    has2 = bool(glob.glob(str(raw / "tess*-s_lc.fits")))
    if not has2:
        from astropy.io import fits
        from astropy.wcs import WCS
        aps = {}
        for f in sorted(pv.cut_dir(inv["ra"], inv["dec"]).glob("*_11x11_astrocut.fits")):
            s = int(re.search(r"tess-s(\d{4})-", f.name).group(1))
            with fits.open(f) as hd:
                h1 = hd[1].header
                w = WCS(hd[2].header)
                x, y = w.world_to_pixel_values(inv["ra"], inv["dec"])
                c0, r0 = int(h1["1CRV4P"]), int(h1["2CRV4P"])
            xi, yi = int(round(float(x))), int(round(float(y)))
            aps[f"s{s}"] = np.array([[c0 + xi + dx, r0 + yi + dy] for dx in (-1, 0, 1) for dy in (-1, 0, 1)], int)
        np.savez(TRI / f"{key}_apertures.npz", **aps)
        print("  FFI apertures (3x3) for", sorted(int(k[1:]) for k in aps))
    return has2


def contrast_curves(key):
    """ExoFOP contrast curves saved by the novelty agent -> TRICERATOPS CSV (arcsec,dmag) + filter."""
    c = CANDS[key]
    d = OUT / "novelty" / "imaging" / str(c["tic"])
    out = []
    for f in sorted(d.glob("*")) if d.exists() else []:
        if f.suffix.lower() not in (".dat", ".txt", ".csv", ".tbl"):
            continue
        rows = []
        for line in open(f, errors="ignore"):
            parts = re.split(r"[,\s]+", line.strip())
            try:
                v = [float(p) for p in parts if p != ""]
            except ValueError:
                continue
            if len(v) == 4:          # Gemini speckle: r_in, r_out, r_centre, dmag
                rows.append([v[2], v[3]])
            elif len(v) >= 2:        # arcsec, dmag[, rms]
                rows.append(v[:2])
        if len(rows) < 5:
            continue
        a = np.array(rows)
        if a[:, 0].max() > 50:          # mas -> arcsec
            a[:, 0] /= 1000.0
        a = a[(a[:, 0] > 0) & np.isfinite(a[:, 1])]
        n = f.name.lower()
        filt = ("K" if any(s in n for s in ("brg", "brgamma", "_k", "-ks", "ks_", "kcont", "2.15")) else
                "J" if "-j_" in n or "_j" in n or "1.24" in n else
                "Vis" if any(s in n for s in ("562", "692", "_v_", "_r_")) else "TESS")
        dst = TRI / f"{key}_cc_{f.stem}.csv"
        np.savetxt(dst, a, delimiter=",", fmt="%.4f")
        out.append((str(dst), filt, float(np.interp(0.5, a[:, 0], a[:, 1]))))
    return out


if __name__ == "__main__":
    cmd, key = sys.argv[1], sys.argv[2]
    if cmd == "prep":
        prep(key)
        for x in contrast_curves(key):
            print("  contrast curve", x)


MOL = Path(__file__).resolve().parent.parent / ".tools"


def run_molusc(key, n=500000):
    """MOLUSC: Gaia contrast + RUWE + the deepest AO/speckle curve (one curve per CLI run)."""
    import subprocess
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = CANDS[key]
    inv = json.load(open(OUT / "inventory.json"))[key]
    mass = inv["tic_row"].get("mass") or c["mstar"]
    sc = SkyCoord(inv["ra"] * u.deg, inv["dec"] * u.deg)
    ra = sc.ra.to_string(unit=u.hour, sep="hms", precision=3)
    dec = sc.dec.to_string(unit=u.deg, sep="dms", precision=2, alwayssign=True)
    ccs = contrast_curves(key)
    pick = None
    pref = {"K": 3, "TESS": 2, "J": 1}
    for f, filt, d05 in ccs:
        if filt in pref and (pick is None or (pref[filt], d05) > (pref[pick[1]], pick[2])):
            pick = (f, filt, d05)
    prefix = str(TRI / f"molusc_{key}")
    opts = ["--gaia", "--ruwe"]
    if pick:
        a = np.loadtxt(pick[0], delimiter=",")
        mas = TRI / f"molusc_{key}_ao_mas.txt"
        with open(mas, "w") as fh:
            fh.write("Sep Contrast\n")
            for s_, d_ in a:
                fh.write(f"{s_ * 1000:.1f} {d_:.4f}\n")
        opts += ["--ao", str(mas), "--filter", {"K": "K", "TESS": "I", "J": "J"}[pick[1]]]
    # "--" so that a negative declination is not parsed as an option
    cmd = [str(MOL / "molenv" / "bin" / "python"), "BinaryStarGUI.py"] + opts + ["--", prefix, ra, dec, str(n), str(mass)]
    print(" ".join(cmd), flush=True)
    with open(TRI / f"molusc_{key}.log", "w") as log:
        r = subprocess.run(cmd, cwd=MOL / "molusc" / "code", stdout=log, stderr=subprocess.STDOUT)
    kept = Path(prefix + "_kept.csv")
    print(key, "MOLUSC exit", r.returncode, "kept file", kept.exists(), "AO used", pick)
    return kept if kept.exists() else None


if __name__ == "__main__" and sys.argv[1] == "molusc":
    run_molusc(sys.argv[2])
