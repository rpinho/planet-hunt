"""Round 8A: independent photons from the ground. For each TOI host, download the public
TFOP SG1 photometry tables on ExoFOP (read-only), find any light curve that happens to cover a
predicted transit of the CANDIDATE (they were all taken for other planets), and measure the
candidate's box depth there.

    .venv/bin/python round8/ground.py fetch      # tables -> cache/r8/exofop/<tic>/
    .r6venv/bin/python round8/ground.py measure  # -> runs/round8_p1/ground.json

Depth: relative target flux, detrended per light curve with a linear (airmass-like) baseline
fitted out of the predicted transit (+-0.75 duration + timing margin), box of the candidate
duration at the predicted time; error = OOT scatter x red-noise beta.
"""
from __future__ import annotations

import io
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from planet_hunt import config  # noqa: E402
from targets import CANDS  # noqa: E402

EXO = config.CACHE / "r8" / "exofop"
OUT = config.RUNS / "round8_p1"
GET = "https://exofop.ipac.caltech.edu/tess/get_file.php?id="


def curl(url, path):
    for k in range(6):
        r = subprocess.run(["curl", "-s", "-f", "-L", "--retry", "3", "--retry-all-errors", "--connect-timeout", "30",
                            "-m", "600", "-o", str(path), url], capture_output=True)
        if r.returncode == 0 and path.exists() and path.stat().st_size > 200:
            return True
        time.sleep(10 * (k + 1))
    return False


def table_files(tic):
    d = json.load(open(EXO / f"target_{tic}.json"))
    tags = {x["tstag"]: x for x in d["time_series"]}
    out = []
    for f in d["files"]:
        if f["ftag"] not in tags:
            continue
        n = f["fname"].lower()
        if (n.endswith(".tbl") or n.endswith(".csv") or (n.endswith(".txt") and "measure" in n) or n.endswith(".dat")) \
                and "neb" not in n and "notes" not in n and "readme" not in n:
            out.append((f, tags[f["ftag"]]))
    return out


def cmd_fetch():
    for key, c in CANDS.items():
        if not (EXO / f"target_{c['tic']}.json").exists():
            continue
        d = EXO / str(c["tic"])
        d.mkdir(exist_ok=True)
        for f, ts in table_files(c["tic"]):
            p = d / f["fname"]
            if not p.exists():
                ok = curl(GET + str(f["fid"]), p)
                print(key, "ok" if ok else "FAIL", f["fname"], flush=True)


def read_table(p):
    """AIJ measurement tables (tab separated, header line) and MuSCAT2-style CSV/whitespace tables."""
    txt = open(p, errors="ignore").read()
    lines = [l for l in txt.splitlines() if l.strip()]
    if not lines:
        return None
    hdr_i = 0
    for i, l in enumerate(lines[:50]):
        ll = l.lower()
        if ("bjd" in ll or "time" in ll or "jd" in ll) and ("flux" in ll or "rel" in ll or "mag" in ll):
            hdr_i = i
            break
    hdr = lines[hdr_i].rstrip("\r\n").lstrip("#")
    sep = "\t" if "\t" in hdr else ("," if "," in hdr else None)
    if sep is None:
        hdr = hdr.strip()
    cols = [c.strip() for c in (hdr.split(sep) if sep else hdr.split())]
    rows = []
    for l in lines[hdr_i + 1:]:
        parts = [x.strip() for x in (l.split(sep) if sep else l.split())]
        if len(parts) < len(cols) - 1:
            continue
        rows.append(parts[:len(cols)])
    if not rows:
        return None
    if len(rows[0]) == len(cols) - 1 and cols and cols[0] == "":
        cols = cols[1:]

    def col(name_opts):
        for n in name_opts:
            for j, cn in enumerate(cols):
                if cn.lower() == n.lower():
                    try:
                        return np.array([float(r[j]) if j < len(r) else np.nan for r in rows])
                    except ValueError:
                        continue
        return None
    t = col(["BJD_TDB", "BJD_TDB_MOBS", "bjd_tdb", "BJD", "time_bjd", "time"])
    if t is None:
        tj = [j for j, cn in enumerate(cols) if "bjd" in cn.lower()]
        if tj:
            t = np.array([float(r[tj[0]]) for r in rows])
    f = col(["rel_flux_T1", "rel_flux_T1_n", "flux", "rel_flux", "flux_target", "Flux"])
    e = col(["rel_flux_err_T1", "rel_flux_err_T1_n", "flux_err", "rel_flux_err", "err", "Err"])
    if f is None and "flux_star_0" in cols:
        # MuSCAT2 raw photometry: target = star 0; comparison = sum of the next 5 brightest stars
        fs = {cn: col([cn]) for cn in cols if cn.startswith("flux_star_")}
        tgt = fs.pop("flux_star_0")
        med = sorted(fs.items(), key=lambda kv: -np.nanmedian(kv[1]))[:5]
        comp = np.sum([v for _, v in med], axis=0)
        f = tgt / comp
        e = None
    if f is None:
        fj = [j for j, cn in enumerate(cols) if "rel_flux_t1" in cn.lower() and "err" not in cn.lower()]
        if fj:
            f = np.array([float(r[fj[0]]) for r in rows])
    if t is None or f is None:
        return {"cols": cols[:40]}
    if np.nanmedian(t) < 2.4e6 and np.nanmedian(t) > 50000:      # MJD-like
        t = t + 2400000.5 if np.nanmedian(t) < 1e5 else t
    ok = np.isfinite(t) & np.isfinite(f)
    return {"t": t[ok], "f": f[ok] / np.nanmedian(f[ok]), "e": e[ok] if e is not None else None, "cols": cols[:40]}


def cmd_measure():
    res = {}
    for key, c in CANDS.items():
        d = EXO / str(c["tic"])
        if not d.exists():
            continue
        P, t0, dur = c["P"], c["t0"] + 2457000, c["dur_h"] / 24
        eph = OUT / f"{key}_ephem.json"
        sig_d = 0.02
        if eph.exists():
            e = json.load(open(eph))
            P, t0 = e["P"], e["T0_bjd"]
        rows = []
        for f, ts in table_files(c["tic"]):
            p = d / f["fname"]
            if not p.exists():
                continue
            tb = read_table(p)
            if not tb or "t" not in tb:
                rows.append({"file": f["fname"], "error": "unparsed", "cols": (tb or {}).get("cols")})
                continue
            t, fl = tb["t"], tb["f"]
            if eph.exists():
                n = np.round((np.median(t) - t0) / P)
                sig_d = float(np.sqrt(e["T0_err"] ** 2 + (n * e["P_err"]) ** 2))
            ph = ((t - t0) / P + 0.5) % 1 - 0.5
            x = ph * P
            n_ep = np.round((t - t0) / P)
            tc = t0 + np.median(n_ep) * P
            cover_in = np.mean(np.abs(t - tc) < dur / 2) if len(t) else 0
            inwin = np.abs(x) < dur / 2
            frac_cov = 0.0
            if len(t) > 5:
                grid = np.linspace(tc - dur / 2, tc + dur / 2, 50)
                cad = np.median(np.diff(np.sort(t)))
                frac_cov = float(np.mean([np.any(np.abs(t - g) < max(cad, 1 / 1440)) for g in grid]))
            row = {"file": f["fname"], "date": ts["tsdate"], "tel": ts["tstel"], "filt": ts["tsfilt"], "toi": ts["tstoi"],
                   "t_start": float(t.min()), "t_end": float(t.max()), "n": int(len(t)),
                   "nearest_pred_tc": float(tc), "pred_minus_start_h": float((tc - t.min()) * 24),
                   "pred_minus_end_h": float((tc - t.max()) * 24), "timing_sigma_min": sig_d * 1440,
                   "in_transit_coverage": frac_cov}
            # measure only if the predicted window (incl. timing 1 sigma) is mostly covered with OOT on a side
            margin = sig_d
            oot = (np.abs(x) > dur / 2 + margin) & (np.abs(x) < 3 * dur + margin)
            if frac_cov > 0.5 and oot.sum() > 10 and inwin.sum() > 5:
                A = np.vstack([np.ones(oot.sum()), t[oot] - tc]).T
                b, *_ = np.linalg.lstsq(A, fl[oot], rcond=None)
                base = b[0] + b[1] * (t - tc)
                r = fl / base - 1
                s = np.std(r[oot])
                depth = -np.mean(r[inwin])
                err = s * np.sqrt(1 / inwin.sum() + 1 / oot.sum())
                # red noise at the duration scale
                cad = np.median(np.diff(np.sort(t)))
                nb = max(int(dur / 2 / cad), 1)
                if oot.sum() > 4 * nb and nb > 1:
                    ro = r[oot][: (oot.sum() // nb) * nb].reshape(-1, nb).mean(1)
                    beta = max(1.0, np.std(ro) / (s / np.sqrt(nb)))
                else:
                    beta = 1.0
                row.update(depth_ppm=depth * 1e6, err_ppm=err * beta * 1e6, beta=beta, expected_ppm=c["depth_ppm"],
                           z=depth / (err * beta), z_vs_expected=(depth - c["depth_ppm"] * 1e-6) / (err * beta))
            rows.append(row)
        res[key] = rows
        cov = [r for r in rows if r.get("in_transit_coverage", 0) > 0]
        print(key, len(rows), "tables;", len(cov), "touch a predicted candidate transit", flush=True)
        for r in cov:
            print("   ", r["date"], r["tel"], r["filt"], f"cov {r['in_transit_coverage']:.2f}",
                  f"depth {r.get('depth_ppm', float('nan')):.0f} +- {r.get('err_ppm', float('nan')):.0f}", r["file"][:60])
    json.dump(res, open(OUT / "ground.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    {"fetch": cmd_fetch, "measure": cmd_measure}[sys.argv[1]]()
