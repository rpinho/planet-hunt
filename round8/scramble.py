"""Round 8A: per-star false-alarm calibration of the round-5 search (look-elsewhere included).

For stars with NO hold-out data, the only honest significance is "how often does the identical
round-5 search, run on this star's own noise, return a signal as strong as the candidate?".

Null realisations: the star's SPOC 2-min sectors with every known-planet AND candidate transit
cadence removed first, then the remaining flux values (and errors) circularly shifted by a random
lag inside each sector over the remaining timestamps (times, gaps and cadence kept). Not inverted:
inverting turns flares into deep "transits" (v1 did this; kept in scramble_v1_inverted/).
Within-sector noise and systematics are preserved; any coherent multi-sector signal is destroyed.
A null signal counts only if it would pass the round-5 automatic triage: >= 3 transits, no EB
flags, season chi2 p > 1e-3, 5%-clip SNR >= 6, not at a TESS-orbit alias. The full round-5 search
(planet_hunt/hostsearch.search_star, 0.5-30 d, same masks and detrending) is then run and the
strongest returned signal's SNR (>= 3 transits) recorded.

    .r6venv/bin/python round8/scramble.py --n 20 --workers 3 [--only 198385543 ...]
Output: runs/round8_p1/scramble/<tic>_<k>.json, summary via --collect.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from planet_hunt import config  # noqa: E402

OUT = config.RUNS / "round8_p1" / "scramble"
STARS = {198385543: ("toi1846", 0.944060, 8.565), 354944123: ("toi4342_d2", 1.957689, 8.41),
         286201103: ("toi6008", 7.220152, 7.93), 179985715: ("toi249", 19.09285, 8.66),
         177254988: ("toi4405", 3.600819, 8.79)}
R5_TARGETS = Path(__file__).resolve().parent.parent / "runs/round5_hosts/targets.csv"


def job(args):
    tic, k = args
    import pandas as pd
    import round5_hosts as r5
    from planet_hunt import hostsearch as hs
    path = OUT / f"{tic}_{k:03d}.json"
    if path.exists():
        return tic, k, "cached"
    t_start = time.time()
    try:
        row = pd.read_csv(R5_TARGETS).set_index("tic").loc[tic]
        secs = r5.load_lc(tic)
        rng = np.random.default_rng(1000 * tic % (2 ** 31) + k)
        lags = []
        known = json.loads(row["known"])
        key, Pc, _ = STARS[tic]
        cand = json.load(open(OUT.parent / f"{key}_lcgate.json"))["cand"]
        windows = [(kk["period"], kk["t0"], kk["duration_h"] / 24) for kk in known
                   if np.isfinite(kk["period"]) and np.isfinite(kk["t0"])]
        windows.append((cand["P"], cand["t0"], cand["dur_h"] / 24))
        if key == "toi4342_d2":
            windows.append((18.908226, 2904.5666, 2.44 / 24))
        for d in secs:
            t = np.asarray(d["time"], float)
            m = np.zeros(len(t), bool)
            for (P, t0, du) in windows:
                m |= np.abs(((t - t0) / P + 0.5) % 1 - 0.5) * P < 0.75 * du + 0.02
            keep = ~m
            for kk in ("time", "flux", "ferr", "cx", "cy"):
                if kk in d and hasattr(d[kk], "__len__") and len(d[kk]) == len(t):
                    d[kk] = np.asarray(d[kk])[keep]
            n = int(keep.sum())
            lag = int(rng.integers(n // 10, n - n // 10))
            lags.append(lag)
            d["flux"] = np.roll(np.asarray(d["flux"], float), lag)
            d["ferr"] = np.roll(np.asarray(d["ferr"], float), lag)
        teff = float(row["teff"]) if np.isfinite(row["teff"]) else 3500.0
        res = hs.search_star(secs, known, float(row["rstar"]), float(row["mstar"]), teff, max_signals=3)

        keep_keys = ("period", "t0", "duration_h", "snr", "depth_ppm", "n_transits", "snr_clip5", "season_chi2_p",
                     "oe_sigma", "eb_flags", "cautions", "tr_chi2_p", "dur_ratio", "ratio_flags", "known_overlap", "sde")
        out = {"tic": int(tic), "k": k, "lags": lags, "elapsed_s": time.time() - t_start,
               "signals": [{kk: s.get(kk) for kk in keep_keys} for s in res.get("signals", [])],
               "star_peaks": res.get("star_peaks", [])[:6], "saturated": res.get("rotation_saturated")}
        path.write_text(json.dumps(out, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
        return tic, k, f"signals {[(round(x['period'], 3), round(x['snr'], 1)) for x in out['signals']]} {out['elapsed_s']:.0f}s"
    except Exception as exc:
        (OUT / f"{tic}_{k:03d}.err").write_text(traceback.format_exc())
        return tic, k, f"error {exc!r}"[:150]


def clip_floor(snr, n, n_sim=4000, seed=0):
    rng = np.random.default_rng(seed)
    k = max(1, int(np.ceil(0.05 * n)))
    z = rng.normal(snr / np.sqrt(n), 1.0, size=(n_sim, n))
    z.sort(axis=1)
    return float(np.percentile(z[:, :-k].sum(1) / np.sqrt(n - k), 1))


def triage(s):
    """The round-5 automatic survivor rules (round5_triage.py) that apply to a single star."""
    if not (s.get("n_transits") or 0) >= 3:
        return False
    if s.get("eb_flags"):
        return False
    if (s.get("known_overlap") or 0) > 0.5:
        return False
    n = int(s["n_transits"])
    p = s.get("tr_chi2_p")
    if n < 6 and p is not None and np.isfinite(p) and p < 1e-4:
        return False
    rf = s.get("ratio_flags") or ""
    if "rotation" in str(rf):
        return False
    sp = s.get("season_chi2_p")
    if sp is not None and np.isfinite(sp) and sp < 1e-3:
        return False
    c5 = s.get("snr_clip5")
    if c5 is not None and np.isfinite(c5) and c5 < clip_floor(s["snr"], n):
        return False
    dr = s.get("dur_ratio")
    if dr is not None and np.isfinite(dr) and dr > 1.6:
        return False
    return True


def collect():
    summ = {}
    for tic, (key, P, snr_c) in STARS.items():
        v, vr = [], []
        short = P < 6.0                       # round-5 search ranges: 0.5-6 d and 6-30 d
        for f in sorted(OUT.glob(f"{tic}_*.json")):
            d = json.load(open(f))
            ok = [x for x in d["signals"] if triage(x)]
            v.append(max([x["snr"] for x in ok], default=0.0))    # nothing passed triage -> 0
            vr.append(max([x["snr"] for x in ok if (x["period"] < 6.0) == short], default=0.0))
        v, vr = np.array(v), np.array(vr)
        if not len(v):
            continue
        # Gumbel fit to the per-realisation maxima for the tail
        from scipy.stats import gumbel_r
        # Gumbel fit to the realisations that produced a triage-passing signal; realisations with
        # none (recorded as 0) count as below the candidate: p = frac_nonzero x Gumbel tail
        nz = v[v > 0]
        loc, sc = gumbel_r.fit(nz) if len(nz) > 3 else (np.nan, np.nan)
        p_g = float(gumbel_r.sf(snr_c, loc, sc) * len(nz) / len(v)) if len(nz) > 3 else None
        summ[key] = {"tic": tic, "candidate_snr_r5": snr_c, "n": int(len(v)), "null_max_snr": v.round(2).tolist(),
                     "median": float(np.median(v)), "max": float(v.max()), "n_ge_candidate": int((v >= snr_c).sum()),
                     "p_emp": float((1 + (v >= snr_c).sum()) / (1 + len(v))), "p_gumbel": p_g,
                     "frac_with_signal": float(np.mean(v > 0)),
                     "expected_in_184_star_survey": p_g * 184 if p_g is not None else None,
                     "same_range": {"range": "0.5-6 d" if short else "6-30 d", "null_max_snr": vr.round(2).tolist(),
                                    "n_ge_candidate": int((vr >= snr_c).sum()),
                                    "p_emp": float((1 + (vr >= snr_c).sum()) / (1 + len(vr))),
                                    "frac_with_signal": float(np.mean(vr > 0)),
                                    "p_gumbel": float(gumbel_r.sf(snr_c, *gumbel_r.fit(vr[vr > 0])) * np.mean(vr > 0)) if (vr > 0).sum() > 3 else None}}
        print(key, summ[key]["n"], "null max SNR median", round(summ[key]["median"], 2), "max", round(v.max(), 2),
              "n>=cand", summ[key]["n_ge_candidate"], "p_gumbel", p_g, "x184", (p_g or 0) * 184,
              "| same range", summ[key]["same_range"]["range"], "n>=cand", summ[key]["same_range"]["n_ge_candidate"],
              "p_gumbel", summ[key]["same_range"]["p_gumbel"], "max", max(summ[key]["same_range"]["null_max_snr"]))
    json.dump(summ, open(OUT.parent / "scramble_summary.json", "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--only", type=int, nargs="*")
    ap.add_argument("--collect", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if a.collect:
        return collect()
    tics = a.only or list(STARS)
    jobs = [(t, k) for k in range(a.n) for t in tics]
    with cf.ProcessPoolExecutor(a.workers) as ex:
        for tic, k, msg in ex.map(job, jobs, chunksize=1):
            print(f"{tic} #{k}: {msg}", flush=True)


if __name__ == "__main__":
    import os
    os.chdir(Path(__file__).resolve().parent.parent)
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    main()
