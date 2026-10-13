"""Round 5 track A: collect search results, injection-recovery, automatic triage.

    .venv/bin/python round5_triage.py collect     # runs/round5_hosts/signals.csv, triage.csv
    .venv/bin/python round5_triage.py inject      # runs/round5_hosts/completeness.json + plot

Triage (signal level), in order:
  known      the period matches a known TOI/planet of the star (1:1) or more
             than half of its in-transit cadences overlap a known transit
             window (imperfect masking), or it is a TOI/CTOI/confirmed planet
             of the star at any simple alias
  threshold  SNR (red-noise) >= 7.5 and >= 3 transits
  systematic same period (0.1%) on >= 2 other stars
  automatic  multisector.vet red flags (odd/even > 3 sigma, secondary,
             sinusoidal, V-shape, few deep events), rotation harmonic,
             momentum-dump excess, per-season depth p < 1e-3, depth vs
             CROWDSAP slope > 3 sigma, robust 5%-clip SNR < 6 or an SNR drop > 2 beyond the noise expectation
Survivors go to deep vetting (pixel localisation, TRICERATOPS).
"""
from __future__ import annotations

import argparse
import glob
import json

import numpy as np
import pandas as pd

from planet_hunt import config
import round5_hosts as r5

OUT = r5.OUT
SNR_MIN, NTR_MIN = 7.5, 3


def load_results():
    rows, stars = [], []
    for f in sorted(glob.glob(str(OUT / "search" / "*.json"))):
        res = json.load(open(f))
        st = res["row"]
        stars.append({**st, "rot_P": res["rotation"].get("P_rot"), "rot_amp_ppm": res["rotation"].get("amp_ppm"),
                      "active": res["rotation"].get("active"), "rotation_saturated": res["rotation_saturated"],
                      "saturation_reasons": "; ".join(res["saturation_reasons"]), "n_sectors": res["n_sectors"],
                      "T_all": res["T_all"], "noise_ppm_2min": res["noise_ppm_2min"],
                      "n_known": len(res["known"]), "n_signals": len(res["signals"])})
        for s in res["signals"]:
            rows.append({"tic": st["tic"], "host_toi": st.get("host_toi"), "tmag": st["tmag"], "teff": st["teff"],
                         "rstar": st["rstar"], "n2min": st["n2min"], "elat": st["elat"],
                         "rotation_saturated": res["rotation_saturated"], "rot_P": res["rotation"].get("P_rot"),
                         **{k: s.get(k) for k in ("iteration", "period", "t0", "duration_h", "depth_ppm",
                                                  "depth_err_ppm", "snr", "snr_white", "snr_clip5", "beta_med",
                                                  "n_transits", "sde", "oe_sigma", "sec_sigma", "sec_any_sigma",
                                                  "sin_ratio", "shape", "tr_chi2_p", "season_chi2_p", "rp_earth",
                                                  "dur_ratio", "known_overlap", "range", "quick_snr")},
                         "eb_flags": "; ".join(s["eb_flags"]), "cautions": "; ".join(s["cautions"]),
                         "ratio_flags": "; ".join(s["ratio_flags"]),
                         "dumps_near": s["dumps"].get("n_near_transit"), "dumps_exp": s["dumps"].get("expected"),
                         "known_periods": "; ".join(f"{k['name']} {k['period']:.5f}" for k in res["known"]),
                         "per_sector": json.dumps(s.get("per_sector", []))})
    return pd.DataFrame(rows), pd.DataFrame(stars)


def crowd_slope(tic, per_sector):
    """Depth vs (1 - CROWDSAP) across sectors: a blend's PDCSAP depth tracks the
    contamination correction. Returns slope significance (sigma)."""
    try:
        secs = r5.load_lc(int(tic))
    except Exception:
        return np.nan
    cr = {d["meta"]["sector"]: d["meta"].get("crowdsap", np.nan) for d in secs}
    ps = [p for p in per_sector if np.isfinite(cr.get(p["sector"], np.nan))]
    if len(ps) < 4:
        return np.nan
    x = np.array([1 - cr[p["sector"]] for p in ps]); d = np.array([p["depth_ppm"] for p in ps])
    e = np.array([p["err_ppm"] for p in ps])
    if np.ptp(x) < 0.01:
        return 0.0
    w = 1 / e ** 2
    A = np.vstack([np.ones_like(x), x - np.average(x, weights=w)]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    coef = cov @ (A.T @ (w * d))
    return float(coef[1] / np.sqrt(cov[1, 1]))


_CAT = {}


def _cats():
    if not _CAT:
        toi, ps = r5.load_catalogs()
        ctoi = pd.read_csv(r5.R5 / "ctoi_r5.csv", comment="#", low_memory=False)
        ctoi = ctoi.rename(columns={"TIC ID": "tic", "Period (days)": "period", "CTOI": "ctoi"})
        _CAT.update(toi=toi, ps=ps, ctoi=ctoi)
    return _CAT["toi"], _CAT["ps"], _CAT["ctoi"]


def catalog_match(tic, P):
    """TOI / CTOI / confirmed planets of this TIC at P or a simple alias."""
    toi, ps, ctoi = _cats()
    hits = []
    fr = (1, 2, 0.5, 3, 1 / 3, 1.5, 2 / 3, 4, 0.25)
    for df, name, pc in ((toi[toi.tic == tic], "TOI", "period"), (ctoi[ctoi.tic == tic], "CTOI", "period"),
                         (ps[ps.tic == tic], "confirmed", "pl_orbper")):
        for rr in df.itertuples():
            Pk = pd.to_numeric(getattr(rr, pc), errors="coerce")
            if not np.isfinite(Pk) or Pk <= 0:
                continue
            for f in fr:
                if abs(P / (Pk * f) - 1) < 0.005:
                    lab = {"TOI": f"TOI-{getattr(rr, 'TOI', '')}", "CTOI": f"CTOI {getattr(rr, 'ctoi', '')}",
                           "confirmed": str(getattr(rr, "pl_name", ""))}[name]
                    hits.append(f"{lab} (P={Pk:.5f}, ratio {f:.3g})")
    return hits


def clip_floor(snr, n, n_sim=4000, seed=0):
    rng = np.random.default_rng(seed)
    k = max(1, int(np.ceil(0.05 * n)))
    z = rng.normal(snr / np.sqrt(n), 1.0, size=(n_sim, n))
    z.sort(axis=1)
    clip = z[:, :-k].sum(1) / np.sqrt(n - k)
    return float(np.percentile(clip, 1))


def cmd_collect(a):
    sig, stars = load_results()
    stars.to_csv(OUT / "stars.csv", index=False)
    sig.to_csv(OUT / "signals.csv", index=False)
    print(f"{len(stars)} stars searched, {len(sig)} signals; rotation-saturated stars: "
          f"{int(stars['rotation_saturated'].sum())}")
    toi, ps = r5.load_catalogs()
    # cross-star shared periods (spacecraft systematics)
    Ps = sig["period"].to_numpy()
    shared = [int(np.sum((np.abs(Ps / p - 1) < 0.001) & (sig["tic"].to_numpy() != t))) for p, t in zip(Ps, sig["tic"])]
    sig["shared_period_stars"] = shared
    status, reasons = [], []
    for k, s in sig.iterrows():
        why = []
        cm = catalog_match(int(s.tic), float(s.period))
        rf = s.ratio_flags if isinstance(s.ratio_flags, str) else ""
        is_known = bool(cm) or (s.known_overlap or 0) > 0.5 or any(
            f.startswith("P = ") and "/" not in f.split("x")[0] and f.split(" ")[2] in ("TOI", "") for f in rf.split("; ") if f)
        if cm or (s.known_overlap or 0) > 0.5:
            status.append("known"); reasons.append("; ".join(cm) or f"overlaps known transits ({s.known_overlap:.0%})")
            continue
        if not (s.snr >= SNR_MIN and s.n_transits >= NTR_MIN):
            status.append("below threshold"); reasons.append(f"SNR {s.snr:.1f}, {s.n_transits} transits")
            continue
        if s.shared_period_stars >= 2:
            why.append(f"same period on {s.shared_period_stars} other stars")
        if isinstance(s.eb_flags, str) and s.eb_flags:
            why.append(s.eb_flags)
        if s.n_transits < 6 and np.isfinite(s.tr_chi2_p) and s.tr_chi2_p < 1e-4:
            why.append(f"few transits with inconsistent depths (p={s.tr_chi2_p:.0e})")
        if "rotation" in rf:
            why.append("rotation harmonic")
        if s.rotation_saturated:
            why.append("star rotation-saturated")
        if (s.dumps_near or 0) > max(3, 3 * (s.dumps_exp or 0)):
            why.append(f"momentum dumps {s.dumps_near} vs {s.dumps_exp:.1f} expected")
        if np.isfinite(s.season_chi2_p) and s.season_chi2_p < 1e-3:
            why.append(f"season depths differ (p={s.season_chi2_p:.0e})")
        # robust clip: drop the deepest 5% of transits (>= 1). Compare with
        # what a constant signal of the same SNR in white noise gives
        # (Monte Carlo over per-transit SNRs); a signal carried by a few
        # events falls below the 1st percentile
        if np.isfinite(s.snr_clip5):
            lo = clip_floor(s.snr, int(s.n_transits))
            sig.loc[k, "clip_floor"] = lo
            if s.snr_clip5 < lo:
                why.append(f"5%-clip SNR {s.snr_clip5:.1f} below the 1% floor {lo:.1f} for a constant signal")
        cs = crowd_slope(s.tic, json.loads(s.per_sector))
        sig.loc[k, "crowd_slope_sigma"] = cs
        if np.isfinite(cs) and abs(cs) > 3:
            why.append(f"depth tracks CROWDSAP ({cs:.1f} sigma)")
        if np.isfinite(s.dur_ratio) and s.dur_ratio > 1.6:
            why.append(f"too long for the star ({s.dur_ratio:.1f}x)")
        if why == ["star rotation-saturated"]:
            status.append("flagged (rotation-saturated star)"); reasons.append("; ".join(why))
            continue
        status.append("rejected" if why else "survivor"); reasons.append("; ".join(why))
    sig["status"], sig["reasons"] = status, reasons
    sig = sig.sort_values(["status", "snr"], ascending=[True, False])
    sig.to_csv(OUT / "triage.csv", index=False)
    print(sig["status"].value_counts().to_dict())
    cols = ["tic", "host_toi", "period", "depth_ppm", "snr", "snr_clip5", "n_transits", "rp_earth", "sde", "reasons"]
    for st in ("survivor", "rejected"):
        sub = sig[sig.status == st]
        print(f"\n== {st} ({len(sub)})")
        with pd.option_context("display.width", 250, "display.max_colwidth", 90):
            print(sub[cols].head(40).to_string(index=False))


def cmd_inject(a):
    rows = []
    for f in glob.glob(str(OUT / "inject" / "*.json")):
        d = json.load(open(f))
        r = d["result"]
        rows.append({"id": d["inj"]["id"], "tic": d["tic"], "tmag": d["tmag"], "rstar": d["rstar"],
                     "n2min": d["n2min"], "P": d["inj"]["period"], "rp": d["inj"]["rp_earth"],
                     "depth_nom": d["inj"].get("depth_ppm_nominal"), "saturated": d.get("saturated"),
                     "snr_true": r.get("snr_true_quick"), "matched": r.get("matched"),
                     "recovered": r.get("recovered"), "recovered_clean": r.get("recovered_clean", False),
                     "snr_found": r.get("snr"), "alias": r.get("alias")})
    df = pd.DataFrame(rows)
    # through triage: no red flags and the star not flagged rotation-saturated
    # (the statistical cuts - clip floor, season and crowding tests - are set
    # at <= 1% false rejection and are not re-applied here)
    df["recovered_triage"] = df.recovered_clean.astype(bool) & (df.saturated != True)
    df.to_csv(OUT / "injections.csv", index=False)
    out = {"n": len(df), "recovered": float(df.recovered.mean()), "recovered_clean": float(df.recovered_clean.mean()),
           "recovered_triage": float(df.recovered_triage.mean()),
           "matched_any_snr": float(df.matched.mean())}
    grid = {}
    for (plo, phi) in ((1, 3), (3, 10), (10, 30)):
        for (rlo, rhi) in ((0.5, 1.0), (1.0, 1.5), (1.5, 2.0)):
            m = (df.P >= plo) & (df.P < phi) & (df.rp >= rlo) & (df.rp < rhi)
            grid[f"P {plo}-{phi} d, Rp {rlo}-{rhi}"] = {"n": int(m.sum()),
                                                       "recovered": float(df.recovered[m].mean()) if m.any() else None,
                                                       "recovered_triage": float(df.recovered_triage[m].mean()) if m.any() else None}
    out["grid"] = grid
    # recovery vs expected SNR (search efficiency independent of planet size)
    bins = [0, 5, 7.5, 10, 15, 1e9]
    eff = {}
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (df.snr_true >= lo) & (df.snr_true < hi)
        eff[f"{lo}-{hi if hi < 1e8 else 'inf'}"] = {"n": int(m.sum()),
                                                   "recovered": float(df.recovered[m].mean()) if m.any() else None,
                                                   "recovered_triage": float(df.recovered_triage[m].mean()) if m.any() else None}
    out["by_true_snr"] = eff
    # by planet radius alone and by Tmag
    out["by_rp"] = {f"{a}-{b}": {"n": int(((df.rp >= a) & (df.rp < b)).sum()),
                                 "recovered": float(df.recovered[(df.rp >= a) & (df.rp < b)].mean()),
                                 "recovered_triage": float(df.recovered_triage[(df.rp >= a) & (df.rp < b)].mean())}
                    for a, b in ((0.5, 0.75), (0.75, 1.0), (1.0, 1.25), (1.25, 1.5), (1.5, 1.75), (1.75, 2.0))}
    out["by_sectors"] = {lab: {"n": int(m.sum()), "recovered": float(df.recovered[m].mean())}
                         for lab, m in (("8-12 sectors", df.n2min <= 12), ("13-29", (df.n2min > 12) & (df.n2min < 30)),
                                        (">=30 (CVZ)", df.n2min >= 30))}
    for lab, m in (("saturated stars", df.saturated == True), ("other stars", df.saturated != True)):
        out[f"recovered_{lab}"] = {"n": int(m.sum()), "frac": float(df.recovered[m].mean()) if m.any() else None}
    (OUT / "completeness.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    plot_injections(df)


def plot_injections(df):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    rec = df.recovered.astype(bool)
    ax[0].scatter(df.P[~rec], df.rp[~rec], s=14, c="#bbbbbb", label="missed")
    ax[0].scatter(df.P[rec], df.rp[rec], s=14, c="#1f6fb4", label="recovered (SNR>=7.5, >=3 tr)")
    ax[0].set_xscale("log"); ax[0].set_xlabel("injected period (d)"); ax[0].set_ylabel("injected Rp (Earth radii)")
    ax[0].legend(fontsize=8, frameon=False); ax[0].set_title(f"{len(df)} injections, {rec.mean():.0%} recovered")
    b = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 20, 30, 60])
    x, y = [], []
    for lo, hi in zip(b[:-1], b[1:]):
        m = (df.snr_true >= lo) & (df.snr_true < hi)
        if m.sum() >= 3:
            x.append(0.5 * (lo + hi)); y.append(rec[m].mean())
    ax[1].plot(x, y, "o-", c="#1f6fb4"); ax[1].axvline(7.5, ls=":", c="k")
    ax[1].set_xlabel("SNR of the injected transit at its true ephemeris"); ax[1].set_ylabel("fraction recovered")
    ax[1].set_ylim(-0.03, 1.03)
    fig.tight_layout()
    (OUT / "round5_hosts").mkdir(exist_ok=True)
    fig.savefig(config.RUNS / "round5_hosts" / "injection_recovery.png", dpi=130)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("collect"); sub.add_parser("inject")
    a = ap.parse_args()
    {"collect": cmd_collect, "inject": cmd_inject}[a.cmd](a)


if __name__ == "__main__":
    main()
