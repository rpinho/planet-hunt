"""Round 5 track A: assemble the vetting evidence for every triage survivor.

    .venv/bin/python round5_assemble.py      # runs/round5_hosts/candidates.csv + candidates.json

Per candidate: search statistics (2-min, triage.csv), pixel localisation of
the candidate and of the star's known planet on the same pixels (deep_vet),
fixed-aperture TESScut depth / odd-even / secondary and its ratio to the 2-min
depth compared with the same ratio for the known planet, SPOC TCE matches,
TRICERATOPS FPP/NFPP, and the hand-written verdict from
vetting/round5_verdicts.json.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from planet_hunt import config
import round5_hosts as r5
from round5_deep import _name, DEEP

OUT = r5.OUT


def _load(name):
    p = DEEP / f"{name}.json"
    return json.load(open(p)) if p.exists() else None


def control_2min_depth(tic, spec):
    """Box depth of the known planet in OUR 2-min light curve (other known
    planets and the candidate masked), for a like-for-like fixed/2-min ratio."""
    from planet_hunt import hostsearch as hs, multisector as ms
    row = r5._target_row(int(tic))
    secs = r5.load_lc(int(tic))
    P, t0, dur = spec["period"], spec["t0"], spec["duration_h"] / 24
    others = [w for w in spec["exclude_windows"]]
    om = lambda t: np.any([ms.in_tr(t, a, b, c / 1.5, 1.5) for a, b, c in others], axis=0) if others else np.zeros(len(t), bool)
    lc = hs.detrend(secs, hs.W_SHORT if P < hs.SPLIT_P else hs.W_LONG,
                    lambda t: om(t) | ms.in_tr(t, P, t0, dur, 1.6), None)
    lc = hs.subset(lc, ~om(lc["t"]))
    st = ms.transit_stats(lc, P, t0, dur, row["rstar"], row["mstar"], row["teff"])
    return st["depth_ppm"], st["depth_err_ppm"]


def main():
    tri = pd.read_csv(OUT / "triage.csv")
    sel = tri[tri.status.isin(["survivor", "flagged (rotation-saturated star)"])].copy()
    specs = json.load(open(config.ROOT / "vetting" / "round5_targets.json"))["targets"]
    tces = json.load(open(OUT / "spoc_tces.json")) if (OUT / "spoc_tces.json").exists() else {}
    vpath = config.ROOT / "vetting" / "round5_verdicts.json"
    verdicts = json.load(open(vpath)) if vpath.exists() else {}
    rows = []
    for s in sel.itertuples():
        name = _name(s.tic, s.period)
        c = _load(name)
        ctl_spec = next((x for x in specs if x.get("control_of") == name), None)
        ctl = _load(ctl_spec["name"]) if ctl_spec else None
        r = {"tic": int(s.tic), "host_toi": None if s.host_toi != s.host_toi else int(s.host_toi),
             "name": name, "period": s.period, "t0": s.t0, "duration_h": s.duration_h,
             "depth_ppm": s.depth_ppm, "depth_err_ppm": s.depth_err_ppm, "snr": s.snr,
             "snr_clip5": s.snr_clip5, "n_transits": int(s.n_transits), "rp_earth": s.rp_earth,
             "sde": s.sde, "oe_sigma": s.oe_sigma, "sec_sigma": s.sec_sigma, "season_p": s.season_chi2_p,
             "crowd_slope_sigma": s.crowd_slope_sigma, "tmag": s.tmag, "teff": s.teff, "rstar": s.rstar,
             "n2min": int(s.n2min), "triage": s.status}
        if c:
            L = c["localization"]
            r.update({"psf_dchi2_target": L["target_dchi2"], "psf_offset_arcsec": L["offset_arcsec"],
                      "psf_offset_err": L["offset_err_arcsec"], "psf_best_sep": L["hypotheses"][0]["sep_arcsec"],
                      "psf_best_tmag": L["hypotheses"][0]["tmag"], "psf_n_sectors": L["n_sectors"],
                      "fixed_depth_ppm": c["fixed_lc"]["depth_ppm"], "fixed_err_ppm": c["fixed_lc"]["err_ppm"],
                      "fixed_oe_sigma": c["fixed_lc"].get("oe_sigma"), "fixed_sec_sigma": c["fixed_lc"].get("sec_sigma"),
                      "rho_fit": (c.get("transit_fit") or {}).get("rho_cgs"),
                      "rho_tic": (c.get("transit_fit") or {}).get("rho_tic_cgs"),
                      "gaia_ruwe": c["gaia_target"].get("ruwe"),
                      "n_gaia_63as": len(c["gaia_neighbours_63as"]),
                      "closest_gaia": (min(c["gaia_neighbours_63as"], key=lambda g: g["sep_arcsec"])
                                       if c["gaia_neighbours_63as"] else None)})
            r["fixed_over_2min"] = r["fixed_depth_ppm"] / s.depth_ppm
            r["fixed_over_2min_err"] = r["fixed_err_ppm"] / s.depth_ppm
        if ctl:
            r.update({"ctl_name": ctl_spec.get("known_name"), "ctl_psf_dchi2": ctl["localization"]["target_dchi2"],
                      "ctl_fixed_depth_ppm": ctl["fixed_lc"]["depth_ppm"], "ctl_fixed_err_ppm": ctl["fixed_lc"]["err_ppm"],
                      "ctl_catalog_depth_ppm": ctl_spec["depth_ppm"]})
            r["ctl_fixed_over_catalog"] = r["ctl_fixed_depth_ppm"] / ctl_spec["depth_ppm"]
            d2, e2 = control_2min_depth(s.tic, ctl_spec)
            r["ctl_2min_depth_ppm"], r["ctl_2min_err_ppm"] = d2, e2
            r["ctl_fixed_over_2min"] = r["ctl_fixed_depth_ppm"] / d2
            r["ctl_fixed_over_2min_err"] = r["ctl_fixed_over_2min"] * np.hypot(r["ctl_fixed_err_ppm"] / r["ctl_fixed_depth_ppm"], e2 / d2)
            if "fixed_over_2min" in r:
                q = r["fixed_over_2min"] / r["ctl_fixed_over_2min"]
                r["ratio_vs_ctl"] = q
                r["ratio_vs_ctl_sigma"] = (q - 1) / (q * np.hypot(r["fixed_over_2min_err"] / r["fixed_over_2min"],
                                                                   r["ctl_fixed_over_2min_err"] / r["ctl_fixed_over_2min"]))
        t = tces.get(str(int(s.tic)), {})
        hits = t.get("candidates", {}).get(f"{s.period:.5f}", [])
        r["spoc_tce"] = "; ".join(f"{h['run']} P={h['period']:.4f} MES {h['mes']:.1f}" for h in hits) or (
            "none" if t else "not checked")
        tp = DEEP / f"triceratops_{name}.json"
        if tp.exists():
            tr = json.load(open(tp))
            r.update({"fpp": tr["FPP_mean"], "fpp_std": tr["FPP_std"], "nfpp": tr["NFPP_mean"]})
        v = verdicts.get(name, {})
        r.update({"lit": v.get("literature", ""), "verdict": v.get("verdict", ""), "class": v.get("class", ""),
                  "notes": v.get("notes", "")})
        rows.append(r)
    df = pd.DataFrame(rows).sort_values("snr", ascending=False)
    df.to_csv(OUT / "candidates.csv", index=False)
    (OUT / "candidates.json").write_text(json.dumps(rows, indent=1, default=str))
    cols = ["name", "host_toi", "period", "depth_ppm", "snr", "rp_earth", "psf_dchi2_target", "ctl_psf_dchi2",
            "fixed_over_2min", "ctl_fixed_over_2min", "ratio_vs_ctl_sigma", "fixed_oe_sigma", "fixed_sec_sigma",
            "spoc_tce", "fpp", "nfpp"]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(df[[c for c in cols if c in df]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
