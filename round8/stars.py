"""Round 8A: host-star parameters, Gaia DR3 RV/binarity checks, and spectral type from the
Mamajek dwarf table (version 2022.04.16) by several independent indicators.

    .venv/bin/python round8/stars.py   -> runs/round8_p1/stars.json

Spectral type: nearest table row in each of Teff, Bp-Rp, M_G, M_Ks and R*, then the median
row index (no single indicator is trusted; round 7 got this wrong once by using Teff alone).

Gaia RV test: DR3 gives per-source summary statistics, not epoch RVs. For a star with
rv_nb_transits epochs, rv_chisq_pvalue < 0.01 (with rv_renormalised_gof > 4) flags variability
(Katz et al. 2023). The robust peak-to-peak (rv_amplitude_robust) bounds the RV semi-amplitude
roughly at K <~ amplitude/2; we compare it with the K an eclipsing stellar companion of
0.08 Msun would induce at the candidate period. Gaia RVs only exist for G_RVS <~ 14.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from planet_hunt import config  # noqa: E402
from targets import CANDS  # noqa: E402

OUT = config.RUNS / "round8_p1"
EEM = config.CACHE / "r8" / "ref" / "EEM_dwarf_UBVIJHK_colors_Teff.txt"


def load_eem():
    rows = []
    started = False
    for line in open(EEM):
        if line.startswith("#SpT"):
            if started:
                break
            started = True
            hdr = line[1:].split()
            continue
        if not started or line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 10:
            continue
        rows.append(parts)
    names = hdr
    tab = {n: [] for n in names}
    for r in rows:
        for n, v in zip(names, r):
            tab[n].append(v)
    def col(n):
        out = []
        for v in tab[n]:
            try:
                out.append(float(v.replace("*", "")) if "." in v or v.lstrip("-").isdigit() else np.nan)
            except ValueError:
                out.append(np.nan)
        return np.array(out)
    return tab["SpT"], {n: col(n) for n in ("Teff", "Bp-Rp", "M_G", "M_Ks", "R_Rsun", "Msun")}


def spt(values):
    names, T = load_eem()
    picks = {}
    for k, v in values.items():
        if v is None or not np.isfinite(v):
            continue
        x = T[k]
        ok = np.isfinite(x)
        i = int(np.nanargmin(np.where(ok, np.abs(x - v), np.inf)))
        picks[k] = (i, names[i])
    idx = int(np.round(np.median([i for i, _ in picks.values()])))
    return {"per_indicator": {k: n for k, (i, n) in picks.items()}, "adopted": names[idx],
            "range": [names[min(i for i, _ in picks.values())], names[max(i for i, _ in picks.values())]]}


def gaia(source_id):
    from astroquery.gaia import Gaia
    Gaia.ROW_LIMIT = 5
    q = f"""SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag, bp_rp,
      ruwe, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude, non_single_star, phot_variable_flag,
      radial_velocity, radial_velocity_error, rv_nb_transits, rv_chisq_pvalue, rv_renormalised_gof,
      rv_amplitude_robust, rv_method_used, grvs_mag, teff_gspphot, duplicated_source
      FROM gaiadr3.gaia_source WHERE source_id = {source_id}"""
    for k in range(5):
        try:
            r = Gaia.launch_job(q).get_results()
            if len(r) == 0:
                return None
            return {c: (r[c][0].item() if hasattr(r[c][0], "item") else r[c][0]) for c in r.colnames}
        except Exception as exc:
            print("gaia retry", repr(exc)[:100]); time.sleep(15 * (k + 1))


def gaia_cone(ra, dec, gmag_hint=None):
    """Fallback when the TIC's Gaia ID is not a DR3 ID: brightest DR3 source within 10 arcsec."""
    from astroquery.gaia import Gaia
    q = f"""SELECT TOP 5 source_id, phot_g_mean_mag FROM gaiadr3.gaia_source
      WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {10 / 3600}))
      ORDER BY phot_g_mean_mag ASC"""
    r = Gaia.launch_job(q).get_results()
    return int(r["source_id"][0]) if len(r) else None


def nss(source_id):
    from astroquery.gaia import Gaia
    out = {}
    for tab in ("nss_two_body_orbit", "nss_acceleration_astro", "nss_non_linear_spectro", "vari_eclipsing_binary"):
        for k in range(3):
            try:
                r = Gaia.launch_job(f"SELECT source_id FROM gaiadr3.{tab} WHERE source_id = {source_id}").get_results()
                out[tab] = len(r)
                break
            except Exception as exc:
                print("nss retry", tab, repr(exc)[:80]); time.sleep(10)
    return out


def k_companion(P_d, m1, m2=0.08):
    G, Ms = 6.674e-11, 1.989e30
    P = P_d * 86400
    return (2 * np.pi * G / P) ** (1 / 3) * m2 * Ms / ((m1 + m2) * Ms) ** (2 / 3) / 1000   # km/s, sin i = 1


def main():
    inv = json.load(open(OUT / "inventory.json"))
    res = {}
    for key, c in CANDS.items():
        t = inv[key]["tic_row"]
        sid = int(t["GAIA"]) if t.get("GAIA") not in (None, "None", "--", "") else None
        g = gaia(sid) if sid else None
        if not g:
            sid = gaia_cone(inv[key]["ra"], inv[key]["dec"])
            g = gaia(sid) if sid else {}
        g = g or {}
        dist = 1000 / g["parallax"] if g and g.get("parallax") else t.get("d")
        mg = g["phot_g_mean_mag"] - 5 * np.log10(dist / 10) if g else None
        mks = t["Kmag"] - 5 * np.log10(dist / 10) if t.get("Kmag") else None
        sp = spt({"Teff": t["Teff"], "Bp-Rp": g.get("bp_rp") if g else None, "M_G": mg, "M_Ks": mks,
                  "R_Rsun": t["rad"]})
        kc = k_companion(c["P"], t["mass"] or c["mstar"])
        rv = {}
        if g and g.get("rv_nb_transits"):
            amp = g.get("rv_amplitude_robust")
            rv = {"n_epochs": g["rv_nb_transits"], "chisq_p": g.get("rv_chisq_pvalue"),
                  "renorm_gof": g.get("rv_renormalised_gof"), "amp_robust_kms": amp,
                  "rv_err_kms": g.get("radial_velocity_error"), "grvs": g.get("grvs_mag"),
                  "K_0.08Msun_companion_kms": kc,
                  "excludes_stellar_EB_on_target": bool(amp is not None and np.isfinite(amp) and amp < kc and (g.get("rv_chisq_pvalue") is None or g["rv_chisq_pvalue"] > 0.01))}
        res[key] = {"tic": c["tic"], "gaia_dr3": sid, "tic_row": t, "gaia": g, "nss": nss(sid) if sid else {},
                    "dist_pc": dist, "M_G": mg, "M_Ks": mks, "spt": sp, "rho_tic_cgs": (t["mass"] or c["mstar"]) / (t["rad"] or c["rstar"]) ** 3 * 1.408,
                    "gaia_rv": rv}
        print(key, sp["adopted"], sp["per_indicator"], "RUWE", g.get("ruwe"), "RV", rv, flush=True)
    json.dump(res, open(OUT / "stars.json", "w"), indent=1, default=lambda o: None if o is np.ma.masked else str(o))


if __name__ == "__main__":
    main()
