"""One HTML report for methods A (M dwarfs in FFIs) and B (HARPS RVs).

    .venv/bin/python mimic_report.py
    .venv/bin/python mimic_report.py --out somewhere/report.html

Reads runs/mdwarf/* and runs/rv/*, plus optional hand-written notes in
runs/mimic_notes.json:
  {"lead": "<html>", "a_notes": {"<tic>_<signal>": "text"},
   "rv_notes": {"<star>_<signal>": "text"}, "hd32564": "<html>",
   "comparison": [["his claim", "what we found"], ...], "a_summary": "<html>",
   "b_summary": "<html>", "c_summary": "<html>"}
Writes reports/mdwarf_rv_<date>.html (self-contained, figures embedded).
"""
from __future__ import annotations

import argparse
import base64
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

from planet_hunt import config
from planet_hunt.plots import fig_to_jpeg, rv_system_figure, stacked_transit_figure
from report import CSS, fmt

esc = html.escape
A = config.RUNS / "mdwarf"
B = config.RUNS / "rv"
THREAD = "https://x.com/mehran__jalali/status/2108339123669111033"


def repo_url() -> str:
    """Public URL of this repository (from the git remote)."""
    import subprocess
    try:
        u = subprocess.run(["git", "-C", str(config.ROOT), "config", "--get", "remote.origin.url"],
                           capture_output=True, text=True, timeout=10).stdout.strip()
    except Exception:
        u = ""
    if u.startswith("git@github.com:"):
        u = "https://github.com/" + u.split(":", 1)[1]
    return u.removesuffix(".git") or "https://github.com"


def img(jpeg: bytes, alt: str) -> str:
    return f'<img loading="lazy" alt="{esc(alt)}" src="data:image/jpeg;base64,{base64.b64encode(jpeg).decode()}">'


def flags(eb: str, caution: str) -> str:
    out = [f'<span class="flag eb">{esc(f)}</span>' for f in str(eb or "").split("; ") if f and f != "nan"]
    out += [f'<span class="flag caution">{esc(f)}</span>' for f in str(caution or "").split("; ") if f and f != "nan"]
    return " ".join(out) or '<span class="flag pass">no flags</span>'


def kv(pairs) -> str:
    return '<div class="kv">' + "".join(f"<div><span>{esc(k)}</span> {v}</div>" for k, v in pairs) + "</div>"


# --------------------------------------------------------------------------
# Method A
# --------------------------------------------------------------------------
CTRL_OFFSETS: list[float] = []   # centroid offsets of re-found known planets (calibration)


def a_card(r, vet: dict | None, note: str, i: int | None = None, fate: tuple | None = None) -> str:
    z = np.load(A / "lc" / f"{int(r.tic)}_{int(r.signal)}.npz")
    row = r._asdict()
    row["season_depths"] = json.loads(r.season_depths) if isinstance(r.season_depths, str) else []
    fig = stacked_transit_figure(z, row, f"TIC {int(r.tic)}" + (f" #{int(r.signal)}" if r.signal > 1 else ""))
    chip = f' <span class="chip {fate[1]}">{esc(fate[0])}</span>' if fate else ""
    parts = [f'<div class="card"><h3>{f"#{i} · " if i else ""}TIC {int(r.tic)}'
             f'{"" if r.signal == 1 else f" signal {int(r.signal)}"} · {esc(str(r.catalog_label))}{chip}</h3>',
             kv([("Period", f"{r.period:.5f} d"), ("Depth", f"{r.depth_ppm:,.0f} ppm"),
                 ("Duration", f"{r.duration_h:.2f} h"), ("SNR", fmt(r.snr)),
                 ("Transits", f"{int(r.n_transits)} in {int(r.n_sectors)} sectors"),
                 ("Radius (if real)", f"{fmt(r.rp_earth)} R⊕"),
                 ("Insolation", f"{fmt(r.insolation, 2)} × Earth"), ("Zone", esc(str(r.hz))),
                 ("Star", f"Tmag {fmt(r.Tmag)}, {fmt(r.teff, 0)} K, {fmt(r.rstar, 2)} R☉, {fmt(r.d, 0)} pc"),
                 ("2-min sectors", f"{int(r.n_spoc)}")]),
             f"<p>{flags(r.eb_flags, r.caution_flags)}</p>"]
    if vet:
        c, nb, sp, g = vet.get("centroid", {}), vet.get("neighbours", {}), vet.get("spoc", {}), vet.get("gaia", {})
        cal = (f" (re-found known planets on these crowded, faint stars: {min(CTRL_OFFSETS):.0f}-"
               f"{max(CTRL_OFFSETS):.0f}″, median {np.median(CTRL_OFFSETS):.0f}″)") if CTRL_OFFSETS else ""
        cline = (f"<b>Light-centroid test</b>: the in-transit shift implies a source ~{c['offset_arcsec']:.0f}″ "
                 f"from the light centre{cal}. On crowded faint stars this test cannot tell target from neighbour."
                 if c.get("ok") else f"<b>Light-centroid test</b>: {esc(c.get('verdict', '–'))}")
        lines = [cline,
                 f"<b>Depth vs crowding</b>: {esc(vet.get('crowding', {}).get('verdict', 'not run'))}",
                 f"<b>Neighbours</b> bright enough to fake {r.depth_ppm:,.0f} ppm within 52″: "
                 f"{nb.get('n_able', '?')} of {nb.get('n_within', '?')} TIC stars",
                 f"<b>SPOC</b>: {esc(sp.get('verdict', '–'))}",
                 f"<b>Difference image</b> (TESScut, 5 recent sectors): "
                 f"{esc(vet.get('diffimage', {}).get('verdict', 'not run'))}",
                 f"<b>Gaia RUWE</b> {fmt(g.get('ruwe'), 2)} (>1.4 hints at an unresolved companion); "
                 f"Gaia eclipsing-binary list: {len(g.get('gaia_eb', []))}; VSX variables within 30″: {len(g.get('vsx', []))}; "
                 f"published FFI candidate lists (DTARPS-S, Olmschenk+21): {len(g.get('dtarps', [])) + len(g.get('olmschenk', []))}"]
        parts.append('<div class="fu">' + "<br>".join(lines) + "</div>")
    if note:
        parts.append(f'<div class="note">{note}</div>')
    parts.append(img(fig_to_jpeg(fig), f"TIC {r.tic}") + "</div>")
    return "\n".join(parts)


def method_a(notes: dict) -> tuple[str, dict]:
    res = pd.read_csv(A / "results.csv")
    tg = pd.read_csv(A / "targets.csv")
    vets = {p.stem: json.loads(p.read_text()) for p in (A / "vet").glob("*.json")} if (A / "vet").exists() else {}
    ok = res[res["error"].isna()] if "error" in res else res
    n_stars = ok["tic"].nunique()
    tce = res[res["tce"]]
    known = tce[tce["catalog_status"].isin(["known_planet", "known_toi", "known_ctoi"])]
    fp = tce[tce["catalog_status"] == "known_fp"]
    unc = tce[tce["catalog_status"].isin(["new", "host_other"])].copy()
    unc["eb_flags"] = unc["eb_flags"].fillna("")
    clean = unc[unc["eb_flags"] == ""]
    surv = notes.get("a_survivors", [])
    for k, v in vets.items():
        m = res[(res["tic"] == v["tic"]) & (res["signal"] == v["signal"])]
        if len(m) and m["catalog_status"].iloc[0] in ("known_planet", "known_toi") and v.get("centroid", {}).get("ok"):
            CTRL_OFFSETS.append(float(v["centroid"]["offset_arcsec"]))

    # positive controls: TOIs on control stars with P <= 30 d
    from planet_hunt.catalogs import load_toi, period_match
    toi = load_toi()
    ctrl_rows = []
    for t in tg[tg["control"]].itertuples():
        for x in toi[(toi["tic"] == t.tic) & (toi["period"] <= 30)].itertuples():
            hit = tce[(tce["tic"] == t.tic) & tce["period"].apply(lambda p: period_match(p, x.period) is not None)]
            disp = x.tfop if isinstance(x.tfop, str) else x.tess_disp
            ctrl_rows.append({"toi": x.TOI, "tic": t.tic, "P": x.period, "depth": x.depth_ppm, "disp": disp,
                              "Tmag": t.Tmag, "found": len(hit) > 0,
                              "snr": float(hit["snr"].iloc[0]) if len(hit) else np.nan})
    ctrl = pd.DataFrame(ctrl_rows)
    real = ctrl[~ctrl["disp"].isin(["FP", "FA"])] if len(ctrl) else ctrl
    stats = {"n_stars": n_stars, "n_tce": len(tce), "n_known": len(known), "n_unc": len(unc),
             "n_clean": len(clean), "n_ctrl": len(real), "n_ctrl_ok": int(real["found"].sum()) if len(real) else 0,
             "med_sectors": float(tg["n_ffi"].median()), "n_never2min": int((tg["n_spoc"] == 0).sum()),
             "n_survive": len(surv)}

    ctrl_html = "".join(
        f"<tr><td>TOI-{c.toi}</td><td>{c.tic}</td><td class='num'>{fmt(c.Tmag)}</td><td class='num'>{fmt(c.P, 3)}</td>"
        f"<td class='num'>{fmt(c.depth, 0)}</td><td>{esc(str(c.disp))}</td>"
        f"<td>{'<span class=\"chip ok\">found</span>' if c.found else '<span class=\"chip bad\">missed</span>'}"
        f"</td><td class='num'>{fmt(c.snr)}</td></tr>" for c in ctrl.sort_values("P").itertuples()) if len(ctrl) else ""

    order = clean.sort_values("snr", ascending=False)
    best = order.head(int(notes.get("a_best", 6)))
    def fate_of(key):
        if key in surv:
            return ("survives every check we can run", "ok")
        f = notes.get("a_fates", {}).get(key)
        return (f, "bad") if f else None
    cards = "\n".join(a_card(r, vets.get(f"{r.tic}_{r.signal}"), notes.get("a_notes", {}).get(f"{r.tic}_{r.signal}", ""), i,
                             fate_of(f"{r.tic}_{r.signal}"))
                      for i, r in enumerate(best.itertuples(), 1))
    extra = [k for k in notes.get("a_extra_cards", [])]
    extra_cards = ""
    for key in extra:
        tic, sig = map(int, key.split("_"))
        m = res[(res["tic"] == tic) & (res["signal"] == sig)]
        if len(m):
            extra_cards += a_card(next(m.itertuples()), vets.get(key), notes.get("a_notes", {}).get(key, ""),
                                  fate=("eclipsing binary (method C)", "warn"))

    def fate(r):
        key = f"{r.tic}_{r.signal}"
        if r.eb_flags:
            return "rejected: " + r.eb_flags, "bad"
        if key in surv:
            return "survives the checks we can run", "ok"
        n = notes.get("a_fates", {}).get(key)
        return (n, "warn") if n else ("clean, see notes", "warn")
    rows = "".join(
        f"<tr><td>{r.tic}{'' if r.signal == 1 else f' #{r.signal}'}</td><td class='num'>{fmt(r.Tmag)}</td>"
        f"<td class='num'>{fmt(r.period, 4)}</td><td class='num'>{fmt(r.depth_ppm, 0)}</td>"
        f"<td class='num'>{fmt(r.duration_h, 2)}</td><td class='num'>{fmt(r.snr)}</td>"
        f"<td class='num'>{fmt(r.rp_earth)}</td><td class='num'>{fmt(r.insolation, 2)}</td>"
        f"<td class='wrap'><span class='chip {fate(r)[1]}'>{esc(fate(r)[0])}</span></td></tr>"
        for r in unc.sort_values(["eb_flags", "snr"], ascending=[True, False]).itertuples())

    # C: eclipsing binaries / stellar companions among the flagged signals
    f = tce["eb_flags"].fillna("")
    eb = tce[(tce["depth_ppm"] > 3000) & (tce["sin_ratio"] < 0.3) & ~f.str.contains("dominated")
             & f.str.contains("second eclipse|secondary|odd/even|too big|V-shaped")]
    eb = eb.sort_values("snr", ascending=False).drop_duplicates("tic")

    def eb_cat(r):
        g = vets.get(f"{r.tic}_{r.signal}", {}).get("gaia", {})
        hits = [f"VSX {v['Name']} ({v['Type']}, P={v['Period']} d)" for v in g.get("vsx", [])]
        hits += ["Gaia DR3 eclipsing binary" for _ in g.get("gaia_eb", [])]
        return "; ".join(hits) or ("none within 30″" if g else "not checked")
    eb_rows = "".join(
        f"<tr><td>{r.tic}</td><td class='num'>{fmt(r.Tmag)}</td><td class='num'>{fmt(r.period, 4)}</td>"
        f"<td class='num'>{fmt(r.depth_ppm / 1e4, 1)}%</td><td class='num'>{fmt(r.duration_h, 2)}</td>"
        f"<td class='wrap'>{esc(eb_cat(r))}</td><td class='wrap'>{esc(str(r.eb_flags))}</td>"
        f"<td class='wrap'>{esc(notes.get('c_notes', {}).get(f'{r.tic}_{r.signal}', ''))}</td></tr>"
        for r in eb.itertuples())

    html_a = f"""
<h2>A · Faint red dwarfs, many years of TESS images stacked</h2>
<p><b>What Mehran's thread did.</b> A red dwarf 1/6 the Sun's width, 3% deep transits lasting 1.4 h
every 11.7 days, found by averaging dozens of transits from TESS full-frame images, 2018 to 2025.</p>
<p><b>What we did.</b> We took {stats['n_stars']} red dwarfs (Teff below 3,900 K, under 0.6 Sun widths,
TESS magnitude 11 to 15, within 100 parsecs) within 12° of the two ecliptic poles, where TESS stares
for a year at a time. The median star has {stats['med_sectors']:.0f} sectors of full-frame-image light curves
(TESS-SPOC, or QLP where TESS-SPOC is missing). {stats['n_never2min']} of the stars were never 2-minute
targets, so NASA's SPOC planet search never ran on them. For each star we stitched every sector,
searched each observing season with BLS, summed the seasons, refined the best periods on all
years at once, then masked the signal and searched again.</p>
{notes.get('a_summary', '')}
<div class="tiles">{''.join(f'<div class="tile"><b>{v}</b><span>{esc(k)}</span></div>' for k, v in (
        ("red dwarfs searched", stats['n_stars']), ("signals over threshold", stats['n_tce']),
        ("known planets / TOIs re-found", stats['n_known']), ("known TOIs recovered (controls)", f"{stats['n_ctrl_ok']}/{stats['n_ctrl']}"),
        ("uncatalogued signals", stats['n_unc']), ("uncatalogued and clean", stats['n_clean']),
        ("survive every check", stats['n_survive'])))}</div>

<h3>Does the search work? The known planets on these stars</h3>
<p>{len(ctrl)} TOIs with periods under 30 days sit on our stars. "Found" means a signal at the TOI's period
(or a simple multiple). False alarms (FA) and false positives (FP) are listed but not counted.</p>
<div class="tablewrap"><table><tr><th>TOI</th><th>TIC</th><th class="num">Tmag</th><th class="num">P (d)</th>
<th class="num">Depth ppm</th><th>Status</th><th>Our search</th><th class="num">Our SNR</th></tr>{ctrl_html}</table></div>

<h3>The best uncatalogued signals</h3>
<p>Each figure: (a) the season-summed periodogram, (b) every transit folded together, (c) the depth in
each observing season (a planet's stays the same), (d) the phase-0.5 window where an eclipsing
binary's second eclipse would be.</p>
{cards or '<p>No uncatalogued signal passed the light-curve checks.</p>'}
{extra_cards}

<h3>Every uncatalogued signal and its fate</h3>
<details><summary>Show all {stats['n_unc']} uncatalogued signals</summary>
<div class="tablewrap"><table><tr><th>TIC</th><th class="num">Tmag</th><th class="num">P (d)</th>
<th class="num">Depth ppm</th><th class="num">Dur h</th><th class="num">SNR</th><th class="num">Rp R⊕</th>
<th class="num">S⊕</th><th>Fate</th></tr>{rows}</table></div></details>

<h3>C · By-product: eclipsing binaries and stellar companions</h3>
{notes.get('c_summary', '')}
<div class="tablewrap"><table><tr><th>TIC</th><th class="num">Tmag</th><th class="num">P (d)</th>
<th class="num">Depth</th><th class="num">Dur h</th><th>Variable-star catalogs</th><th>Why it is not a planet</th><th>Notes</th></tr>
{eb_rows or '<tr><td colspan=8>none</td></tr>'}</table></div>
"""
    return html_a, stats


# --------------------------------------------------------------------------
# Method B
# --------------------------------------------------------------------------
def method_b(notes: dict) -> tuple[str, dict]:
    res = pd.read_csv(B / "results_all.csv")
    meta = json.loads((B / "stars_all.json").read_text())
    from planet_hunt import rv
    known = rv.load_known_planets()
    ids = rv.simbad_ids([m["name"] for m in meta])
    rec = rv.recovery_table(res, meta, known, ids, rv.load_star_list())
    res["flags"] = res["flags"].fillna("")
    res["cautions"] = res["cautions"].fillna("") if "cautions" in res else ""
    n_stars = sum(1 for m in meta if not m.get("error"))
    unk = res[res["status"] != "known_planet"]
    clean = unk[unk["flags"] == ""]
    stats = {"n_stars": n_stars, "n_sig": len(res), "n_known_sig": int((res["status"] == "known_planet").sum()),
             "n_rec": int(rec["recovered"].sum()), "n_rec_all": len(rec), "n_unk": len(unk),
             "n_clean": len(clean)}
    rec_bins = "".join(
        f"<tr><td>{lab}</td><td class='num'>{len(s)}</td><td class='num'>{int(s['recovered'].sum())}</td></tr>"
        for lab, s in (("1-10 d", rec[rec["P"] < 10]), ("10-100 d", rec[(rec["P"] >= 10) & (rec["P"] < 100)]),
                       ("100-1000 d", rec[(rec["P"] >= 100) & (rec["P"] < 1000)]), ("over 1000 d", rec[rec["P"] >= 1000])))

    # HD 32564
    hd = res[res["star"] == "HD32564"]
    z = np.load(B / "star" / "HD32564.npz")
    labs = ["(b?)", "(c)", "(d)", "(e?)", "(not in thread)"]
    hd_fig = img(fig_to_jpeg(rv_system_figure(z, "HD 32564: HARPS RVBank, 193 nights 2009-2016", labs)), "HD 32564")
    his = [(2.6, 4.4), (11.1, 3.7), (23.6, 3.1), (49.1, 1.7)]
    hd_rows = ""
    for k, r in enumerate(hd.itertuples()):
        h = his[k] if k < len(his) else (None, None)
        hd_rows += (f"<tr><td>{r.signal}</td><td class='num'>{'–' if h[0] is None else h[0]}</td>"
                    f"<td class='num'>{'–' if h[1] is None else h[1]}</td><td class='num'>{r.period:.4f} ± {r.period_err:.4f}</td>"
                    f"<td class='num'>{r.K:.2f} ± {r.K_err:.2f}</td><td class='num'>{r.msini_earth:.1f}</td>"
                    f"<td class='num'>{r.fap:.0e}</td><td class='wrap'>{esc((r.flags + '; ' if r.flags else '') + str(r.cautions))}</td></tr>")

    # survivors
    show = clean.sort_values("fap")
    surv_rows = "".join(
        f"<tr><td>{esc(r.star)}</td><td>{esc(str(r.SpType))}</td><td class='num'>{r.n_nights}</td>"
        f"<td class='num'>{r.period:.3f}</td><td class='num'>{r.K:.2f}</td><td class='num'>{r.msini_earth:.1f}</td>"
        f"<td class='num'>{r.fap:.0e}</td><td>{esc(r.status)}</td>"
        f"<td class='wrap'>{notes.get('rv_notes', {}).get(f'{r.star}_{r.signal}', '') or esc(str(r.cautions))}</td></tr>"
        for r in show.itertuples())
    figs = ""
    for star in notes.get("rv_figures", []):
        p = B / "star" / f"{star}.npz"
        if p.exists():
            figs += f'<div class="card"><h3>{esc(star)}</h3>{notes.get("rv_fig_notes", {}).get(star, "")}' + \
                img(fig_to_jpeg(rv_system_figure(np.load(p), star)), star) + "</div>"

    html_b = f"""
<h2>B · Hidden planets in archival HARPS radial velocities</h2>
<p><b>What Mehran's thread did.</b> HD 32564, a Sun-like star 158 light-years away: 191 nights of
HARPS velocities from 2009 to 2016 hide planets at 2.6, 11.1 and 23.6 days (and a 49.1-day
candidate), each 9 to 13 Earth masses.</p>
<p><b>What we did.</b> The HARPS RVBank (Trifonov et al. 2020) re-reduced every public HARPS spectrum
up to 2018. For each of the {n_stars} stars with at least 50 nights we averaged each night,
treated data before and after the 2015 fibre upgrade as two instruments, and ran a generalised
Lomb-Scargle periodogram. If the strongest peak had a false-alarm probability below 0.1% we fitted a
Keplerian orbit there (also trying its 1-day and 1-year aliases), refitted everything together, and
looked again. Each signal was then checked against nine activity indicators (Ca II H&K from
Perdelwitz et al. 2024, H-alpha, sodium D, line width, chromatic index, FWHM, contrast, bisector),
the aliases of stronger signals, and whether it has the same strength in both halves of the data.</p>

<h3>Validation: HD 32564</h3>
{notes.get('hd32564', '')}
<div class="tablewrap"><table><tr><th>#</th><th class="num">Thread P (d)</th><th class="num">Thread K (m/s)</th>
<th class="num">Our P (d)</th><th class="num">Our K (m/s)</th><th class="num">m sin i (M⊕)</th>
<th class="num">FAP</th><th>Checks</th></tr>{hd_rows}</table></div>
{hd_fig}

<h3>The other {n_stars - 1} stars</h3>
{notes.get('b_summary', '')}
<div class="tiles">{''.join(f'<div class="tile"><b>{v}</b><span>{esc(k)}</span></div>' for k, v in (
        ("stars with 50+ nights", n_stars), ("significant signals", stats['n_sig']),
        ("matched a known planet", stats['n_known_sig']),
        ("known planets recovered (K ≥ 2 m/s)", f"{stats['n_rec']}/{stats['n_rec_all']}"),
        ("not a known planet", stats['n_unk']), ("…and pass the activity / alias / coherence checks", stats['n_clean'])))}</div>
<p>Calibration: of the {stats['n_rec_all']} confirmed planets on these stars with a published semi-amplitude
of at least 2 m/s, the blind search re-found {stats['n_rec']}:</p>
<div class="tablewrap"><table><tr><th>Period</th><th class="num">Known planets</th><th class="num">Re-found</th></tr>{rec_bins}</table></div>

<h3>Signals that are not known planets and pass the automatic checks</h3>
<p>"Not known" means no planet in the NASA Exoplanet Archive or exoplanet.eu at that period. Ambers in the
last column are cautions, not rejections.</p>
<div class="tablewrap"><table><tr><th>Star</th><th>Type</th><th class="num">Nights</th><th class="num">P (d)</th>
<th class="num">K (m/s)</th><th class="num">m sin i (M⊕)</th><th class="num">FAP</th><th>Catalog</th><th>Read</th></tr>{surv_rows}</table></div>
{figs}
"""
    return html_b, stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out")
    args = ap.parse_args()
    notes_path = config.RUNS / "mimic_notes.json"
    notes = json.loads(notes_path.read_text()) if notes_path.exists() else {}
    ha, sa = method_a(notes)
    hb, sb = method_b(notes)
    comp = "".join(f"<tr><td class='wrap'>{a}</td><td class='wrap'>{b}</td></tr>" for a, b in notes.get("comparison", []))
    title = "Planet Hunt: M Dwarfs and HARPS"
    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>{CSS}
.answer p,.answer li{{font-size:1.02rem}}</style></head><body><main>
<h1>{title}</h1>
<p class="sub">Two searches of public data modelled on <a href="{THREAD}">Mehran Jalali's thread</a>,
where Claude found a habitable-zone planet candidate around a red dwarf in TESS images and a hidden
three-planet system in archival HARPS velocities. Generated {pd.Timestamp.now():%Y-%m-%d}. Nothing has been
submitted anywhere.</p>
<div class="card tldr answer"><b>The answer.</b> {notes.get('lead', '').replace('{CTRL}', f"{sa['n_ctrl_ok']} of {sa['n_ctrl']}")}</div>
{ha}
{hb}
<h2>His claims next to ours</h2>
<div class="tablewrap"><table><tr><th>Thread</th><th>This run</th></tr>{comp}</table></div>
<h2>What it would take to believe one</h2>
{notes.get('believe', '')}
<footer>Code: <a href="{repo_url()}">{repo_url().replace("https://", "")}</a>
(<code>mdwarf_hunt.py</code>, <code>rv_hunt.py</code>, <code>mimic_report.py</code>). Data: NASA TESS via MAST
(TESS-SPOC and QLP high-level science products), HARPS RVBank (Trifonov et al. 2020, A&amp;A 636, A74),
Ca II H&amp;K (Perdelwitz et al. 2024, A&amp;A 683, A125), TIC v8.2, Gaia DR3, ExoFOP TOI/CTOI lists,
NASA Exoplanet Archive, exoplanet.eu, SIMBAD.</footer>
</main></body></html>"""
    out = Path(args.out) if args.out else config.REPORTS / f"mdwarf_rv_{pd.Timestamp.now():%Y-%m-%d}.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(doc)
    print(f"[report] {out} ({out.stat().st_size / 1e6:.1f} MB); A: {sa}; B: {sb}")


if __name__ == "__main__":
    main()
