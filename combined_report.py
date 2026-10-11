"""Step 4b: ONE report across several sectors, leading with the answer.

    .venv/bin/python combined_report.py --sectors 104 105 106 107
    .venv/bin/python combined_report.py --sectors 104 105 106 107 --out newest.html \
        --lead runs/combined_lead.html

Reads runs/sector_NNNN/{results.csv, summary.json, followup/, vet/, notes.json}
for every sector (run hunt.py, followup.py and vet.py on each first) and writes
one self-contained HTML page:
  * the answer first: stars, signals, how many survive, the best candidates
  * why these sectors (TOI-catalog coverage per sector)
  * the funnel per sector and every uncatalogued signal with its fate
  * candidate cards with figures, multi-sector and "on this star?" checks
  * what it would take to believe one, and where to send it.
`--lead` replaces the auto-generated opening paragraph with your own HTML.
"""
from __future__ import annotations

import argparse
import datetime
import json
from pathlib import Path

import pandas as pd

from planet_hunt import config
from planet_hunt.catalogs import toi_coverage
from planet_hunt.search import followup_targets
from report import CSS, esc, fmt, signal_card

KNOWN = ["known_planet", "known_toi", "known_ctoi"]


def mjd_to_date(mjd: float) -> str:
    d = datetime.date(1858, 11, 17) + datetime.timedelta(days=float(mjd))
    return f"{d:%b} {d.day}, {d.year}"


def fate(fu: dict | None, vt: dict | None, caution: str = "") -> tuple[str, str, int]:
    """(short label, css class, rank) for one uncatalogued clean signal.
    Lower rank = better candidate."""
    cv = ((vt or {}).get("centroid") or {}).get("verdict", "")
    v = (fu or {}).get("verdict", "not checked")
    if cv.startswith("SHIFT"):
        return "rejected: centroid moves (neighbour)", "bad", 9
    tr = ((vt or {}).get("transits") or {})
    if tr.get("ok") and not tr.get("consistent", True):
        return "rejected: individual transits disagree", "bad", 7
    if v.startswith(("FADES", "NOT SEEN")):
        return "rejected: fades in other sectors", "bad", 8
    marginal = cv.startswith("MARGINAL")
    spoc = ((vt or {}).get("spoc") or {})
    if spoc.get("match"):   # NASA's pipeline raised a TCE at this period
        return ("survives: SPOC TCE at the same period" + (" (centroid marginal)" if marginal else ""),
                "ok", 1 + marginal)
    if v.startswith(("PERSISTS", "REAPPEARS")):
        return ("survives: seen in other sectors" + (" (centroid marginal)" if marginal else ""),
                "ok" if not marginal else "warn", 1 + marginal)
    # Soft problems: not disqualifying, but each one pushes a signal down the list.
    soft = []
    if marginal:
        soft.append("centroid marginal")
    if tr.get("ok") and tr.get("p", 1) < 1e-3:
        soft.append("one dip dominates")
    caution = str(caution or "")
    if "very short" in caution:
        soft.append("very short for the star")
    if "TESS orbit" in caution:
        soft.append("period near TESS orbit")
    if spoc.get("verdict", "").startswith("no SPOC"):
        soft.append("SPOC did not flag it")
    rank = 3 + sum(x != "SPOC did not flag it" for x in soft)
    where = "only one sector exists" if v.startswith("only one") else "other sectors too noisy"
    return (f"open: {where}" + (f" ({'; '.join(soft)})" if soft else ""), "warn", min(rank, 6))


def load_sector(sector: int) -> dict:
    run = config.RUNS / f"sector_{sector:04d}"
    res = pd.read_csv(run / "results.csv")
    res["eb_flags"] = res["eb_flags"].fillna("")
    res["caution_flags"] = res["caution_flags"].fillna("")
    res["sector"] = sector
    obs = pd.read_csv(config.TABLE_DIR / f"sector_{sector:04d}_spoc120.csv",
                      usecols=["t_min", "t_max"])
    rd = lambda d: {p.stem: json.loads(p.read_text()) for p in (run / d).glob("*.json")}  # noqa: E731
    notes = json.loads((run / "notes.json").read_text()) if (run / "notes.json").exists() else {}
    return {"sector": sector, "run": run, "res": res,
            "summ": json.loads((run / "summary.json").read_text()),
            "followups": rd("followup"), "vets": rd("vet"), "notes": notes,
            "start": mjd_to_date(obs["t_min"].min()), "end": mjd_to_date(obs["t_max"].max())}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sectors", type=int, nargs="+", required=True)
    ap.add_argument("--best", type=int, default=3, help="candidate cards up top")
    ap.add_argument("--lead", help="HTML file with your own opening paragraph")
    ap.add_argument("--title", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    secs = [load_sector(s) for s in sorted(args.sectors)]
    cov = toi_coverage()
    span = f"{secs[0]['sector']}–{secs[-1]['sector']}" if len(secs) > 1 else str(secs[0]["sector"])
    title = args.title or f"Planet Hunt: TESS Sectors {span}"

    # ---- pool every sector ------------------------------------------------
    n_stars = sum(d["summ"]["n_searched"] for d in secs)
    n_ctrl = sum(d["summ"]["n_controls"] for d in secs)
    n_ctrl_ok = sum(d["summ"]["controls_recovered"] for d in secs)
    allres = pd.concat([d["res"] for d in secs], ignore_index=True)
    tce = allres[allres["tce"]]
    n_tce = len(tce)
    n_known = int(tce["catalog_status"].isin(KNOWN).sum())
    n_fp = int((tce["catalog_status"] == "known_fp").sum())
    unc = tce[tce["catalog_status"].isin(["new", "host_other"])]
    unc_clean = unc[unc["eb_flags"] == ""]

    cands = []
    for d in secs:
        for _, r in followup_targets(d["res"], borderline=True).iterrows():
            key = f"{r.tic}_{r.signal}"
            fu, vt = d["followups"].get(key), d["vets"].get(key)
            edge = r.eb_flags != ""
            if edge and fu is None:      # borderline signal nobody re-checked
                continue
            label, cls, rank = fate(fu, vt, r.caution_flags)
            if edge:
                label = f"second look ({r.eb_flags}): " + label
            cands.append({"d": d, "r": r, "key": key, "fu": fu, "vt": vt, "edge": edge,
                          "label": label, "cls": cls, "rank": rank})
    cands.sort(key=lambda c: (c["rank"], -c["r"].snr))
    survivors = [c for c in cands if c["rank"] < 7]
    n_second = sum(c["edge"] for c in cands)
    multi = [c for c in survivors if c["rank"] <= 2]
    single = [c for c in survivors if c["rank"] > 2]
    best = cands[:args.best]

    # ---- the answer -------------------------------------------------------
    if args.lead:
        lead = Path(args.lead).read_text()
    else:
        def one(c):
            r = c["r"]
            return (f"TIC {r.tic} (sector {c['d']['sector']}, P = {r.period:.2f} d, "
                    f"{r.depth_ppm:,.0f} ppm, ~{fmt(r.rp_earth)} R⊕): {c['label']}")
        lead = (f"<p><b>{n_stars} stars searched in {len(secs)} sectors. {n_tce} signals over threshold. "
                f"{len(unc_clean)} are uncatalogued and pass the basic checks, plus {n_second} borderline ones given a second look; "
                f"{len(survivors)} of those {len(cands)} survive "
                f"every test we can run ({len(multi)} confirmed in other sectors, {len(single)} not testable yet).</b></p>"
                + ("<ul>" + "".join(f"<li>{esc(one(c))}</li>" for c in best) + "</ul>" if best else ""))

    tiles = "".join(f'<div class="tile"><b>{v}</b><span>{esc(k)}</span></div>' for k, v in (
        ("stars searched", n_stars), ("signals over threshold", n_tce),
        ("known planets / TOIs re-found", n_known), ("uncatalogued signals re-checked", len(cands)),
        ("survive every test run", len(survivors)), ("positive controls recovered", f"{n_ctrl_ok}/{n_ctrl}")))

    # ---- why these sectors -----------------------------------------------
    cov_rows = "".join(
        f"<tr><td>S{d['sector']}</td><td>{d['start']} – {d['end']}</td>"
        f"<td class='num'>{cov.get(d['sector'], 0)}</td><td class='num'>{d['summ']['n_searched']}</td>"
        f"<td class='num'>{d['summ']['n_tce']}</td><td class='num'>"
        f"{int(((d['res']['tce']) & d['res']['catalog_status'].isin(KNOWN)).sum())}</td>"
        f"<td class='num'>{int((d['res']['tce'] & (d['res']['eb_flags'] == '') & d['res']['catalog_status'].isin(['new', 'host_other'])).sum())}</td>"
        f"<td class='num'>{sum(1 for c in survivors if c['d'] is d)}</td>"
        f"<td class='num'>{d['summ']['controls_recovered']}/{d['summ']['n_controls']}</td></tr>"
        for d in secs)
    ref = {s: cov.get(s, 0) for s in (85, 90, 95)}
    ref_txt = ", ".join(f"sector {s}: {n}" for s, n in ref.items())

    # ---- every uncatalogued signal ---------------------------------------
    def unc_row(sector, r, fate_label, cls):
        return (f"<tr><td>S{sector}</td><td>{r.tic}{'' if r.signal == 1 else f' #{r.signal}'}</td>"
                f"<td class='num'>{fmt(r.Tmag, 1)}</td><td class='num'>{fmt(r.period, 3)}</td>"
                f"<td class='num'>{fmt(r.depth_ppm, 0)}</td><td class='num'>{fmt(r.snr)}</td>"
                f"<td class='num'>{fmt(r.rp_earth)}</td><td class='wrap'><span class='chip {cls}'>"
                f"{esc(fate_label)}</span></td></tr>")
    rows = [unc_row(c["d"]["sector"], c["r"], c["label"], c["cls"]) for c in cands]
    seen = {(c["d"]["sector"], c["r"].tic, c["r"].signal) for c in cands}
    flagged = unc[(unc["eb_flags"] != "") & ~pd.Series(
        [(r.sector, r.tic, r.signal) in seen for r in unc.itertuples()], index=unc.index)]
    flagged = flagged.sort_values("snr", ascending=False)
    rows += [unc_row(r.sector, r, "rejected by light-curve checks: " + r.eb_flags, "bad")
             for r in flagged.itertuples()]

    # ---- cards --------------------------------------------------------------
    def card(c, i=None):
        d = c["d"]
        head = f"#{i} · " if i else ""
        return (f'<div class="fate {c["cls"]}">{head}Sector {d["sector"]} · {esc(c["label"])}</div>'
                + signal_card(d["run"], c["r"], d["notes"], d["followups"], d["vets"],
                              extra_title=f"S{d['sector']} · "))
    best_html = "\n".join(card(c, i) for i, c in enumerate(best, 1)) or \
        "<p>No uncatalogued signal passed the basic checks.</p>"
    rest_html = "\n".join(card(c) for c in cands[args.best:])

    css = CSS + """
.fate{margin:22px 0 -8px;font-weight:700;font-size:.92rem}.fate.ok{color:var(--ok)}
.fate.bad{color:var(--bad)}.fate.warn{color:var(--warn)}
.answer p,.answer li{font-size:1.02rem}"""
    html_doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(title)}</title>
<style>{css}</style></head><body><main>
<h1>{esc(title)}</h1>
<p class="sub">A blind transit search of the newest TESS data ({secs[0]['start']} – {secs[-1]['end']}),
vetted, cross-matched against the TOI, community-TOI and confirmed-planet catalogs, re-checked in each
star's other sectors and tested for "is the dip on this star?". Generated {pd.Timestamp.now():%Y-%m-%d}.</p>

<div class="card tldr answer"><b>The answer.</b> {lead}</div>
<div class="tiles">{tiles}</div>

<h2>1 · Why these sectors</h2>
<p>The TESS team turns signals into TOIs months after the data arrive. For mature sectors the TOI list
is thick ({ref_txt}), so an uncatalogued signal there is almost always junk. These sectors are fresh:
few or no TOIs yet. That is the only place a real, uncatalogued planet is plausible. It also means
the “new” label is cheap here: the TESS team simply has not looked yet.</p>
<div class="tablewrap"><table><tr><th>Sector</th><th>Observed</th><th class="num">TOIs so far</th>
<th class="num">Stars</th><th class="num">Signals</th><th class="num">Known re-found</th>
<th class="num">New &amp; clean</th><th class="num">Survive</th><th class="num">Controls</th></tr>{cov_rows}</table></div>

<h2>2 · The best {len(best)} candidates</h2>
<p>Ranked: seen again in other sectors first, then not-yet-testable, then rejected (shown only if
fewer than {args.best} survive). Each figure: (a) raw light curve and trend, (b) periodogram,
(c) folded transit, (d) odd vs even transits, (e) the phase-0.5 window where an eclipsing binary's
second eclipse would sit.</p>
{best_html}

<h2>3 · Every uncatalogued signal and its fate</h2>
<p>Rp assumes the dip is on the target and undiluted. “Rejected by light-curve checks” = red flags
from the automated vetting (odd/even, secondary eclipse, variability, gap edges, shared timing).</p>
<div class="tablewrap"><table><tr><th>Sector</th><th>TIC</th><th class="num">Tmag</th>
<th class="num">P (d)</th><th class="num">Depth ppm</th><th class="num">SNR</th><th class="num">Rp R⊕</th>
<th>Fate</th></tr>{''.join(rows)}</table></div>
<p>Also over threshold: {n_known} re-discoveries of known planets/TOIs and {n_fp} known false positive(s).
Positive controls (stars chosen because they host a findable TOI): {n_ctrl_ok} of {n_ctrl} recovered.</p>

{f'<h2>4 · The other uncatalogued signals that passed the basic checks</h2>{rest_html}' if rest_html else ''}

<h2>{'5' if rest_html else '4'} · What it would take to believe one</h2>
<ol>
<li><b>Other sectors.</b> A real planet comes back at the same period and depth. <code>followup.py</code>
does this; for the newest sectors many stars have only older sectors years away, or none.</li>
<li><b>On this star?</b> TESS pixels are 21″ wide. <code>vet.py</code> runs a flux-weighted centroid test
and lists neighbours bright enough to fake the depth. It is coarse: calibrated on sector 90 it flagged
TOI-409.01 (a known centroid false positive) only as “marginal”, and one of four known planets also came
out “marginal”. The real test is a difference image from the target pixel file
(<code>lk.search_targetpixelfile</code>), e.g. with
<a href="https://github.com/stevengiacalone/triceratops">TRICERATOPS</a> for a false-positive probability.</li>
<li><b>The official pipeline.</b> Check the SPOC data-validation report and TCE list for the sector on MAST,
and ExoFOP for the TIC, before anything else. They often already contain the signal.</li>
<li><b>Then</b> high-resolution imaging and radial velocities. None of that is a citizen-scientist job.</li>
</ol>
<p><b>Where a survivor would go:</b> post the TIC and ephemeris on
<a href="https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess">Planet Hunters TESS</a>
Talk, or, after the checks above, an <a href="https://exofop.ipac.caltech.edu/tess/">ExoFOP</a> CTOI.
Nothing has been submitted anywhere.</p>

<details><summary>Method</summary>
<p>Per sector: 150 SPOC 2-minute targets (Tmag &lt; 11, FGK dwarfs, low contamination), 15 of them
TOI-host positive controls. PDCSAP flux, 1-day Savitzky-Golay flatten, BLS from 0.5 d to half the
baseline, two signals per star, SNR ≥ 7.1 (lower of white and red noise), ≥ 2 transits, SDE ≥ 7.
Light-curve red flags: odd/even, secondary eclipse, V-shape, sinusoidal variability, persistence after
masking, gap-edge transits, transit times shared with other stars. Cross-match by TIC and period
(with aliases) to the ExoFOP TOI and CTOI lists and the NASA Exoplanet Archive. Follow-up: narrow BLS in
each other sector, plus stitching of sectors within about a year. Code:
the planet-hunt starter kit (<code>combined_report.py</code>).</p>
<pre>for s in {' '.join(str(d['sector']) for d in secs)}; do
  .venv/bin/python hunt.py --sector $s --n 150
  .venv/bin/python followup.py --sector $s
  .venv/bin/python vet.py --sector $s
done
.venv/bin/python combined_report.py --sectors {' '.join(str(d['sector']) for d in secs)}</pre>
</details>
<footer>Data: NASA TESS mission via MAST (SPOC light curves), ExoFOP-TESS, NASA Exoplanet Archive.
Nothing here has been vetted beyond the automated checks described.</footer>
</main></body></html>"""
    tag = "_".join(str(d["sector"]) for d in secs)
    out = Path(args.out) if args.out else config.REPORTS / f"planet_hunt_sectors_{tag}.html"
    out.write_text(html_doc, encoding="utf-8")
    print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB; {len(cands)} uncatalogued re-checked, "
          f"{len(survivors)} survive, best {len(best)})")
    for c in cands:
        r = c["r"]
        print(f"  S{c['d']['sector']} TIC {r.tic} #{r.signal} P={r.period:.4f} d={r.depth_ppm:.0f} "
              f"SNR={r.snr:.1f} Rp={r.rp_earth:.1f}: {c['label']}")


if __name__ == "__main__":
    main()
