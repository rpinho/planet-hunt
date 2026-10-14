"""Step 4: build ONE self-contained HTML report for a hunt.

    .venv/bin/python report.py --sector 90
    .venv/bin/python report.py --sector 90 --top 20 --out my_report.html

Reads runs/validation/* and runs/sector_NNNN/*, redraws the diagnostic
figures from the saved arrays, embeds them as JPEG data URIs, and writes
reports/planet_hunt_sector_NNNN.html (phone-friendly, dark-mode aware).

Optional hand-written notes: runs/sector_NNNN/notes.json, a dict
{"<tic>_<signal>": "your read of this signal"} -- shown on the cards.
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
from planet_hunt.catalogs import toi_coverage
from planet_hunt.data import CleanLC
from planet_hunt.plots import candidate_figure, fig_to_jpeg

STATUS = {  # catalog status -> (chip text, css class)
    "known_planet": ("Known planet", "ok"),
    "known_toi": ("Known TOI", "toi"),
    "known_fp": ("Known false positive", "bad"),
    "known_ctoi": ("Known community TOI", "toi"),
    "host_other": ("Known host, other period", "warn"),
    "new": ("NEW: not in catalogs", "new"),
}

esc = html.escape


def img_tag(jpeg: bytes, alt: str) -> str:
    b64 = base64.b64encode(jpeg).decode()
    return f'<img loading="lazy" alt="{esc(alt)}" src="data:image/jpeg;base64,{b64}">'


def load_signal(run: Path, row) -> tuple[CleanLC, dict]:
    z = np.load(run / "lc" / f"{int(row.tic)}_{int(row.signal)}.npz")
    lc = CleanLC(tic=int(row.tic), sector=0, time=z["time"], flux=z["flux"],
                 flux_err=np.zeros_like(z["flux"]), trend=z["trend"], flat=z["flat"])
    return lc, {"periods": z["periods"], "power": z["power"]}


def fmt(v, nd=1, dash="–"):
    try:
        if v is None or not np.isfinite(float(v)):
            return dash
    except (TypeError, ValueError):
        return esc(str(v))
    return f"{float(v):,.{nd}f}"


def flags_html(eb, caution) -> str:
    out = []
    for f in str(eb or "").split("; "):
        if f and f != "nan":
            out.append(f'<span class="flag eb" title="likely not a planet">✕ {esc(f)}</span>')
    for f in str(caution or "").split("; "):
        if f and f != "nan":
            out.append(f'<span class="flag caution" title="caution">! {esc(f)}</span>')
    return " ".join(out) or '<span class="flag pass">✓ passes basic checks</span>'


def status_chip(status: str, label: str) -> str:
    text, cls = STATUS.get(status, (status, "warn"))
    return f'<span class="chip {cls}" title="{esc(label)}">{esc(text)}</span>'


def followup_html(fu: dict) -> str:
    """Multi-sector re-check box (from followup.py)."""
    v = fu["verdict"]
    cls = "ok" if v.startswith(("PERSISTS", "REAPPEARS")) else (
        "bad" if v.startswith(("FADES", "NOT SEEN")) else "warn")
    rows = "".join(
        f"<tr><td>S{o['sector']}</td><td class='num'>{fmt(o.get('period'), 4)}</td>"
        f"<td class='num'>{fmt(o.get('depth_ppm'), 0)}</td><td class='num'>{fmt(o.get('snr'))}</td>"
        f"<td class='num'>{fmt(o.get('expected_snr'), 0)}</td></tr>"
        for o in fu["other_sectors"])
    st = fu.get("stitched")
    if st:
        rows += (f"<tr><td><b>stitched {st['n_sectors']}</b></td><td class='num'>{fmt(st['period'], 5)}</td>"
                 f"<td class='num'>{fmt(st['depth_ppm'], 0)}</td><td class='num'><b>{fmt(st['snr'])}</b></td>"
                 f"<td class='num'>{fmt(st['expected_snr'], 0)}</td></tr>")
    return (f'<div class="fu {cls}"><b>Other-sector check:</b> {esc(v)}'
            f'<details><summary>per-sector numbers</summary><div class="tablewrap"><table>'
            f'<tr><th>Sector</th><th class="num">best P (±1.5%)</th><th class="num">Depth ppm</th>'
            f'<th class="num">SNR</th><th class="num">SNR if real</th></tr>{rows}</table></div>'
            f'<p><small>Per-sector rows: the best peak within ±1.5% of the period, so noise always finds '
            f'<i>something</i> at SNR ~4-5; a real planet sits at the same period every time. '
            f'“SNR if real” is the SNR the candidate depth would have there, with that sector’s noise (lower of the white- and red-noise estimates, as for detection).</small></p>'
            f'</details></div>')


def vet_html(vt: dict) -> str:
    """'Is the dip on this star?' box (from vet.py)."""
    c, nb = vt.get("centroid", {}), vt.get("neighbours", {})
    v = c.get("verdict", "not measured")
    cls = "bad" if v.startswith("SHIFT") else ("warn" if v.startswith(("MARGINAL", "inconclusive", "not"))
                                               else "ok")
    if "n_able" in nb:
        names = ", ".join(f"TIC {a['tic']} (T={a['Tmag']:.1f}, {a['sep_arcsec']:.0f}\")"
                          for a in nb["able"][:4])
        nbt = (f"{nb['n_able']} of {nb['n_within']} TIC neighbours within {nb['radius_arcsec']:.0f}\" are bright "
               f"enough (Δmag &lt; {nb['dmag_limit']:.1f}) to produce this depth if fully eclipsed"
               + (f": {esc(names)}." if names else ".") + f" Target share of the light: {100 * nb['dilution']:.0f}%.")
        if nb["n_able"] == 0:
            cls = "ok" if cls == "ok" or v.startswith("inconclusive") else cls
    else:
        nbt = "neighbour query failed."
    tv = vt.get("transits", {}).get("verdict")
    if tv and tv.startswith("INCONSISTENT"):
        cls = "bad"
    trt = f'<br><b>Transit by transit:</b> {esc(tv)}.' if tv else ""
    sp = vt.get("spoc", {}).get("verdict")
    spt = f'<br><b>NASA pipeline:</b> {esc(sp)}.' if sp else ""
    return (f'<div class="fu {cls}"><b>On this star?</b> Centroid: {esc(v)}.<br>'
            f'<small>Neighbours: {nbt}</small>{trt}{spt}</div>')


def auto_read(r) -> str:
    """A plain-language, deliberately sceptical one-paragraph read."""
    bits = []
    rp = r.rp_earth
    if np.isfinite(rp):
        size = ("Earth/super-Earth" if rp < 1.8 else "sub-Neptune" if rp < 4
                else "Neptune-to-Saturn" if rp < 8 else "Jupiter-size" if rp < 22
                else "stellar (too big for a planet)")
        bits.append(f"Depth {r.depth_ppm:,.0f} ppm on a {fmt(r.star_rad, 2)} R☉ star implies "
                    f"Rp ≈ {rp:.1f} R⊕ ({size}), if the dip is on this star and undiluted.")
    bits.append(f"{int(r.n_transits)} transits with data; SNR {r.snr:.1f} "
                f"(white-noise {r.snr_white:.1f}, red-noise {r.snr_red:.1f}); SDE {r.sde:.1f}.")
    if np.isfinite(r.dur_ratio):
        bits.append(f"Duration is {r.dur_ratio:.2f}x the maximum for a central transit of this star.")
    status = r.catalog_status
    if status == "new":
        bits.append("Not in the TOI list or confirmed-planet table. Default assumption: noise, "
                    "a systematic, stellar variability or a blended eclipsing binary until proven otherwise.")
    elif status == "host_other":
        bits.append("The star already has catalogued signal(s) at another period: this peak is often "
                    "a residual or alias of the known signal, or a systematic.")
    elif status == "known_fp":
        bits.append("TFOPWG has already ruled this a false positive, even though it passes our simple "
                    "checks: a reminder that pixel-level and follow-up vetting catch what light-curve stats cannot.")
    else:
        bits.append("Re-discovery of a catalogued signal: the pipeline is doing its job.")
    return " ".join(bits)


def signal_card(run: Path, r, notes: dict, followups: dict, vets: dict | None = None,
                extra_title: str = "") -> str:
    """One signal card: numbers, flags, auto read, follow-up/vet boxes, figure."""
    lc, bls = load_signal(run, r)
    st = r.to_dict()
    title = (f"TIC {r.tic}" + (f" signal {r.signal}" if r.signal > 1 else "")
             + f"  P={r.period:.4f} d  depth={r.depth_ppm:,.0f} ppm  [{STATUS[r.catalog_status][0]}]")
    jpeg = fig_to_jpeg(candidate_figure(lc, bls, st, title))
    note = notes.get(f"{r.tic}_{r.signal}")
    fu = followups.get(f"{r.tic}_{r.signal}")
    vt = (vets or {}).get(f"{r.tic}_{r.signal}")
    links = (f'<a href="https://exofop.ipac.caltech.edu/tess/target.php?id={r.tic}">ExoFOP</a> · '
             f'<a href="https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery=TIC%20{r.tic}">MAST</a>')
    return f"""
<div class="card" id="tic{r.tic}_{r.signal}">
  <h3>{extra_title}TIC {r.tic}{f' · signal {r.signal}' if r.signal > 1 else ''}{' · positive control' if r.control else ''}
  {status_chip(r.catalog_status, r.catalog_label)}</h3>
  <div>{flags_html(r.eb_flags, r.caution_flags)}</div>
  <div class="kv">
    <div><span>Period</span> {r.period:.4f} d</div><div><span>Depth</span> {r.depth_ppm:,.0f} ± {r.depth_err_ppm:,.0f} ppm</div>
    <div><span>Duration</span> {r.duration_h:.2f} h</div><div><span>SNR / SDE</span> {r.snr:.1f} / {r.sde:.1f}</div>
    <div><span>Tmag</span> {fmt(r.Tmag, 2)}</div><div><span>Star</span> {fmt(r.Teff, 0)} K, {fmt(r.star_rad, 2)} R☉</div>
    <div><span>Odd/even</span> {fmt(r.oe_sigma)}σ</div><div><span>Secondary</span> {fmt(r.depth_sec_ppm, 0)} ppm ({fmt(r.sec_sigma)}σ)</div>
  </div>
  <p>{'' if r.catalog_status == 'new' else esc(r.catalog_label) + '. '}{esc(auto_read(r))}</p>
  {followup_html(fu) if fu else ''}
  {vet_html(vt) if vt else ''}
  {f'<div class="note"><b>Analyst note:</b> {esc(note)}</div>' if note else ''}
  <p><small>{links}</small></p>
  {img_tag(jpeg, f'diagnostics for TIC {r.tic}')}
</div>"""


CSS = """
:root{color-scheme:light;--bg:#f6f5f1;--card:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#7a7972;
--line:#e3e2dc;--accent:#2a78d6;--ok:#1a7f37;--okbg:#e5f4ea;--toi:#2a78d6;--toibg:#e6f0fb;--bad:#c4312f;
--badbg:#fbe9e8;--warn:#8a5a00;--warnbg:#fdf2dc;--new:#6250d6;--newbg:#eeebfb;--code:#efeee9}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--bg:#121211;--card:#1a1a19;
--ink:#f0efec;--ink2:#c3c2b7;--muted:#9a9990;--line:#2e2e2b;--accent:#3987e5;--ok:#5cc27a;--okbg:#17301f;
--toi:#6aa8f0;--toibg:#16263a;--bad:#ef7b79;--badbg:#3a1b1a;--warn:#e3b04f;--warnbg:#33280f;--new:#a99ff0;
--newbg:#25213d;--code:#24241f}}
:root[data-theme="dark"]{color-scheme:dark;--bg:#121211;--card:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--muted:#9a9990;
--line:#2e2e2b;--accent:#3987e5;--ok:#5cc27a;--okbg:#17301f;--toi:#6aa8f0;--toibg:#16263a;--bad:#ef7b79;
--badbg:#3a1b1a;--warn:#e3b04f;--warnbg:#33280f;--new:#a99ff0;--newbg:#25213d;--code:#24241f}
*{box-sizing:border-box}html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Inter,Roboto,sans-serif}
main{max-width:980px;margin:0 auto;padding:20px 16px 64px}
h1{font-size:1.65rem;line-height:1.2;margin:.2em 0 .1em}h2{font-size:1.2rem;margin:2.2em 0 .6em}
h3{font-size:1rem;margin:0 0 .4em}p{margin:.5em 0}.sub{color:var(--ink2);margin:0 0 1em}
a{color:var(--accent)}code{background:var(--code);padding:1px 5px;border-radius:4px;font-size:.88em}
pre{background:var(--code);padding:10px 12px;border-radius:8px;overflow-x:auto;font-size:.85em}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:16px;margin:14px 0}
.tldr{border-left:4px solid var(--accent)}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:10px;margin:12px 0}
.tile{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 12px}
.tile b{display:block;font-size:1.6rem;line-height:1.1}.tile span{color:var(--ink2);font-size:.82rem}
.tablewrap{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:var(--card)}
table{border-collapse:collapse;width:100%;font-size:.82rem}
th,td{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top;white-space:nowrap}
td.wrap{white-space:normal;min-width:180px}th{position:sticky;top:0;background:var(--card);color:var(--ink2);font-weight:600}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
.chip{display:inline-block;padding:1px 8px;border-radius:99px;font-size:.78rem;font-weight:600;white-space:nowrap}
.chip.ok{color:var(--ok);background:var(--okbg)}.chip.toi{color:var(--toi);background:var(--toibg)}
.chip.bad{color:var(--bad);background:var(--badbg)}.chip.warn{color:var(--warn);background:var(--warnbg)}
.chip.new{color:var(--new);background:var(--newbg)}
.flag{display:inline-block;font-size:.76rem;padding:0 6px;border-radius:5px;margin:1px 0;white-space:nowrap}
.flag.eb{color:var(--bad);background:var(--badbg)}.flag.caution{color:var(--warn);background:var(--warnbg)}
.flag.pass{color:var(--ok);background:var(--okbg)}
img{width:100%;height:auto;border-radius:8px;border:1px solid var(--line);margin-top:8px;display:block}
.kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:4px 14px;font-size:.88rem;margin:6px 0}
.kv div span{color:var(--muted)}.note{background:var(--warnbg);border-radius:8px;padding:8px 10px;margin-top:8px;font-size:.9rem}
.fu{border-radius:8px;padding:8px 10px;margin-top:8px;font-size:.9rem;border:1px solid var(--line)}
.fu.ok{background:var(--okbg)}.fu.bad{background:var(--badbg)}.fu.warn{background:var(--warnbg)}
.verdict{font-size:1.15rem;font-weight:700}.verdict.ok{color:var(--ok)}.verdict.bad{color:var(--bad)}
details summary{cursor:pointer;color:var(--accent);font-weight:600}
ul,ol{padding-left:1.25em}li{margin:.25em 0}footer{color:var(--muted);font-size:.8rem;margin-top:3em}
"""


def validation_section() -> str:
    parts = []
    for res_path in sorted((config.RUNS / "validation").glob("*/result.json")):
        v = json.loads(res_path.read_text())
        f, p, c = v["found"], v["published"], v["checks"]
        ok = v["recovered"]
        jpeg = (res_path.parent / "figure.jpg").read_bytes()
        rows = "".join(
            f"<tr><td>{lab}</td><td class='num'>{fmt(f[k], nd)}</td><td class='num'>{fmt(p.get(k), nd)}</td>"
            f"<td class='num'>{fmt(100 * c[ck], 1)}%</td></tr>"
            for lab, k, ck, nd in (("Period (d)", "period", "period_rel_err", 4),
                                   ("Depth (ppm)", "depth_ppm", "depth_rel_err", 0),
                                   ("Duration (h)", "duration_h", "duration_rel_err", 2)))
        parts.append(f"""
<div class="card">
  <div class="verdict {'ok' if ok else 'bad'}">{'✓ Recovered' if ok else '✕ Not recovered'}: {esc(v['planet'])}
  (TIC {v['tic']}), TESS sector {v['sector']}</div>
  <p>Blind search: the pipeline was not told the period. BLS's strongest peak is compared with the
  published values ({esc(p.get('source', ''))}). SNR {fmt(f['snr'])}, SDE {fmt(f['sde'])},
  {f['n_transits']} transits, implied radius {fmt(f['rp_earth'], 2)} R⊕ (published {fmt(p.get('rp_earth'), 2)} R⊕).</p>
  <div class="tablewrap"><table><tr><th></th><th class="num">Found</th><th class="num">Published</th>
  <th class="num">Rel. error</th></tr>{rows}</table></div>
  {img_tag(jpeg, 'validation figure')}
</div>""")
    return "\n".join(parts) or "<p>No validation run found: run <code>validate.py</code> first.</p>"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sector", type=int, default=90)
    ap.add_argument("--top", type=int, default=15, help="number of candidate cards")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    run = config.RUNS / f"sector_{args.sector:04d}"
    res = pd.read_csv(run / "results.csv")
    summ = json.loads((run / "summary.json").read_text())
    notes = json.loads((run / "notes.json").read_text()) if (run / "notes.json").exists() else {}
    followups = {p.stem: json.loads(p.read_text()) for p in (run / "followup").glob("*.json")}
    vets = {p.stem: json.loads(p.read_text()) for p in (run / "vet").glob("*.json")}
    res["eb_flags"] = res["eb_flags"].fillna("")
    res["caution_flags"] = res["caution_flags"].fillna("")
    tce = res[res["tce"]].copy()
    clean = tce[tce["eb_flags"] == ""]
    n_new = int((tce["catalog_status"] == "new").sum())
    n_new_clean = int((clean["catalog_status"] == "new").sum())
    n_known = int(tce["catalog_status"].isin(["known_planet", "known_toi", "known_ctoi"]).sum())
    cand_keys = [f"{r.tic}_{r.signal}" for r in clean.itertuples()
                 if r.catalog_status in ("new", "host_other")]
    fu_done = [k for k in cand_keys if k in followups]
    fu_survive = [k for k in fu_done if followups[k]["verdict"].startswith(("PERSISTS", "REAPPEARS"))]
    n_host_clean = int((clean["catalog_status"] == "host_other").sum())

    # --- positive controls --------------------------------------------------
    ctrl_rows = []
    for tic, g in res[res["control"] == True].groupby("tic", sort=False):  # noqa: E712
        g = g.sort_values("signal")
        hit = g[g["tce"] & g["catalog_status"].isin(["known_planet", "known_toi", "known_fp", "known_ctoi"])]
        r = (hit if len(hit) else g).iloc[0]
        ok = len(hit) > 0
        ctrl_rows.append(
            f"<tr><td>{tic}</td><td class='num'>{fmt(r.Tmag, 2)}</td>"
            f"<td>{'✓ recovered' if ok else '✕ missed'}</td><td class='num'>{fmt(r.period, 3)}</td>"
            f"<td class='wrap'>{esc(r.catalog_label)}</td></tr>")
    n_ctrl_ok = summ["controls_recovered"]

    # --- summary table ------------------------------------------------------
    trs = []
    for r in tce.itertuples():
        trs.append(
            f"<tr><td>{r.tic}{'' if r.signal == 1 else f' #{r.signal}'}{' ⋆' if r.control else ''}</td>"
            f"<td class='num'>{fmt(r.Tmag, 2)}</td><td class='num'>{fmt(r.period, 4)}</td>"
            f"<td class='num'>{fmt(r.depth_ppm, 0)}</td><td class='num'>{fmt(r.duration_h, 2)}</td>"
            f"<td class='num'>{fmt(r.snr)}</td><td class='num'>{fmt(r.rp_earth)}</td>"
            f"<td class='wrap'>{flags_html(r.eb_flags, r.caution_flags)}</td>"
            f"<td class='wrap'>{status_chip(r.catalog_status, r.catalog_label)}<br>"
            f"<small>{esc(r.catalog_label)}</small></td></tr>")

    # --- candidate cards ----------------------------------------------------
    order = {"new": 0, "host_other": 1, "known_ctoi": 2, "known_toi": 3, "known_fp": 4,
             "known_planet": 5}
    cards_df = clean.assign(_o=clean["catalog_status"].map(order)).sort_values(
        ["_o", "snr"], ascending=[True, False]).head(args.top)
    rejects = []
    for key in ("persists", "gap edges", "deep", "systematic", "secondary", "odd/even"):
        m = tce[tce["eb_flags"].str.contains(key, regex=False)
                & ~tce.index.isin([r.name for r in rejects])]
        if len(m) and len(rejects) < 3:
            rejects.append(m.iloc[0])

    cards_html = "\n".join(signal_card(run, r, notes, followups, vets) for _, r in cards_df.iterrows())
    rejects_html = "\n".join(signal_card(run, r, notes, followups, vets) for r in rejects)

    s = summ
    n_toi = toi_coverage().get(args.sector, 0)
    if n_toi >= 30:
        context = (f"That is the expected answer: sector {args.sector} has already been combed by the TESS team "
                   f"({n_toi} TOIs list it), and an uncatalogued BLS hit is almost always noise, a systematic, "
                   "stellar variability or a blended eclipsing binary.")
    else:
        context = (f"Sector {args.sector} is fresh: only {n_toi} TOIs list it so far, so “not in the catalogs” "
                   "partly means “not looked at yet”. A survivor here still needs the checks in section 6.")
    html_doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Planet Hunt: TESS Sector {args.sector}</title>
<style>{CSS}</style></head><body><main>
<h1>Planet Hunt: TESS Sector {args.sector}</h1>
<p class="sub">A blind transit search of {s['n_searched']} bright dwarf stars with BLS, vetted and cross-matched
against the ExoFOP TOI list and the NASA Exoplanet Archive. Generated {s['date']}.</p>

<div class="card tldr"><b>Bottom line.</b> The pipeline works: it blindly recovered the known planet
below and {n_ctrl_ok} of {s['n_controls']} known-TOI positive controls. Across {s['n_searched']} stars it raised
{s['n_tce']} threshold-crossing signals (SNR ≥ 7.1, ≥ 2 transits): {n_known} are catalogued planets/TOIs,
{n_new} are not in any catalog, and only {n_new_clean} of those (plus {n_host_clean} extra signal(s) on known
hosts) pass the automated eclipsing-binary / systematics checks. Re-checked against each star's other TESS
sectors, <b>{len(fu_survive)} of {len(fu_done)}</b> survive{'s' if len(fu_survive) == 1 else ''}. {context}</div>

<h2>1 · Validation: can it find a known planet?</h2>
{validation_section()}

<h2>2 · The hunt in numbers</h2>
<div class="tiles">
 <div class="tile"><b>{s['n_searched']}</b><span>stars searched (Tmag &lt; {s['settings']['tmag']})</span></div>
 <div class="tile"><b>{s['n_tce']}</b><span>signals over threshold</span></div>
 <div class="tile"><b>{n_known}</b><span>known planets / TOIs</span></div>
 <div class="tile"><b>{len(tce) - len(clean)}</b><span>flagged EB / systematic / variable</span></div>
 <div class="tile"><b>{n_new_clean}</b><span>uncatalogued &amp; pass basic checks</span></div>
</div>
<p>Data: {s['download_mb']:.0f} MB of SPOC 2-minute PDCSAP light curves from MAST; up to
{s['settings'].get('max_signals', 1)} signals searched per star (find, mask, search again). Runtime
{s['runtime_min']} min on the Mac Mini.</p>

<h3>Positive controls: {n_ctrl_ok} / {s['n_controls']} recovered</h3>
<p>{s['n_controls']} of the stars were chosen <i>because</i> they host a TOI that one sector should reveal
(period 0.5–12.5 d, TOI SNR ≥ 10). Misses are informative: shallow signals that need several sectors,
or a stronger systematic/variability winning the BLS peak.</p>
<div class="tablewrap"><table><tr><th>TIC</th><th class="num">Tmag</th><th>Result</th>
<th class="num">BLS P (d)</th><th>Catalog</th></tr>{''.join(ctrl_rows)}</table></div>

<h2>3 · Every signal over threshold</h2>
<p>⋆ = positive control. “#2” = second signal on the same star. Flags in red mean “probably not a
planet”; amber means “be suspicious”. Rp assumes the dip is on the target star and undiluted.</p>
<div class="tablewrap"><table>
<tr><th>TIC</th><th class="num">Tmag</th><th class="num">P (d)</th><th class="num">Depth ppm</th>
<th class="num">Dur h</th><th class="num">SNR</th><th class="num">Rp R⊕</th><th>Vetting</th><th>Catalog</th></tr>
{''.join(trs)}</table></div>

<h2>4 · Top {len(cards_df)} signals that pass the basic checks</h2>
<p>Uncatalogued first, then signals on known hosts, then re-discoveries. Each figure: (a) the raw light
curve and the trend removed, (b) the BLS periodogram, (c) the phase-folded transit, (d) odd vs even
transits, (e) the phase-0.5 window where an eclipsing binary's secondary eclipse would sit.</p>
{cards_html}

<h2>5 · What rejects look like</h2>
<p>The most common impostors, for calibration.</p>
{rejects_html}

<h2>6 · How to check a candidate yourself</h2>
<p>Step 5 below is automated in <code>followup.py</code> and its verdicts are on the cards above. It was
positive-controlled on two known TOIs (TOI-2810.01 and TOI-733 b), which both <i>persist</i> when sectors are stitched.</p>
<ol>
<li><b>Look at every transit individually.</b> Plot each event (panel a, zoomed). Real transits look alike;
one dominant event plus noise is a red flag. Check none sits on a data-gap edge or a momentum dump
(the quality flags in the FITS file).</li>
<li><b>Odd/even and secondary.</b> Different odd/even depths or any eclipse at phase 0.5 points to an
eclipsing binary at twice the period.</li>
<li><b>Shape and duration.</b> Planet transits are U-shaped with short ingress/egress; the duration must be
possible for the star's density. V-shapes and over-long durations suggest a grazing EB or a giant star.</li>
<li><b>Is it on this star?</b> TESS pixels are 21″ wide, so a faint eclipsing binary nearby can mimic a
shallow transit. Do a centroid / difference-image test on the target pixel file
(<code>lk.search_targetpixelfile</code>), or use <a href="https://github.com/stevengiacalone/triceratops">TRICERATOPS</a>
to compute the false-positive probability from neighbours.</li>
<li><b>Does it repeat in other sectors?</b> A real planet shows up at the same period and depth every
time TESS looks. <code>followup.py --tic N</code> checks each other sector and stitches nearby ones; a real
signal's SNR grows roughly as √(number of sectors), noise does not.</li>
<li><b>Check ExoFOP</b> for the TIC: existing TOIs, CTOIs, imaging and spectroscopy notes, and the
SPOC/QLP data-validation reports, which run these tests (and more) professionally.</li>
</ol>

<h2>7 · What to do next</h2>
<ul>
<li><b>Planet Hunters TESS</b> (<a href="https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess">Zooniverse</a>)
is the realistic route for a citizen scientist (check the project page for its current status): volunteers flag
events by eye (especially single transits, which BLS misses) and the science team vets and publishes them;
volunteers have been credited on planet discoveries.</li>
<li><b>ExoFOP CTOI submission</b>: a Community TOI can be submitted by a registered ExoFOP user once a signal has
survived the checks above, ideally including multi-sector recovery and a centroid test. Submitting unvetted BLS
hits wastes follow-up observers' time, so the bar is high.</li>
<li><b>Before anyone should believe a candidate</b>: multi-sector confirmation, pixel-level centroid analysis,
a statistical-validation false-positive probability (e.g. TRICERATOPS &lt; 1.5%), high-resolution imaging to rule
out close companions, and ultimately radial velocities (mass) or ground-based transit photometry.</li>
<li><b>Ways to extend this kit</b>: hunt the newest sector (signals not yet alerted as TOIs), fainter stars
(<code>--tmag 12</code>), M dwarfs (small planets give deeper transits), Transit Least Squares
(installed: <code>transitleastsquares</code>) for better sensitivity to small planets, multi-sector stitching
for longer periods, and a single-transit search for long-period planets.</li>
</ul>

<details><summary>Method details</summary>
<p><b>Targets.</b> SPOC 2-minute targets in sector {args.sector} with Tmag &lt; {s['settings']['tmag']},
log g &gt; 4.1, 3800 K &lt; Teff &lt; 6500 K, TIC contamination ratio &lt; 0.1: a seeded random draw plus
{s['n_controls']} TOI-host positive controls.
<b>Cleaning.</b> PDCSAP flux, default quality mask, 0.25 d trimmed after each gap &gt; 0.5 d, 1-day
Savitzky-Golay flatten with iterative sigma masking, upward 4σ clip; after detection the star is re-flattened
with its transits masked so the filter does not eat the depth.
<b>Search.</b> astropy BoxLeastSquares on 10-min bins, 0.5 d to half the baseline, frequency-uniform grid,
durations 0.7–8 h limited to ≤ 2× a central transit across a Sun-like star at each period.
<b>Threshold.</b> SNR ≥ 7.1 using the lower of white-noise and red-noise (transit-timescale binned scatter)
estimates, ≥ 2 transits with data, and periodogram SDE ≥ 7. Each star is searched again after masking its first
signal.
<b>Vetting flags.</b> depth &gt; 3% or implied Rp &gt; 2 RJup; odd/even depth difference &gt; 3σ; secondary at
phase 0.5 &gt; 3σ; V-shape; signal persisting after its transits are masked (variability); transit times shared
with ≥ 2 other stars (systematic); half or more of the transits within 0.25 d (+ half a duration) of a data gap; duration implausible for the star; period near TESS's 13.7-day orbit;
modest periodogram peak (SDE &lt; 9).
<b>Cross-match.</b> by TIC ID and period (including 2×, ½×, 3×, ⅓×, ¼× aliases) against the ExoFOP TOI list,
the ExoFOP community-TOI (CTOI) list and the NASA Exoplanet Archive pscomppars table.
<b>Follow-up.</b> followup.py: narrow BLS (±1.5% of P) in each other SPOC 2-min sector, then stitch sectors
within ~1 year, re-flatten with the transits masked and fold; verdict PERSISTS if stitched SNR ≥ 10 and ≥ 60% of
the SNR a real signal would have, FADES if below that while a real signal would have had SNR ≥ 12.</p>
<p><b>Re-run</b> (on the Mac Mini, in <code>.</code>):</p>
<pre>.venv/bin/python validate.py
.venv/bin/python hunt.py --sector {args.sector} --n 150
.venv/bin/python followup.py --sector {args.sector}
.venv/bin/python report.py --sector {args.sector}</pre>
</details>
<footer>Data: NASA TESS mission via MAST (SPOC light curves), ExoFOP-TESS TOI catalog, NASA Exoplanet Archive.
Code: . (planet_hunt package). Nothing in this report has been vetted beyond the automated checks described.</footer>
</main></body></html>"""

    out = Path(args.out) if args.out else config.REPORTS / f"planet_hunt_sector_{args.sector:04d}.html"
    out.write_text(html_doc, encoding="utf-8")
    print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB, {len(cards_df)} cards, {len(rejects)} rejects)")


if __name__ == "__main__":
    main()
