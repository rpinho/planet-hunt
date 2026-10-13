"""Round 7: assemble the ranked table (CSV + markdown) from revet.csv, targets.csv,
the novelty checks and the deep-vet / TRICERATOPS outputs.

    .venv/bin/python round7/table.py   -> runs/round7_eumemic_weak/round7_table.csv, round7_table.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from planet_hunt import config  # noqa: E402

OUT = config.RUNS / "round7_eumemic_weak"

# hand verdicts after diagnostics / pixels (see the report for the reasoning)
DEEP = {
    144439305: {"pixel": "on target (best source 7/8 sectors; neighbours dchi2 >= 25 worse)",
                "verdict": "CANDIDATE (strong)"},
    119876352: {"pixel": "not run (rejected on light curve)",
                "verdict": "REJECTED: dip sits in a +-3,000 ppm ~1.7 h oscillation; depth varies 1.7-8.4 ppt by sector"},
    166422847: {"pixel": "not run (known neighbour signal)",
                "verdict": "KNOWN: TOI-3282.01 (APC) 33\" away, same 2.94783 d period"},
    270188035: {"pixel": "OFF target: Gaia T=15.8 star at 41.5\" (target dchi2 3,653)",
                "verdict": "REJECTED: blended EB on a neighbour"},
    231706739: {"pixel": "unresolved: target vs T=19.0 star at 22\" (dchi2 56 vs 60)",
                "verdict": "WEAK CANDIDATE (pixel-unresolved, SNR 8)"},
}


def tri(tic):
    p = OUT / "tri" / f"triceratops_{tic}.json"
    if not p.exists():
        return ""
    d = json.load(open(p))
    s = f"{d['FPP_mean']:.3f} / {d['NFPP_mean']:.4f}"
    q = OUT / "tri" / f"triceratops_{tic}_molusc.json"
    if q.exists():
        e = json.load(open(q))
        s += f" (MOLUSC+cleared: {e['FPP_mean']:.3f} / {e['NFPP_mean']:.4f})"
    return s


def main():
    r = pd.read_csv(OUT / "revet.csv")
    tg = pd.read_csv(OUT / "targets.csv")[["tic", "period", "depth", "tmag", "teff", "rstar", "d", "ra", "dec", "n_sectors", "E_hold"]]
    nov = pd.read_csv(OUT / "novelty_live_toi_ctoi.csv")
    r = r.merge(tg, left_on=["tic", "P_his"], right_on=["tic", "period"], how="left")
    r = r.merge(nov, on=["tic", "period"], how="left")
    r["depth_eum_ppm"] = r.depth * 1e6
    r["pixel"] = r.tic.map(lambda t: DEEP.get(t, {}).get("pixel", "-"))
    r["fpp_nfpp"] = r.tic.map(tri)
    def verdict(x):
        if x.tic in DEEP:
            return DEEP[x.tic]["verdict"]
        return {"REFUTED": "not reproduced (hold-out inconsistent > 3 sigma)",
                "INCONCLUSIVE": "inconclusive (hold-out too weak or flags)",
                "SURVIVOR": "survivor"}[x.verdict]
    r["final"] = r.apply(verdict, axis=1)
    order = {"CANDIDATE (strong)": 0, "WEAK CANDIDATE (pixel-unresolved, SNR 8)": 1}
    r["rank_key"] = r.final.map(lambda v: order.get(v, 2 if v.startswith(("KNOWN", "REJECTED")) else
                                         (3 if v.startswith("inconclusive") else 4)))
    r = r.sort_values(["rank_key", "z_hold"], ascending=[True, False]).reset_index(drop=True)
    r.index += 1
    r.to_csv(OUT / "round7_table.csv", index_label="rank")
    lines = ["| # | TIC | P (d) | depth ppm (his / ours) | SNR his / ours-his-era / ours-all | z hold-out (expected if real) | Rp (Re) | transits | pixel | FPP / NFPP | TOI/CTOI live | verdict |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, x in r.iterrows():
        lines.append(f"| {i} | {x.tic} | {x.P_full:.5f} | {x.depth_eum_ppm:.0f} / {x.depth_full_ppm:.0f} | "
                     f"{x.snr_eum:.1f} / {x.snr_his_ours:.1f} / {x.snr_full:.1f} | {x.z_hold:.1f} ({x.z_hold_expected:.1f}) | "
                     f"{x.rp_full:.2f} | {int(x.n_tr_full)} | {x.pixel} | {x.fpp_nfpp or '-'} | {x.toi_ctoi_live} | {x.final} |")
    (OUT / "round7_table.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:14]))
    print(r.final.value_counts().to_string())


if __name__ == "__main__":
    main()
