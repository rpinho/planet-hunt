TRICERATOPS rerun for TOI-6000 c (2026-10-30), with the public Keck/NIRC2 contrast curves from ExoFOP.

Input: 33 SPOC 2-minute sectors, 2-minute phase bins, 120 s exposure, SPOC apertures from 12 sectors, Gaia background, default priors, no MOLUSC. Five seeded runs of 200,000 draws per configuration.

| Configuration | 14 neighbours kept: FPP, NFPP | Removed: FPP, NFPP |
|---|---|---|
| No imaging | 0.094 +/- 0.012, 0.048 | 0.045 +/- 0.005, 0 |
| Keck Hcont | 0.059 +/- 0.005, 0.059 | 0.0002 +/- 0.0001, 0 |
| Keck Hcont + Jcont | same as Hcont, same seeds | same as Hcont |

`tri6000.py` needs a TRICERATOPS install with list-valued contrast curves, the round-4 aperture and light-curve products, and `net.py` from the job directory. It records provenance; it is not a turnkey script.
