# planet-hunt: a TESS transit-search starter kit

A small, readable pipeline that searches public NASA TESS light curves for
transiting planets, vets what it finds, checks it against the catalogs of known
objects and writes one self-contained HTML report.

It is built to be read and extended, not to compete with the TESS team's
pipelines (SPOC, QLP). Every step is a plain function in `planet_hunt/`.

## Setup

```bash
git clone https://github.com/rpinho/planet-hunt && cd planet-hunt
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements.txt
```

All downloads and library caches go to `./cache` (`planet_hunt/config.py`
redirects astropy, astroquery and lightkurve there). Delete `./cache` to free
the space. One sector of 150 stars is about 300 MB.

## The steps

```bash
.venv/bin/python validate.py                  # 1. recover a known planet (pi Men c, sector 1)
.venv/bin/python hunt.py --sector 90 --n 150  # 2+3. search a sector, vet, cross-match
.venv/bin/python followup.py --sector 90      # 3b. re-check survivors in the star's other sectors
.venv/bin/python vet.py --sector 90           # 3c. on this star? transits agree? did SPOC see it?
.venv/bin/python report.py --sector 90        # 4. reports/planet_hunt_sector_0090.html
.venv/bin/python combined_report.py --sectors 104 105 106 107   # 4b. one report, many sectors
```

| Script | What it does | Output |
|---|---|---|
| `validate.py` | Blind BLS search of pi Mensae c (TOI-144.01). Compares period, depth and duration with the NASA Exoplanet Archive. `--tic/--planet/--sector` to try another planet. | `runs/validation/<planet>/` |
| `hunt.py` | Picks bright FGK dwarfs (Tmag < 11) with SPOC 2-min data in the sector, plus 15 known-TOI hosts as positive controls. Downloads the PDCSAP light curves, detrends them, runs BLS (twice per star: find, mask, search again), measures SNR, applies the vetting flags and labels every signal against the catalogs. | `runs/sector_NNNN/results.csv`, `summary.json`, `lc/*.npz` |
| `followup.py` | For clean uncatalogued signals, checks the star's other TESS sectors and stitches the nearby ones. A real planet's SNR grows about as the square root of the number of sectors; noise does not. | `runs/sector_NNNN/followup/*.json` |
| `vet.py` | For the same signals: flux-weighted centroid shift in transit and TIC neighbours bright enough to fake the depth (`planet_hunt/centroid.py`), a transit-by-transit depth consistency test, and a lookup of SPOC's own threshold-crossing events in the public DV reports (`planet_hunt/spoc.py`). `--borderline` also re-checks signals whose only red flag is one transit at a gap edge. | `runs/sector_NNNN/vet/*.json` |
| `combined_report.py` | One report across several sectors that leads with the answer: stars, signals, survivors, the best candidates, then every uncatalogued signal and its fate. `--lead` for your own opening paragraph. | `reports/planet_hunt_sectors_*.html` |
| `report.py` | Builds one HTML report (validation, numbers, controls, the full signal table, candidate cards with figures, rejects, next steps). Optional `runs/sector_NNNN/notes.json` adds your own notes to the cards. | `reports/planet_hunt_sector_NNNN.html` |

## Re-running on another sector

```bash
.venv/bin/python hunt.py --sector 92 --n 150          # any sector with SPOC 2-min data
.venv/bin/python followup.py --sector 92
.venv/bin/python report.py --sector 92
```

The first MAST query for a new sector takes about 4 minutes; after that it is
cached. Useful options: `--tmag 12` (fainter stars), `--n 300`,
`--controls 0` (no positive controls), `--max-signals 3`, `--seed 7` (a
different random draw).

**Which sector to choose.** Sector 90 (March 2025) was picked because the TOI
catalog is mature there: if something real were in it, the TESS team would
very likely have alerted it already, so the "known vs new" label means
something. The newest sectors (about 100 to 107 as of October 2026) have
barely been through the TOI process, so a real planet could still be
uncatalogued there. The "new" bucket will also be much larger and noisier.
Check catalog coverage with `cache/catalogs/toi.csv` (`Sectors` column).

## How the search and vetting work

- **Cleaning** (`data.py`): PDCSAP flux with the default quality mask, normalised. 0.25 d is trimmed after every gap longer than 0.5 d (thermal ramps). A 1-day Savitzky-Golay flatten removes variability, and outliers more than 4 sigma *above* the trend are clipped. Once a signal is found, the star is re-flattened with its transits masked so the filter does not eat the depth.
- **Search** (`search.py`): astropy `BoxLeastSquares` on 10-min bins, with periods from 0.5 d to half the baseline on a frequency-uniform grid. At each period only physically plausible durations are tried (up to 2x a central transit across a Sun-like star).
- **Threshold**: SNR >= 7.1, taking the lower of the white-noise and red-noise estimates, plus >= 2 transits and periodogram SDE >= 7.
- **Red flags** (probably not a planet): depth > 3%; implied radius > 2 RJup; odd/even depth difference > 3 sigma; secondary eclipse at phase 0.5; V-shape; mostly sinusoidal signal (spots, pulsation, contact binary); the same period coming back after masking; transits at data-gap edges; transit times shared with other stars (spacecraft systematic).
- **Amber flags**: duration implausible for the star, period near TESS's 13.7-d orbit, a modest periodogram peak.
- **Catalogs** (`catalogs.py`): ExoFOP TOI list, ExoFOP CTOI list and NASA Exoplanet Archive `pscomppars`. Matching uses the TIC ID plus the period, including 2x, 1/2, 3x, 1/3 and 1/4 aliases. All are public and cached for 3 days.

## Results so far

| Run | Stars | Signals | Known re-found | Uncatalogued, fully checked | Survive |
|---|---|---|---|---|---|
| Sector 90 (mature TOI list) | 150 | 33 | 10 | 4 | 0 |
| Sectors 104-107 (May-Sep 2026, no TOIs yet) | 600 | 135 | 63 (48/60 controls) | 18 | 4, all at SNR 7-8 |

In the newest sectors the best survivor (TIC 11528365, 8.17 d, 814 ppm) is a
watch-list item, not a candidate: it sits at the detection threshold, SPOC did
not raise a TCE for it, and the star's earlier sector does not clearly show it.
The lesson: "no TOI yet" does not mean "unexamined". SPOC's search and DV
reports for new sectors are public within weeks and catch nearly every real
signal on bright stars, so the better hunting ground is fainter stars, M
dwarfs and single transits. Also: the gap-edge rule is harsh (it rejected the
known planet TOI-262 b), hence `--borderline`.

The checks are calibrated on known objects: the transit-by-transit test flags
none of 71 re-found known planets/TOIs (except one caught at half its period);
the centroid test calls TOI-409.01, a known centroid false positive, only
"marginal", so treat it as coarse.

## Be honest about "NEW"

A BLS peak that is in no catalog is almost always noise, a systematic, stellar
variability or a blended eclipsing binary. On the first run (sector 90, 150
stars), the 4 signals that passed every automated check all faded in other
sectors. Before anyone should believe a candidate, it needs:
multi-sector recovery, a pixel-level centroid test (is the dip on this star?),
a statistical false-positive probability (e.g. TRICERATOPS), high-resolution
imaging and, ultimately, radial velocities. The real routes for a citizen
scientist are Planet Hunters TESS (Zooniverse) and, for a well-vetted
signal, an ExoFOP CTOI submission.

## Ideas for extending it (with Claude)

- Use Transit Least Squares (`transitleastsquares` is installed) for better small-planet sensitivity.
- Stitch several sectors before searching, to reach longer periods.
- Add a single-transit search for long-period planets.
- Add centroid / difference-image vetting from target pixel files (`lk.search_targetpixelfile`).
- Search M dwarfs (`select_targets` filters): small planets give deeper transits there.
- Run on the newest sector, ahead of the TOI alerts.
