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

# Method A: faint red dwarfs in full-frame images (see below)
.venv/bin/python mdwarf_hunt.py select --n 300 && .venv/bin/python mdwarf_hunt.py fetch
.venv/bin/python mdwarf_hunt.py search && .venv/bin/python mdwarf_hunt.py vet
# Method B: archival HARPS radial velocities (see below)
.venv/bin/python rv_hunt.py --validate        # HD 32564 first
.venv/bin/python rv_hunt.py --min-nights 50   # every RVBank star with >= 50 nights
.venv/bin/python mimic_report.py              # one report for A and B
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
| Method A: 300 faint M dwarfs near the ecliptic poles, ~39 FFI sectors each | 300 | 270 | 28 (28/29 TOI controls) | 4 | 1 good (TIC 167723381, 2.09 d, 1.7 R⊕ if on target; caveats) + 1 marginal |
| Method B: HARPS RVBank, 333 stars with >= 50 nights | 333 | 301 | 185 signals (128/190 planets with K >= 2 m/s) | 42 pass automatic checks | HD 32564 reproduced; a few to look at |

**Methods A and B (modelled on Mehran Jalali's thread).** Method B reproduces his
HD 32564 system blind: 2.56 d (or its 1-day alias 1.63 d), 11.06, 23.61 and 49.2 d,
K = 4.35 / 3.48 / 2.84 / 1.71 m/s, 9-12 Earth masses, and no catalog lists it. Elsewhere
in the RVBank, most uncatalogued signals are stellar activity; HD 71334 (89 d),
MCC 150 (622 d), HD 207129 (1,398 d) and HD 125881 (41 d) are worth a closer look.
Method A did not find a habitable-zone planet. Of 240 uncatalogued transit signals
on 300 faint M dwarfs, 4 passed the light-curve checks. One is a blend (its depth tracks
TESS-SPOC's crowding correction, a test every known planet passes); one is doubtful (SPOC's
2-minute search of the star saw nothing). TIC 167723381 (2.09 d) survives with caveats and
TIC 298780765 (1.17 d) is marginal. Calibration lesson: on these crowded faint stars the
light-centroid shift is useless as an on-target test (known planets show 2-60"). By-products:
uncatalogued M-dwarf eclipsing binaries such as TIC 55452500 (17.7 d, eclipses up to 16% deep).

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

## Method A: faint red dwarfs in full-frame images

Idea credited to [Mehran Jalali's public thread](https://x.com/mehran__jalali/status/2108339123669111033),
where Claude found a habitable-zone planet candidate around a red dwarf by averaging dozens of TESS
transits from the full-frame images (FFIs), 2018 to 2025.

Why it is a gap: SPOC's planet search only runs on 2-minute targets, and most nearby M dwarfs
fainter than Tmag ~12 never were. Their light curves exist anyway, extracted from the FFIs by
TESS-SPOC and by MIT's Quick-Look Pipeline (QLP); QLP's own search mostly covers Tmag < 13.5. A
small planet gives a deep transit on a small star, and near the ecliptic poles TESS stares for a
year at a time, so a star there has 30-70 sectors to stack.

| Step | What it does |
|---|---|
| `mdwarf_hunt.py select` | TIC dwarfs with Teff 2700-3900 K, R 0.08-0.6 R☉, Tmag 11-15, d < 100 pc, within 12° of either ecliptic pole (`planet_hunt/ffi.py`). Lists every SPOC, TESS-SPOC and QLP light curve per star (one MAST query per 150 stars). Picks stars with no 2-minute data at all, plus TOI hosts as positive controls. |
| `mdwarf_hunt.py fetch` | Downloads one FFI light curve per sector (TESS-SPOC first, QLP otherwise), keeps time/flux/centroid as a small `.npz` and deletes the FITS. Stops if `./cache` passes 8 GB. |
| `mdwarf_hunt.py search` | Per star (`planet_hunt/multisector.py`): running-median detrending per sector with flare clipping; **semi-coherent BLS**: one BLS per observing season on a shared, star-aware period grid (0.5-30 d), log-likelihood powers summed over seasons; the top peaks are then refined **coherently on all years**. Statistics on the unbinned multi-cadence data with per-sector weights and red-noise inflation: SNR, odd/even, secondary eclipse, V-shape, sinusoid test, per-transit and per-season depth consistency, radius, insolation and habitable-zone status (Kopparapu et al. 2014). Masks the signal and searches again (up to 4 signals). Cross-match with TOI, CTOI and confirmed planets; the same period on several stars is flagged as a systematic. |
| `mdwarf_hunt.py vet` | For clean uncatalogued signals (and known planets, for calibration): centroid shift in transit combined over sectors, TIC neighbours bright enough to fake the depth, SPOC TCEs for any 2-minute data, Gaia DR3 RUWE, Gaia eclipsing-binary list and VSX. |

A full search of one star (stitch 40 sectors, ~10^5 trial periods, 4 rounds) takes 20-60 s.
300 stars x ~37 sectors is about 11,000 light curves and 2.5 GB after shrinking.

## Method B: hidden planets in archival HARPS radial velocities

Also from [the same thread](https://x.com/mehran__jalali/status/2108339123669111033): three planets of
9-13 Earth masses around HD 32564 in 191 nights of public HARPS data (2009-2016).

Data: the HARPS RVBank (Trifonov et al. 2020, A&A 636, A74; CDS J/A+A/636/A74), SERVAL velocities
of every public HARPS spectrum up to 2018 corrected for nightly zero points, with activity indicators;
plus per-spectrum Ca II H&K R'HK from Perdelwitz et al. (2024, A&A 683, A125), fetched through VizieR
TAP. `rv_hunt.py` (`planet_hunt/rv.py`), for every star with at least 50 nights:

1. nightly weighted means; data before and after the 2015 fibre upgrade are two instruments
   (separate zero points and jitter);
2. generalised Lomb-Scargle of the residuals; if the best peak's Baluev false-alarm probability
   is below 0.1%, add a Keplerian there. The peak and its 1-day and 1-year aliases are each fitted
   and the best likelihood wins (the losers are reported: HD 32564's innermost signal is
   2.56 d or 1.63 d, the data cannot tell). Everything (offsets, a linear trend, all orbits, jitter)
   is refitted jointly; eccentric orbits are kept when the BIC improves by more than 10;
3. repeat until nothing significant is left (max 6 signals);
4. vet each signal: a significant peak at the same period in any of nine activity indicators, or a
   correlation with one, rejects it; so do 1-day / 1-year / lunar aliases of a stronger signal,
   periods at 1 yr or the lunar month, periods longer than ~0.6x the baseline, a different
   amplitude in the two halves of the data, and unconstrained eccentric fits. Peaks at 2x or 1/2x
   the period in an indicator, harmonics of other signals and periods inside the typical rotation
   range are cautions;
5. cross-match by every SIMBAD identifier and by position against the NASA Exoplanet Archive
   (`pscomppars` and every period in `ps`) and exoplanet.eu.

Calibration: the blind search re-finds the planets of HD 40307, HD 69830, GJ 163, GJ 876, HD 51608
and many more; of the confirmed planets on these stars with a published K >= 2 m/s it recovers
about two thirds. The misses are mostly long-period eccentric giants with few nights and small
planets that were found with transit priors, more data or activity modelling.

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
- Add a single-transit search for long-period planets.
- Add centroid / difference-image vetting from target pixel files (`lk.search_targetpixelfile`).
- Model stellar activity in the RVs with a Gaussian process instead of rejecting by indicator.
- Add newer public HARPS / ESPRESSO data (ESO archive, after 2018) to the RV stars.
- Run on the newest sector, ahead of the TOI alerts.
