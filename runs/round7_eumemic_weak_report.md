# Round 7: mining Eumemic's weak signals

Branch `round7-eumemic-weak` (worktree `.`), 2026-10-13. Every process ran at `nice -n 20`. Nothing was pushed, posted, emailed, uploaded or submitted. GitHub, ExoFOP, MAST, the NASA Exoplanet Archive, arXiv, Zenodo, VizieR and X were only read; there are no new accounts, issues, stars or follows.

## Bottom line

- **Mined:** all 264,468 public BLS signals in Eumemic's seven `all_signals*.csv.gz` tables (K/M dwarfs, TOI hosts, faint M dwarfs T 11.5-13.5, the TCE re-check and the 20-100 d passes).
  - 35,407 passed the table-level cuts (SNR 6-9, at least 3 transits, small star, clean shape, not his and not catalogued).
  - We inventoried the products for 3,898 of them and selected **120 signals on 110 stars** with the most extra data.
  - **Extra data, median:** 3 sectors he did not have on top of his 3. Almost all of it is QLP sectors his reader dropped.
- **One strong new candidate: TIC 144439305**, the best result of this round.
  - **Planet:** P = 0.9398094 ± 0.0000023 d, 569 ± 58 ppm, Rp = 1.28 (1.20-1.36) R⊕, Teq ≈ 875 K.
  - **Star:** M3.5 dwarf at 62 pc, T = 12.15, 3,505 K, 0.46 R☉.
  - **He had it:** it is his rank-1 signal at SNR 8.35, below his SNR 9 pre-filter. It was never vetted.
  - **Hold-out:** it reproduces in data he never had: QLP S96/103/104 give z = 5.2, with depth 468 ± 91 ppm against 590 ± 94 in his-era data.
  - **Fresh data:** it reproduces again in TESScut S105/106, which no pipeline has processed: 909 ± 137 ppm, z = 6.6.
  - **Pixels:** on target.
  - **Binary checks:** clean. Gaia's epoch RVs exclude an eclipsing star on the target.
  - **FPP:** TRICERATOPS gives FPP 0.078 / NFPP 0.002 on its own, and **0.008 ± 0.002 / 0.00004** with Gaia/MOLUSC limits and pixel-excluded neighbours cleared. That is formally below the validation bar, but with no high-resolution imaging.
  - **Novelty:** clean in every source checked.
- **One weak candidate:** TIC 231706739, 16.566 d, 1.8 R⊕. Hold-out z = 4.3, full SNR 7.4, the pixels can't tell the target from a T = 19 star 22″ away, and FPP is 0.13. It isn't claimable.
- **Three strong hold-out signals that are not planets:**
  - TIC 166422847 is **TOI-3282.01**: a neighbour 33″ away with the same period, leaking into the aperture.
  - TIC 270188035 is a **blended EB on a Gaia star 41.5″ away**.
  - TIC 119876352 is a **short-period oscillation, not a transit**.
- **Everything else is noise.** Excluding those five, the other 113 hold-out z-values are N(0,1): mean 0.06, std 1.04, KS p = 0.79. A real signal of the claimed depth would have given a median z of 6.2. 85 of 120 are not reproduced at > 3σ.
- **Expected noise survivors** of the pre-registered test:
  - random ephemerides on the same stars predict **1.0**;
  - running the whole procedure on inverted light curves gives **0 of 120**.

  We got 2 survivors. One is the candidate (no random ephemeris out of 300 reached its z) and one is an astrophysical non-planet. **The noise test controls noise, not blends;** the pixels and diagnostics did that part.
- **Completeness (120 injections of his own depths at other periods):** 46% pass, rising to 60% where the hold-out data should give z ≥ 6. So "not reproduced" means "not confirmed", not "proven false": about 30% of real signals at his strength would also land there.

## 1. His data, thresholds and known gaps

- **Source:** `github.com/eumemic/exoplanets`, a shallow sparse clone of commit `da722613` (2026-10-13 07:39 UTC). No later push as of this report.
- **Tables mined:** `results_public/all_signals.csv.gz` (K/M dwarfs T ≤ 11.5), `_toi_hosts`, `_faint_t1` (11.5-12.5), `_faint_t2` (12.5-13.5), `_tce`, `_long`, `_long_toi_hosts`. That's 264,468 signals on about 132k stars.
  - The full-sky search (4.04 M stars) publishes only its 2,426 vetted signals (`candidates_allsky.csv`), not every signal.
- **His promotion path:**
  1. Pre-filter at **SNR ≥ 9**, ≥ 3 transits, duration and shape cuts.
  2. Then LEO-Vetter (17 tests), **MES ≥ 7.1**, per-sector χ²/dof ≤ 5, centroid within 15″ and no capable TIC neighbour, then TRICERATOPS.
  - For the full sky only, a classifier also pulls some SNR 7.5-9 signals into tier L1.
  - **So in the earlier searches the SNR < 9 band was never vetted at all.** TOI-6000 c sat there at SNR 7.05.
- **Pipeline gaps we used or checked:**
  - **QLP files with empty uncertainty columns were silently dropped** in S94-97 and S101-104 (24-82% of QLP light curves there, about 900k overall; PR #6 fixed the reader, but the re-search is still "planned").
    - Our extractor keeps them.
    - On our 120 signals, the his-era sector count we infer by excluding exactly those files reproduces his `n_sectors` exactly for 113. In 6 of the other 7 we label more sectors "his" than he counted, which is conservative. In one (TIC 257567390) we label fewer, so its hold-out may include data he had; it is in the not-reproduced pile anyway. That is the "data he did not have" below.
  - **QLP timing error** from the second half of S99 to S104 (MAST notice; minutes). Recorded per target.
    - For TIC 144439305 the S103/S104 O−C values are −10 and −3 min, within the per-sector scatter (±2-14 min) and small against a 52-min transit.
  - **Two more, from his issue tracker:**
    - QLP S74-79 are missing from his full-sky index (#9);
    - the tess-point pseudo-sector "1751" (#14). We hit it too and filtered it.

## 2. Selection (`round7/select.py`)

| Step | Signals |
|---|---|
| All signals, 7 tables | 264,468 |
| His SNR 6-9 (below his pre-filter) | 134,403 |
| ≥ 3 good transits; P ≥ 0.5 d | 107,841 |
| Not flagged repeat / rotation alias / TOI alias | 98,107 |
| Small star (ρ★ > 2.5 g/cc; the TIC cut comes later) | 78,000 |
| Odd/even < 3σ, secondary < 3σ | 71,385 |
| No single dominant event (max_ses_frac < 1/√n + 0.35) | 61,990 |
| Duration consistent with the star (0.25-1.6 × expected) | 54,476 |
| SDE ≥ 8 | 35,560 |
| Deduplicated across his tables; not in any of his lists or tiers; no TOI/CTOI/confirmed planet on the star at an alias (10-12 tables) | **35,407** on 31,566 stars |
| MAST product inventory: every TOI-host, long-period and K/M signal + the best 900 each of faint t1/t2 | 3,898 on 3,295 stars |
| Extra data exists (E_hold ≥ 2) and TIC Teff < 4200 K or R★ < 0.7 | 935 |
| **Top 120 by expected power with all sectors** (E_full + 0.5 E_hold) | **120 signals, 110 stars** |

- **Expected power:** E_full = SNR_his × √(N_all / N_his) and E_hold = SNR_his × √(N_new / N_his).
- **Sector sources:** N_all comes from MAST (2-min > TESS-SPOC > QLP per sector, through S107). N_his is his `n_sectors`.
- **What was picked:**
  - 64 bright K/M dwarfs and 56 faint M dwarfs;
  - median Teff 3,690 K, R★ 0.54;
  - T 8.3-13.5, P 0.53-29.4 d.
- **No TOI-host signal made the cut.** Those stars already have many sectors, so new data adds little, and round 5 searched the 184 best of them.

## 3. The test (`round7/revet.py`), fixed before any result was seen

**Data.** Each star gets all its sectors: SPOC 2-min where it exists, else TESS-SPOC, else QLP (including the empty-error QLP files). That is 766 files: 81 2-min, 172 TESS-SPOC and 513 QLP. Each sector is labelled:
- **his:** a sector ≤ 105 that is not an empty-error QLP file;
- **new:** everything else.

**Pre-registered procedure per signal:**
1. Run a narrow BLS (his period ± 0.5%) on the his-era data alone to get the ephemeris.
2. **Hold-out:** measure the depth at that ephemeris in the new data alone, with red-noise errors, giving z_hold.
3. Run a narrow BLS on all data, plus odd/even, secondary and per-season tests.
4. **Survivor** = all of these:
   - z_hold ≥ 2.5;
   - hold-out depth within 3σ of the his-era depth;
   - SNR_full ≥ 8 within 0.1% of his period;
   - |odd/even| < 3, secondary < 3σ, season χ² p > 0.001.

**Noise.**
- **Random ephemerides:** for each star, z_hold at 300 random ephemerides (P within ±5% but > 0.5% from his, random T0) in the same new data gives a per-star null pass rate.
- **Inverted flux:** the whole procedure is also run on every star's inverted light curve.

**Completeness.** Each target gets a transit with his own depth injected at another period (×0.7-1.3), and the identical procedure is run.

## 4. Results

| Verdict | n |
|---|---|
| Survivor (pre-registered) | 2 (TIC 144439305, TIC 119876352) |
| Not reproduced: hold-out depth > 3σ below his-era depth | 85 |
| Inconclusive: hold-out too weak, or flags | 33 |

**Noise.**
- Sum of the per-star null pass rates: **0.98 expected survivors from noise** (Gaussian 0.6% × 118 = 0.73).
- Inverted light curves: **0 of 120** pass.
- Neither survivor's z_hold was reached by any of its 300 random ephemerides (z = 5.2 and 7.0; Gaussian p ~10⁻⁷ and 10⁻¹²).

**Population.**
- All 118 hold-out z: mean 0.32, std 1.72.
- Without the five strong or known signals: mean 0.06, std 1.04, KS p = 0.79 against N(0,1).
- Inverted: mean 0.09, std 1.06.
- Real signals at the claimed depths would give a median z of 6.2 (81 of 120 have an expected z ≥ 5).
- **Reading:** even the best-ranked SNR 6-9 signals in his tables are about 95% noise. Our full-data SNR is 0.84× his on median, which is what noise does when you add data.

**Completeness (120 injections).**
- 46% pass overall: 0% where the expected hold-out z is < 4, 39% at 4-6 and 60% at ≥ 6.
- Losses:
  - the his-era narrow BLS locks onto the wrong peak 26% of the time;
  - winner's-curse depth bias (his-era depth ×1.09, hold-out ×0.93 of the injected depth);
  - season and product flags.
- **30% of injected real signals are labelled "not reproduced".** Read that row as "not confirmed".
- Implied, roughly: the 120 contain about 2-10 astrophysical signals, and we found 5.

### Ranked table (top 12; all 120 rows are in the appendix and in `runs/round7_eumemic_weak/round7_table.csv`)

Depth: his BLS depth / ours (all data). SNR: his / ours on his-era data / ours on all data. z hold-out: measured (expected if his depth were real). FPP/NFPP: TRICERATOPS. "TOI/CTOI live": ExoFOP tables fetched 2026-10-13 20:11 UTC, same star or within 3′ at the same period or an alias.

| # | TIC | P (d) | Depth ppm | SNR | z hold-out | Rp (R⊕) | Transits | Pixel | FPP / NFPP | TOI/CTOI live | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **144439305** | 0.939809 | 667 / 569 | 8.4 / 6.3 / 9.9 | **5.2** (6.5); S105/106 cutouts z 6.6 | 1.28 (1.20-1.36) | 132 | **on target**: best source in 7/8 sectors; neighbours Δχ² ≥ 25 worse | 0.078 / 0.0020 (MOLUSC + cleared: **0.008 / 0.00004**) | none | **STRONG CANDIDATE, new** |
| 2 | 231706739 | 16.56563 | 1069 / 749 | 6.6 / 6.5 / 7.4 | 4.3 (7.4) | 1.83 | 8 | unresolved: target vs a T = 19.0 star at 22″ (Δχ² 56 vs 60) | 0.13 / 0.0000 | none | weak candidate, new; not claimable |
| 3 | 166422847 | 2.94783 | 1661 / 2728 | 8.0 / 6.3 / 14.4 | 11.4 (3.1) | 2.22 | 36 | – | – | **TOI-3282.01 (APC) 33″, same P** | KNOWN: neighbour's signal |
| 4 | 119876352 | 0.58703 | 3053 / 3016 | 8.5 / 8.1 / 11.0 | 7.0 (8.8) | 3.27 | 173 | – | – | none | REJECTED: ±3,000 ppm ~1.7 h oscillation; depth varies 1.7-8.4 ppt by sector |
| 5 | 270188035 | 19.23305 | 5340 / 6259 | 8.9 / 10.1 / 14.1 | 3.5 (5.1) | 3.43 | 6 | **off target**: Gaia T = 15.8 star at 41.5″ (target Δχ² 3,653) | – | none | REJECTED: blended EB |
| 6 | 395108648 | 13.01769 | 56526 / 82230 | 8.7 / 6.9 / 7.7 | 2.6 (2.2) | 15.2 | 7 | – | – | none | inconclusive (EB-size, odd/even) |
| 7 | 49623957 | 7.15440 | 25150 / 40883 | 8.2 / 4.2 / 12.7 | 1.9 (1.1) | 8.7 | 9 | – | – | none | inconclusive (EB-size) |
| 8 | 111227963 | 0.71619 | 6759 / 3821 | 7.6 / 5.9 / 8.3 | 1.7 (2.7) | 4.0 | 137 | – | – | none | inconclusive |
| 9 | 201294141 | 20.23497 | 1159 / 517 | 7.6 / 3.6 / 6.0 | 1.4 (5.6) | 1.5 | 8 | – | – | none | inconclusive |
| 10 | 281591952 | 13.00514 | 1421 / 2310 | 7.5 / 5.0 / 16.4 | 1.2 (3.8) | 3.0 | 13 | – | – | none | inconclusive |
| 11 | 321981039 | 19.41724 | 1133 / 1512 | 7.6 / 4.9 / 12.2 | 0.9 (4.3) | 2.4 | 5 | – | – | none | inconclusive |
| 12 | 387132154 | 6.65310 | 837 / 622 | 8.2 / 3.1 / 4.5 | 0.7 (3.3) | 1.4 | 36 | – | – | none | inconclusive |

**Novelty of the five deep-looked signals, per source** (all read 2026-10-13):

| Source | 144439305 | 231706739 | 270188035 | 119876352 | 166422847 |
|---|---|---|---|---|---|
| ExoFOP TOI + CTOI (live), star and 5′ | none | none | none | none | **TOI-3282.01 at 33″, same P** |
| ExoFOP target page (imaging / spectra / notes) | empty | empty | – | – | – |
| NASA Exoplanet Archive (`ps`, `toi`) | none | none | none | none | – |
| SPOC TCEs (DV reports) | no 2-min data, no DV | 2 sectors of 2-min, no DV | no DV | no DV | – |
| arXiv API ("TIC n", metadata only) | 0 | 0 | 0 | 0 | 0 |
| Eumemic (all lists, tiers, follow-up, vetted, validation tables, issues #1-17) | **raw signal only**: rank 1, SNR 8.35, `all_signals_faint_t1` | raw signal only (`all_signals`) | raw only | raw only | raw only |
| Marchand registry (1,250 scanned, candidates, preregistered, external, known TICs) | none | none | none | none | none |
| Paras Chopra, Aiden Huang, Avijit Thawani target and signal lists; Samson Fraser repos (code search) | none | none | none | none | none |
| GitHub code search (all public) | none relevant | 2 star lists (tess-quiet `boring_tic_ids.txt`, Astro_Transit_GPU), no claim | none | none | – |
| Zenodo search | 0 | 0 | 0 | 0 | – |
| X recent search (`xurl`, 4 queries, 0 posts read) | 0 | 0 | 0 | 0 | – |
| VSX, TESS EB catalogue, Gaia DR3 NSS / variability, LEO-Vetter M-dwarf list | none | – | – | – | – |

Limits of these checks:
- arXiv's API searches metadata, not tables, and ADS full text was not available.
- X recent search covers 7 days only.
- GitHub code search skips files over 384 kB and gzipped tables (the positive control, TOI-4342's TIC in Eumemic's repo, returned 13 hits).

## 5. TIC 144439305: what we know

**Star.**
- TIC 144439305 = Gaia DR3 6532568669359184000.
- RA 358.86748, Dec −44.12371. T = 12.15, G = 13.29.
- ϖ = 16.07 mas (62.2 pc). Teff 3,505 K, R★ 0.465, M★ 0.462 (TIC 8.2). RUWE 1.07.
- 2-min data: none. TESS-SPOC S2, S29, S69; QLP S96, S103, S104; TESScut S105, S106. **No further TESS sectors planned through about S134** (tess-point 0.9.5).

**Signal.**
- **Ephemeris:** P = 0.9398094 ± 0.0000023 d, T0 = BJD_TDB 2460155.4571 ± 0.0027, from per-sector times with errors scaled to χ²/dof = 1. Predicted transits stay within ±6 min through mid-2027.
- **Depth and size:** 569 ± 58 ppm, SNR 9.9 with all light curves. The limb-darkened fit gives Rp/R★ = 0.0252 (0.0237-0.0267), **Rp = 1.28 (1.20-1.36) R⊕**, b ≈ 0.33 (0.10-0.60).
- **Density:** ρ★ from the fit is 9.9 (6.0-12.8) g/cc against 6.5 from the TIC (0.9σ).
- **Orbit and climate:** a/R★ = 6.7, a = 0.0145 au, S ≈ 139 S⊕, Teq ≈ 875 K.
- **Mass and RV:** a rocky planet would be about 2 M⊕, giving K ≈ 2.4 m/s.

**Data the signal survives:**

| Data set | Who has searched it | Depth (ppm) | z |
|---|---|---|---|
| His era: TESS-SPOC S2, 29, 69 | Eumemic (SNR 8.35) | 590 ± 94 | 6.3 |
| Dropped by his reader: QLP S96, 103, 104 | nobody | 468 ± 91 | **5.2** |
| TESScut S105 + S106, fixed aperture (no light-curve products exist yet) | nobody | 909 ± 137 | **6.6** |
| Fixed aperture, all 8 sectors | – | 658 ± 78 (χ² p = 0.29) | 8.4 |

**Vetting.**
- **Light-curve checks:**
  - odd/even 0.7σ;
  - phase 0.5: 2 ± 65 ppm (P/2 excluded; equal odd/even rules out 2P);
  - secondary at any phase 1.9σ;
  - season χ² p = 0.57, per-transit χ² p = 0.31;
  - TESS-SPOC 589 ± 93 vs QLP 479 ± 92 ppm (product-independent);
  - sinusoid ratio 0.38.
- **Pixels** (`deep_vet.py`, 8 TESScut sectors):
  - the best Gaia source is the target, best in 7 of 8 sectors;
  - the free-fit centroid is 9.9 ± 6.0″ from the target, consistent with it;
  - every neighbour fits worse: T = 18.8 at 24.7″ by Δχ² 25 (and it would need a 37-47% eclipse), T = 14.2 and 15.4 at 48″ by Δχ² 61-66;
  - the ring test favours the target (χ² 53 vs 95 or more).
- **Binary checks:**
  - **Gaia DR3 epoch RVs:** 24 epochs, χ² p = 0.95, per-epoch σ ≈ 6 km/s. They exclude an eclipsing star on the target (a 0.08 M☉ companion at this period would give K ≈ 26 km/s).
  - No VSX, TESS-EB or Gaia NSS/variability entry.
- **TRICERATOPS** (5 runs × 10⁶ draws, 3×3-pixel FFI apertures):
  - Default: FPP 0.078 ± 0.002, NFPP 0.0020. The residual is STP (a planet on an unresolved bound companion, 0.076).
  - With MOLUSC (Gaia contrast + RUWE; 235,272 of 500,000 companions survive) and the neighbours beyond 30″ cleared (pixel Δχ² ≥ 61): **FPP 0.0082 ± 0.0019, NFPP 0.00004**. The 24.7″ star is kept, as in round 6.
  - These meet FPP < 0.015 and NFPP < 0.001, **but there is no high-resolution imaging.** The STP and background-star scenarios within about 1-2″ are constrained only by Gaia and priors.
- **Archival plates: not run.** The star has moved only about 4.9″ since the 1977 SERC-J epoch, roughly 1.5 plate FWHM. A background star would need ΔT ≤ 7.4 to mimic 569 ppm, which the plates can't exclude that close to a T = 12 star.

**Context.**
- Small ultra-short-period planets around M dwarfs are a known population (LHS 3844 b, GJ 367 b, TOI-1685 b, TOI-1442 b).
- A 1.3 R⊕, 875 K planet at 62 pc is a reasonable thermal-emission target, but we haven't computed an ESM.
- The signal sat below Eumemic's SNR 9 pre-filter, and half of its data was dropped by the QLP bug. **Both gaps this round set out to exploit are what hid it.**

## 6. What's worth an RNAAS or CTOI

1. **TIC 144439305: yes, and quickly.** It's the cleanest candidate we've had.
   - Three independent data sets agree, two of which nobody has searched.
   - It's on target in the pixels, binary-clean, and formally below the FPP bar with Gaia limits.
   - **Credit Eumemic:** his public raw table lists the signal (rank 1, SNR 8.35, unvetted). Any note should say so, the same way AAS82279 v2 does.
   - **Suggested order (Ricardo's call; nothing sent):**
     1. CTOI upload: the fastest timestamp, and it needs his ExoFOP account.
     2. A short RNAAS: "A 1.3 R⊕ ultra-short-period planet candidate around a nearby M3.5 dwarf, recovered from data a public search dropped."
     3. Optionally a preregistered ephemeris on Zenodo.
   - **Follow-up:** 569 ppm on a T = 12, I ≈ 11 M dwarf is within reach of 1-2 m telescopes in the red over several transits (it transits nightly). CHEOPS could too, and speckle imaging would remove the STP residual.
2. **TIC 231706739: no.** SNR 7.4, the pixels can't separate the target from a faint neighbour, and FPP is 0.13. Revisit if new sectors appear (none are planned).
3. **Nothing else.** The other "inconclusive" signals are noise-consistent or EB-sized.

## 7. Caveats

- **Selection:** the 120 were chosen for the most new data. Stars whose weak signal can't be tested with new data (most TOI hosts and CVZ stars) weren't tested at all, so this is not a census of his weak band.
- **Completeness:** 46% overall. "Not reproduced" is not "false": by the injections, about 30% of real signals at his strength would land there.
- **What each test controls:**
  - the noise tests control noise and systematics, not blends;
  - the pixels and diagnostics, applied to 4 of the 5 signals with strong hold-outs, did that;
  - a pre-registered survivor (TIC 119876352) failed on sight, so the statistical pass alone is not a candidate.
- **Precision of the analysis:**
  - his-era vs new labelling is per sector and data-driven (empty-error flags), and it matches his `n_sectors` for 113 of 120 signals;
  - products are mixed (TESS-SPOC/QLP/2-min) with per-sector normalisation, and the product depths of the candidate agree;
  - no ground-based or imaging data exist for either candidate.

## 8. Resources and files

- **Disk:** about 0.7 GB new (TESScut +0.50 GB, FFI npz +0.06, 2-min npz +0.02, Eumemic clone and novelty tables 0.13, run outputs 0.02). Free space stayed at 46-50 GB.
- **Network:** an estimated 1.5-2 GB transferred: 766 light-curve files, 24 TESScut cutouts, a 112 MB sparse repo, MAST/Gaia/ExoFOP queries. The Mini's en0 counter moved 4.5 GB over the session, but that includes all other traffic on the machine.
- **X:** 4 search queries, 0 posts read.
- **Code** (`round7/`): `select.py` (stages 1-2), `revet.py` (fetch, hold-out test, nulls, injections, collect), `diag.py`, `ephem.py`, `trifold.py`, `table.py`, `deepvet_targets.json`. Pixel vetting reuses `deep_vet.py`; TRICERATOPS reuses `run_triceratops.py` and `round6/run_tri.py` (MOLUSC, clearing); MOLUSC comes from round 6's `.tools`.
- **Outputs** (`runs/round7_eumemic_weak/`):
  - selection: `stage1*.csv`, `stage2_*.csv`, `targets.csv`;
  - per-signal JSON: `revet/`, `inject/`;
  - results: `revet.csv`, `round7_table.csv/.md`;
  - diagnostics: `diag/*.png` (folds), `deepvet/*.json/png` (pixels), `tri/` (TRICERATOPS, MOLUSC);
  - novelty: `novelty_live_toi_ctoi.csv`, `spoc_tces.json`.
  - The novelty source files are in `<data>/r7/novelty/`.

## Appendix: all 120 signals

In this table, Rp and depth are box-model values from the narrow re-search on all data (the limb-darkened Rp for TIC 144439305 is 1.28 R⊕). SNR is his / ours on his-era data / ours on all data. Rows are ordered by verdict class, then by hold-out z.

| # | TIC | P (d) | depth ppm (his / ours) | SNR his / ours-his-era / ours-all | z hold-out (expected if real) | Rp (Re) | transits | pixel | FPP / NFPP | TOI/CTOI live | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 144439305 | 0.93981 | 667 / 549 | 8.4 / 6.3 / 9.5 | 5.2 (6.5) | 1.19 | 132 | on target (best source 7/8 sectors; neighbours dchi2 >= 25 worse) | 0.078 / 0.0020 (MOLUSC+cleared: 0.008 / 0.0000) | none | CANDIDATE (strong) |
| 2 | 231706739 | 16.56563 | 1069 / 749 | 6.6 / 6.5 / 7.4 | 4.3 (7.4) | 1.83 | 8 | unresolved: target vs T=19.0 star at 22" (dchi2 56 vs 60) | 0.130 / 0.0000 | none | WEAK CANDIDATE (pixel-unresolved, SNR 8) |
| 3 | 166422847 | 2.94783 | 1661 / 2728 | 8.0 / 6.3 / 14.4 | 11.4 (3.1) | 2.22 | 36 | not run (known neighbour signal) | - | TOI 3282.01 at 33" P=2.9478 x1 | KNOWN: TOI-3282.01 (APC) 33" away, same 2.94783 d period |
| 4 | 119876352 | 0.58703 | 3053 / 3016 | 8.5 / 8.1 / 11.0 | 7.0 (8.8) | 3.27 | 173 | not run (rejected on light curve) | - | none | REJECTED: dip sits in a +-3,000 ppm ~1.7 h oscillation; depth varies 1.7-8.4 ppt by sector |
| 5 | 270188035 | 19.23305 | 5340 / 6259 | 8.9 / 10.1 / 14.1 | 3.5 (5.1) | 3.43 | 6 | OFF target: Gaia T=15.8 star at 41.5" (target dchi2 3,653) | - | none | REJECTED: blended EB on a neighbour |
| 6 | 395108648 | 13.01769 | 56526 / 82230 | 8.7 / 6.9 / 7.7 | 2.6 (2.2) | 15.24 | 7 | - | - | none | inconclusive (hold-out too weak or flags) |
| 7 | 49623957 | 7.15440 | 25150 / 40883 | 8.2 / 4.2 / 12.7 | 1.9 (1.1) | 8.68 | 9 | - | - | none | inconclusive (hold-out too weak or flags) |
| 8 | 111227963 | 0.71619 | 6759 / 3821 | 7.6 / 5.9 / 8.3 | 1.7 (2.7) | 4.03 | 137 | - | - | none | inconclusive (hold-out too weak or flags) |
| 9 | 201294141 | 20.23497 | 1159 / 517 | 7.6 / 3.6 / 6.0 | 1.4 (5.6) | 1.53 | 8 | - | - | none | inconclusive (hold-out too weak or flags) |
| 10 | 281591952 | 13.00514 | 1421 / 2310 | 7.5 / 5.0 / 16.4 | 1.2 (3.8) | 3.00 | 13 | - | - | none | inconclusive (hold-out too weak or flags) |
| 11 | 321981039 | 19.41724 | 1133 / 1512 | 7.6 / 4.9 / 12.2 | 0.9 (4.3) | 2.44 | 5 | - | - | none | inconclusive (hold-out too weak or flags) |
| 12 | 387132154 | 6.65310 | 837 / 622 | 8.2 / 3.1 / 4.5 | 0.7 (3.3) | 1.38 | 36 | - | - | none | inconclusive (hold-out too weak or flags) |
| 13 | 450827167 | 28.82747 | 32621 / 37719 | 8.9 / 6.3 / 6.6 | 0.7 (2.7) | 11.52 | 3 | - | - | none | inconclusive (hold-out too weak or flags) |
| 14 | 293568983 | 3.60907 | 1489 / 969 | 8.0 / 5.4 / 5.8 | 0.6 (4.3) | 1.12 | 45 | - | - | none | inconclusive (hold-out too weak or flags) |
| 15 | 295047206 | 7.79527 | 15320 / 9569 | 8.2 / 5.3 / 8.9 | 0.5 (0.8) | 4.93 | 23 | - | - | none | inconclusive (hold-out too weak or flags) |
| 16 | 464620974 | 1.39042 | 690 / 416 | 8.7 / 5.1 / 6.9 | 0.5 (4.5) | 1.33 | 100 | - | - | none | inconclusive (hold-out too weak or flags) |
| 17 | 422880113 | 8.25943 | 18248 / 12425 | 7.6 / 4.0 / 4.9 | 0.4 (0.5) | 4.57 | 14 | - | - | none | inconclusive (hold-out too weak or flags) |
| 18 | 299754340 | 6.92186 | 642 / 312 | 7.1 / 3.9 / 5.6 | 0.4 (4.8) | 1.27 | 13 | - | - | none | inconclusive (hold-out too weak or flags) |
| 19 | 220488995 | 17.40104 | 348 / 617 | 6.8 / 4.4 / 10.7 | 0.3 (4.2) | 1.80 | 5 | - | - | none | inconclusive (hold-out too weak or flags) |
| 20 | 234298028 | 8.45075 | 2639 / 4331 | 7.3 / 6.5 / 10.6 | 0.2 (0.6) | 3.60 | 16 | - | - | none | inconclusive (hold-out too weak or flags) |
| 21 | 111227963 | 5.19475 | 7482 / 17806 | 7.3 / 4.2 / 19.0 | 0.2 (1.5) | 8.69 | 19 | - | - | none | inconclusive (hold-out too weak or flags) |
| 22 | 396169984 | 24.51249 | 11164 / 30061 | 8.7 / 8.9 / 10.6 | 0.1 (1.0) | 11.15 | 2 | - | - | none | inconclusive (hold-out too weak or flags) |
| 23 | 355638067 | 2.79524 | 1837 / 1201 | 7.3 / 5.0 / 10.8 | 0.1 (2.4) | 2.17 | 67 | - | - | none | inconclusive (hold-out too weak or flags) |
| 24 | 111227963 | 17.88016 | 11996 / 122471 | 8.9 / 4.8 / 20.6 | -0.1 (1.9) | 22.79 | 4 | - | - | none | inconclusive (hold-out too weak or flags) |
| 25 | 403236288 | 23.09722 | 728 / 436 | 8.4 / 2.5 / 6.8 | -0.1 (3.4) | 1.49 | 9 | - | - | none | inconclusive (hold-out too weak or flags) |
| 26 | 155941095 | 16.23853 | 4842 / 1785 | 8.5 / 3.6 / 6.3 | -0.1 (4.0) | 2.32 | 10 | - | - | none | inconclusive (hold-out too weak or flags) |
| 27 | 115347568 | 5.68576 | 2891 / 3769 | 8.1 / 7.1 / 13.3 | -0.5 (1.5) | 4.37 | 15 | - | - | none | inconclusive (hold-out too weak or flags) |
| 28 | 387266782 | 1.47208 | 706 / 388 | 7.7 / 2.3 / 4.4 | -0.5 (3.7) | 1.18 | 85 | - | - | none | inconclusive (hold-out too weak or flags) |
| 29 | 389149025 | 0.56645 | 418 / 397 | 6.6 / 4.3 / 5.3 | -0.8 (0.0) | 1.06 | 146 | - | - | none | inconclusive (hold-out too weak or flags) |
| 30 | 419791995 | 5.19295 | 795 / 737 | 6.8 / 5.7 / 9.6 | -0.8 (0.8) | 1.25 | 33 | - | - | none | inconclusive (hold-out too weak or flags) |
| 31 | 450827167 | 13.72844 | 14699 / 13749 | 8.4 / 4.9 / 4.8 | -0.9 (0.8) | 6.96 | 8 | - | - | none | inconclusive (hold-out too weak or flags) |
| 32 | 140069648 | 12.19933 | 132016 / 156490 | 8.3 / 3.8 / 14.0 | -1.2 (1.4) | 9.26 | 11 | - | - | none | inconclusive (hold-out too weak or flags) |
| 33 | 342151941 | 0.68369 | 1669 / 1848 | 8.3 / 8.6 / 8.7 | -1.4 (0.1) | 2.03 | 140 | - | - | none | inconclusive (hold-out too weak or flags) |
| 34 | 231024400 | 26.88262 | 5150 / 14951 | 8.7 / 4.8 / 12.1 | nan (nan) | 4.37 | 4 | - | - | none | inconclusive (hold-out too weak or flags) |
| 35 | 369579378 | 13.43895 | 2079 / 849 | 7.6 / 4.8 / 6.1 | nan (nan) | 1.53 | 8 | - | - | none | inconclusive (hold-out too weak or flags) |
| 36 | 120586821 | 26.97016 | 3374 / 4878 | 8.8 / 4.6 / 7.2 | 2.1 (13.6) | 4.16 | 3 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 37 | 360740748 | 1.15077 | 264 / 176 | 6.7 / 5.7 / 5.3 | 1.9 (7.1) | 1.00 | 126 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 38 | 259783240 | 12.35242 | 1081 / 548 | 8.6 / 4.2 / 6.7 | 1.9 (9.8) | 0.96 | 12 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 39 | 139200209 | 21.61155 | 1345 / 1091 | 7.9 / 7.2 / 5.9 | 1.9 (9.6) | 1.54 | 5 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 40 | 161439851 | 5.52962 | 213 / 131 | 7.5 / 6.1 / 5.9 | 1.7 (6.1) | 0.69 | 27 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 41 | 403236288 | 28.26001 | 906 / 959 | 8.2 / 5.2 / 10.2 | 1.6 (7.8) | 2.21 | 5 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 42 | 49450836 | 1.41869 | 2592 / 1559 | 8.5 / 6.9 / 5.9 | 1.5 (7.3) | 2.33 | 77 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 43 | 70553159 | 14.26135 | 340 / 381 | 7.6 / 3.9 / 8.2 | 1.3 (8.8) | 1.09 | 7 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 44 | 257500473 | 1.70315 | 2237 / 1620 | 8.7 / 5.2 / 6.9 | 1.3 (8.8) | 3.00 | 47 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 45 | 231865504 | 1.08649 | 444 / 221 | 7.4 / 5.3 / 5.0 | 1.1 (6.0) | 0.70 | 153 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 46 | 237313698 | 4.72667 | 342 / 165 | 6.9 / 5.2 / 4.9 | 1.1 (7.1) | 0.89 | 28 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 47 | 201289095 | 26.98043 | 49979 / 196134 | 8.8 / 9.8 / 16.4 | 1.1 (8.1) | 25.49 | 3 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 48 | 80262344 | 5.02449 | 970 / 342 | 6.2 / 5.0 / 4.5 | 1.1 (6.1) | 1.04 | 19 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 49 | 419609686 | 2.22266 | 576 / 137 | 6.2 / 5.0 / 4.7 | 1.0 (7.9) | 0.82 | 56 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 50 | 277803433 | 1.58214 | 815 / 704 | 8.0 / 7.5 / 9.0 | 1.0 (4.6) | 1.20 | 115 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 51 | 115356289 | 10.54647 | 870 / 328 | 7.9 / 6.0 / 5.0 | 1.0 (8.5) | 1.20 | 9 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 52 | 234299972 | 23.28426 | 1516 / 907 | 6.6 / 5.1 / 7.1 | 0.9 (10.3) | 1.28 | 7 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 53 | 161458236 | 13.88106 | 548 / 776 | 7.3 / 6.2 / 12.5 | 0.9 (88.7) | 1.74 | 6 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 54 | 355687587 | 0.82202 | 156 / 70 | 6.8 / 5.7 / 5.8 | 0.9 (6.5) | 0.59 | 218 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 55 | 100066059 | 22.46795 | 4274 / 1776 | 8.4 / 4.4 / 6.6 | 0.9 (6.9) | 2.75 | 7 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 56 | 201251602 | 1.99856 | 844 / 465 | 7.3 / 5.4 / 5.0 | 0.8 (5.6) | 1.31 | 110 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 57 | 290369580 | 20.59415 | 6089 / 13634 | 8.3 / 4.7 / 19.0 | 0.8 (12.1) | 5.71 | 7 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 58 | 320659505 | 6.17271 | 480 / 221 | 6.2 / 5.5 / 4.9 | 0.8 (7.1) | 0.95 | 21 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 59 | 352405596 | 13.30950 | 3393 / 13805 | 8.5 / 6.9 / 38.0 | 0.8 (6.1) | 7.12 | 7 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 60 | 466306858 | 1.68641 | 416 / 211 | 6.5 / 5.0 / 5.2 | 0.7 (6.6) | 0.94 | 84 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 61 | 257463364 | 28.67642 | 2213 / 980 | 8.2 / 6.5 / 6.6 | 0.7 (6.4) | 1.93 | 9 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 62 | 161376559 | 11.01958 | 1025 / 296 | 6.9 / 5.7 / 5.1 | 0.7 (5.8) | 0.89 | 12 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 63 | 341307991 | 4.52091 | 353 / 149 | 7.6 / 5.4 / 4.8 | 0.6 (5.0) | 0.78 | 31 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 64 | 257567390 | 10.76347 | 3190 / 2363 | 8.9 / 3.8 / 5.3 | 0.5 (7.2) | 1.64 | 9 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 65 | 235010431 | 0.53261 | 495 / 234 | 7.6 / 5.6 / 5.4 | 0.5 (7.9) | 0.97 | 265 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 66 | 408037269 | 25.79416 | 556 / 2074 | 8.8 / 5.1 / 10.5 | 0.4 (5.5) | 2.75 | 1 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 67 | 320204496 | 0.67915 | 243 / 88 | 6.2 / 4.4 / 4.9 | 0.4 (5.3) | 0.63 | 209 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 68 | 161208104 | 12.20804 | 675 / 297 | 6.6 / 5.1 / 5.2 | 0.4 (9.0) | 0.99 | 7 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 69 | 206354321 | 0.95566 | 121 / 60 | 7.8 / 6.1 / 5.6 | 0.3 (5.5) | 0.44 | 234 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 70 | 281499025 | 1.94702 | 653 / 285 | 7.3 / 6.7 / 5.3 | 0.3 (6.5) | 1.00 | 115 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 71 | 410371059 | 21.58303 | 4465 / 3287 | 8.3 / 4.9 / 14.8 | 0.2 (7.5) | 3.58 | 8 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 72 | 277803435 | 25.48647 | 1900 / 4551 | 8.1 / 5.2 / 15.8 | 0.2 (5.3) | 3.27 | 5 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 73 | 320204496 | 1.23833 | 359 / 267 | 6.8 / 6.1 / 6.0 | 0.2 (7.6) | 1.09 | 115 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 74 | 222188854 | 1.10723 | 9441 / 7081 | 8.6 / 9.0 / 11.2 | 0.2 (11.4) | 5.16 | 47 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 75 | 219309872 | 1.93365 | 1006 / 515 | 7.7 / 5.1 / 4.7 | 0.1 (4.1) | 1.19 | 72 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 76 | 403209095 | 13.94278 | 1352 / 785 | 6.2 / 6.1 / 5.1 | 0.0 (5.8) | 1.80 | 5 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 77 | 91514466 | 5.36315 | 552 / 231 | 7.1 / 5.4 / 6.2 | 0.0 (7.5) | 1.04 | 22 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 78 | 281499025 | 6.92094 | 1376 / 605 | 7.5 / 4.9 / 5.6 | -0.0 (5.6) | 1.45 | 31 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 79 | 79180679 | 7.57793 | 1051 / 336 | 7.2 / 4.8 / 4.5 | -0.0 (6.3) | 0.92 | 14 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 80 | 77111228 | 18.97224 | 607 / 289 | 7.2 / 5.3 / 6.7 | -0.1 (9.2) | 0.79 | 6 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 81 | 355687587 | 17.86975 | 520 / 252 | 7.1 / 5.1 / 5.1 | -0.1 (5.8) | 1.12 | 10 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 82 | 355545170 | 1.14091 | 975 / 470 | 7.1 / 5.4 / 5.9 | -0.1 (5.3) | 1.29 | 146 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 83 | 407734245 | 19.82299 | 9094 / 2444 | 8.6 / 6.0 / 5.8 | -0.2 (5.6) | 3.06 | 7 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 84 | 201287627 | 26.15103 | 1001 / 664 | 7.6 / 6.2 / 5.6 | -0.2 (14.3) | 1.87 | 8 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 85 | 7143366 | 24.25159 | 1335 / 1079 | 7.9 / 8.7 / 7.2 | -0.3 (13.8) | 1.94 | 4 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 86 | 139458319 | 1.81823 | 988 / 482 | 7.5 / 5.3 / 5.0 | -0.3 (5.9) | 1.31 | 47 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 87 | 2027133756 | 4.55326 | 10579 / 8888 | 6.4 / 5.4 / 6.3 | -0.3 (9.6) | 6.33 | 15 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 88 | 325385847 | 2.89962 | 9391 / 1957 | 8.5 / 4.3 / 4.4 | -0.3 (8.3) | 2.70 | 26 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 89 | 409858636 | 25.29284 | 11927 / 5072 | 8.6 / 10.7 / 10.2 | -0.3 (78.8) | 4.33 | 6 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 90 | 410443108 | 18.08708 | 4732 / 1654 | 8.7 / 6.2 / 7.6 | -0.3 (11.7) | 2.36 | 10 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 91 | 234298748 | 1.25764 | 466 / 202 | 6.4 / 4.8 / 7.8 | -0.4 (6.4) | 1.05 | 124 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 92 | 116660700 | 15.57057 | 2080 / 1452 | 7.5 / 6.0 / 6.7 | -0.4 (10.4) | 2.44 | 7 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 93 | 404452683 | 14.33965 | 1245 / 3732 | 9.0 / 5.5 / 12.8 | -0.4 (4.7) | 2.89 | 2 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 94 | 234298748 | 20.54759 | 908 / 5205 | 6.3 / 3.8 / 34.7 | -0.4 (5.8) | 5.34 | 5 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 95 | 219338378 | 11.61013 | 63470 / 53185 | 8.6 / 5.7 / 8.7 | -0.5 (8.9) | 10.58 | 9 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 96 | 355738167 | 19.18222 | 5769 / 18270 | 8.1 / 6.4 / 22.4 | -0.5 (8.7) | 3.99 | 6 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 97 | 213945195 | 12.03035 | 3259 / 1834 | 8.3 / 4.7 / 7.2 | -0.6 (4.9) | 2.14 | 9 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 98 | 319246937 | 13.90706 | 603 / 432 | 7.2 / 6.7 / 8.1 | -0.6 (12.2) | 1.29 | 9 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 99 | 122847583 | 10.06221 | 1084 / 400 | 7.7 / 6.5 / 4.8 | -0.7 (8.9) | 1.16 | 11 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 100 | 351339274 | 13.31529 | 429 / 649 | 8.4 / 7.4 / 10.5 | -0.9 (6.5) | 1.03 | 5 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 101 | 406277365 | 1.06096 | 197 / 182 | 7.0 / 6.7 / 4.9 | -0.9 (7.9) | 0.96 | 209 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 102 | 219121412 | 1.04081 | 273 / 272 | 6.4 / 6.3 / 5.5 | -0.9 (5.8) | 0.77 | 129 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 103 | 403208541 | 12.88082 | 721 / 556 | 6.9 / 5.9 / 7.3 | -0.9 (11.6) | 1.67 | 13 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 104 | 238813286 | 29.36765 | 3100 / 2659 | 8.2 / 7.5 / 11.0 | -1.1 (7.5) | 3.01 | 9 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 105 | 469874861 | 23.93180 | 5170 / 2550 | 8.7 / 6.1 / 5.8 | -1.1 (9.5) | 1.89 | 5 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 106 | 381923220 | 4.06671 | 587 / 284 | 6.2 / 5.0 / 7.9 | -1.1 (6.8) | 1.04 | 41 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 107 | 404452683 | 11.89583 | 1150 / 1611 | 8.5 / 5.0 / 7.4 | -1.1 (3.6) | 1.90 | 6 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 108 | 220378788 | 6.71582 | 712 / 232 | 6.4 / 4.8 / 5.3 | -1.1 (5.8) | 1.05 | 21 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 109 | 292330430 | 0.85293 | 983 / 472 | 7.8 / 5.9 / 5.5 | -1.2 (4.6) | 1.34 | 227 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 110 | 281656667 | 29.45000 | 4803 / 4559 | 8.3 / 5.2 / 10.7 | -1.2 (8.0) | 2.94 | 4 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 111 | 80431622 | 10.96278 | 44428 / 82940 | 8.2 / 6.3 / 7.7 | -1.3 (2.3) | 11.71 | 6 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 112 | 277775567 | 16.93540 | 812 / 810 | 8.1 / 4.9 / 5.7 | -1.3 (3.1) | 1.90 | 11 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 113 | 161384122 | 1.27387 | 275 / 128 | 6.2 / 5.5 / 4.5 | -1.4 (6.5) | 0.66 | 106 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 114 | 231268304 | 8.37350 | 5397 / 1264 | 7.8 / 5.8 / 5.0 | -1.4 (9.8) | 1.93 | 19 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 115 | 118586898 | 19.56828 | 2627 / 1807 | 7.9 / 4.8 / 5.2 | -1.7 (1.7) | 2.48 | 6 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 116 | 404289580 | 16.28371 | 1652 / 1074 | 8.1 / 6.8 / 8.5 | -1.7 (16.4) | 1.64 | 6 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 117 | 206354321 | 23.44207 | 465 / 212 | 7.4 / 6.3 / 6.9 | -2.0 (6.7) | 0.83 | 10 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 118 | 372106864 | 0.90214 | 635 / 298 | 7.1 / 6.3 / 4.8 | -2.4 (9.4) | 1.02 | 155 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 119 | 144003231 | 15.30868 | 1246 / 564 | 7.7 / 5.9 / 5.6 | -2.5 (11.1) | 1.37 | 10 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |
| 120 | 387263115 | 2.44871 | 492 / 344 | 6.7 / 5.3 / 4.9 | -2.8 (6.6) | 1.29 | 60 | - | - | none | not reproduced (hold-out inconsistent > 3 sigma) |

## 9. Novelty re-check for the RNAAS draft (2026-10-13, 21:28-22:27 UTC, live, read-only)

TIC 144439305 was checked by TIC ID, by Gaia DR3 ID, by 2′ position and by period aliases. Nothing claims or lists a transit signal on the star except Eumemic's raw, unvetted table.

- **Sources with no hit:**
  - ExoFOP TOI/CTOI lists and the target page;
  - the NASA Exoplanet Archive (`toi`, `ps`, `pscomppars`);
  - arXiv (metadata only);
  - **T16** (full MRT tables, 11,143 candidates; T16 does have a Sector 2 light curve of this star at MAST and did not list it);
  - **RAVEN** (all Zenodo tables);
  - the Th0rgal registry, brahle/exoplanet-db and 14 other community repos;
  - GitHub search, Zenodo, and PHT Talk (search plus all 1,590 "Planets!" threads);
  - SIMBAD/VizieR (no variability or EB entries; KELT lists it as non-variable);
  - MAST: no DV products, and no match in 307,754 SPOC TCEs (S1-107).
- **Eumemic** (commit da722613): still only three raw rows in `all_signals_faint_t1` (0.93981 d SNR 8.35; 1.9588 d SNR 6.92; 2.085 d SNR 4.79). The signal is in no vetted, candidate, tier, follow-up or ephemeris list.
- **Not covered:** arXiv/ADS full text, and the full per-file contents of brahle/exoplanet-candidates (clone incomplete).

**Corrections to this report, found while checking the draft:**
- **Spectral type:** about M2-M2.5V from BP−RP, M_G, M_Ks, Teff and R★ against the Mamajek table, not M3.5.
- **Gaia RVs:** Gaia DR3 gives RV summary statistics (24 epochs, constancy p = 0.95, robust peak-to-peak 9.3 km/s), not epoch RVs.
- **Teq:** 875 K assumes albedo 0.3. It is 956 K at zero albedo.
- **Eumemic "rank 1":** this rank is per star, not table-wide.
- **Star count:** the 7 tables cover 129,536 stars.
- **S105/106:** the cutouts are not sky-subtracted (FLUX_BKG = 0). Simple background-subtracted photometry gives 683 ± 138 ppm (z = 5.0), against 909 ± 137 (z = 6.6) from the PSF-scene aperture.
- **Radius:** an independent batman/emcee fit (HLSP S2-S104 plus TESScut S105/106) gives Rp = 1.23 (+0.12/-0.09) R⊕.

All of these are documented in the draft's checks.md.
