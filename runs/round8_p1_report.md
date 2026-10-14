# Round 8A: taking the P1 low-SNR candidates to a verdict

Branch `round8-p1` (worktree `.`, started from `round7-rnaas`), 2026-10-13, 19:00-22:00 EDT. Heavy compute ran at `nice -n 15`. Nothing was pushed, posted, emailed, uploaded or submitted. MAST, ExoFOP, Gaia, the NASA Exoplanet Archive, arXiv, GitHub, Zenodo and PHT Talk were only read.

## Bottom line

- **One candidate clears every gate: TOI-4405 at 3.6008 d, 1.35 R⊕.**
  - It is **reproduced in independent photons** with the ephemeris fixed in advance. TESS Year-1 FFI sectors S1-S13 were never part of the round-5 search, and they give 357 ± 75 ppm (z = 4.8) against 316 ± 43 ppm in the discovery data. That is 0.5σ apart, and 4.2 was expected if the signal is real. None of 2,000 random ephemerides in the same hold-out data got above z = 3.9.
  - **Pixels:** on target in 42 TESScut sectors (Δχ² 3.3, offset 3.8 ± 3.1″). Every neighbour that could mimic it is at Δχ² ≥ 44.
  - **FPP:** TRICERATOPS gives **0.014 ± 0.003, NFPP 0**, with the SOAR speckle curve (digitised from ExoFOP's PDF), MOLUSC (Gaia + RUWE + SOAR) and pixel-cleared neighbours. That is just under the 0.015 bar, with run-to-run spread 0.010-0.017.
  - **Dynamics:** a 1.09 / 3.60 / 5.92 d three-planet system is stable (SPOCK p = 0.96).
  - **Not strictly novel.** The signal is in Eumemic's public raw table (rank 1, SNR 9.51), and his vetting **failed** it ("inconsistent transit SNRs"). We do not reproduce a real inconsistency (per-epoch χ² p = 0.025, all six TESS years positive), but any note must credit him and say so.
- **Watch, not claimable yet:** TOI-1846 0.944 d (0.72 R⊕) and TOI-6008 7.220 d (0.92 R⊕).
  - Neither star has any TESS data that the discovery search did not use, so the hold-out test cannot be run.
  - The look-elsewhere test (the full round-5 search repeated on 60 scrambled copies of each star) puts the per-star false-alarm probability at about 0.2-3%. None of the 60 nulls reached either candidate's SNR. Across the 184 stars searched, that means **~0.4-6 noise signals this strong are expected**, depending on the tail fit.
  - **TOI-1846** is clean otherwise: on target, and FPP 0.011 / NFPP 0 with three contrast curves + MOLUSC.
  - **TOI-6008**'s pixels lean **off target**: TESScut localisation Δχ² 23 (control 4), robust to dropping any sector, toward a crowded patch 20-35″ west. Its TRICERATOPS NFPP stays at about 0.008 even with the AO curve and MOLUSC (FPP ≈ 0.04).
  - **TESS returns to both:** TOI-1846 in S117-S132 and TOI-6008 in S120-S128 (2027). Those are the real tests, so preregister the ephemerides now.
- **TOI-4342 1.958 d: watch, and not ours to claim.** Nothing new exists for this star (no TESS data after S103 and none planned through S134). The round-6 verdict stands: FPP 0.0065 / NFPP 0.00097 with MOLUSC, loosely on target. The signal is in Eumemic's raw table (SNR 8.90) and his vetted table, where it **failed** ("inconsistent transit SNRs").
- **Rejected:**
  - **TIC 11528365** (8.17 d) is refuted by the hold-out: −29 ± 70 ppm in S3/S30/S97/S108 against 786 ± 125 ppm in its discovery sector, 5.7σ apart.
  - **TIC 350522259** (0.739 d) fades: the second half of the data gives z = 1.2 against 6.9 in the first, the per-epoch χ² p is 0.002, the hold-out is 2.7σ low, and the pixels are inconclusive.
  - **TOI-249** (19.09 d) is consistent with noise: scrambled nulls beat its SNR 2 times in 20, it fits as grazing (b ≈ 0.92), FPP is 0.18 even with every constraint, and **Eumemic published it on 10-08**.
- **Spectral types re-derived** (median of the nearest Mamajek row in Teff, Bp−Rp, M_G, M_Ks and R★):
  - TOI-1846 is M2.5V, not M3.
  - TOI-6008 is M4.5V, not M5.
  - TOI-4342 is K9V-M0V.
  - TOI-249 is K8V.
  - TOI-4405 is K6V.
  - TIC 350522259 is M3.5V.
  - TIC 11528365 is G8V (range G2-G9: its TIC Teff is hotter than its absolute magnitudes).

## Ranked: worth an RNAAS?

1. **TOI-4405 3.601 d: yes, as a candidate note crediting Eumemic. Ricardo's call.**
   - It is the only P1 signal reproduced in photons the discovery never touched.
   - The pixels are clean, FPP is at the bar, and the system is stable.
   - **What weakens it:**
     - Eumemic's public raw signal and failed vetting have to be acknowledged.
     - The SOAR curve is our digitisation of a PDF, not the observer's numbers. Ask the observer (Ziegler) for the file, or say so.
     - FPP 0.014 sits right at the bar and depends on MOLUSC's RUWE limit: without MOLUSC it is 0.045.
     - No ground transit exists: none of the 5 public SG1 light curves covers a predicted transit, and 330 ppm is too shallow for 1-m class telescopes anyway.
   - **Before submitting:** a CTOI-quality ephemeris and a preregistered transit list.
2. **TOI-1846 0.944 d: not yet.**
   - Formally it validates (FPP 0.011, NFPP 0), but the signal's reality rests on one dataset, and the survey look-elsewhere expects about 2-6 signals this strong from noise.
   - Preregister now (Zenodo), then test against S117-S132 in 2027.
   - The note would also need Eumemic credit: it is in his raw tables at SNR 8.25.
3. **TOI-6008 7.220 d: no.**
   - It has the same no-hold-out problem, plus a pixel test that disfavours the target.
   - It is also in Eumemic's raw tables, and in a public BLS output (JorgeAcebes).
   - Revisit after S120-S121.
4. **TOI-4342 1.958 d: no** (nothing has changed since round 6). It would at most be a sentence in a TOI-4342 18.9-d note, credited to Eumemic.
5. **TOI-249, TIC 350522259, TIC 11528365: no.**

## Method (pre-registered before any new data was loaded)

- **Ephemerides** are the discovery values in `round8/targets.py`, committed before any download (`3274f53`). Nothing was re-fitted before the hold-out depth was measured.
- **Data:** every TESS product on MAST as of 2026-10-13 (`round8/inventory.py`; tables in `runs/round8_p1/inventory.json`):
  - SPOC 2-min PDCSAP, TESS-SPOC, QLP;
  - TGLC / eleanor-lite where nothing else exists;
  - 11×11 TESScut cutouts for every sector (own fixed-aperture photometry with per-cadence sky subtraction, since S105+ cutouts have none).
  - **No light-curve product or cutout after S108 exists for any of the seven stars.**
  - **Hold-out data** = sectors the discovery search never saw:
    - round-5 hosts: FFI-only sectors and anything after S107;
    - TIC 350522259: everything except the 27 TESS-SPOC sectors round 5 used;
    - TIC 11528365: everything except S107.
  - FFI products of sectors that also have 2-min data share photons, so they were used only as product checks.
- **Estimator** (`round8/lcgate.py`):
  - per sector and product, a wotan biweight with candidate and known transits masked;
  - the depth is the weighted least-squares amplitude of a box whose per-cadence value is the exposure fraction inside transit, so 2-min, 200-s, 10-min and 30-min data measure the same quantity;
  - errors are robust scatter × a red-noise β at the transit duration.
- **Hold-out verdicts:**
  - "reproduced" = z_hold ≥ 2.5 and within 3σ of the discovery depth;
  - "refuted" = more than 3σ below it;
  - "no power" = expected z < 2.5.
- **Checks:**
  - odd/even; phase 0.5 (P/2); phases 1/3 and 2/3 (P/3); 3P residue classes;
  - an any-phase secondary scan;
  - a fine period scan (is the fixed P the local optimum?);
  - per-epoch, per-sector, per-year and per-half consistency;
  - 2,000 random-ephemeris nulls, on all the data and on the hold-out;
  - batman + emcee fits, free and with a TIC-density prior (R★ error propagated);
  - product depth ratios against the host's known planet as an on-target control.
- **Look-elsewhere for stars with no hold-out** (`round8/scramble.py`):
  - The identical round-5 search (0.5-30 d) is run on null copies of each star's 2-min data.
  - Every known and candidate transit cadence is removed, then the remaining flux is circularly shifted by a random lag inside each sector, which keeps the noise and systematics and destroys coherence.
  - The strongest signal that passes the round-5 automatic triage is recorded.
  - A first version inverted the flux, but that turns flares into "transits"; it is kept separately and not used.
- **Pixels:**
  - `deep_vet.py` (TESScut PSF-scene difference-image localisation against Gaia DR3), every sector, with the host's known planet as control.
  - A 2-min-TPF template-shift test (`round8/tpfloc.py`, the round-6 method generalised) works on TOI-4342 (T = 11.0). At T ≈ 12.3 it fails its own control (the known planets' difference images reach only SNR 1-3 per sector), so it was **not** used for TOI-6008 or TOI-4405.
- **FPP:**
  - TRICERATOPS (`round6/run_tri.py`, 3 runs × 10⁶ draws each);
  - SPOC 2-min apertures;
  - every machine-readable ExoFOP contrast curve (and the digitised SOAR curve for TOI-4405);
  - MOLUSC (Gaia contrast + RUWE + the deepest AO curve, 5×10⁵ companions);
  - neighbours excluded by the pixels set to zero depth.
  - A bug fix: MOLUSC's CLI parsed a negative declination as an option; it is now passed after `--`.
- **Other checks:**
  - Gaia DR3 RV summary statistics against the K that a 0.08 M☉ eclipsing companion would induce;
  - SPOCK for multis;
  - every public SG1 ground light curve on ExoFOP, scanned for chance coverage of a candidate transit (`round8/ground.py`).
- **Novelty:** a separate read-only pass (section 3).

## 1. Summary table

Rp is from the limb-darkened fit with the TIC density prior (16/50/84%, R★ error included). FPP/NFPP is the most constrained configuration. "Null p" is the per-star probability that the full round-5 search returns a triage-passing signal at least as strong from scrambled noise (Gumbel tail, 60 nulls; 20 for TOI-249).

| Candidate | P (d) | Depth (ppm, all) | Rp (R⊕) | Hold-out | Pixels | FPP / NFPP | Null p (per star) | Novelty | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| **TOI-4405** | 3.600797 | 326 ± 37 | **1.35** (1.21-1.50) | **REPRODUCED** z = 4.8 (exp. 4.2) | on target (Δχ² 3.3; neighbours ≥ 44) | **0.014 ± 0.003 / 0** | 0/60; ~0.006-0.008 (the hold-out makes it moot) | Eumemic raw (failed his vetting) | **STRONG CANDIDATE** |
| TOI-1846 | 0.944057 | 264 ± 32 | 0.72 (0.66-0.79) | no data exists | on target (Δχ² 0.1; neighbours ≥ 39) | 0.011 ± 0.001 / 0 | 0/60; ~0.01-0.03 | Eumemic raw | WATCH (no independent data) |
| TOI-6008 | 7.21999 | 1297 ± 178 | 0.91 (0.84-0.99) | no data exists | **disfavoured** (Δχ² 23 vs control 4) | ≈ 0.04 / 0.008 (AO + MOLUSC, no clearing) | 0/60; ~0.002-0.017 | Eumemic raw + JorgeAcebes BLS output | WATCH, leaning FP |
| TOI-4342 (1.958 d) | 1.957687 | 294 ± 37 | 1.10 (1.01-1.20) | no data exists | loose (round 6/8) | 0.0065 / 0.00097 (round 6) | 1/60 (at 27 d); ~0.02-0.03 | Eumemic raw + vetted (failed) | WATCH, not novel |
| TIC 350522259.01 | 0.738772 | 1081 ± 161 | 1.04 (0.95-1.12) | 2.7σ low (z 0.9, exp. 4.1) | inconclusive | 0.27 / 0.093 (round 5) | – | novel | REJECT (fades) |
| TIC 11528365 | 8.1715 | – | – | **REFUTED** (5.7σ) | dip absent in pixels | – | – | novel | REJECT |
| TOI-249 (19.09 d) | 19.09278 | 1118 ± 134 | grazing (b 0.92) | no data exists | on target | 0.18 / 0 | 2/20; ~0.08 | **Eumemic candidate (10-08)** | REJECT / not ours |

## 2. Per-candidate verdicts

### 2.1 TOI-4405 (TIC 177254988), 3.6008 d: STRONG CANDIDATE

**Star.** K6V (K5-K6 across indicators). Teff 4,361 ± 128 K, R★ 0.687 ± 0.059 R☉, M★ 0.68 ± 0.08 M☉, d = 156 pc, T = 12.38. RUWE 0.95, no Gaia NSS or EB entry, TIC contamination 0.011.

**Known system.** TOI-4405 b = .01 at 5.921 d (confirmed in `pscomppars`). TOI-4405.02 at 1.093 d (PC).

**Data.**
- 31 SPOC 2-min sectors (S27-S98) are the discovery set.
- 11 FFI-only Year-1 sectors (S1, S3-S11, S13; TESS-SPOC 30-min) are the hold-out.
- S12 and S32 exist only as TESScut, and both cutouts are unusable.
- No data after S98 exists, and none is planned through S134.

**Signal.**
- **Hold-out** (pre-registered): 357 ± 75 ppm, z = 4.8 (4.2 expected) against 316 ± 43 ppm in the discovery data. Consistent at 0.5σ.
  - Random-ephemeris null on the hold-out: max z 3.9 in 2,000 trials.
  - 10 of 11 hold-out sectors are positive.
  - Plot: `runs/round8_p1/toi4405_holdout_fold.png`.
- **All data:** 326 ± 37 ppm, z = 8.8, 242 transits.
- **Ephemeris:** P = 3.600797 ± 0.000018 d, T0 = BJD_TDB 2459053.1815 ± 0.0035.
- **Fit:** Rp = 1.35 (1.21-1.50) R⊕, b = 0.53 (0.28-0.67), T14 = 1.87 h. Free-density fit ρ = 3.1 (1.0-5.4) against 2.95 g/cc from the TIC.
- **Derived:** a = 0.0405 au, a/R★ = 12.7, S ≈ 94 S⊕, Teq ≈ 870 K (A = 0).

**Gates.**
- **Light-curve checks:**
  - odd/even 0.4σ; phase 0.5: 19 ± 37 ppm;
  - P/3 excluded at 7.0σ; all three 3P residue classes positive (z 4.9-5.5);
  - secondary scan maximum z 2.6 (scan std 1.03);
  - P is the local optimum, and the next peak is at z 3.7.
- **Consistency:**
  - per-year χ² p = 0.43, with all six TESS years positive;
  - halves z 7.5 and 4.8;
  - per-sector p = 0.055 (37 of 42 positive);
  - per-epoch χ² p = 0.025 (mild overdispersion; the strongest single transit is z 3.9, and dropping the top 5% leaves z 6.8).
  - This is the only hint of the "inconsistent transit SNRs" that failed Eumemic's vetting. It is not significant here.
- **Products:**
  - TESS-SPOC/2-min ratio 1.06 ± 0.24 (control TOI-4405 b: 1.02);
  - TESScut 0.99 ± 0.23 (control 0.95);
  - QLP 0.79 ± 0.19 (control 0.95).
  - Round 5's "fixed aperture only 0.58×" was a small-sample effect: with 42 sectors the fixed aperture gives 296 ± 51 ppm.
- **Pixels** (deep_vet, 42 sectors):
  - offset 3.8 ± 3.1″, target Δχ² 3.3 (control b 13.9);
  - the T = 19.1 star at 11.7″ (needs a 9% eclipse) is at Δχ² 44.5, and the T = 16.1 star at 38.3″ (needs 0.3%) at 63.
- **Gaia RV:** 17 epochs, constancy p = 0.53, robust peak-to-peak 12.1 km/s. A 0.08 M☉ eclipsing companion would give K = 13.3 km/s. With per-epoch σ of ~4 km/s that would have failed the constancy test, so an on-target stellar EB is excluded.
- **TRICERATOPS:**

  | Configuration | FPP | NFPP |
  |---|---|---|
  | base | 0.115 | 0.007 |
  | MOLUSC only (Gaia + RUWE; no curve, no clearing) | 0.030 | 0.009 |
  | SOAR curve + stars beyond 10″ cleared | 0.045 | 0 |
  | + MOLUSC (Gaia + RUWE + SOAR) | **0.0138 ± 0.0028** (runs 0.014, 0.010, 0.017) | 0 |

  The residual is DTP 0.008, STP 0.007 and BEB 0.006.
- **Stability:** SPOCK p(stable) 0.96 (16-84%: 0.95-0.97, 151 draws). Mutual Hill spacings are 39 (.02-cand) and 13.5 (cand-.01). The period ratios are 3.29 and 1.644, so it is not in resonance.
- **Ground:** none of the 14 public SG1 tables (5 nights, all on .01/.02) covers a predicted candidate transit.

**Verdict.** It passes every gate. FPP sits at the bar, and that rests on two inputs: MOLUSC's RUWE limit and a digitised speckle curve. **It is the only P1 candidate worth an RNAAS.**

### 2.2 TOI-1846 (TIC 198385543), 0.9441 d: WATCH

**Star.** M2.5V by every indicator (round 5 said M3). Teff 3,512 K, R★ 0.411 R☉, M★ 0.404 M☉, 47.2 pc. RUWE 1.30.

**Data.**
- 29 SPOC 2-min sectors (S17-S86) are all of the data. Every FFI sector is the same photons.
- TESScut lists S16, but the cutout is empty: the star falls off the detector.
- **So no hold-out exists.** TESS returns in S117-S119, S122-S125 and S128-S132 (tess-point).

**Signal.** 264 ± 32 ppm, z = 8.2, 630 transits.
- **Light-curve checks:**
  - odd/even 0.6σ; phase 0.5 −1 ± 32 ppm; P/3 excluded at 6.8σ; 3P classes all positive;
  - per-year p = 0.72; per-sector p = 0.30 (27 of 29 positive);
  - halves z 6.1 and 5.5.
- **Ephemeris:** P = 0.9440574 ± 0.0000023 d, T0 = 2459649.7122 ± 0.0014.
- **Fit:** Rp = 0.72 (0.66-0.79) R⊕, b = 0.78 with the density prior. The free fit gives ρ = 21 (5-34) against 8.2: short and grazing-ish, as round 5 noted.
- **Products:** TESS-SPOC 1.00 (control 1.03), TESScut 0.99 (control 1.06).

**Pixels.** On target: offset 0.6 ± 1.8″, Δχ² 0.1. The nearest mimic-capable star (T = 17.3 at 18.4″) is at Δχ² 39.

**Gaia RV.** It excludes an on-target stellar EB: peak-to-peak 11.8 km/s against K(0.08 M☉) = 28 km/s.

**TRICERATOPS.**

| Configuration | FPP | NFPP |
|---|---|---|
| base | 0.167 | 0.006 |
| 'Alopeke 832 nm + PHARO Brγ + 'Alopeke 562 nm | 0.024 | 0.009 |
| + stars beyond 15″ cleared | 0.015 | 0 |
| + MOLUSC | **0.011 ± 0.001** | 0 |

**Ground.** One MuSCAT2 night (2022-05-03) covers a transit, but its precision is ±1,500-3,000 ppm against 264 ppm needed. Useless.

**Look-elsewhere.** In 60 scrambled nulls the strongest triage-passing signal reached SNR 7.97 against the candidate's 8.57 (none exceeded it). The Gumbel tail gives 0.9% (same period range) to 3.2% (any period) per star, or **~1.7-6 expected noise signals at this strength across the 184-star survey**.

**Verdict.** Statistically "validated" by FPP but **not established as real**, since there are no independent photons. Watch, and preregister for 2027.

### 2.3 TOI-6008 (TIC 286201103), 7.220 d: WATCH, leaning false positive

**Star.** M4.5V (round 5 said M5). Teff 3,098 K, R★ 0.235 R☉, M★ 0.205 M☉, 23 pc, T = 12.29. RUWE 1.28.
- Crowded galactic-plane field: TIC contamination ratio 0.79.
- Palomar AO reports a companion at 7.18″, ΔBrγ = 4.98.

**Data.** 9 SPOC 2-min sectors are all of the data, so no hold-out exists. TESS returns in S120, S121, S127 and S128.

**Signal.** 1297 ± 178 ppm, z = 7.3, 33 transits.
- **Light-curve checks:**
  - every sector positive (p = 0.86); per-year p = 0.59; per-epoch p = 0.65;
  - odd/even 1.1σ; P/3 excluded at 5.8σ.
- **Fit:** Rp = 0.91 (0.84-0.99) R⊕, b = 0.56. The free-fit ρ is 27 against 22 from the TIC.
- **Ephemeris:** P = 7.219997 ± 0.000109 d.

**Pixels: the problem.**
- TESScut (all 9 sectors): target Δχ² 23.4 against 3.8 for the control TOI-6008 b. The best position is 20.9 ± 10.1″ west, and every jackknife (drop one sector) gives −17″ to −25″.
- Candidate sources there:
  - two T = 18.0-18.1 stars at 31-33″ (Δχ² 3-4; they would need 25-27% eclipses);
  - a T = 14.6 star at 35.6″ (Δχ² 8; would need a 1% eclipse).
- The TPF test can't decide: its own control fails at this brightness.
- **Products:** TESScut/2-min 0.51 ± 0.15 against 0.61 for the control (consistent). The TPF aperture series (core → 3.5 px) rises from 0.74 to 1.88 relative to the control, which is what an off-target source does, but only at < 2σ.

**Ground.** 41 public SG1 tables from 10 nights. They target TOI-6008 b and TFOP's internal "TOI-6008.22", a 16.36-d SHERLOCK signal, not ours. None covers a predicted 7.22-d transit.

**TRICERATOPS.** We did not clear neighbours, because the pixels do not exclude them.

| Configuration | FPP | NFPP |
|---|---|---|
| base (mean of runs completed so far: 0.106, 0.097) | ≈ 0.10 | 0.007 |
| PHARO Brγ + MOLUSC | ≈ 0.04 (run 1: 0.044) | 0.008 |

Runs were stopped early to free the machine: 2 of 3 runs (base) and 1 of 3 (AO + MOLUSC); the per-run values are in the `tri/tri_toi6008_*.log` files. Run-to-run scatter on this star is about 0.01. **NFPP stays near 0.008**, above the 0.001 bar: the neighbours carry real probability.

**Look-elsewhere.** Null maximum SNR 7.76 in 60 nulls against the candidate's 7.93 (none exceeded it): 0.2% (same range) to 1.7% (any period) per star, or ~0.4-3 expected across the survey.

**Verdict.** It has no independent confirmation, and the pixels point off target in a crowded field. Watch until S120-S121.

### 2.4 TOI-4342 (TIC 354944123), 1.958 d: WATCH, not novel

**Nothing new exists.** The 9 SPOC sectors (S13-S103) are all of the data, and TESS does not return through S134.

**Our re-measurement agrees with round 6.**
- 294 ± 37 ppm, z = 8.0; Rp = 1.10 (1.01-1.20) R⊕.
- Odd/even 0.4σ, phase 0.5 at 1.6σ, per-sector p = 0.89 (9 of 9 positive).
- P = 1.9576868 ± 0.0000081 d.

**Pixels.** The round-8 TPF template test on the same pixels finds it near the target (Δχ² 4.8, or 8.7 at the control planet's position). But mimic-capable neighbours sit at only Δχ² 12-17, as round 6 found (24.7). "Loosely on target."

**FPP.** Round 6: 0.0065 / 0.00097 with MOLUSC. MOLUSC's RUWE module again refused this star, so round 8 adds nothing here.

**Look-elsewhere.** 1 of 60 nulls exceeded the candidate, with a 27-d signal; none did in the same 0.5-6 d range. That is about 2-3% per star.

**Novelty.** In Eumemic's `all_signals_toi_hosts` (SNR 8.90) and `vetted_signals_toi_hosts`, where it **failed** ("inconsistent transit SNRs"). Our per-epoch χ² p is 0.052.

**Verdict.** Unchanged from round 6. Watch; credit Eumemic if it is ever mentioned.

### 2.5 TIC 350522259.01, 0.739 d: REJECT

**Star.** M3.5V (M3.5-M4), T = 14.46, 90 pc, RUWE 1.05.

**Hold-out.** QLP S31/68/95/97/98, TGLC S1/4/11 and TESScut S90 give 305 ± 326 ppm against 1333 ± 186 in the discovery data. That is 2.7σ low, where 4.1σ was expected.
- The TESScut-only sectors S27, 66, 87-89, 93, 94 and 96 are too noisy at T = 14.5 for the 3×3 aperture: scatter above 2% per cadence, so they were dropped.

**Signal behaviour.**
- The first half of the data gives z 6.9 and the second half z 1.2; the per-epoch χ² p is 0.002.
- The random-ephemeris null on all data has std 1.54, with a maximum of 22.6 (heavy-tailed noise), so the full-data z of 6.7 is only empirical p = 0.003.

**Pixels.** The fit again prefers an empty position 118″ away (Δχ² 78 at the target).

**Verdict.** A fading signal in non-Gaussian noise. Reject.

### 2.6 TIC 11528365, 8.17 d: REJECT

G8V (G2-G9), T = 10.70.

**Hold-out.** S3 (TESS-SPOC), S30 and S97 (2-min) and S108 (TESScut) give −29 ± 70 ppm. That is 5.7σ below the 786 ± 125 ppm of the S107 discovery, where z = 11 was expected. Per-sector χ² p = 10⁻⁷.

**Pixels.** All-sector fixed-aperture depth 178 ppm.

**Verdict.** A one-sector event. Refuted, as round 4 already suspected.

### 2.7 TOI-249 (TIC 179985715), 19.09 d: REJECT (and not ours)

**Star.** K8V (K7-K9), RUWE 2.03: an unresolved companion is likely.

**Data.** Every 2-min sector, including S105-S107, was already in round 5, so no hold-out exists and TESS does not return.

**Signal.** 1118 ± 134 ppm in 9 transits.
- With the density prior the fit wants b = 0.92 (grazing).
- The free fit leaves Rp unconstrained (2.1-13 R⊕).

**Look-elsewhere.** Scrambled nulls reached SNR 10.1 and beat the candidate's 8.66 in 2 of 20 runs: about 8% per star, ~15 expected across the survey.

**TRICERATOPS.** Zorro 562 + 832 nm curves, stars beyond 25″ cleared, plus MOLUSC: FPP 0.18 ± 0.02 (SEB 0.06, EBx2P 0.04).

**Novelty.** Eumemic published it as a candidate on 2026-10-08 (commit `b6d84c88`); it is now his tier C ("grazing").

**Verdict.** Reject.

## 3. Novelty gate (live, read-only, 2026-10-13 23:09 UTC to 10-14 ~00:30 UTC)

Full per-source records are in `runs/round8_p1/novelty/novelty.json` (with file and line references) and `src_*.json`.

**Sources read:**
- ExoFOP TOI (8,151 rows) and CTOI (5,138), plus the 7 target pages;
- NASA Exoplanet Archive (`ps`, `pscomppars`, `toi`);
- arXiv API (24 queries) and the LaTeX source of 6 relevant papers;
- all 132 SPOC DV tcestats CSVs (341,351 TCEs). The latest multi-sector run is still S1-S96; single sectors reach S107;
- Eumemic `da722613` (all tables, including the gz raw ones, and issues #1-17);
- brahle/exoplanet-candidates `702e6fa9` and exoplanet-db `3cd5afa1` (including its public API);
- Th0rgal `19e5e9e8`;
- 10 other community repos, T16 (all ApJS MRT tables), RAVEN (Zenodo 19661443);
- PHT Talk (29 queries), Zenodo (30), GitHub code search (34).

| Candidate | Verdict | Decisive hit |
|---|---|---|
| TOI-4405 3.601 d | CREDIT-NEEDED | Eumemic `all_signals_toi_hosts` L4600 (rank 1, P 3.60082524, SNR 9.51); `vetted_signals_toi_hosts` L578 **failed** ("inconsistent transit SNRs", MES 10.07). RAVEN's other signals on the star (1.974 d, 1.093 d) are not aliases. |
| TOI-1846 0.944 d | CREDIT-NEEDED | Eumemic `all_signals_faint_t1` L17644 (SNR 8.25) and `all_signals_toi_hosts` L6953 (SNR 8.26); raw only. |
| TOI-6008 7.220 d | CREDIT-NEEDED | Eumemic raw (SNR 7.46 / 6.79). JorgeAcebes/tess-exoplanet-detection has a BLS result at 7.22016 d (2026-01-08), mislabelled "TOI-6008 b" and not claimed. |
| TOI-4342 1.958 d | CREDIT-NEEDED | Eumemic raw (SNR 8.90) and vetted (**failed**). |
| TOI-249 19.09 d | NOT NOVEL | Eumemic candidate since 2026-10-08, now tier C. **Correction for the KB:** it is listed but not on his *counted* list. |
| TIC 350522259 | NOVEL | nothing anywhere |
| TIC 11528365 | NOVEL | nothing anywhere (his full-sky raw signals are not public) |

**Limits:**
- arXiv/ADS full text was not searched beyond six papers.
- QLP has no public TCE list.
- GitHub code search is unreliable for numbers.
- Zenodo and PHT Talk were searched by metadata and keyword only.
- X, Reddit and Bluesky were not checked.

## 4. If TOI-4405 goes to an RNAAS: the method papers it would have to cite

ExoFOP accepts a community candidate from an RNAAS only if every detection and vetting method is published in a refereed journal and cited. This is the list for what we actually used:

| What | Reference |
|---|---|
| TESS mission | Ricker et al. 2015, JATIS 1, 014003 |
| SPOC pipeline, 2-min light curves | Jenkins et al. 2016, SPIE 9913, 99133E |
| PDC detrending | Smith et al. 2012, PASP 124, 1000; Stumpe et al. 2012, PASP 124, 985; Stumpe et al. 2014, PASP 126, 100 |
| TESS-SPOC FFI light curves (the hold-out data) | Caldwell et al. 2020, RNAAS 4, 201 |
| QLP (product check) | Huang et al. 2020, RNAAS 4, 204 and 206; Kunimoto et al. 2021, RNAAS 5, 234 |
| TESScut cutouts | Brasseur et al. 2019, ASCL 1905.007 (software record, not refereed). TESScut is used only for pixel vetting. |
| Detection: BLS | Kovács, Zucker & Mazeh 2002, A&A 391, 369 |
| Detrending: wotan | Hippke et al. 2019, AJ 158, 143 |
| Difference-image localisation | Bryson et al. 2013, PASP 125, 889 (plus Twicken et al. 2018, PASP 130, 064502 for the DV centroid tests) |
| Statistical validation: TRICERATOPS | Giacalone et al. 2021, AJ 161, 24 (also defines the FPP < 0.015 / NFPP < 0.001 bar) |
| Bound-companion limits: MOLUSC | Wood et al. 2021, AJ 162, 128 |
| Speckle imaging used (SOAR HRCam) | Tokovinin 2018, PASP 130, 035002; Ziegler et al. 2020, AJ 159, 19 (TESS SOAR survey). Ask the observer before publishing with the curve. |
| Transit model and sampler | Kreidberg 2015, PASP 127, 1161 (batman); Foreman-Mackey et al. 2013, PASP 125, 306 (emcee) |
| Limb darkening | Claret 2017, A&A 600, A30. **We used approximate coefficients:** replace them with Claret's tables before quoting Rp in a paper. |
| Stellar parameters | Stassun et al. 2019, AJ 158, 138 (TIC 8); Gaia Collaboration, Vallenari et al. 2023, A&A 674, A1 (DR3); Pecaut & Mamajek 2013, ApJS 208, 9 (spectral type) |
| Gaia RV constancy | Katz et al. 2023, A&A 674, A5 |
| Stability | Tamayo et al. 2020, PNAS 117, 18194 (SPOCK); Rein & Liu 2012, A&A 537, A128 (REBOUND) |
| Mass-radius relation (stability only) | Chen & Kipping 2017, ApJ 834, 17 |
| Software | Astropy Collaboration 2013, 2018, 2022; Ginsburg et al. 2019, AJ 157, 98 (astroquery) |
| Prior public signal | Eumemic, github.com/eumemic/exoplanets (commit `da722613`, `all_signals_toi_hosts`), not refereed. Credit it in the text. |

**Gaps:**
- The **hold-out test and the scrambled-search null have no refereed method paper.** They are standard statistical practice, but they cannot be presented as the detection or vetting method of record. Keep BLS + TRICERATOPS + difference imaging as the methods of record and describe the hold-out as supporting evidence.
- If the multiplicity boost is invoked, cite Lissauer et al. 2012, ApJ 750, 112. We have not applied it.

## 5. Caveats

- **"Independent" means different photons, not a different instrument.** TOI-4405's hold-out is TESS Year-1 FFI data, and no ground transit exists.
- **For TOI-1846, TOI-6008, TOI-4342 and TOI-249 the requested hold-out test is impossible.** No TESS photon outside the discovery set exists today. The scrambled-search null substitutes for it as a look-elsewhere calibration. It controls noise and systematics, not astrophysical false positives.
- **Survey-level look-elsewhere:** at per-star false-alarm rates of 0.2-3%, the 184-star round-5 survey expects about 0.4-6 noise signals as strong as these four. Statistically, one or two of TOI-1846, TOI-6008, TOI-4342 d2 and TOI-4405 could be noise. TOI-4405 is the exception, because it passed a fixed-ephemeris test in new data.
- **TRICERATOPS FPPs assume the signal is astrophysical.** Three inputs are constructed rather than measured:
  - the pixel-clearing radii come from our TESS-resolution localisation;
  - the SOAR curve for TOI-4405 is our digitisation;
  - MOLUSC's RUWE limit refused TOI-4342, as in round 6.
- **The TPF template-shift method fails its own control at T ≈ 12.3** (the known planets' difference images reach only SNR 1-3 per sector). Pixel verdicts for TOI-6008 and TOI-4405 rest on the calibrated TESScut PSF-scene method.
- **Light-curve fits used approximate quadratic limb darkening by Teff band.** Radii could shift by a few percent with proper coefficients.
- **A process note.** The scrambled-search null went through two discarded versions. The first inverted the flux, which turns flares into transits. In the second, the stop left orphaned worker processes running, which wrote stale files. Both sets were moved aside: `scramble_v1_inverted/`, `scramble_v2a/`, `scramble_stale/`. The reported numbers come from a clean run.

## 6. Resources and files

**Downloads:**
- about 1.5 GB into `cache/r8`: light curves, ExoFOP SG1 tables, and 20 TPFs for TOI-6008/TOI-4405 (the TPFs were deleted afterwards; 0.6 GB remains);
- about 1 GB of new TESScut cutouts in the shared `cache/tesscut`;
- 1.3 GB of novelty sources in `<data>/r8/novelty/` (0.7 GB after cleanup).

**Disk:** free space on `/` fell from 44 GB to **24 GB**, below the 25 GB floor, about 21:00. About 4 GB of that was ours; the rest is other activity on the machine (the shared `cache/` is 23 GB and was about 18 GB at the start). We freed 1.6 GB of re-downloadable session data: the 20 TPFs of the failed TPF method, the SPOC TCE CSVs and the arXiv sources. That brought it back to 30 GB. No further downloads are planned. **Ricardo may want to check what else is filling the disk.** The bandwidth cap was lifted mid-session; at most 4 download workers were used.

**Code (`round8/`):**

| File | What |
|---|---|
| `targets.py` | pre-registered candidates |
| `inventory.py`, `fetch.py` | product inventory and downloads |
| `lcgate.py` | light-curve gate and hold-out test |
| `scramble.py` | look-elsewhere nulls |
| `tpfloc.py` | TPF localisation |
| `tri.py`, `run_tri.sh` | TRICERATOPS and MOLUSC |
| `stars.py` | stellar parameters, Gaia RV, Mamajek spectral types |
| `ground.py` | SG1 ground light curves |
| `stability.py` | SPOCK |
| `plot_fold.py` | hold-out fold plot |

Pixel specs are in `vetting/round8_targets.json`.

**Outputs (`runs/round8_p1/`):**

| Path | What |
|---|---|
| `<key>_lcgate.json/.log`, `<key>_fold.npz`, `<key>_ephem.json` | light-curve gate, folds, ephemerides |
| `pixels/` | deep_vet results |
| `tpf/` | TPF localisation |
| `tri/` | TRICERATOPS, MOLUSC, contrast curves |
| `scramble/`, `scramble_summary.json` | look-elsewhere nulls |
| `stars.json` | stellar parameters |
| `ground.json` | ground-data coverage |
| `toi4405_stability.json` | SPOCK |
| `novelty/` | the novelty pass |
| `toi4405_holdout_fold.png` | hold-out fold plot |
