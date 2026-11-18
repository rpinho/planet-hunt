![planet-hunt banner](assets/banner.png)

# planet-hunt: the record

**AI for Astronomy.**

The TESS planet candidates I've worked on with AI agents, in one place. Every section names the star, what we found, how it was checked, the verdict, and where it stands now. The signals remain **candidates**, not confirmed planets.

What was easy, what was hard, and what is still open: [RETRO.md](RETRO.md). Current status: [RESEARCH_STATUS.md](RESEARCH_STATUS.md).

Sister hubs: [AI for biology](https://github.com/rpinho/ai-biology) · [math](https://github.com/rpinho/math-certificates) · [physics](https://github.com/rpinho/ai-physics)

This version leaves out the TOI-4342 section, which is in review at RNAAS. It also leaves out the round-8 novelty table and caveats under TOI-4405, because they cover TOI-4342 too. The full file is [PAPER.md](PAPER.md).

## What this is, and what it isn't

None of our notes has been peer reviewed. Four Research Notes were sent to RNAAS in October 2026. RNAAS declined three of them under a [new policy](https://doi.org/10.3847/2515-5172/aeba86) for notes on transits in a single source or a small number of sources; the TOI-4342 note is still in review. The corrected versions on Zenodo are the stable record.

An uncatalogued signal or low statistical false-positive probability alone does not establish novelty, on-target origin or planetary status.

This file puts the record of every candidate in one place, for readers and reviewers who want it in a single document. Each section below is copied word for word from the public branch where the work was done: the claim, the method, the vetting checks, the verdict and the caveats. Only the heading levels are changed. The branch reports were written before the October 30 corrections, so where a number differs, the corrected Zenodo version in the status line is the one to use. Verdict labels (strong candidate, WATCH, REJECT) are the reports' own words.

Analysis and publication materials were prepared by AI agents under my direction.

---

# TOI-6000 c (TIC 259233660), 4.21 d

Branch: [round4-tic259233660](https://github.com/rpinho/planet-hunt/tree/round4-tic259233660)

## Status

Candidate. First reported publicly by Eumemic. RNAAS declined the note on November 14, 2026 under the new policy. Corrected version: [Zenodo 23895484](https://doi.org/10.5281/zenodo.23895484). On November 16 an [ExoFOP observing note](https://exofop.ipac.caltech.edu/tess/edit_obsnotes.php?id=259233660) on TOI-6000 pointed TESS follow-up observers to the corrected version.

## Verdict

_Source: [round4-tic259233660/vetting/round4_tic259233660_verdict.json](https://github.com/rpinho/planet-hunt/blob/round4-tic259233660/vetting/round4_tic259233660_verdict.json)_

```json
{
 "target": "TIC 259233660 (TOI-6000), candidate second signal",
 "verdict": "strong candidate - worth drafting a CTOI for Ricardo to decide; not a validated planet",
 "recommendation": "draft CTOI (decision is Ricardo's; nothing submitted)",
 "summary": "A 4.21344-d, ~790 ppm, 1.15 +/- 0.09 Earth-radius transit next to the USP planet TOI-6000 b. Detected in SPOC 2-min (SNR 8.8, null-calibrated), TESS-SPOC FFI and fixed-aperture TESScut photometry, localised on the target to 2.0 +/- 2.0 arcsec with every Gaia neighbour within 63 arcsec excluded (dchi2 >= 43), and its transit shape implies the same stellar density as TOI-6000 b (0.02 sigma). TRICERATOPS with the public Keck/NIRC2 contrast curves and the pixel-cleared neighbours gives FPP 0.0008, NFPP 0. It is not validated because the statistical tests assume the signal is real, all TESS photons were used by the discovery search, SNR is below 10 and the depth drifts down ~10%/yr at ~2 sigma.",
 "for": [
  "SPOC 2-min (33 sectors, TOI-6000 b masked): 781 +/- 88 ppm, SNR 8.8 with red-noise errors; 3000 random ephemerides give SNR std 1.006 and max 3.2.",
  "Robust: 105/149 transits positive (sign test z = 4.9 calibrated), Student-t depth 798 +/- 112 ppm (SNR 6.8 calibrated); look-elsewhere 1%-FA threshold ~6.4 sigma.",
  "Each half separately: Cycle 2 915 +/- 134 ppm (SNR 6.8), Cycles 4-7 ~668 +/- 115 ppm (SNR 5.8).",
  "Pixel localisation, 39 sectors: 2.0 +/- 2.0 arcsec (TOI-6000 b on the same pixels: 2.3 +/- 1.8 arcsec); TIC 1884419644 (13 arcsec) now disfavoured at dchi2 44.6 vs 1.6 for the target.",
  "Shape: density 7.1 +/- 4.8 g/cc free vs TOI-6000 b 7.0 (+2.3/-4.0) and TIC 9.6 +/- 0.9; shared-density joint fit costs dchi2 1.8; odd/even 0.15 sigma; secondary -168 +/- 88 ppm.",
  "Not an alias or b harmonic (P/Pb = 9.385); not momentum dumps (1 of 169 transits within one duration, 4 expected); no depth-contamination correlation (0.1 sigma).",
  "TRICERATOPS: FPP 0.089 (no constraints) -> 0.059 (+ Keck H/J contrast curves) -> 0.0008, NFPP 0 (+ pixel-cleared neighbours); 0.0014 on the older, smeared fold."
 ],
 "against": [
  "No independent photons: the discovery search used the TESS-SPOC FFI light curves of 39 sectors plus all 33 2-min sectors; Sector 108 does NOT observe the star (tess-point); next TESS coverage is Sectors 117-132 (2027).",
  "Depth declines ~10%/yr: 2-min -80 +/- 46 ppm/yr (1.8 sigma; null-calibrated p ~ 0.09 for the -1.3 sigma version), fixed aperture -154 +/- 56 ppm/yr (2.7 sigma); TOI-6000 b is flat in both products. Found after looking at season depths (a posteriori).",
  "Per-sector depth scatter chi2/dof 1.76 vs 1.05 +/- 0.28 at random ephemerides (p ~ 0.017; S40 high, S47 and S75 low).",
  "Fixed-aperture/2-min depth ratio 1.35 +/- 0.13 vs 1.09 +/- 0.06 for TOI-6000 b (1.8 sigma); direction not blend-like for a target-centred aperture.",
  "SPOC never raised a second TCE (expected MES ~6 < 7.1 threshold, so not surprising) and no published work searched for more planets.",
  "TRICERATOPS assumes the dip is astrophysical; neighbour clearing comes from our own TESS pixel fit, not seeing-limited ground photometry."
 ],
 "bars": {
  "ctoi_fpp_lt_0p5": "met under every configuration (0.0008-0.089; 0.27 in the earlier unconstrained run)",
  "validation_fpp_lt_0p015_nfpp_lt_0p001": "met formally only with the pixel-based neighbour clearing (FPP 0.0008, NFPP 0); without it NFPP = 0.058. Do not claim validation: signal reality is not yet independently established."
 },
 "literature": "Gomez Barrientos et al. 2025 (AJ 170, 148; arXiv 2508.02782) validate TOI-6000 b only; no RVs, TTVs, residual transit search or injection-recovery; no mention of a second signal. Exoplanet Archive and ExoFOP list one planet; SPOC TCEs: only planet 1 in every multi-sector run through s14-s86.",
 "next": [
  "Ground-based transit at the posterior ephemeris (P = 4.213442 d, T0 = BTJD 2634.7373; ~10-15 min 1-sigma in 2026-27): ~0.8 ppt, 1.5 h; Palomar/WIRC reached 172 ppm per 10 min on this star in J, so one transit would give ~10 sigma.",
  "TESS Sectors 117-132 (2027): an independent check; predicted transit near BTJD 4303.26 +/- 0.008 d.",
  "If Ricardo approves: CTOI on ExoFOP with these numbers and the caveats above."
 ]
}
```

## False-positive probability rerun

_Source: [round4-tic259233660/vetting/fpp_keck_2026-10-30/README.md](https://github.com/rpinho/planet-hunt/blob/round4-tic259233660/vetting/fpp_keck_2026-10-30/README.md)_

TRICERATOPS rerun for TOI-6000 c (2026-10-30), with the public Keck/NIRC2 contrast curves from ExoFOP.

Input: 33 SPOC 2-minute sectors, 2-minute phase bins, 120 s exposure, SPOC apertures from 12 sectors, Gaia background, default priors, no MOLUSC. Five seeded runs of 200,000 draws per configuration.

| Configuration | 14 neighbours kept: FPP, NFPP | Removed: FPP, NFPP |
|---|---|---|
| No imaging | 0.094 +/- 0.012, 0.048 | 0.045 +/- 0.005, 0 |
| Keck Hcont | 0.059 +/- 0.005, 0.059 | 0.0002 +/- 0.0001, 0 |
| Keck Hcont + Jcont | same as Hcont, same seeds | same as Hcont |

`tri6000.py` needs a TRICERATOPS install with list-valued contrast curves, the round-4 aperture and light-curve products, and `net.py` from the job directory. It records provenance; it is not a turnkey script.

---

# TIC 144439305, 0.94 d

Branch: [round7-rnaas](https://github.com/rpinho/planet-hunt/tree/round7-rnaas)

## Status

Candidate. RNAAS declined the note on November 15, 2026 under the new policy. Corrected version: [Zenodo 23836131](https://doi.org/10.5281/zenodo.23836131).

## Claim, method and verdict

_Source: [round7-rnaas/runs/round7_eumemic_weak_report.md](https://github.com/rpinho/planet-hunt/blob/round7-rnaas/runs/round7_eumemic_weak_report.md)_

#### Bottom line

- **Mined:** all 264,468 public BLS signals in Eumemic's seven `all_signals*.csv.gz` tables (K/M dwarfs, TOI hosts, faint M dwarfs T 11.5-13.5, the TCE re-check and the 20-100 d passes).
  - 35,407 passed the table-level cuts (SNR 6-9, at least 3 transits, small star, clean shape, not his and not catalogued).
  - We inventoried the products for 3,898 of them and selected **120 signals on 110 stars** with the most extra data.
  - **Extra data, median:** 3 sectors he did not have on top of his 3. Almost all of it is QLP sectors his reader dropped.
- **One strong new candidate: TIC 144439305**, the best result of this round.
  - **Planet:** P = 0.9398094 ± 0.0000023 d, 569 ± 58 ppm, Rp = 1.28 (1.20-1.36) R⊕, Teq ≈ 875 K.
  - **Star:** early M dwarf (about M2–M2.5; first estimated as M3.5, see the corrections at the end) at 62 pc, T = 12.15, 3,505 K, 0.46 R☉.
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

#### 3. The test (`round7/revet.py`), fixed before any result was seen

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

#### 5. TIC 144439305: what we know

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

#### 7. Caveats

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

#### 9. Novelty re-check for the RNAAS draft (2026-10-13, 21:28-22:27 UTC, live, read-only)

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

## False-positive probability rerun

_Source: [round7-rnaas/vetting/fpp_2026-10-30/README.md](https://github.com/rpinho/planet-hunt/blob/round7-rnaas/vetting/fpp_2026-10-30/README.md)_

Reruns for the TIC 144439305 correction (2026-10-30).

FPP (TRICERATOPS): selected 200-second cadence at its true exposure, original apertures and cleared-neighbour list, five runs per configuration. MOLUSC runs use five companion samples at N=1,000,000. `verified_padding_probes.csv` shows that N changes the result: TRICERATOPS zero-pads the surviving MOLUSC companions to N, so the draw count sets the companion prior. Values in the note are means and sample standard deviations of each `tri_*.json`.

Hold-out: `tic144_holdout.json` has the QLP comparison with the public orbit-edge mask (first 1.0 d and last 0.5 d of each orbit, Sectors 94-104), z calibrated against 2,000 random ephemerides.

## Hold-out rerun

_Source: [round7-rnaas/vetting/holdout_2026-10-31/README.md](https://github.com/rpinho/planet-hunt/blob/round7-rnaas/vetting/holdout_2026-10-31/README.md)_

Hold-out rerun for the 2026-10-31 correction.

Each hold-out here uses the ephemeris and duration fitted before its hold-out sectors entered any fit: mask-A null z, 2,000 random ephemerides. `*_pre.json` is the value in the note. `*_note.json` is the earlier value at the printed ephemeris, which was fitted with hold-out sectors, kept for comparison. Printed ephemerides are for predicting transits.

---

# TOI-4405 (TIC 177254988), 3.60 d

Branch: [round8-p1](https://github.com/rpinho/planet-hunt/tree/round8-p1)

## Status

Candidate, in the combined note. RNAAS declined the combined note on November 16, 2026 under the new policy. Corrected version: [Zenodo 23895440](https://doi.org/10.5281/zenodo.23895440).

## Claim, method and verdict

_Source: [round8-p1/runs/round8_p1_report.md](https://github.com/rpinho/planet-hunt/blob/round8-p1/runs/round8_p1_report.md)_

#### Method (pre-registered before any new data was loaded)

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

##### 2.1 TOI-4405 (TIC 177254988), 3.6008 d: STRONG CANDIDATE

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

---

# TIC 101653309 and TIC 262843927

Branch: [round8-crowd](https://github.com/rpinho/planet-hunt/tree/round8-crowd)

## Status

Candidates, in the combined note. The 14.52-day signal on TIC 262843927 was first reported by Kunimoto et al. (2025). RNAAS declined the combined note on November 16, 2026 under the new policy. Corrected version: [Zenodo 23895440](https://doi.org/10.5281/zenodo.23895440).

## Claim, method and verdict

_Source: [round8-crowd/runs/round8_crowd_report.md](https://github.com/rpinho/planet-hunt/blob/round8-crowd/runs/round8_crowd_report.md)_

#### Bottom line

- **Tested:** 4,026 unvetted or below-threshold public signals on 3,500 small stars, from Eumemic plus seven other community searches. Each was tested pre-registered on TESS data its searcher never used, with noise nulls, inverted light curves and 4,026 injections. This is the round-7 method at 33× the scale.
- **Survivors:** 17 of the 4,021 signals with usable data pass. **Noise alone predicts about 3:** the whole procedure run on inverted light curves passes 3 of 4,021.
- **The deep vet sorts the 17:**
  - 5 are stellar variability (oscillations or spots), not transits.
  - 4 are blended eclipsing binaries on a neighbouring star (the pixels localise the dip 14-70″ away).
  - 2 (two periods on one star) are V-shaped 1.5-2% dips on a star that supplies only 12% of the light in its pixels; they can't be localised.
  - 1 is an on-target grazing EB. TRICERATOPS gives FPP 0.97, the star has RUWE 2.4 and a Gaia acceleration solution, and Eumemic already lists it as tier C.
  - 2 are planet-like but the pixels can't separate them from faint neighbours.
  - **3 are planet-like, on target and, as far as every source we checked shows, unclaimed.**
- **Best lead: TIC 101653309 b?**, 0.7034 d, 1.30 (1.21-1.51) R⊕, around an M1.5-M2 dwarf (T = 11.94) at 60 pc.
  - It is the twin of round 7's TIC 144439305. The signal sat in Eumemic's raw table at SNR 8.1, below his SNR 9 pre-filter.
  - **Hold-out:** it reproduces in a QLP sector his reader dropped (S94: z = 4.9, null p < 1/300) and is consistent in TESScut S105, which no pipeline has processed (337 ± 165 ppm).
  - **Pixels:** on target. The best Gaia source is the target in 5 of 6 sectors; the fixed-aperture depth is 542 ± 54 ppm, constant across 6 sectors (p = 0.36).
  - **Binary checks:** Gaia RVs are constant (26 epochs, p = 0.40), which excludes an on-target eclipsing star.
  - **FPP:** TRICERATOPS gives FPP 0.080 / NFPP 0.040 by default, and 0.041 / 0.034 with Gaia/MOLUSC limits and the neighbours beyond 10″ cleared. The residual is one T = 18.0 star 5.9″ away, which the pixels disfavour at Δχ² 19 (about 4.3σ); clearing it too gives 0.0071 ± 0.0006 / 0. Validation therefore hinges on resolving that one star, which any seeing-limited ground photometry can do.
  - **Novelty:** clean everywhere, except that the raw signal sits unvetted in Eumemic's table.
- **Second lead: TIC 262843927**, 24.109 d, 2.1 R⊕, warm (S ≈ 3 S⊕, Teq ≈ 340 K), around an M1 dwarf (T = 11.98) at 83 pc. On target, FPP 0.027 ± 0.002 and NFPP 0.
  - The same star carries a **second signal at 14.517 d (1.9 R⊕)** in Eumemic's follow-up list (not in any tier). It also passes our hold-out (z = 3.7).
  - So this is probably a two-planet system near 5:3. The 24.1 d signal is in Eumemic's raw table only (SNR 8.5).
- **Third lead: TIC 410418450**, 6.820 d, 3.5 R⊕, around an M dwarf at 136 pc (T = 13.2). On target; FPP 0.013 ± 0.001, NFPP 0.0035.
- **The rest of the crowd's weak band is noise.**
  - The 3,458 Eumemic SNR 6-9 hold-out z-values have mean 0.06 and std 1.10.
  - If his claimed depths were real, they would give a median z of 4-5.5, and 50% would pass by the injections.
  - Correcting for completeness, about 99% of the tested discards are noise. Roughly 0.7% are real astrophysical signals, mostly variables and blends.
- **Needs Ricardo:** whether TIC 101653309 (and TIC 262843927) go the TIC 144439305 way: CTOI, RNAAS and preregistered times. No note is drafted. The citation list is in §8. Credit Eumemic in each case.

#### 3. The test (stages 3-5, `round8/holdout.py`), pre-registered

**Pre-registration.** The labelling rules, modes and survivor criteria were written in the docstring and committed (`bd1c893`, 19:54 EDT) before any round-8 hold-out result was computed.

**Labelling.** Each sector is marked as the source's data ("theirs") or new, by a per-source rule:
- Eumemic's early searches: identical to round 7.
- Eumemic's full-sky catalogue: also excludes empty-error QLP files and QLP-only S74-79 (his issue #9).
- Paras/Avijit: lightkurve TESS-SPOC > QLP with non-finite errors dropped; conservative when unsure.
- Huang, Marchand, Valdez, ExoNit, Rahle: their listed sectors.

**FIXED mode** (their data span ≥ 150 d; 3,975 signals). This is the round-7 procedure unchanged:
1. Narrow BLS on their data gives the ephemeris.
2. The depth at that ephemeris in new data alone gives z_hold.
3. **Survivor** = all of these:
   - z_hold ≥ 2.5, and the hold-out depth within 3σ of their-data depth;
   - SNR_full ≥ 8 within 0.1% of the source period;
   - |odd/even| < 3, secondary < 3σ, season χ² p > 0.001.
- **Noise:** 300 random ephemerides per signal.

**WINDOW mode** (single-sector sources, whose period can't be propagated over years; 46 signals with data):
- The period comes from their sector.
- A narrow BLS on new data within ±0.5%, with T0 free, gives z_hold.
- **Survivor** = z_hold ≥ 4 and above all 60 null windows of the same width, plus the same consistency and full-data cuts.

**Both modes:**
- The whole procedure is also run on the inverted light curve.
- **Injections** put each target's source depth at another period (×0.7-1.3) and run the identical procedure.

**Positive control.** Round 7's TIC 144439305 run through the new engine as an E1 row gives z_hold 5.16, depth 468 ± 91 ppm (identical to round 7), and passes.

##### TIC 101653309 (best)

**Star.**
- TIC 101653309 = Gaia DR3 6667277333597987456. RA 302.6893, Dec −49.2954.
- T = 11.94, G = 13.01, BP−RP 2.18.
- ϖ = 16.55 mas (60.3 pc). Teff 3,627 K, R★ 0.459, M★ 0.456 (TIC 8.2): about M1.5-M2V by Teff and radius. RUWE 1.20.
- Gaia variability class SOLAR_LIKE (6.7 mmag RMS spot modulation); TESS shows a ~7 d rotation, 580 ppm.
- Data: TESS-SPOC S13, S27, S67; QLP S94 (empty errors, dropped by Eumemic's reader) and S104; TESScut S105 (no light curve exists).

**Signal.**
- **Ephemeris:** P = 0.703442 ± 0.000005 d, T0 = BTJD 1659.684 ± 0.012 (per-sector mid-times, χ²/dof 0.45).
- **Depth:** 499-542 ppm.
- **Fit:** Rp/R★ 0.0260 (0.0242-0.0301), **Rp = 1.30 (1.21-1.51) R⊕**, b 0.55 (0.23-0.93), ρ★ 3.6 (0.4-6.3) g/cc against the TIC's 6.6 (1.1σ).
- **Orbit and climate:** a = 0.0119 au, a/R★ 5.6, S ≈ 230 S⊕, Teq ≈ 990 K (albedo 0.3) or 1,090 K (albedo 0).
- **Mass and RV:** a rocky planet of ~2 M⊕ gives K ≈ 2.4 m/s.

**What the signal survives:**

| Data | Who searched it | Depth (ppm) | z |
|---|---|---|---|
| Their data: TESS-SPOC S13/27/67 + QLP S104 | Eumemic (SNR 8.08, rank 2 on the star, `all_signals_faint_t1`) | 483 ± 67 | 7.2 |
| QLP S94, dropped by his reader | nobody | 546 ± 112 | **4.9** (0 of 300 random ephemerides reach it) |
| TESScut S105 (fixed aperture, unprocessed) | nobody | 337 ± 165 | 2.0 |
| Fixed TESScut aperture, all 6 sectors | – | 542 ± 54 (χ² p = 0.36) | 10.0 |

**Vetting.**
- Light curve:
  - odd/even 0.3σ;
  - phase 0.5: −30 ± 59 ppm;
  - any-phase secondary 1.4σ;
  - season χ² p = 0.47; per-transit χ² p = 0.009 (marginal);
  - TESS-SPOC 561 ± 89 vs QLP 453 ± 77 ppm;
  - BLS power at P/2 and 2P lower than at P; sinusoid ratio 0.37.
- **Pixels:** the target is the best Gaia source (free-fit offset 5.4 ± 1.8″, pulled away from every neighbour), best in 5 of 6 sectors. The ring test is flat.
  - **Neighbours:**
    - the T = 18.0 star at 5.9″ needs an 11% eclipse and fits worse by Δχ² 19;
    - the T = 14.6 star at 20.3″ needs 0.9% and fits worse by Δχ² 51;
    - every other Gaia star within 63″ fits worse by Δχ² ≥ 34.
- **Gaia:** 26 RV epochs with constancy p = 0.40 and per-epoch σ ≈ 3.6 km/s. An eclipsing 0.08 M☉ companion at this period would give K ≈ 29 km/s, so an on-target EB is excluded. No NSS solution; not a Gaia EB.
- **TRICERATOPS:**
  - default (5 runs × 10⁶ draws): FPP 0.080 ± 0.001, NFPP 0.040;
  - with MOLUSC (Gaia contrast + RUWE; 238,029 of 500,000 companions survive) and the neighbours beyond 10″ cleared by the pixels, keeping the 5.9″ star: **FPP 0.041 ± 0.001, NFPP 0.034**. Almost all of the residual is an eclipsing binary on that T = 18.0 star (TIC 1990153847, NEB 0.033);
  - clearing it as well (pixel Δχ² 19, about 4.3σ): **FPP 0.0071 ± 0.0006, NFPP 0**, mostly TP 0.90 and PTP 0.08.
  - So the honest status is FPP 0.04 and NFPP 0.03: just above the validation bar of FPP < 0.015 and NFPP < 0.001. A ground-based seeing-limited transit, which resolves 6″ easily, would settle the neighbour. The STP scenario (0.003) is small because MOLUSC and Gaia already constrain companions.
  - **No high-resolution imaging exists.** The ExoFOP target page is empty.
- **Archival plates:** not run.

**Novelty (live, 2026-10-14 02:20 UTC):**
- **No hit in:**
  - ExoFOP TOI, CTOI and target page;
  - the NASA Exoplanet Archive;
  - SPOC TCE S1-107, and no DV products on MAST;
  - arXiv (TIC and Gaia ID);
  - Zenodo;
  - PHT Talk;
  - T16 and RAVEN;
  - all 14 community clones and GitHub code search (9 non-astronomy hits).
- **Eumemic:** only the raw row (0.70344 d, SNR 8.08) in `all_signals_faint_t1`; it appears in no list, tier or follow-up.

**Context.** It is round 7's TIC 144439305 again: an M-dwarf USP rocky planet below Eumemic's SNR 9 pre-filter, with half its new data in the QLP files his reader drops.

##### TIC 262843927

**Star.**
- TIC 262843927 = Gaia DR3 4673750877131752832.
- T = 11.98, 83 pc, Teff 3,725 K, R★ 0.564, RUWE 1.31.
- Gaia RVs constant (23 epochs, p = 0.98).
- Data: SPOC 2-min S1-3, TESS-SPOC S4, 28-30, 34, 61, 64, 68, 69 (his era); QLP S95-97, empty-error files his reader dropped (new).

**Signal.**
- P = 24.10896 ± 0.00036 d; depth 957-1,201 ppm; Rp 2.11 (1.93-2.31) R⊕; b 0.40.
- ρ★ 4.0 (1.6-7.6) against the TIC's 4.4.
- S ≈ 3 S⊕, Teq ≈ 340 K. A warm sub-Neptune.
- Hold-out: 720 ± 219 ppm, z 3.3, in QLP S95-97 (dropped by his reader); SNR_full 10.1; 13 transits.
- No odd/even (1.4σ) and no secondary.

**Pixels.** On target (offset 0.0 ± 4.1″).

**TRICERATOPS.** FPP 0.027 ± 0.002 (TP 0.82, PTP 0.15, STP 0.03), NFPP 0.

**Second signal.** Eumemic's `followup_faint_t1.csv` lists the same star at 14.51683 d (1,089 ppm, no flags), and it is in no tier.
- Run through our identical procedure, it passes too: hold-out z 3.65, null p 0.003, SNR_full 12.4, Rp 1.9 R⊕.
- A two-transiting-planet system would earn the usual multiplicity boost (Lissauer+ 2012; Guerrero+ 2021). That would bring both signals' FPP well under the bar, but we haven't computed it formally.
- **Credit:** the 14.5 d signal is Eumemic's (follow-up list). The 24.1 d one is in his raw table only (SNR 8.53, rank 2).

**Novelty.** Clean, except Eumemic's tables and ExoNit's lot-dossier veto list, which names the star only as a "deep-systematic" host, not a candidate.

###### 14.52 d signal: pixel vet and FPP

Using the saved full-data ephemeris (P = 14.51686896 d, T0 = 1325.657972 BTJD, duration 2.444 h), `deep_vet.py` reused exactly the eight cached cutout sectors of the 24.1 d vet, with six yielding usable transits (S1, 3, 29, 61, 68, 97; S31 failed the existing quality gate and S95 had no usable transit). **The on-target origin is not established:** the joint difference-image centroid is 21.4 ± 20.8″ from the target (jackknife uncertainty), with target Δχ² 61.3; the target is preferred in five of six sectors, but S97 favours the T = 18.65 Gaia source 4673750838475299968 at 42.0″. Removing S97 in the built-in jackknife brings the free centroid to 3.7″, so this is an unresolved, sector-sensitive localisation rather than a secure blended-EB identification. Fixed-aperture depths average 753 ± 208 ppm and fail consistency (χ² = 22.79, 5 dof, p = 0.000371); S97 is −1,056 ± 470 ppm. `round6/run_tri.py`, using the original 3×3 apertures unchanged, 4-minute bins and three repeats of 10⁶ draws, gives **FPP 0.05125 ± 0.00162, NFPP 0** (error is the population standard deviation across repeats, as in the reference runner). No neighbours were cleared, no contrast curve or multiplicity boost was applied, and no MOLUSC variant was run because the 24.1 d reference had none. **Verdict: unresolved planet-like signal, not statistically validated**; its hold-out pass alone does not establish a second on-target planet. All light curves and cutouts were reused from cache, with zero sector-data downloads. Full results are in `runs/round8_crowd/tic262843927_1452_vet.json`, `deepvet/tic262843927_1452.json` and `tri/triceratops_tic262843927_1452.json`.

#### 9. Caveats

- **Selection and scope:**
  - 3,500 of the 7,092 small stars with any testable signal were tested, the highest-power half.
  - Stars with no data outside the source's own (most CVZ stars and anything Eumemic searched with all sectors) can't be tested this way at all.
  - **This is a census of the testable crowd discards, not of all of them.**
- **Window mode is thin:** single-sector sources (ExoNit, Rahle, Valdez) contribute 46 tested signals, because their stars are mostly FGK. The 0.01 d period precision of ExoNit forced the window test.
- **The noise controls control noise, not astrophysics.** 14 of 17 pre-registered survivors are stellar variability or blends.
  - The pixel test was decisive four times (blends 14-70″ away), as Avijit Thawani's retractions predicted.
  - Without the deep vet, this round would have "found" 17 planets.
- **TESScut S105-107** were used only in the deep vet. A scaled hold-out on unprocessed sectors (our own photometry) would add data to every star, but costs about 5-20 MB per sector per star.
- **No high-resolution imaging or archival plates** exist for any lead. TRICERATOPS's STP and background scenarios within about 1-2″ are constrained only by Gaia and priors.
- **Per-transit χ² for TIC 101653309** is marginal (p = 0.009), but the fixed-aperture per-sector depths agree (p = 0.36). It's probably QLP noise underestimation; worth a look in a note.

## Hold-out rerun

_Source: [round8-crowd/vetting/holdout_2026-10-31/README.md](https://github.com/rpinho/planet-hunt/blob/round8-crowd/vetting/holdout_2026-10-31/README.md)_

Hold-out rerun for the 2026-10-31 correction.

Each hold-out here uses the ephemeris and duration fitted before its hold-out sectors entered any fit: mask-A null z, 2,000 random ephemerides. `*_pre.json` is the value in the note. `*_note.json` is the earlier value at the printed ephemeris, which was fitted with hold-out sectors, kept for comparison. Printed ephemerides are for predicting transits.
