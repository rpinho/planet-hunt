![planet-hunt banner](assets/banner.png)

# planet-hunt: the record

**AI for Astronomy.**

The TESS planet candidates I've worked on with AI agents, in one place. Every section names the star, what we found, how it was checked, the verdict, and where it stands now. The signals remain **candidates**, not confirmed planets.

What was easy, what was hard, and what is still open: [RETRO.md](RETRO.md). Current status: [RESEARCH_STATUS.md](RESEARCH_STATUS.md).

Sister hubs: [AI for biology](https://github.com/rpinho/ai-biology) · [math](https://github.com/rpinho/math-certificates) · [physics](https://github.com/rpinho/ai-physics)

[PAPER-openversions.md](PAPER-openversions.md) is the same file without the TOI-4342 section, which is in review at RNAAS.

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

# TOI-4342 (TIC 354944123), 18.9 d and 1.96 d

Branch: [round6-toi4342](https://github.com/rpinho/planet-hunt/tree/round6-toi4342)

## Status

Candidate (18.9 d) and lower-confidence candidate (1.96 d). The 18.9-day signal was first reported publicly by Eumemic. The note is in review at RNAAS. Corrected version: [Zenodo 23824790](https://doi.org/10.5281/zenodo.23824790). This section is left out of [PAPER-openversions.md](PAPER-openversions.md) while the note is in review. One correction to the report below: its caveats say Claude Opus 5.5 was the only AI tool. That's wrong. OpenAI Codex was used too, and the corrected note's AI disclosure names both.

## Deep vet: claim, method, vetting checks and verdict

_Source: [round6-toi4342/runs/round6_toi4342_report.md](https://github.com/rpinho/planet-hunt/blob/round6-toi4342/runs/round6_toi4342_report.md)_

#### Verdict

**18.908-d candidate ("d18"): STRONG CANDIDATE — not validated.**
- The signal is real and periodic. Every one of the 10 TESS windows shows the transit at the predicted time. The S101 "null" from round 5 was a pipeline artifact, not missing signal.
- It is on the target star to within about 11″. No Gaia neighbour can produce it.
- The shape matches the star's density and is coplanar with planets b and c. The system is dynamically stable.
- TRICERATOPS with the Gemini/Zorro contrast curve and the pixel-excluded neighbours removed gives **FPP = 0.018 ± 0.002** and **NFPP = 0**. NFPP meets the bar (< 0.001), but FPP sits just above the 0.015 bar, and it stays there in every configuration: 0.016 with the contrast curve alone, and 0.016 ± 0.011 when MOLUSC companion limits are added. **So the candidate is not validated.** Almost all of the residual is the TRICERATOPS "STP" scenario: a planet transiting an unresolved bound companion.
- It supports an RNAAS candidate note. A draft is in `runs/round6_toi4342_rnaas/`.

**1.958-d candidate ("d2"): WEAK-TO-MODERATE CANDIDATE.**
- Its period is clean and it is consistent across seasons and between the two halves of the data.
- But it is a 290-ppm, SNR 8.8 signal, its pixel localisation is loose (a few neighbours are only Δχ² ≈ 22-25 away), and the archival plates cannot exclude a background star bright enough to make 284 ppm.
- TRICERATOPS gives FPP = 0.013 (contrast curve and neighbours beyond 30″ cleared), or 0.0065 ± 0.0014 with MOLUSC and NFPP = 0.0010 (right at the bar).
- With MOLUSC the numbers formally reach the bar (FPP 0.0065, NFPP 0.00097). But NFPP sits on the threshold, and TRICERATOPS assumes the signal is astrophysical. At SNR 8.8, after ~10⁵ trials per star, that is exactly what is not yet established. **Not validated.**
- Mention it in the note as a secondary signal at most. Don't make it the headline.

#### Key numbers (d18)

Fit: 2-min PDCSAP, 9 sectors, circular orbit, stellar-density prior ρ★ = 3.86 ± 0.36 g cm⁻³ (Parc et al. 2026).

| Quantity | Value |
|---|---|
| P | 18.908226 ± 0.000042 d |
| T0 | 2459904.5666 ± 0.0021 BJD_TDB (BTJD 2904.5666) |
| Depth | 1182 ± 107 ppm |
| Rp/R★ | 0.0344 ± 0.0016 |
| Rp | 2.25 ± 0.13 R⊕ (R★ = 0.598 ± 0.018 R☉) |
| b, i | 0.76 ± 0.03, 88.96 ± 0.07° |
| a/R★ | 41.7 ± 1.3 (a = 0.116 au) |
| T14 | 2.44 ± 0.11 h. Round 5's 1.89 h was the BLS box width, not T14. |
| S, Teq | 5.6 S⊕, ~420 K (A = 0) |
| Fixed-ephemeris SNR | 11.7 (2-min, red-noise scaled). SAP 10.7, QLP 8.8, TESS-SPOC FFI 6.7, 20-s 7.4 |
| Free-density fit | ρ_circ = 8.8 [2.8, 13.5] g cm⁻³ (16/50/84%), consistent with ρ★ = 3.86 at ~1σ. Round 5's "8.8 vs 3.9" tension is not significant. |
| Turtelboom et al. 2025 prediction | 18.5 (+2.7/−2.0) d, 2.24 R⊕, i = 88.4°, ~1170 ppm. Our values are 18.908 d, 2.25 R⊕, i = 88.96°, 1182 ppm. |

#### What each test showed

##### 1. Every TESS photon (`round6/fetch.py`, `inventory.py`, `ttv_known.py`)

- **Products downloaded (487 MB):**
  - SPOC 2-min light curves and TPFs for all 9 sectors that ever observed the star (S13, 27, 66, 67, 93, 94, 101, 102, 103)
  - 20-s light curves for S93-S103
  - QLP (8 sectors), TESS-SPOC FFI (S13-S67) and TARS
- **tess-point (pointing table through Sector 134, ~2028): the star is not observed again.** FFIs add no sector that 2-min lacks; they add only cadences inside the same sectors.
- **Predicted windows.** There are 14 windows in the sector spans. 10 have data, the other 4 fall in data gaps.
  - The 10 with data include one window round 5 never used (S13 E−66, in cadences SPOC flags 8192, "planet-search exclude") and a partial S103 window (42% coverage in 2-min, 77% in 20-s).
- **Per-window fit.** Each window was fit with the shape fixed, a local linear baseline, the depth scale free, and the other three planets' models subtracted.

| Epoch | Sector | Tc (BTJD) | 2-min depth scale | O−C (min) | SAP | 20-s | QLP |
|---|---|---|---|---|---|---|---|
| −66 | 13 (flagged cadences) | 1656.624 | 0.77 ± 0.31 | (−20 ± 105) | – | – | 0.78 ± 0.27 |
| −65 | 13 | 1675.532 | 1.16 ± 0.23 | +2 ± 8 | 1.06 | – | 0.67 ± 0.30 |
| −45 | 27 | 2053.697 | 1.10 ± 0.23 | 0 ± 6 | 0.86 | – | 0.69 ± 0.26 |
| 11 | 66 | 3112.557 | 1.41 ± 0.23 | −7 ± 4 | 1.37 | – | (no QLP) |
| 12 | 67 | 3131.465 | 0.98 ± 0.22 | −1 ± 7 | 0.93 | – | 1.04 ± 0.26 |
| 13 | 67 | 3150.374 | 0.62 ± 0.23 | +10 ± 33 | 0.64 | – | 0.92 ± 0.29 |
| 49 | 93 | 3831.069 | 0.81 ± 0.22 | +14 ± 14 | 0.73 | 0.94 ± 0.21 | 0.76 ± 0.24 |
| 64 | 101 | 4114.693 | 1.18 ± 0.22 | −18 ± 14 | 1.16 | 1.17 ± 0.26 | 0.94 ± 0.21 |
| 65 | 102 | 4133.601 | 1.08 ± 0.22 | +16 ± 10 | 0.88 | 0.97 ± 0.21 | 0.68 ± 0.18 |
| 67 | 103 (partial) | 4171.417 | 1.12 ± 0.40 | +5 ± 39 | 0.77 | 0.77 ± 0.26 | – |

- **All 10 windows are positive.** The depth-consistency χ² is 8.7 for 9 dof. Combined: 1.03 ± 0.08 (13.6σ from summing the windows).
- **Dropping windows.** Without the 3 most significant windows the signal is still 9.9σ. Without S101 it is 12.5σ. Without S101, the partial S103 and the flagged S13 window it is 12.0σ.
- **Season split.** The four seasons give 0.97 ± 0.17, 0.99 ± 0.14, 0.60 ± 0.22 and 1.02 ± 0.15 (χ² 2.8 / 3 dof). The two halves of the data give 0.98 ± 0.11 and 0.89 ± 0.12.
- **Noise null.** The same comb statistic was run on the data with the real windows removed, at 3,000 random periods between 12 and 30 d, taking the maximum over phase. Median 4.6, 99.9th percentile 8.3, maximum 8.9. The observed value is 13.4. Not once in 3,000 periods did noise reach it.
- Plots: `round6_toi4342/plots/d18_every_window.png`, `folds_d18_d2.png`.

##### 2. The S101 null: explained (pipeline artifact)

- **The cause.** The S101 2-min data have a 2.4-d gap (BTJD 4112.04-4114.47). The transit starts about 0.17 d after the data resume.
  - Round 5's detrender trims 0.25 d after every gap of more than 0.5 d. That threw away the ingress and most of the transit: only 19 of ~56 in-transit cadences were left, and they gave −122 ± 404 ppm.
  - With all the cadences, the window shows 1.18 ± 0.22 (2-min), 1.16 (SAP), 1.17 ± 0.26 (20-s) and 0.94 ± 0.21 (QLP).
- **Caveat.** The event sits near the start of an orbit, where thermal ramps happen. It is weighted like any other window, but the result does not depend on it (12.5σ without it). Plot: `plots/d18_s101_detail.png`.
- **Period and cycle count.** The coherent comb peaks at 18.9082 d. The alias tests all point the same way:

| Test | Depth |
|---|---|
| P/2 (in-between epochs) | 53 ± 78 ppm (a planet at P/2 would give ~1150) |
| P/3 | −88 ± 82 and 130 ± 117 ppm |
| P/4 | −42 and −78 ppm |
| 2P (odd vs even epochs) | 1179 ± 105 vs 1106 ± 154 ppm |
| 3P (three residue classes) | 1096 ± 129, 1426 ± 179, 1036 ± 157 ppm |

  - So P/2, P/3, P/4, 2P and 3P are all excluded.
  - In-season triplets (S66-S67: E11, 12, 13; S101-S103: E64, 65, 67) fix the period locally far better than the ±0.33 d needed to change the cycle count between seasons. No nearby-period alias survives (`plots/d18_period_scan.png`).
  - The longer-period peaks (37.8 and 64.6 d, SNR ~11.8) are subsets of the same events.
- **Supported:** P = 18.90823 d, linear cycle count E = −66 … 67.
- **TTV check.** Timing posteriors were taken from the per-transit likelihood. The d18 O−C give χ² = 7.9 for 6 dof against a linear ephemeris, so no TTV is needed (rms 10.7 min, median error 8.9 min). The S101 and S102 epochs scatter by −18 and +16 min, which is within errors.

##### 3. Dynamics (`ttv_known.py`, `ttv_sin.py`, `dynamics.py`, `stability.py`, `stability_lit.py`)

- **b and c TTVs from TESS.**
  - c: χ² = 10.8 / 17 against linear.
  - b: χ² = 104.7 / 35, driven by a few outlier epochs (e.g. E221 at −30 min). The per-transit errors on b look about 1.7× too small.
  - A sinusoid at the b-c near-2:1 super-period (152.6 d) improves b by Δχ² = 16.6 (amplitude 6.8 min), and c by 1.9 (5.9 min). With b's errors inflated to χ²_ν = 1, that is not significant.
  - Bottom line: TESS cannot measure TTVs below ~7 min here. Parc et al. 2026 also report no TTVs.
- **N-body predictions** (rebound/IAS15, circular and Rayleigh e ~ 0.03, M-R or literature masses).
  - Half-amplitudes predicted over the TESS baseline: b 2-3 min, c 1-7 min, d18 1.4-3 min.
  - P_d18/P_c = 1.7690, 1.1% wide of 7:4. That is a third-order resonance, so its TTV term goes as e³ and is negligible.
  - Adding d18 changes the predicted b/c TTVs by less than 2 min.
  - χ² against the TESS timings: b+c model 108 (linear 115.5); b+c+d18 model 117.5 (linear 123.4).
  - **The TESS timings are consistent with d18 and cannot confirm it dynamically.** Plot: `plots/ttv.png`.
- **Stability (SPOCK FeatureClassifier, 10⁹ inner orbits):**

| Configuration | p(stable), nominal | Median of draws | Draws with p > 0.5 |
|---|---|---|---|
| b+c+d18 (M-R masses) | 0.96 | 0.90 [16-84%: 0.76-0.93] | 98% |
| b+c+d18 (Parc masses: b 7.3, c 4.8 M⊕) | 0.95 | 0.89 | 95% |
| d2+b+c | 0.97 | 0.96 | 100% |
| d2+b+c+d18 (M-R masses) | 0.97 | 0.93 [0.77-0.95] | 97% |
| b+c+d18 + 47.5-d RV candidate (Parc masses) | 0.94 | 0.89 [0.77-0.92] | 97% |
| d2+b+c+d18 + 47.5-d candidate | 0.97 | 0.93 [0.73-0.95] | 97% |

  - Mutual Hill spacings: d2-b 28.7, b-c 15.9, c-d18 14.1, d18-(47.5 d) 18.2 (M-R masses). All are above the ~10 usually needed for long-term stability.. MEGNO: the 4- and 5-planet configurations were integrated with WHFast for 2×10⁵ orbits of d2 (~1,070 yr), nominal plus 8 and 4 random draws. Every run gave ⟨Y⟩ = 1.96-2.004, i.e. quasi-periodic, with no escapes..
  - The 47.5-d RV candidate of Parc et al. 2026 (m sin i 17.8 M⊕) was included as a coplanar outer planet where noted.
  - **Both candidates are dynamically allowed.**

##### 4. ExoFOP (read only)

- **High-resolution imaging (two speckle observations, no companions detected):**

| Instrument | Date | Filter | Contrast | File |
|---|---|---|---|---|
| Gemini-S/Zorro | 2021-07-23 | 832 nm | Δm = 4.7 at 0.1″, 5.9 at 0.5″, 6.1 at 1.15″ | `exofop/imaging/TOI4342I-ef20210723-832_sensitivity.dat` |
| SOAR/HRCam | 2021-10-01 | I | — | PDF only. ExoFOP's "6.5 mag at 1″" summary is reached only at ~3″ on the plot. |

  - The Zorro curve was converted to `tri/zorro832.csv` and passed to TRICERATOPS with `filt=TESS`, since 832 nm sits inside the TESS band.
- **SG1 ground-based light curves.** There are 8 LCO 1-m i′ observations (2021-05-30 to 2021-09-14), all on b or c. Data are in `exofop/tseries/` (51 files).
  - None comes within 2 d of a predicted d18 transit.
  - For d2, the 2021-08-31 CTIO light curve ends 0.2 h before a predicted transit, so it offers no constraint.
- **Lists.** The live TOI and CTOI lists have only .01/.02 on this star. No CTOI exists near 18.9 or 1.96 d or their aliases.

##### 5. Contamination (`contam/contamination.json`, `contam/archival_imaging.png`)

- **Gaia DR3 target.** No companion signature: RUWE 0.98, ipd_frac_multi_peak 0, ipd_gof_harmonic_amplitude 0.011, non_single_star 0, no NSS solution. The RV scatter is consistent with noise (p = 0.23). TIC contamination ratio 0.0095.
- **Neighbours within 150″.** 57 Gaia sources, the closest at 16.8″. 32 of them could in principle mimic 1150 ppm, if fully in the aperture:
  - G = 18.4 at 18.9″ (needs a 59% eclipse)
  - 17.7 at 37.6″ (22%)
  - 15.0 at 60.9″ (2.6%)
  - 15.0 at 94.5″ (2.7%)
  - 14.0 at 140.5″ (1.2%)
- **No bound companion in Gaia.** No co-moving source; the one similar-motion star at 95″ is three times farther away.
- **Archival plates.** The SERC-J plate of 1976.66 is the oldest. Pan-STARRS does not reach Dec −78.
  - The star moved 7.4″ by 2026, which is 2.2× the plate's FWHM.
  - Today's position shows no source down to ΔBP ≈ 7.25. That means ΔT ≲ 6.9-7.25 for a source of the same colour, and ≲ 7.6-7.9 for typical bluer background colours.
- **Verdict:**
  - For **1150 ppm** (needs ΔT ≤ 7.4), a hidden background star is largely excluded: a same-colour one would need a ≥ 69% eclipse, and typical colours are excluded outright.
  - For **284 ppm** (needs ΔT ≤ 8.9), stars with ΔT ≈ 7.3-8.9 are not excluded.

##### 6. A pixel test that works here (`tpfvet.py`, `tpfvet2.py`, `loccal.py`)

- **Why round 5 failed.** The PSF test used TESScut cutouts and a parametric PSF, and planet c "failed" it (Δχ² 263). Most of that came from one sector, S103, where the scene fit was off (target-in-aperture fraction 1.29). A per-sector Gaussian centroid on the 2-min TPFs also gave known planets 11-14″ "offsets" with underestimated errors.
- **Method.** For each sector and each planet, a stacked difference image (in-transit minus a linear per-pixel out-of-transit model). The noise image comes from 40 fake epochs.
  - Each difference image is fit with the same sector's **out-of-transit image used as an empirical PRF**, shifted over a 2″ sky grid (±110″) using the TPF WCS.
  - χ² surfaces are scaled by the reduced χ² and summed across all 9 sectors.
  - Because the template is the star itself, no PSF model is needed.
- **Calibration on the same pixels:**

| Signal | Best offset | Δχ²(target) | Closest mimic-capable neighbour |
|---|---|---|---|
| b | (0, +4)″ | 9.5 | Δχ² ≥ 131 |
| c | (−2, 0)″ | 3.5 | Δχ² ≥ 73 |
| **d18** | (−2, +4)″ | **11.4** (8.0 calibrated by b+c per sector) | Δχ² ≥ **79** |
| d2 | (−10, +14)″ | 28.0 (9.7 calibrated) | Δχ² **24.7** (24″, T = 18.1) |

  - So on-target known planets give Δχ²(target) of 3.5-9.5 and off-target hypotheses ≥ 73.
  - For d18, the 3σ region (Δχ² < 11.8) reaches at most 11″ from the target. Every Gaia star that could mimic it is at Δχ² ≥ 79: the 24″ star at 79, the 18.9″ star at 87, the others at > 90.
  - Stars beyond the 11×11 TPF (≥ 105″, e.g. G = 14.0 at 140″) put less than ~0.1% of their flux in the SPOC aperture and would need eclipses deeper than 100%.
- **Depth versus aperture**, as (depth in aperture / depth in SPOC aperture) for the candidate divided by the same ratio for the known planet:

| Aperture | d18 vs b | d18 vs c |
|---|---|---|
| Core 1-2 px | 1.05 ± 0.13 | 1.19 ± 0.16 |
| r ≤ 1.5 px | 1.06 ± 0.08 | 1.13 ± 0.10 |
| r ≤ 2.5 px | 0.99 ± 0.09 | 1.04 ± 0.10 |
| r ≤ 3.5 px | 0.92 ± 0.11 | 0.99 ± 0.13 |
| Annulus | 0.58 ± 0.39 | — |

  - A source off the target would show a shallower core and a deeper outer aperture. The candidate behaves like b and c. d2 behaves the same within its larger errors.
- Plot: `plots/pixel_localisation.png`.

##### 7. Shape, density and coplanarity (`mcmc.py`, `mcmc_prior.py`)

- **Free circular density.** The d18 fit gives 8.8 [2.8, 13.5] g cm⁻³. The same method gives 4.3 [2.6, 5.0] for b and 3.6 [1.7, 4.6] for c, against ρ★ = 3.86.
  - The photo-eccentric ratio g = (ρ_circ/ρ★)^(1/3) is 1.31 [0.90, 1.52]. That is consistent with a circular orbit, or with modest eccentricity near periastron.
  - With the density prior the shape is fully consistent with the host: b = 0.76, T14 = 2.44 h.
- **Coplanarity.** The inclinations are b 88.83 ± 0.22° (Parc et al. 2026), c 89.61 ± 0.26° and d18 88.96 ± 0.07°. The three agree within ~0.7°, as expected for a compact coplanar system.
  - At i ≈ 89° and a/R★ ≈ 42, d18 has b ≈ 0.76, which is why its transit is shorter than a central one would be.
- **Other light-curve checks.**
  - Odd/even differ by 0.5σ.
  - Phase-0.5 secondary: 0.03 ± 0.07 of the depth.
  - The highest any-phase dip is 3.3σ against a scan SD of 1.06, which is chance level.

##### 8. Future coverage (`future.py`, `data/future_windows.json`)

- **TESS:** no sectors through S134 (~2028).
- **Ground:** the star is circumpolar from Chile, South Africa and Australia, but never higher than about 42°.
- **Predicted d18 transits through 2027 that are fully observable** (transit ± [30 min + 1σ], Sun < −12°, altitude > 25°). Timing 1σ is ±11-12 min, including a 10-min TTV allowance.

| Date and time (UTC) | Site | Max airmass |
|---|---|---|
| 2026-10-28 02:10 | CTIO | 1.73 |
| 2026-12-04 21:46 | SAAO | 2.36 |
| 2027-05-24 01:56 | SAAO | 1.65 |
| 2027-06-11 23:44 | SAAO | 1.77 |
| 2027-06-30 21:32 | SAAO | 1.92 |
| 2027-07-19 19:20 | SAAO | 2.1 |
| 2027-08-07 17:08 | SSO | 1.72 |
| 2027-08-26 14:55 | SSO | 1.61 |
| 2027-09-14 12:43 | SSO | 1.53 |

  - A 1.2-ppt, 2.4-h transit on a T = 11.0, I ≈ 10.4 star is within reach of 1-m class telescopes (LCO, as SG1 did for b/c at 1.4 ppt).
  - Mid-transit times: Tn = 2459904.5666 + 18.908226 n BJD_TDB.
- **d2:** 104 of 227 windows are observable, but 290 ppm in 1.4 h is not realistic for 1-m class ground photometry. It needs CHEOPS or a larger aperture.

##### 9. Literature (live, 2026-10-13; `lit/literature.json`)

- **Sources checked:**
  - Searches: arXiv full text ("TOI-4342", "TOI 4342", 354944123), Google Scholar, Semantic Scholar/OpenAlex citations
  - Exoplanet Archive, the ExoFOP TOI/CTOI lists, and the SPOC TCEs for s1-s96 and s101-s103
  - Papers read: Tey et al. 2023, Turtelboom et al. 2025, Parc et al. 2026 (A&A 708, A81) and two Sep-2026 arXiv catalog papers
- **Result:** nobody reports an 18.9-d or 1.96-d signal, or any alias of either. No SPOC TCE matches.
- **What we could not access:** ADS full text (needs a token). The A&A page for Parc et al. returned 403, so the arXiv v2 text was used.
- **Published system parameters:** Parc et al. 2026 measure b = 7.3 ± 1.3 M⊕ and c = 4.8 ± 1.4 M⊕, circular orbits, no TTVs, a 47.5-d RV candidate, and P_rot = 14.7 d.

##### Statistical validation (TRICERATOPS; `round6/run_tri.py`, `tri/`)

- **Setup.** SPOC optimal apertures from all 9 sectors, a Gaia background population, and the fold of the 2-min data with the other planets removed (4-min bins). Each configuration was run 2-3 times with 10⁶ draws. The bar is FPP < 0.015 and NFPP < 0.001.

| Signal | Configuration | FPP | NFPP |
|---|---|---|---|
| d18 | no contrast curve (2 runs; round 5 got 0.12 / 0.004) | 0.099 | 0.0003 |
| d18 | Zorro 832-nm contrast curve | 0.0161 ± 0.0002 | 0.0005 |
| d18 | contrast curve + every star beyond 12″ cleared (pixel test) | **0.0176 ± 0.0022** | **0** |
| d18 | as above + MOLUSC companion population | 0.016 ± 0.011 (runs 0.004, 0.012, 0.031) | 0 |
| d2 | contrast curve | 0.0149 | 0.0025 |
| d2 | contrast curve + stars beyond 30″ cleared | 0.0130 | 0.0010 |
| d2 | as above + MOLUSC | 0.0065 ± 0.0014 | 0.00097 |

- **Where the d18 FPP comes from.** The surviving false-positive probability is almost all STP (planet around an unresolved bound companion, 0.015-0.02). Planet-on-primary-with-companion (PTP, 0.02-0.08) counts as a planet. EB scenarios are below 10⁻³.
- **MOLUSC** (Wood et al. 2021; `.tools/molusc`) used the Zorro curve (I band) and the Gaia contrast test.
  - Its RUWE module refused this star ("magnitude or colour outside bounds"), so Gaia's RUWE of 0.98 is **not** in the companion limits.
  - The SOAR curve exists only as a PDF and was not used.
  - The Parc et al. ESPRESSO RVs are not public (no CDS table), so MOLUSC got no RV constraint either.
  - Large run-to-run scatter: MOLUSC kept 226,866 companions, which TRICERATOPS pads out to 10⁶ draws.
- **TRICERATOPS ignores the multiplicity boost.** Third candidates in multi-transiting systems are false positives far less often than single candidates (Lissauer et al. 2012). Many validation papers multiply the prior by that factor, which would put d18 well below 0.015. I have not applied it, and whether to use that argument is Ricardo's call.
- **Technical notes.**
  - pandas 3 makes `.values` read-only, which crashes TRICERATOPS's MOLUSC path. `run_tri.py` hands it a writable copy.
  - MOLUSC needed a one-line CLI fix (`star_jitter`) and scipy < 1.14. Neither is committed; both are under `.tools/`.

#### d2 (1.95769 d) in brief

- **Fit** (ρ★ prior): depth 287 ± 40 ppm, Rp = 1.11 ± 0.08 R⊕, b = 0.56, T14 = 1.38 h, T0 = BTJD 2900.1282 ± 0.0035.
- **Significance.** SNR 8.8 (2-min) and 8.9 (SAP). QLP gives 6.9 and the 20-s data 6.9. TESS-SPOC FFI gives only 3.4.
- **Consistency.**
  - Seasons: χ² 2.5 / 3.
  - Halves: 0.93 ± 0.16 and 0.95 ± 0.14.
  - Odd/even: 0.5σ.
  - Phase-0.5: 0.21 ± 0.11 (1.9σ).
- **Aliases.** P/2, P/3 and P/4 are excluded (phase 1/3, 1/2, 3/4 amplitudes of 0.12, 0.21 and 0.08 ± 0.11). 2P and 3P are excluded because all residue classes show the transit.
- **Free-density fit:** 4.7 [2.0, 7.4] g cm⁻³, consistent with the star.
- **Pixels.** The source is near the target, but loosely: 3σ radius ~28″ uncalibrated, Δχ²(target) 9.7 after calibration. The neighbours at 24″ and 18.9″ are only Δχ² 24.7 and 49.8 away.
- **Background stars.** The archival plates cannot exclude one with ΔT 7.3-8.9.
- **Statistics.** FPP 0.013 (contrast curve and neighbours beyond 30″ cleared), or 0.0065 ± 0.0014 with MOLUSC, NFPP 0.0010 (right at the bar). Stable in every configuration tested.
- **Caveat.** Round 5's own noise model predicts 1-4 pure-noise signals above SNR 8.5 across the survey.

#### What is still needed for a statistical validation of d18

FPP has to come down from ~0.017 to below 0.015. NFPP already meets the bar. Any one of these should do it:

1. **One ground-based transit** with an LCO 1-m in i′ (ideally simultaneous g′ or z′). Windows are listed in section 8; the next two are 2026-10-28 02:10 UTC (CTIO) and 2026-12-04 21:46 UTC (SAAO).
   - It is independent photons, so it confirms the signal is real. That is the real gap; everything above shares TESS photons.
   - It localises the event at about 1″ seeing, which clears neighbours properly instead of through TESS pixels.
   - An achromatic depth argues against a planet on a redder companion (STP) and against EBs. That is exactly the remaining FPP.
   - This is what TFOP SG1 would do, and it is the strongest single next step.
2. **Bound-companion limits that include Gaia RUWE and RVs.** Fix MOLUSC's RUWE module for this colour/magnitude (or apply the RUWE limit by hand), digitise the SOAR curve, and use the ESPRESSO RVs if Parc et al. share them.
   - Those RVs would also test d18 directly. The expected semi-amplitude is K ≈ 2 m/s for ~6 M⊕, which matches Turtelboom et al.'s prediction, under an active star (P_rot = 14.7 d, which is not a harmonic of 18.9 d).
3. **The multiplicity-boost argument**, applied explicitly and stated as such.

Without new data, the honest label is **strong candidate (statistically at, but not below, the validation threshold)**. That is suitable for an RNAAS candidate note and a CTOI, not a validation claim.

#### Caveats

1. **Same photons.** Every dataset here shares photons with the discovery data. QLP, TESS-SPOC and 20-s are different reductions of the same TESS exposures.
   - The genuinely new inputs are: the two windows round 5 could not use (flagged S13 E−66, and the trimmed part of S101 E64), the 20-s cadence, and the speckle contrast curve.
   - An independent ground transit is the real test.
2. **QLP depths run about 20% shallower** than SPOC 2-min (0.80 ± 0.09 of the model). This is typical of QLP's KSP detrending on short transits and its different dilution model. The SPOC SAP depth (0.95 ± 0.08) shows the signal is not a PDC artifact.
3. **Pixel localisation is TESS-resolution.** "On target" means within about 11″ at 3σ. Companions closer than that are constrained only by Gaia (which resolves ≳ 1″; RUWE 0.98), the Zorro speckle curve (0.1-1.2″, Δm 4.7-6.1 at 832 nm) and the 1976 plate.
4. **TRICERATOPS assumes the signal is astrophysical.** The pixel-clearing radius (12″ for d18, 30″ for d2) comes from our own localisation, not from seeing-limited photometry. TFOP would normally require the latter before calling neighbours "cleared".
5. **Timing.** Timing errors are posterior widths from fixed-shape fits. b's TESS timings scatter about 1.7× more than those errors predict, so the TTV amplitude limits are approximate.
6. **Using someone else's speckle data.** The Zorro contrast curve belongs to the TFOP observers (E. Furlan, per ExoFOP). ExoFOP's data-use etiquette is to contact the observer before publishing with it.
   - The RNAAS draft carries a TODO to confirm the Gemini program ID and ask permission.
   - The draft's AI disclosure names Claude Opus 5.5 (the only AI tool used; no Codex or other model was used) and deliberately makes **no** claim that Ricardo verified anything. He must read and check it himself before submission.
7. **One window per event at most.** No single d18 transit is individually above ~6σ. The detection rests on the coherent set of 10.

## The 1.96-day signal, later verdict

_Source: [round8-p1/runs/round8_p1_report.md](https://github.com/rpinho/planet-hunt/blob/round8-p1/runs/round8_p1_report.md)_

##### 2.4 TOI-4342 (TIC 354944123), 1.958 d: WATCH, not novel

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

## False-positive probability and radial-velocity reruns

_Source: [round6-toi4342/vetting/fpp_rv_2026-10-30/README.md](https://github.com/rpinho/planet-hunt/blob/round6-toi4342/vetting/fpp_rv_2026-10-30/README.md)_

Reruns for the TOI-4342 correction (2026-10-30).

FPP (TRICERATOPS): the original nine-sector 2-minute fold, 2-minute bins, 120 s modeled exposure, five runs per configuration. MOLUSC configurations use five companion samples at N=1,000,000. `verified_padding_probes.csv` shows that N changes those results: TRICERATOPS zero-pads the surviving MOLUSC companions to N, so the draw count sets the companion prior. The values in the note are the means and sample standard deviations of each `tri_*.json`.

Radial velocities: `rv_fit.py` fits a circular orbit at P=18.908226 d, T0=59904.5666 RJD to the 55 public ESPRESSO LBL velocities of Parc et al. 2026 (DACE), with b and c at published ephemerides and the published activity kernel at fixed hyperparameters. `rv.json` holds the results. The limits are conditional: GP hyperparameters are not marginalized.

1.96 d: `short4342.json` holds the depth with and without 0.5 d around orbit boundaries, the orbit-locking controls, and the free-search period check, which fails for the literal period.

`rnaas_fig.py` is the figure script with the panel (a) label corrected to eight sectors.

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

#### 3. Novelty gate (live, read-only, 2026-10-13 23:09 UTC to 10-14 ~00:30 UTC)

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

#### 5. Caveats

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
