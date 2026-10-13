# Round 6: deep vet of TOI-4342 (TIC 354944123) — the 18.9-d and 1.96-d candidates

Branch `round6-toi4342` (worktree `~/.worktrees/planet-hunt-r6-toi4342`), 2026-10-13. Every process ran at `nice -n 20`. Nothing was submitted, uploaded, emailed, pushed or posted. ExoFOP, MAST, Gaia, the STScI/SkyView plate archives, SkyMapper, arXiv and the journals were only read.

## Verdict

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

## Key numbers (d18)

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

## What each test showed

### 1. Every TESS photon (`round6/fetch.py`, `inventory.py`, `ttv_known.py`)

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

### 2. The S101 null: explained (pipeline artifact)

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

### 3. Dynamics (`ttv_known.py`, `ttv_sin.py`, `dynamics.py`, `stability.py`, `stability_lit.py`)

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

### 4. ExoFOP (read only)

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

### 5. Contamination (`contam/contamination.json`, `contam/archival_imaging.png`)

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

### 6. A pixel test that works here (`tpfvet.py`, `tpfvet2.py`, `loccal.py`)

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

### 7. Shape, density and coplanarity (`mcmc.py`, `mcmc_prior.py`)

- **Free circular density.** The d18 fit gives 8.8 [2.8, 13.5] g cm⁻³. The same method gives 4.3 [2.6, 5.0] for b and 3.6 [1.7, 4.6] for c, against ρ★ = 3.86.
  - The photo-eccentric ratio g = (ρ_circ/ρ★)^(1/3) is 1.31 [0.90, 1.52]. That is consistent with a circular orbit, or with modest eccentricity near periastron.
  - With the density prior the shape is fully consistent with the host: b = 0.76, T14 = 2.44 h.
- **Coplanarity.** The inclinations are b 88.83 ± 0.22° (Parc et al. 2026), c 89.61 ± 0.26° and d18 88.96 ± 0.07°. The three agree within ~0.7°, as expected for a compact coplanar system.
  - At i ≈ 89° and a/R★ ≈ 42, d18 has b ≈ 0.76, which is why its transit is shorter than a central one would be.
- **Other light-curve checks.**
  - Odd/even differ by 0.5σ.
  - Phase-0.5 secondary: 0.03 ± 0.07 of the depth.
  - The highest any-phase dip is 3.3σ against a scan SD of 1.06, which is chance level.

### 8. Future coverage (`future.py`, `data/future_windows.json`)

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

### 9. Literature (live, 2026-10-13; `lit/literature.json`)

- **Sources checked:**
  - Searches: arXiv full text ("TOI-4342", "TOI 4342", 354944123), Google Scholar, Semantic Scholar/OpenAlex citations
  - Exoplanet Archive, the ExoFOP TOI/CTOI lists, and the SPOC TCEs for s1-s96 and s101-s103
  - Papers read: Tey et al. 2023, Turtelboom et al. 2025, Parc et al. 2026 (A&A 708, A81) and two Sep-2026 arXiv catalog papers
- **Result:** nobody reports an 18.9-d or 1.96-d signal, or any alias of either. No SPOC TCE matches.
- **What we could not access:** ADS full text (needs a token). The A&A page for Parc et al. returned 403, so the arXiv v2 text was used.
- **Published system parameters:** Parc et al. 2026 measure b = 7.3 ± 1.3 M⊕ and c = 4.8 ± 1.4 M⊕, circular orbits, no TTVs, a 47.5-d RV candidate, and P_rot = 14.7 d.

### Statistical validation (TRICERATOPS; `round6/run_tri.py`, `tri/`)

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

## d2 (1.95769 d) in brief

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

## What is still needed for a statistical validation of d18

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


## Caveats

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

## Files

| Path | What |
|---|---|
| `round6/*.py` (branch `round6-toi4342`) | All code: `fetch`, `inventory`, `prep`, `fit_known`, `ttv_known`, `alias`, `nulltest`, `mcmc`, `mcmc_prior`, `vetlc`, `d2alias`, `tpfvet`, `tpfvet2`, `loccal`, `dynamics`, `ttv_sin`, `stability`, `stability_lit`, `future`, `plots`, `rnaas_fig`, `run_tri` |
| `runs/round6_toi4342/data/` | Per-transit fits, ephemerides, MCMC posteriors, alias and null tests, pixel maps |
| `runs/round6_toi4342/tri/` | TRICERATOPS inputs, outputs and logs |
| `runs/round6_toi4342/dyn/` | N-body TTVs and stability |
| `runs/round6_toi4342/contam/`, `exofop/`, `lit/` | Gaia/archival imaging, ExoFOP reads and downloads, literature |
| `runs/round6_toi4342/plots/` | `d18_every_window.png`, `d18_s101_detail.png`, `d18_period_scan.png`, `folds_d18_d2.png`, `ttv.png`, `pixel_localisation.png` |
| `runs/round6_toi4342_rnaas/` | RNAAS draft (tex, bib, fig1, compiled PDF) |
