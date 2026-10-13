# Round 5, track A: hidden extra planets around small-star planet hosts

Branch `round5-hosts` (worktree `~/.worktrees/planet-hunt-r5-hosts`), 2026-10-12. Data: SPOC 2-min PDCSAP light curves (3,814 files, S1-S107), TESScut 11x11 cutouts for the vetted stars, Gaia DR3, the ExoFOP TOI and CTOI tables and NASA Exoplanet Archive `pscomppars` (all downloaded 2026-10-12), SPOC DV reports. Nothing was posted, uploaded or sent anywhere.

## Bottom line

- **184 stars searched.** These are every TOI host (TFOPWG disposition not FP/FA) and every confirmed-planet host with Teff < 4200 K or R* < 0.7 Rsun, Tmag < 14 and at least 8 sectors of SPOC 2-min data. That is 140 TOI hosts and 44 confirmed hosts. The median star has 12.5 sectors; 60 stars have 30 or more (the CVZ).
- **Completeness (300 injections, 0.5-2 Re, 1-30 d, full pipeline):**
  - 52% recovered overall at SNR >= 7.5 with >= 3 transits; 50% also pass the triage flags.
  - Recovery is 91-97% once the injected transit's true SNR exceeds 7.5, so the search itself loses little above threshold. Completeness is set by SNR.
  - By planet size: 14% at 0.5-0.75 Re, 32% at 0.75-1.0, 53-56% at 1.0-1.5, and 73-81% at 1.5-2.0 Re.
  - Stars flagged as rotation-saturated are effectively blind: 7% recovered there.
- **The method control works.** TOI-6000.02, the round-4 candidate now in the RNAAS, was re-found blind at 4.21349 d, SNR 8.8, 891 ppm.
- **No validated planet, and nothing that clears the published validation bar** (TRICERATOPS FPP < 0.015 with NFPP < 0.001). Of the 16 automatic survivors:
  - **6 are new candidates.** None is in the TOI or CTOI lists, the SPOC TCEs or the literature, and each passed the light-curve and pixel tests we can run. Ranked by evidence: TOI-4342 at 18.91 d, TOI-1846 at 0.944 d, TOI-4342 at 1.958 d, TOI-6008 at 7.22 d, TOI-249 at 19.09 d and TOI-4405 at 3.60 d. All have SNR 7.9-9.6. The last two are weak.
  - **3 are real-looking but not new.** TOI-1449 at 4.867 d and TOI-4572 at 8.683 d are unpromoted SPOC TCEs. TOI-4572 at 4.899 d is twice the period of a candidate in a July 2026 preprint.
  - **1 is ours**: TOI-6000.02.
  - **2 are marginal or doubtful**: TOI-1215 at 4.717 d and TOI-1215 at 23.05 d.
  - **4 are rejected**, including a 0.52 Re signal on LHS 1815 that is a blended eclipsing binary in the pixels.
- **Be careful with the word "new".** From the shape of the SNR distribution, roughly 5-15 noise signals are expected across 184 stars at SNR 7.5-8.5, and about 1-4 above 8.5. So any single candidate below SNR ~9 here has a real chance of being noise. A second, independent dataset is the test that matters.

## 1. Targets

| | |
|---|---|
| Catalogs | ExoFOP TOI (8,151 rows) and CTOI (5,138), and NASA `pscomppars` (6,445), all fetched 2026-10-12 |
| 2-min coverage | MAST per-sector download scripts for S1-S107 give 1.74 M TIC-sector light curves. No per-star MAST query is needed. |
| Cuts | (Teff < 4200 K or R* < 0.7 Rsun) and Tmag < 14 and >= 8 2-min sectors; hosts whose every TOI is FP/FA are dropped |
| Result | **184 stars**: 140 TOI hosts and 44 confirmed hosts. 146 have Teff < 4200 K and 175 have R* < 0.7. 61 lie within 12 deg of an ecliptic pole. Sectors: 8-44, median 12.5, 3,814 light curves in total. |
| Not included | 11 TOI hosts with neither Teff nor R* in ExoFOP. Hosts with fewer than 8 2-min sectors. |

The brief expected "a few hundred" targets. The honest count under these cuts is 184. Most small-star TOIs have fewer than 8 sectors of 2-min data, because many were found in the FFIs only.

## 2. Search

`planet_hunt/hostsearch.py`, driven by `round5_hosts.py search`:

1. **Known planets removed.** Every cadence within 0.75 duration of any TOI (all dispositions) or confirmed transiting planet is dropped, which is a window 1.5x the duration. RV-only planets are not masked; their periods are flagged.
2. **Rotation.** Lomb-Scargle per sector on 30-min bins over 0.1-12 d. On active stars (semi-amplitude > max(1000 ppm, 5x the 30-min noise)) a 3-harmonic series at P_rot, with amplitudes linear in time, is subtracted per sector, and the filter window is shortened to P_rot/6 (minimum 0.25 d).
   - A star is flagged **rotation-saturated** in either of two cases: 2 or more of its 6 strongest periodogram peaks sit at rotation harmonics, or its strongest signal is mostly sinusoidal (sine/box > 0.5).
   - The same applies if red noise has beta > 2.
   - **41 of 184 stars are flagged.** Their signals are listed separately and not reported as candidates.
3. **Detrending.** wotan biweight per sector, with the transit masks applied and the trend interpolated under them. Gaps over 0.5 d are split and 0.25 d is trimmed after each one; flares above 4 sigma are clipped. Two windows are used: 0.5 d feeds the 0.5-6 d search and 1.0 d feeds the 6-30 d search.
4. **BLS.** Fully coherent BLS on 10-min bins of all sectors. The frequency grid keeps the phase drift over the whole baseline (up to ~2,600 d) below 1/4 of a central transit duration, which gives up to ~5x10^5 trial periods per star; durations are set from the star.
   - The 25 strongest distinct peaks in each range are ranked by a quick red-noise SNR.
   - The winner is re-measured on the un-binned 2-min light curve, re-detrended with its own transits masked (`multisector.transit_stats`). Per-sector red-noise beta inflates the errors, and the step also gives odd/even, secondary at phase 0.5 and at any phase, the sinusoid test, per-transit and per-season consistency, and density-implied duration.
   - Up to three signals per star: mask the signal, drop its 1/3-3x aliases from the pool, take the next best.
5. Total CPU: 9.4 h for 184 stars on a shared, heavily loaded Mac mini.

**Data products.** SPOC 2-min PDCSAP only, for homogeneity. I did **not** stack TESS-SPOC FFI light curves into the search, deliberately: round 4's critique found that mixing products creates depth jumps. The FFI pixels come back in vetting instead, as fixed-aperture TESScut photometry of every available sector, which covers the sectors without 2-min data.

## 3. Completeness: injection-recovery

There were 300 injections into the raw normalised 2-min flux before any detrending. Each used a limb-darkened transit (with the TIC density) on a random target: Rp uniform in 0.5-2.0 Re, P log-uniform in 1-30 d, b uniform in 0-0.8, random phase. Each injected star then went through the identical search. "Recovered" means a returned signal within 1% of the injected period, in phase, with SNR >= 7.5 and >= 3 transits. "Triage" also requires no red flags and a star not flagged rotation-saturated. The statistical triage cuts (5% clip, season, crowding) are set at <= 1% false rejection and were not re-applied.

| | 1-3 d | 3-10 d | 10-30 d |
|---|---|---|---|
| 0.5-1.0 Re | 28% (n=39) | 22% (27) | 18% (34) |
| 1.0-1.5 Re | 70% (23) | 41% (39) | 61% (33) |
| 1.5-2.0 Re | 97% / 84% triage (32) | 73% / 68% (37) | 64% / 61% (36) |

- **Overall:** 52.0% recovered, 50.7% with no red flag, 49.7% through triage. Any-SNR period match: 57%.
- **By the injection's true SNR:**

  | True SNR | Recovered |
  |---|---|
  | < 5 | 0% (n=102) |
  | 5-7.5 | 31% (36) |
  | 7.5-10 | 91% (23) |
  | 10-15 | 97% (36) |
  | > 15 | 97% (92) |

  Of the 6 misses above true SNR 7.5, 3 were found at twice the injected period (counted as misses), 2 sat on rotation-saturated stars, and 1 is unexplained.
- **By number of sectors:** 8-12 sectors 42%, 13-29 sectors 55%, CVZ (>= 30) 66%.
- **Saturated stars:** 7% recovered (n=45 injections) against 60% elsewhere.
- **Reading:** a 1.5-2 Re planet out to 30 d would have been found about 2/3 of the time. An Earth-size planet only about 1/3 to 1/2 of the time. Sub-Earths mostly not. The 10-30 d band is genuinely searched this round (61-64% for > 1 Re), which answers round 4's critique.
- Plot: `runs/round5_hosts/injection_recovery.png`.

## 4. Funnel

| Step | Count |
|---|---|
| Signals returned (up to 3 per star) | 496 |
| SNR >= 7.5 and >= 3 transits | 89 |
| Known: matches a TOI, CTOI or confirmed planet of the star at a simple alias, or > 50% overlap with known transits | 22. These include L 98-59's CTOI .04 at 1.048 d, TOI-2084.02, Kepler-138 c and residuals of several USP TOIs at 2P. |
| Rejected automatically | 64 (eclipsing-binary flags, sinusoids, few deep events, crowding-tracking depths, season inconsistency, robust-clip failures, shared periods across stars, or a rotation-saturated star with other flags) |
| Flagged: rotation-saturated star only | 1 (TOI-807) |
| **Survivors to deep vetting** | **16**, plus the 1 flagged |

The robust clip drops the deepest 5% of transits (at least 1) and compares the SNR left over with what a constant signal of the same SNR would retain in white noise (Monte Carlo). A candidate is rejected only below that floor's 1st percentile. This replaces round 4's "few deep events" rule, which killed real long-period planets.

## 5. Deep vetting of the survivors

For every survivor (`round5_deep.py`, `deep_vet.py`, `round5_tces.py`, `run_triceratops.py`):

- **Pixel PSF localisation** on up to 24 TESScut sectors. The Gaia DR3 scene is fitted per sector, and difference images are taken at the candidate ephemeris with the known planets' transits excluded. A joint point-source fit then gives the chi^2 of "dip on the target", Δχ² relative to the free best fit.
  - **The star's own known planet is run through the same pixels as the on-target control.** Round-4 calibration says on-target planets give Δχ² <= 13. This round's controls give 0-12.6 on most stars, but **263 on TOI-4342 and 65 on Kepler-968**. There the test is uninformative, and I say so.
  - Re-running the same control with different windows excluded moved TOI-4572 b from 2.7 to 19.3. Treat Δχ² differences under ~10 as noise.
- **Fixed-aperture FFI depth vs 2-min depth, compared with the same ratio for the known planet on the same pixels.** A blend's depth changes with the aperture; an on-target planet's ratio matches its sibling's.
- Odd/even and secondary eclipse in the fixed aperture, and a transit fit with free stellar density.
- SPOC TCEs from every DV run (single- and multi-sector).
- **TRICERATOPS FPP/NFPP.** 2 runs of 5x10^5 draws, with Gaia neighbours and SPOC optimal apertures, no contrast curves and no pixel-cleared neighbours. That makes these conservative, upper-end FPPs.
- **Literature check** of the TOI/CTOI tables, SPOC TCEs, the NASA archive, arXiv and journal pages. ADS full text was not available.

### Ranked candidates

PSF Δχ²: candidate / known-planet control on the same pixels (on target if <= ~13). Ratio σ: fixed-aperture/2-min depth ratio compared with the control planet's (0 = same as an on-target planet). FPP/NFPP: TRICERATOPS mean of 2-3 runs. OE = odd/even (σ), sec = phase-0.5 depth (σ), both from the 2-min data.

| # | TIC (host) | P (d) | Depth (ppm) | SNR | Rp (Re) | Transits | OE / sec | PSF Δχ² | Ratio σ | FPP / NFPP | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 354944123 (TOI-4342) | 18.9083 | 1150 ± 119 | 9.6 | 2.21 | 6 (7) | 2.0 / 0.8 | 94 / **263 (uninformative)** | -0.7 | 0.12 / 0.004 | **NEW** (dynamically predicted) |
| 2 | 198385543 (TOI-1846) | 0.94406 | 257 ± 30 | 8.6 | 0.72 | 604 | 0.0 / -0.3 | 1.1 / 11.4 | +2.0 | 0.15 / 0.013 | **NEW** |
| 3 | 354944123 (TOI-4342) | 1.95769 | 284 ± 34 | 8.4 | 1.10 | 91 | 0.0 / 1.5 | 1.5 / 263 (uninf.) | +0.8 | 0.09 / 0.001 | **NEW** |
| 4 | 286201103 (TOI-6008) | 7.22015 | 1298 ± 164 | 7.9 | 0.92 | 31 | 1.3 / 0.1 | 11.3 / 4.1 | +0.8 | 0.09 / 0.008 | **NEW** (lowest SNR of the six) |
| 5 | 179985715 (TOI-249) | 19.0928 | 1121 ± 129 | 8.7 | 2.21 | 8 | 0.7 / -1.2 | 7.3 / 0.0 | -1.3 | 0.28 / 0.048 | **NEW, weak** |
| 6 | 177254988 (TOI-4405) | 3.60082 | 345 ± 39 | 8.8 | 1.39 | 170 | 0.2 / 1.0 | 12.5 / 1.4 | -2.4 | 0.10 / 0.010 | **NEW, weak** |
| 7 | 356158613 (TOI-1449) | 4.86736 | 461 ± 42 | 10.9 | 1.42 | 45 | 0.2 / -0.7 | 6.1 / 12.6 | +0.8 | 0.08 / 0.008 | known to SPOC (TCE, not a TOI) |
| 8 | 154940895 (TOI-4572) | 4.89923 | 579 ± 62 | 9.3 | 1.23 | 54 | 1.3 / 2.7 | 1.7 / 19.3 | +1.6 | 0.07 / 0.0000 | in a preprint at P/2 (Tschudi 2026) |
| 9 | 154940895 (TOI-4572) | 8.68309 | 751 ± 95 | 7.9 | 1.41 | 31 | 0.6 / -0.6 | 2.5 / 19.3 | -0.6 | 0.23 / 0.012 | known to SPOC (TCE, not a TOI) |
| 10 | 259233660 (TOI-6000) | 4.21349 | 891 ± 101 | 8.8 | 1.23 | 149 | 0.0 / -0.2 | 5.2 / 1.7 | +0.6 | 0.19 / 0.11 | ours (round 4, RNAAS): control |
| 11 | 453260209 (TOI-1215) | 4.71715 | 412 ± 45 | 9.1 | 1.15 | 47 | 1.0 / 0.6 | **16.0** / 2.6 | +0.6 | **0.67 / 0.56** | marginal (1.8″ twin star) |
| 12 | 453260209 (TOI-1215) | 23.0515 | 736 ± 90 | 8.2 | 1.53 | 7 | 1.3 / -0.1 | 1.4 / 2.6 | **-4.3** | - | doubtful: SPOC TCE, absent in FFI |
| - | 260004324 (TOI-704 / LHS 1815) | 2.47341 | 91 ± 10 | 9.2 | 0.52 | 364 | 1.2 / -0.1 | **62.9** / 11.5 | **+4.2** | - | REJECTED: blended EB |
| - | 321688498 (TOI-2290) | 2.52520 | 347 ± 45 | 7.8 | 1.15 | 87 | 1.1 / -0.8 | **20.0** / 12.2 | **-4.5** | - | REJECTED |
| - | 172370679 (TOI-1899) | 26.7989 | 2185 ± 285 | 7.7 | 3.15 | 5 | 0.3 / 0.1 | **53.1** / 10.8 | -1.2 | - | REJECTED (FFI odd/even 5.8σ) |
| - | 350814607 (Kepler-968) | 2.41507 | 885 ± 117 | 7.6 | 2.17 | 78 | 0.7 / -1.0 | 3.8 / 65 (uninf.) | 0.0 | - | REJECTED: Kepler would have seen it |
| - | 30853470 (TOI-807) | 19.9824 | 240 ± 31 | 7.7 | 1.17 | 45 | 0.7 / 1.3 | **27.5** / 9.7 | -0.4 | - | REJECTED: rotation-saturated star, off target, too long |

T0 values (BTJD) and every statistic are in `runs/round5_hosts/candidates.csv`.

| # | T0 (BTJD) | Duration |
|---|---|---|
| 1 | 1675.535 | 1.89 h |
| 2 | 1765.128 | 0.59 h |
| 3 | 1658.952 | 1.26 h |
| 4 | 1717.795 | 1.11 h |
| 5 | 1357.295 | 1.16 h |
| 6 | 2038.777 | 1.76 h |
| 7 | 1685.066 | 1.67 h |

### Per-candidate notes

**1. TOI-4342, 18.908 d, 2.2 Re (NEW; predicted).** Tey et al. 2023 found b (5.54 d) and c (10.69 d) in 2 QLP sectors; we have 9 sectors of 2-min data.
- **Prediction match.** Turtelboom et al. 2025 (AJ, DYNAMITE) predicted an additional planet at 18.5 (+2.7/-2.0) d, 2.24 Re, 1170 ppm. That is essentially what we see: 18.908 d, 2.21 Re, 1150 ppm.
- **For:** season-consistent (p = 0.72); the fixed-aperture depth matches the 2-min depth (ratio -0.7σ against c); FPP 0.12 and NFPP 0.004 with no contrast curve.
- **Weaknesses:**
  - Only 7 transit windows. One of them (S101) shows about zero depth. That could be a TTV, plausible near 7:4 with c, or a sign that the cycle count is wrong.
  - The robust-clip SNR (6.3) sits just above its floor (6.0).
  - The pixel test is uninformative on this star, because planet c itself "fails" it.
  - Transit-fit density 8.8 [4.2, 12.5] g/cc against 3.9 from the TIC (mild).
- The Parc et al. 2026 RV candidate at 47.5 d is unrelated. **Most interesting result of the round.**

**2. TOI-1846, 0.944 d, 0.72 Re (NEW).**
- **For:**
  - A sub-Earth USP candidate around an M3 dwarf: 604 transits in 23 sectors.
  - Every check is clean: season p 0.86, odd/even 0.0σ, no secondary.
  - On target: Δχ² 1.1 against 11.4 for TOI-1846 b on the same pixels.
  - The 5%-clip SNR is well above its floor.
  - Soubkiou et al. 2025 searched 0.5-40 d and found nothing. By their own injections, ~1 Re planets were below their sensitivity, so a 0.7 Re signal in 2 more years of data is not a contradiction.
- **Against:**
  - SNR 8.6 is near the noise tail.
  - The fixed aperture shows a 1.41x deeper transit (+2.0σ against b).
  - The 0.59-h duration is short: 0.58 of a central transit, and the fit density 32 [16, 54] g/cc is above the TIC 8.2. That points to a high impact parameter, or a wrong host if it is a blend.
  - FPP 0.15.

**3. TOI-4342, 1.958 d, 1.1 Re (NEW).** 91 transits; consistent across sectors; fixed-aperture depth consistent (+0.8σ); FPP 0.09, NFPP 0.0008. The pixel test is uninformative (same star as #1). If both #1 and #3 are real, TOI-4342 would be a 4-planet system: 1.96, 5.54, 10.69 and 18.9 d. Period ratios are 2.83, 1.93 and 1.77; none is a low-order alias of another.

**4. TOI-6008, 7.22 d, 0.92 Re (NEW).**
- **For:**
  - An Earth-size candidate around an M5 dwarf (0.24 Rsun), 31 transits, with every one of the 9 sectors positive (season p 0.66).
  - Odd/even 1.3σ, no secondary.
  - The fixed-aperture depth matches TOI-6008 b's ratio (+0.8σ).
  - Transit density 25.5 [12.4, 40.8] g/cc, consistent with the TIC 22.2.
  - FPP 0.09.
  - Barkaoui et al. 2024 searched fewer sectors with SHERLOCK and reported only a 16.35-d signal, later ruled out. 7.22 d is not an alias of it.
- **Against:**
  - SNR 7.9 is the lowest of the six: the noise tail predicts several such events in this survey.
  - The PSF Δχ² of 11.3 is at the top of the on-target range: the free fit wanders 22″, and the control gives 4.1.
  - The TIC contamination ratio is 0.79.

**5. TOI-249, 19.09 d, 2.2 Re (NEW, weak).**
- 8 transits in 5 sectors over 8 years, so cycle-count aliases are possible.
- 4 of the transits come from S104/S106, which no SPOC multi-sector run has searched yet, so SPOC may still alert it.
- **Concerns:** season p 0.04; fixed-aperture odd/even 3.1σ; Gaia RUWE 2.0 (an unresolved companion is likely); the duration is only 0.33 of a central transit. FPP 0.28, NFPP 0.05.

**6. TOI-4405, 3.60 d, 1.4 Re (NEW, weak).** 170 transits in 31 sectors. But the PSF Δχ² is 12.5 (control 1.4), and the fixed aperture gives only 0.58x the 2-min depth (-2.4σ against TOI-4405 b). FPP 0.10.

**7. TOI-1449, 4.867 d (not new).** This is the strongest signal of the round: SNR 10.9, clean, on target (Δχ² 6.1), FPP 0.08. But it is SPOC TCE 356158613-05 in the S14-S86 run (MES 7.18), never promoted to a TOI. The host is a K7/M0 dwarf; the 5780 K in ExoFOP is a placeholder because the TIC lacks Teff. A CTOI would rest on SPOC's own detection.

**8-9. TOI-4572.**
- **4.899 d:** on target (Δχ² 1.7), FPP 0.07.
  - Tschudi 2026 (arXiv 2607.23781) lists a "TOI-4572.02" at 2.4496 d, which is exactly half our period.
  - Our data prefer 4.899 d. The in-between epochs give 164 ± 64 ppm against 580 ± 62 ppm, so equal alternate events are rejected at 4.6σ.
  - The interpretation is open: either a 4.90-d planet, with the 164 ppm a 2.6σ fluctuation, or an eclipsing binary with a weak secondary.
- **8.683 d:** SPOC TCE 154940895-02 (MES 7.38), not a TOI, FPP 0.23. Not new.

**10. TOI-6000.02** was re-found at the round-4 values. Unconstrained TRICERATOPS gives 0.19 / 0.11, in line with round 4's unconstrained 0.09-0.27. With contrast curves and pixel clearing, round 4 got 0.0008.

**11-12. TOI-1215.**
- **4.717 d:** narrowly fails the on-target threshold (Δχ² 16.0 against control 2.6; per-sector target best 1/10).
  - The star is a 1.8″ pair of near-equal M dwarfs (same Gaia parallax) that no TESS pixel can split.
  - TRICERATOPS puts most of the probability on the companion and nearby stars: FPP 0.67, NFPP 0.56.
  - Marginal: it may orbit the companion, which would change the radius but not necessarily the reality.
- **23.05 d:** SPOC TCE 453260209-02 (MES 8.5). But the fixed-aperture FFI depth is only 265 ± 155 ppm against 736 ± 90 in 2-min (4.3σ low relative to the known planet's ratio), and the FFI aperture shows a 3.8σ dip at phase 0.5. Doubtful.

### Rejections worth knowing

- **LHS 1815 (TOI-704), 2.473 d, 91 ppm, 0.52 Re.** It would have been a tiny planet. In 2-min data it is clean: SNR 9.2, 364 transits, season p 0.92, odd/even 1.2σ, the strongest peak in the periodogram. The pixels disagree:
  - PSF Δχ² 62.9 against 11.5 for LHS 1815 b; the best fit is 20″ off target.
  - The fixed-aperture depth is 2.1x the 2-min depth (+4.2σ against b).
  - In the larger FFI aperture the signal shows odd/even at 7.6σ and a phase-0.5 eclipse at 19.6σ.
  - Verdict: a diluted eclipsing binary (true period probably 4.947 d) leaking into the SPOC aperture. **The fixed-aperture product comparison is what caught it.**
- **TOI-2290 2.525 d, TOI-1899 26.80 d, TOI-807 19.98 d:** off target in the pixels and/or the depth does not survive the FFI aperture.
- **Kepler-968 2.415 d:** an 885-ppm, 2.4-d transit would be unmissable in 4 years of Kepler data, and no such KOI exists. Probably a residual of the Kepler planets masked with decade-old ephemerides.

## 6. Caveats

1. **False-alarm rate.** About 10^5 effective trials per star across 184 stars. The SNR histogram of un-flagged signals falls by ~3x per 0.5 between 6.5 and 7.5. Extrapolating, I expect roughly 5-15 pure-noise signals at SNR 7.5-8.5 and about 1-4 above 8.5 in the whole survey. The six new candidates sit at 7.9-9.6. **Some of them are probably noise.** The next 2-min sectors (TOI-4342, TOI-6008, TOI-1846, TOI-249) are the independent test: SNR must grow as √N.
2. **No independent photons.** Every candidate was found and measured in the same 2-min data. The FFI fixed-aperture photometry shares photons with the 2-min data in overlapping sectors. It tests aperture dependence, not reality.
3. **Pixel localisation is calibration-limited.** The known-planet control fails outright on 2 of 12 stars. On the rest, control Δχ² ranges from 0 to 12.6, and the same control moves by about 15 depending on which windows are excluded. "On target" here means consistent with the target at TESS resolution (21″ pixels): companions inside ~10″ are not excluded. The candidate host stars include a 1.8″ binary (TOI-1215) and a RUWE 2.0 star (TOI-249).
4. **TRICERATOPS** was run without contrast curves and without pixel-cleared neighbours: 2 runs of 5x10^5, with run-to-run scatter of ~0.01-0.03. All FPPs except TOI-1215's are below the CTOI-level 0.5, and none is near the 0.015 validation bar. If ExoFOP has high-resolution imaging for these hosts (not checked star by star), the FPPs would fall, as TOI-6000's did from 0.09 to 0.0008 once contrast curves and pixel clearing were added. That was not done here.
5. **Completeness caveats.**
   - Injections went through the BLS, the SNR/transit threshold, the red-flag vetting and the saturated-star flag. They did not go through the statistical triage cuts or pixel vetting.
   - Injected b is capped at 0.8; grazing planets are not represented.
   - Rotation-saturated stars (22% of the sample) are essentially unsearched.
   - The 0.5-1 d range was searched but not injected.
6. **Search choices.**
   - One coherent BLS pass per star, then iterative masking of the peak pool. A weak second planet hiding under a first signal's sidelobes could be missed; the injections measure this, since they used the same procedure.
   - Up to 3 signals per star.
   - Known-planet masks use catalog ephemerides. Badly propagated ones (Kepler-era epochs) leave residuals, which is what killed the Kepler-968 signal.
7. **Literature.** The check covered the TOI/CTOI tables and SPOC TCEs (systematic), plus arXiv, journal and web pages. There was no ADS full-text access, so a passing mention inside a large catalog paper could be missed.
8. **Data volume.** 837 MB of shrunk light curves in `cache/r5/lc` and about 3 GB of TESScut cutouts were added to `cache/tesscut`. Free disk never dropped below 55 GB. Fetches used 6-8 workers.

## 7. What would move a candidate forward

- **TOI-4342 18.9 d:** fit TTVs jointly with c, check the S101 epoch for gaps, and add the TESS sectors after S103. A ground transit at ~1.2 ppt and 1.9 h is feasible for 1-m class telescopes on a Tmag 11 star.
- **TOI-6008 7.22 d and TOI-1846 0.944 d:** rerun TRICERATOPS with the ExoFOP contrast curves, and wait for new 2-min sectors. TOI-6008's 1.3 ppt depth on a Tmag 12.3 M5 star is reachable from the ground (MuSCAT-class).
- **TOI-1449 4.867 d and TOI-4572 8.683 d:** these are already SPOC TCEs. A CTOI would be a reasonable community contribution, but they are not discoveries.
- **House rule from round 4:** CTOI only after independent confirmation of signal reality. Nothing here meets that yet. The decision is Ricardo's.

## Files

| Path | What |
|---|---|
| `round5_hosts.py` | select / fetch / search / inject |
| `planet_hunt/hostsearch.py` | the per-star search |
| `round5_triage.py` | collect, triage, completeness |
| `round5_deep.py` | pixel specs, folds, plots |
| `round5_tces.py` | SPOC TCE scan |
| `round5_assemble.py` | evidence table |
| `vetting/round5_targets.json` | pixel-localisation specs |
| `vetting/round5_verdicts.json` | verdicts |
| `runs/round5_hosts/` | `targets.csv`, `stars.csv`, `signals.csv`, `triage.csv`, `candidates.csv`, `injections.csv`, `completeness.json`, `spoc_tces.json`, `search/` (one JSON per star), `inject/`, `deep/` (localisation and TRICERATOPS JSON), plots |
| Committed copy | `vetting/round5/` on branch `round5-hosts` holds the tables, completeness and plots (`runs/` is gitignored) |
| Plots | `tic354944123_18p908.png`, `tic198385543_0p944.png`, `tic354944123_1p958.png`, `tic286201103_7p220.png`, `tic179985715_19p093.png` (top 5 new), plus `tic356158613_4p867.png`, `tic154940895_4p899.png`, `tic260004324_2p473.png` (the LHS 1815 blend) and `injection_recovery.png` |
