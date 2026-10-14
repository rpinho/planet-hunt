# Round 8: scaling the hold-out test to the crowd's discards

Branch `round8-crowd` (worktree `.`, from `round7-rnaas`), 2026-10-13. Commits are local, author rpinho; nothing was pushed. Compute ran at `nice -n 15`. GitHub, MAST, ExoFOP, the NASA Exoplanet Archive, arXiv, Zenodo, Gaia and Zooniverse Talk were only read: no posts, issues, uploads or ExoFOP/journal actions.

## Bottom line

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

## 1. Sources (stage 1, `round8/sources.py`)

Every table is a read-only shallow clone in `<data>/r8crowd/src/` (commit in the author column). Named, preregistered, retracted or showcased candidates were never re-tested as leads: Paras Chopra's 3, Avijit Thawani's 3, Aiden Huang's candidate and weak-candidate rows, Marchand's 10 certificate/watch-list stars, moosesteve21's 4, Valdez's "candidate" rows, and those named on GitHub, Zenodo or PHT in the community map.

| Stream | Source | What we took | Rows (stars) |
|---|---|---|---|
| E1 | Eumemic `all_signals*.csv.gz` (7 tables, @da722613, still HEAD on 10-13) | the SNR 6-9 band after the round-7 table cuts, minus the 120 round 7 tested | 35,287 (31,484) |
| E2 | Eumemic `vetted_signals*.csv.gz` | entered his LEO-Vetter stage (SNR ≥ 9) and failed at most one soft test (stellar-density or noise-like), never secondary, odd/even, V-shape or size | 1,435 (1,380) |
| E3 | Eumemic `candidates_allsky.csv` | rows he lists but does not count: tier C, tier L, "likely false alarm" | 760 (758) |
| C1 | Aiden Huang, Aidenn8/tess-planet-search @9257049 | "not vetted" and "below threshold" signals; explicit 2-min sectors per star | 443 |
| C2 | Paras Chopra, paraschopra/tess-planet-candidates @74d5624 | phase-2 signals with no hard flag | 118 |
| C3 | Avijit Thawani, avijit-thawani/tess-niyati-avijit @76ee534 | same pipeline, same cut | 396 |
| C4 | Thomas Marchand, Th0rgal/open-exoplanet-hunt @19e5e9e | best signal of each rejected scanned star, if only noise-like gates failed (SNR, leave-one-epoch-out, coherence, baseline) | 402 of 1,240 |
| C5 | Alejandro Valdez, alexvaldex/nasa-tess-planet-search @71e2553 | "low_score" rows without hard flags; explicit 1-3 sectors | 194 |
| C6 | Marc Balcells, marc252/exonit-results @3c47877 | single-sector strong candidates (L0, "do not submit"); P given to 0.01 d | 924 |
| C7 | Bruno Rahle, brahle/exoplanet-candidates @702e6fa | 185 unclassified Sector-1 dossiers, discovery hypothesis (S1 only) | 185 |

**Not used as sources:**
- **T16 and RAVEN.** T16's 10,091 candidates are published in ApJS, so they are claimed. RAVEN vets TOIs and SPOC TCEs, which are excluded by rule. Both were used in the novelty gate. T16's 2026 single-transit table has no periods, so it doesn't fit a fixed-ephemeris test (a duotransit job for later).
- **Comdex4, CptSyphilis/shadowcatch and moosesteve21** publish no table of unvetted signals.
- **brahle/exoplanet-db** is a catalogue of others' lists (novelty only).

## 2. Funnel and power ranking (stage 2, `round8/rank.py`, `round8/inventory.py`, `round8/catalogs.py`)

Exclusion used live snapshots fetched 2026-10-13 ~23:40 UTC: ExoFOP TOI (8,151) and CTOI (5,138), NEA pscomppars (6,138), 321,600 SPOC TCEs (S1-107), every Eumemic list and tier, and community claims. A match is the same star at the same period or an alias (×1, 2, 3, ½, ⅓, 3/2, 2/3, 4, ¼; ±0.3%).

| Step | Signals | Stars |
|---|---|---|
| All source rows | 40,144 | 35,877 |
| 0.5 ≤ P ≤ 30 d | 39,678 | 35,610 |
| − TOI 2, CTOI 6, confirmed 1, **SPOC TCE 434**, Eumemic list 290 (his promoted E2 rows), community claim 0 | 38,932 | 34,976 |
| Deduplicated across sources | 38,887 | 34,976 |
| MAST inventory (2-min ≤ S107, TESS-SPOC, QLP ≤ S104) for all of them | 38,887 | 34,976 |
| Some independent new data under the source's own labelling rule | 9,346 | – |
| Expected hold-out power E_hold ≥ 2 | 9,173 | 8,375 |
| Small star (TIC Teff < 4200 K or R★ < 0.7, R★ < 1) | 7,879 | 7,092 |
| **Tested: top 3,500 stars by score, all their signals** | **4,026** | **3,500** |
| … tranche 1 (ranks 1-1,500) / tranche 2 (1,501-3,500) | 1,738 / 2,288 | 1,500 / 2,000 |

- **Power ranking.** As in round 7: E_hold = SNR_src·√(N_new/N_theirs) and E_full = SNR_src·√(N_all/N_theirs).
- **Score (the 2×2):** E_full + ½E_hold, ×1.15 for R★ < 0.45 and ×1.1 for P < 3 d. Eumemic-listed E3 rows get ×0.8 because someone already lists them.
  - SNR is capped at 15 for ranking. Big signals that someone rejected (mostly variability and EBs) still enter if they have new data; they just don't crowd the top.
  - The cap and the hard-gate filters for C4/C5 were added before any hold-out result existed.
- **Who made the cut:** E1 3,713 signals, E2 129, C4 74, C1 66, C6 27, C5 7, C2 5, C7 3, E3 2.
  - Most community streams shrink to a handful because their stars are FGK (ExoNit, Valdez, Rahle) or already used every sector (Paras, Avijit, Aiden). There is little new data to test them with, so their weakness is not the crowd's blind spot.
- **Tested stars:** median Teff 3,550 K, R★ 0.45, T 8.4-14.6, P 0.50-30 d.
- **Not tested:** 3,592 lower-ranked stars of the small-star pool (median expected hold-out z below 3, where the injections give ≤ 22% completeness).
- **Where the new data came from:** 83% of the hold-out sectors (8,307 files) are QLP files with empty error columns, which Eumemic's reader and the Paras/Avijit pipeline both drop. The rest are SPOC 2-min S105-107 and sectors the explicit-sector sources didn't use. No TESS-SPOC or QLP light curves exist after S104 yet (checked on MAST 10-13).

## 3. The test (stages 3-5, `round8/holdout.py`), pre-registered

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

## 4. Results

| Verdict | Signals |
|---|---|
| **Survivor (pre-registered)** | **17** |
| Not reproduced: hold-out depth > 3σ below their depth | 2,190 |
| Inconclusive: hold-out too weak, or a flag | 1,814 |
| No usable data (their-era files missing) | 5 |

By stream: E1 9 of 3,713, E2 6 of 129, E3 1 of 2, C4 1 of 74. C1, C2, C5 and C6 (105 signals) gave none. C7's 3 had no usable Sector-1 data.

### How many pass by chance

- **Inverted light curves, the whole procedure:** 3 of 3,975 fixed-mode and 0 of 46 window-mode signals pass. That predicts **about 3 noise survivors among the 4,021** (Poisson 1σ range 1.3-5.9).
- **Random ephemerides:** the per-signal null rate of z_hold ≥ 2.5 with a consistent depth sums to 26.8 (window: 0.7).
  - That is a loose upper bound, because it skips the SNR_full, odd/even, secondary and season cuts that the inverted test includes.
  - Round 7's equivalent was 0.98 for 118 signals; the ratio per signal is the same.
- **Population:**
  - The 3,458 Eumemic E1 hold-out z have mean 0.06 and std 1.10 (KS against N(0,1): p = 0.03; 0.05 without the survivors).
  - 33 have z ≥ 3 against 11 expected for a σ = 1.10 Gaussian. So there is a small excess, a few dozen signals that are astrophysically real, sitting inside a noise population.
  - Inverted z: mean 0.03, std 1.12.
  - The real signals at their claimed depths would have a median expected z of 4.0 (tranche 2) to 5.5 (tranche 1).
- **Reading:** the 17 survivors are not mostly noise. All 17 turned out to be identifiable astrophysical signals (§5), so the ~3 noise passes are either among the "unresolved" or were never there.
  - The noise test controls noise, not astrophysics. 14 of 17 survivors failed on variability, blends or EB shape in the deep vet.
  - **A statistical pass alone is not a candidate.** Round 7 showed this too.

### Completeness (4,026 injections, `runs/round8_crowd/completeness.csv`)

| Expected hold-out z of the injected signal | n | Pass |
|---|---|---|
| < 2 | 311 | 2.6% |
| 2-4 | 842 | 22% |
| 4-6 | 1,028 | 56% |
| 6-10 | 1,236 | 71% |
| ≥ 10 | 335 | 76% |
| **All** | **3,752** | **50.5%** |

- **Losses:** z_hold too low 1,349, SNR_full 1,225, depth inconsistency 716, season 375, odd/even 290.
- **Ephemeris lock:** the their-data BLS locks on the injected period 83% of the time.
- **What it means:**
  - Had the crowd's signals been real at their claimed depths, about 2,000 would pass. 17 do, of which 3 are on-target planet-like.
  - Correcting for completeness, roughly 0.7% of the tested crowd discards are astrophysically real, and about 0.15-0.25% are plausible on-target planets.
  - "Not reproduced" means "not confirmed", not "proven false": at these strengths about 30-45% of real signals would be missed.

## 5. The 17 survivors, deep-vetted (stage 6)

Tools:
- `round8/diag.py`: folds; odd/even; phase 0.5 and any-phase secondary; aliases at P/2, 2P and 3P; per-sector and per-product depths; per-transit χ²; ephemeris from per-sector mid-times.
- `deep_vet.py`: TESScut 11×11 pixels; joint PSF localisation over every Gaia star within 63″; ring test; flux budget; fixed-aperture depths, including TESScut sectors no pipeline has processed.
- TRICERATOPS (`round6/run_tri.py`, 3×3-px FFI apertures, MOLUSC on the best).
- Gaia DR3 RV, RUWE, NSS and variability tables.

| # | Signal (stream) | P (d) | Depth (ppm) | Rp (R⊕) | z hold (exp.) | SNR full | Pixels | FPP / NFPP | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| 1 | **TIC 101653309** (E1) | 0.703442 | 542 ± 54 (fixed ap.) | **1.30** (1.21-1.51) | **4.9** (4.3); TESScut S105 z 2.0 | 9.7 | **on target**: best source in 5/6 sectors; nearest capable star 5.9″ T = 18.0 (needs 11% eclipse) fits worse by Δχ² 19, all others by ≥ 34 | 0.080 / 0.040 default; MOLUSC + cleared > 10″: 0.041 / 0.034 (5.9″ star); also clearing it: **0.0071 / 0** | **STRONG CANDIDATE, new** |
| 2 | **TIC 262843927** (E1) | 24.10896 | 1,201 ± 267 (fixed ap.) | **2.11** (1.93-2.31) | 3.3 (5.2) | 10.1 | **on target**: offset 0.0 ± 4.1″, best in 3/5 sectors | **0.027 ± 0.002 / 0.000** | **CANDIDATE, new**; second planet-like signal at 14.517 d (Eumemic follow-up list) passes too |
| 3 | **TIC 410418450** (E1) | 6.82022 | 3,135 ± 746 (fixed ap.) | **3.54** (3.32-3.78) | 3.2 (3.6) | 9.3 | **on target**: offset 1.2 ± 1.7″, best in 4/8 sectors | **0.013 ± 0.001 / 0.0035** | **CANDIDATE, new** (T = 13.2, faint) |
| 4 | TIC 317300777 (E1) | 2.82108 | 543 ± 92 | 1.46 | 3.2 (2.8) | 8.8 | **unresolved**: a T = 15.1 star 7.3″ away fits better (Δχ² 11.6) and would need a 1.0% eclipse | 0.078 / 0.012 (neighbour not cleared) | weak; not claimable |
| 5 | TIC 348760393 (E1) | 12.38849 | 673 ± 233 | 1.48 | 3.5 (3.9) | 8.2 | **unresolved**: offset 11.9 ± 10.9″; a T = 18.7 star at 12.6″ would need a 24% eclipse | 0.47 ± 0.20 / 0.057 | weak; not claimable. One PHT volunteer tagged a single "#transit" dip in 2018, no period |
| 6 | TIC 198787224 (E3) | 0.646728 | 13,500 | 10-15 (grazing) | 15.9 (16.3) | 47 | on target | **0.97 ± 0.03** / 0 | REJECTED: grazing EB (Rp/R★ 0.17, b ≈ 0.9); RUWE 2.38 and Gaia DR3 Acceleration7 solution (hierarchical). Eumemic tier C (RUWE) |
| 7 | TIC 421279377 (E2) | 2.05752 | 17,500 | 5.8 | 8.4 (7.2) | 14.6 | **unlocalisable**: the target is 12% of the scene's flux (contamination ratio 1.8); every sector failed the quality gate | – | REJECTED: V-shaped 1.5-2% dips, two unrelated periods |
| 8 | TIC 421279377 (E2) | 7.22689 | 20,600 | 6.3 | 3.5 (5.8) | 11.6 | same star | – | REJECTED: as above |
| 9 | TIC 447558713 (E1+E2) | 5.55760 | 1,840 | 2.2 | 6.3 (4.8) | 9.9 | **off target**: T = 15.7 star 14″ away, target worse by Δχ² 251, best 0/8 sectors; ring depth rises outward; fit ρ★ 0.07 vs 6.2 | – | REJECTED: blended EB |
| 10 | TIC 261117566 (E2) | 0.62318 | 1,400 | 2.4 | 3.0 (4.8) | 11.0 | **off target**: source 55-70″ away, target worse by Δχ² 22,263 | – | REJECTED: blended EB |
| 11 | TIC 201752790 (E1) | 2.51266 | 1,130 | 1.9 | 3.9 (5.7) | 8.2 | **off target**: T = 17.9 star 40″ away, Δχ² 783 | – | REJECTED: blended EB |
| 12 | TIC 423840582 (E1) | 16.16735 | 1,200 | 1.1 | 4.2 (4.2) | 9.1 | **off target**: T = 15.6 star 30″ away, Δχ² 628; fit ρ★ 0.4 vs 14.8 | – | REJECTED: blended EB |
| 13 | TIC 363119909 (E2) | 0.55376 | 3,400 | 3.5 | 7.4 (11.0) | 17.0 | – | – | REJECTED: ±7,000 ppm ~1.85 h oscillation; phase-0.5 brightening of 3,100 ± 216 ppm; SPOC TCE at another period |
| 14 | TIC 47393760 (C4) | 1.83171 | 320 | 1.2 | 3.6 (3.4) | 8.8 | – | – | REJECTED: ~3.7 h oscillation; phase-0.5 brightening of −302 ± 39 ppm |
| 15 | TIC 189827896 (E1) | 0.99623 | 2,450 | 3.1 | 3.7 (6.1) | 8.2 | – | – | REJECTED: 0.2 d, 3,700 ppm variability; phase 0.5 −2,074 ± 334 ppm |
| 16 | TIC 267546678 (E2) | 2.23255 | 3,800 | 4.2 | 3.2 (1.9) | 13.0 | – | – | REJECTED: oscillatory structure; new data is one noisy QLP sector |
| 17 | TIC 254643343 (E2) | 0.69971 | 384 | 0.8 | 3.4 (5.1) | 11.0 | – | – | REJECTED: broad dip with sinusoid ratio 1.7 (a sinusoid fits better than a box); duration 1.5× the star's |

z hold (exp.) is the measured hold-out z, with the value expected at the source's own depth in brackets. Rp for the on-target candidates is from the limb-darkened pixel-photometry fit (`deep_vet.py`); the others are box estimates.

## 6. The three leads

### TIC 101653309 (best)

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

### TIC 262843927

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

#### 14.52 d signal: pixel vet and FPP

Using the saved full-data ephemeris (P = 14.51686896 d, T0 = 1325.657972 BTJD, duration 2.444 h), `deep_vet.py` reused exactly the eight cached cutout sectors of the 24.1 d vet, with six yielding usable transits (S1, 3, 29, 61, 68, 97; S31 failed the existing quality gate and S95 had no usable transit). **The on-target origin is not established:** the joint difference-image centroid is 21.4 ± 20.8″ from the target (jackknife uncertainty), with target Δχ² 61.3; the target is preferred in five of six sectors, but S97 favours the T = 18.65 Gaia source 4673750838475299968 at 42.0″. Removing S97 in the built-in jackknife brings the free centroid to 3.7″, so this is an unresolved, sector-sensitive localisation rather than a secure blended-EB identification. Fixed-aperture depths average 753 ± 208 ppm and fail consistency (χ² = 22.79, 5 dof, p = 0.000371); S97 is −1,056 ± 470 ppm. `round6/run_tri.py`, using the original 3×3 apertures unchanged, 4-minute bins and three repeats of 10⁶ draws, gives **FPP 0.05125 ± 0.00162, NFPP 0** (error is the population standard deviation across repeats, as in the reference runner). No neighbours were cleared, no contrast curve or multiplicity boost was applied, and no MOLUSC variant was run because the 24.1 d reference had none. **Verdict: unresolved planet-like signal, not statistically validated**; its hold-out pass alone does not establish a second on-target planet. All light curves and cutouts were reused from cache, with zero sector-data downloads. Full results are in `runs/round8_crowd/tic262843927_1452_vet.json`, `deepvet/tic262843927_1452.json` and `tri/triceratops_tic262843927_1452.json`.

### TIC 410418450

- T = 13.23, 136 pc, Teff 3,482 K, R★ 0.583. Gaia has no RV statistics (too faint). RUWE 1.11.
- Hold-out: QLP S94 (empty errors, dropped), 2,021 ± 635 ppm, z 3.2, against 2,259 in his-era QLP S12-93.
- P = 6.82022 ± 0.00007 d; Rp 3.54 (3.32-3.78) R⊕. QLP dilutes it: 2,126 ppm against 3,135 in the fixed aperture.
- ρ★ 1.7 (0.8-2.6) against the TIC's 4.1, a 2σ tension (longer than expected; possibly eccentric, or a slightly larger star).
- Pixels: on target (offset 1.2 ± 1.7″).
- TRICERATOPS: FPP 0.013 ± 0.001, NFPP 0.0035 (TP 0.85, PTP 0.14). Formally at the FPP bar, but NFPP is above 0.001 and there is no MOLUSC run, imaging or RV.
- Novelty: clean, except Eumemic's raw table (`all_signals_faint_t2`, SNR 8.93).
- A sub-Neptune on a faint star; lower priority.

## 7. What's worth an RNAAS or CTOI (Ricardo's call; nothing drafted or sent)

1. **TIC 101653309: yes, on the same footing as TIC 144439305.** Two independent data sets nobody searched agree with Eumemic-era data. It is on target and binary-clean. Its default FPP (0.08) matches TIC 144439305's; its NFPP is higher only because of one faint neighbour 5.9″ away, which one seeing-limited ground transit would clear.
   - **Suggested order:** CTOI → short RNAAS (or a combined note with TIC 144439305) → preregistered ephemeris on Zenodo.
   - **Credit Eumemic** (raw table, SNR 8.08, unvetted).
2. **TIC 262843927: a good second.** Possibly a two-planet system: 14.5 d is Eumemic's follow-up signal, 24.1 d is our hold-out-confirmed raw signal. FPP 0.027 alone, before any multiplicity boost.
   - A joint note would have to credit Eumemic for the 14.5 d signal, and he may promote it himself any day. Lower novelty for 14.5 d; 24.1 d is ours to report.
3. **TIC 410418450: hold.** It is faint, the stellar-density tension is 2σ, and Gaia gives no RV. Revisit if new sectors land.
4. **Not claimable:** TIC 317300777 and TIC 348760393 (pixel-unresolved), and the 12 rejects.

## 8. Method papers a note would have to cite (ExoFOP CTOI rule)

Detection and vetting used only these published methods:
- **BLS:** Kovács, Zucker & Mazeh 2002, A&A 391, 369 (astropy implementation).
- **wotan** biweight detrending: Hippke et al. 2019, AJ 158, 143.
- **TRICERATOPS:** Giacalone et al. 2021, AJ 161, 24 (and Giacalone & Dressing 2020, ascl:2002.004).
- **MOLUSC:** Wood et al. 2021, AJ 162, 128.
- **Difference-image localisation:** Bryson et al. 2013, PASP 125, 889 (our PSF-scene version follows the same principle).
- **New this round (if a note uses them):**
  - Gaia DR3: Gaia Collaboration, Vallenari et al. 2023, A&A 674, A1; RVs from Katz et al. 2023, A&A 674, A5; NSS from Gaia Collaboration, Arenou et al. 2023, A&A 674, A34; variability classifier from Eyer et al. 2023, A&A 674, A13.
  - Multiplicity boost, if used for TIC 262843927: Lissauer et al. 2012, ApJ 750, 112; Guerrero et al. 2021, ApJS 254, 39.
  - Limb darkening: Claret 2017, A&A 600, A30.
- **Data products:**
  - TESS-SPOC: Caldwell et al. 2020, RNAAS 4, 201;
  - QLP: Huang et al. 2020, RNAAS 4, 204, and Kunimoto et al. 2021, RNAAS 5, 234;
  - TESScut/astrocut: Brasseur et al. 2019, ascl:1905.007;
  - TIC 8.2: Stassun et al. 2019, AJ 158, 138;
  - Astropy Collaboration 2013/2018/2022.
- **Credit (not a refereed method):** Eumemic, github.com/eumemic/exoplanets (`all_signals_faint_t1.csv.gz`, commit da722613), for TIC 101653309's and TIC 262843927's raw signals and the 14.5 d follow-up signal.
- **What doesn't meet the rule:** the hold-out test itself (pre-registered fixed-ephemeris depth in independent sectors) is ours and not published in a refereed journal. It is a simple box-depth measurement, but a note should present it as an analysis, not a published method, and lean on BLS, TRICERATOPS and pixels for the formal vetting.

## 9. Caveats

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

## 10. Resources and files

- **Disk:**
  - The shared caches we add to (`ffi/`, `r5/lc/`, `tesscut/`, `mast/`) grew 5.0 GB during the session; that includes the other agent's additions.
  - Our own inventory and TIC tables (`cache/r8crowd/`) are 50 MB, and the clones and catalogue snapshots (`<data>/r8crowd/`) 0.55 GB.
  - Free space stayed at 30-44 GB.
- **Network:**
  - 24,427 light-curve sector files were planned for the 3,500 stars (2-min files are fetched per star); a few hundred were already in the shared cache from earlier rounds.
  - Plus TESScut cutouts (up to 8 sectors) for the 11 deep-vetted stars, the MAST inventory for ~35,000 stars, the TIC rows and the clones.
  - Speed was 0.8-10 MB/s across the Starlink→Spectrum switch. Every download retries and resumes.
  - We did not re-download the ~500 MB of SPOC TCE bulk tables: the other agent had already pulled them, and we read its copies (`<data>/r8/novelty/tce/csv`) read-only. T16 and RAVEN were also read from its copies.
- **Code** (`round8/`):
  - `common.py`: paths, a resumable downloader, guards;
  - `sources.py`: stage 1;
  - `catalogs.py`: live exclusion tables;
  - `inventory.py`: MAST products;
  - `rank.py`: stage 2;
  - `holdout.py`: stages 3-5, pre-registered;
  - `completeness.py`;
  - `diag.py`, `trifold.py`, `run_tri_batch.sh`: stage 6;
  - `novelty.py`: stage 7.
  - Pixels reuse `deep_vet.py`; TRICERATOPS reuses `round6/run_tri.py`; MOLUSC comes from round 6's `.tools`.
- **Outputs** (`runs/round8_crowd/`):
  - selection: `src_*.csv`, `signals_all.csv`, `signals_excluded.csv`, `funnel_*.csv`, `power_all.csv`, `rank_pool.csv`, `targets.csv` (`targets_tranche1.csv`);
  - test: `holdout/*.json`, `holdout.csv`, `inject/*.json`, `completeness.csv`;
  - vetting: `diag/`, `deepvet/`, `tri/`, `gaia_dr3_*.csv`;
  - novelty: `novelty/*.json`.
- **Committed vs local:** the report, summary CSVs, survivor hold-out JSONs, diagnostics, pixel vets, TRICERATOPS and novelty records are committed. The 4,026 per-signal hold-out JSONs and 4,026 injection JSONs (36 MB) stay local in `runs/round8_crowd/holdout/` and `inject/`, summarised in `holdout.csv` and `completeness.csv`.
