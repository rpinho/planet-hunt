Every numeric occurrence in NOTE.md is listed in reading order, including the header, repeated numbers, citation years, identifiers, release/processing labels and the spelled-out counts. A range has a row for each endpoint. Exact means the same numerical value, including unit conversion; rounding is to the precision used in the unchanged note.

File aliases: R27 = [planets_r27_1102](planets_r27_1102/), R28 = [planets_r28_1102](planets_r28_1102/), R29 = [planets_r29_1102](planets_r29_1102/), R30 = [planets_r30_1102](planets_r30_1102/), R31 = [planets_r31_1102](planets_r31_1102/), R32 = [planets_r32_1102](planets_r32_1102/), R36 = [planets_r36_1102](planets_r36_1102/), R37 = [planets_r37_1102](planets_r37_1102/), R38 = [planets_r38_1102](planets_r38_1102/). JSON locations use JSON Pointer after #; L identifies a text line.

DERIVED.json is recomputed by verify_bundle.py from the copied inputs. Timing endpoints use the first/last cadences, while midpoint values use the cycle nearest each sector midpoint. The sector-40 uncertainty is the nominal epoch/period quadrature without covariance. The 200 EB-injection denominator includes one unavailable fit; none of the 199 valid fits reaches the real negative delta. Processing labels and citation years are documentary, not measurements.

| # | Number as written | Meaning | File and key/line | Match |
|---:|---|---|---|---|
| 1 | 275488151.01 | Header CTOI identifier | `R27/PLAN.json#/targets/12/ctoi` | Exact |
| 2 | 2026 | Header release year | `RELEASE.json#/year` | Exact packaging date |
| 3 | 11 | Header release month | `RELEASE.json#/month` | Exact packaging date |
| 4 | 02 | Header release day | `RELEASE.json#/day` | Exact packaging date |
| 5 | 275488151.01 | CTOI identifier | `R27/PLAN.json#/targets/12/ctoi` | Exact |
| 6 | 2022 | Montalto citation year | `R29/REPORT.md#L9` | Exact |
| 7 | 10014830 | EB KIC identifier | `R29/RESULT.json#/targets/9/eb_within_R/0/kic` | Exact |
| 8 | 275488133 | EB TIC identifier | `R31/REPORT.md#L94` | Exact |
| 9 | 3.0305267 | EB period (d) | `R36/PLAN.json#/EB/period` | Exact |
| 10 | 2011 | Prsa citation year | `R29/REPORT.md#L15` | Exact |
| 11 | 2011 | Slawson citation year | `R29/REPORT.md#L15` | Exact |
| 12 | 2016 | Kirk citation year | `R29/REPORT.md#L15` | Exact |
| 13 | 83 | Separation (arcsec) | `R29/RESULT.json#/targets/9/eb_within_R/0/sep_arcsec` | Rounding: 82.7 |
| 14 | 4 | Separation (TESS pixels) | `R36/PLAN.json#/sectors/0/separation` | Rounding: 3.95960 |
| 15 | 2.6 | Tmag difference | `DERIVED.json#/tmag_difference` | Rounding: 2.6225 |
| 16 | 13.8 | EB Tmag | `R31/outputs/metadata.json#/eb_tmag` | Rounding: 13.838 |
| 17 | 11.2 | Target Tmag | `R31/outputs/metadata.json#/target_tmag` | Rounding: 11.2155 |
| 18 | 43 | Kepler primary (%) | `DERIVED.json#/kepler_primary_percent` | Rounding: 42.712322 |
| 19 | 0.126 | Kepler half-window (d) | `R30/PROTOCOL.md#L45` | Exact |
| 20 | 11 | Catalog primary width (h) | `DERIVED.json#/catalogue_primary_width_hours`; `R38/REPORT.md#L55` | Rounding: 11.128094 |
| 21 | 0.210 | CTOI duration (d) | `DERIVED.json#/ctoi_duration_days` | Exact numeric value (trailing zero) |
| 22 | 5 | CTOI duration (h) | `R27/PLAN.json#/targets/12/duration_hours` | Rounding: 5.04 |
| 23 | 0.21 | ExoFOP literal duration (h) | `R27/REPORT.md#L59` | Exact numeric value (0.210) |
| 24 | 10014894 | Target KIC | `R27/PLAN.json#/targets/12/kic` | Exact |
| 25 | 2009 | Kepler start year | `DERIVED.json#/kepler_start_year`; `R32/outputs/apertures.json#/quarters/0/table_header/DATE-OBS` | Exact |
| 26 | 2013 | Kepler end year | `DERIVED.json#/kepler_end_year`; `R32/outputs/apertures.json#/quarters/17/table_header/DATE-END` | Exact |
| 27 | 18 | Quarters with EB outside | `R32/outputs/verification.json#/geometry_network/eb_outside` | Exact |
| 28 | 0.0 | Target primary depth (ppm) | `R32/outputs/completed_real.json#/results/primary/stack/depth` | Rounding: 0.008334754 |
| 29 | 2.6 | Target primary error (ppm) | `R32/outputs/completed_real.json#/results/primary/stack/sigma` | Rounding: 2.570569 |
| 30 | 424 | Target primary events | `R32/outputs/completed_real.json#/results/primary/stack/n` | Exact |
| 31 | 99.7 | Nominal Gaussian limit (%) | `R32/PROTOCOL.md#L7` | Exact convention |
| 32 | 7 | Target upper limit (ppm) | `R32/outputs/completed_real.json#/results/primary/stack/upper` | Rounding: 7.071696 |
| 33 | 0.56 | TIC contamination ratio | `R27/PLAN.json#/targets/12/tic_contratio` | Rounding: 0.564471543 |
| 34 | 873 | Expected on-target Kepler depth (ppm) | `R27/PLAN.json#/targets/12/expected_kepler_aperture_depth_ppm` | Rounding: 872.975121 |
| 35 | 873 | Injected box depth (ppm) | `R32/PROTOCOL.md#L26` | Exact; specified injection depth |
| 36 | 200 | Random injected phases | `R32/outputs/completed_gate.json#/denominators/injection` | Exact |
| 37 | 200 | Uninjected phases | `R32/outputs/completed_gate.json#/denominators/null` | Exact |
| 38 | 6 | Null upper-limit min (ppm) | `R32/outputs/summary.json#/gate_quantiles/null/upper/min` | Rounding: 6.161228 |
| 39 | 12 | Null upper-limit max (ppm) | `R32/outputs/summary.json#/gate_quantiles/null/upper/max` | Rounding: 11.662037 |
| 40 | 8 | Null upper-limit median (ppm) | `R32/outputs/summary.json#/gate_quantiles/null/upper/median` | Rounding: 8.083434 |
| 41 | 25 | DR25 release label | `R27/REPORT.md#L83` | Exact release label |
| 42 | 0.1 | Planned period tolerance (%) | `DERIVED.json#/planned_period_tolerance_percent` | Exact |
| 43 | 0.165 | CTOI period difference (%) | `DERIVED.json#/period_difference_percent` | Rounding: 0.164879 |
| 44 | 4.2 | Period mismatch (sigma) | `DERIVED.json#/period_difference_sigma` | Rounding: 4.198908 |
| 45 | 14 | EB timing sector, start | `R36/outputs/EB_TESS.json#/sectors/0/sector` | Exact |
| 46 | 15 | EB timing sector, end | `R36/outputs/EB_TESS.json#/sectors/1/sector` | Exact |
| 47 | 0.038 | EB timing advance (d) | `R36/outputs/EB_TESS.json#/sectors/0/timing/best_offset_days` | Exact magnitude; offset is -0.038 |
| 48 | 1.2 | TESS/Kepler depth ratio | `DERIVED.json#/tess_to_kepler_primary_depth_ratio` | Rounding: 1.192281 |
| 49 | 14 | Ephemeris drift sector, start | `R36/PLAN.json#/sectors/0/sector` | Exact |
| 50 | 15 | Ephemeris drift sector, end | `R36/PLAN.json#/sectors/1/sector` | Exact |
| 51 | 0.9 | Boundary timing after EB (h) | `DERIVED.json#/boundary_timing_hours/start` | Rounding: +0.885095; first cadence |
| 52 | 1.3 | Boundary timing before EB (h) | `DERIVED.json#/boundary_timing_hours/end` | Rounding: -1.256089; last cadence |
| 53 | +0.3 | S14 midpoint-cycle timing (h) | `DERIVED.json#/midpoint_cycle_timing_hours/14` | Rounding: +0.317863 |
| 54 | -0.8 | S15 midpoint-cycle timing (h) | `DERIVED.json#/midpoint_cycle_timing_hours/15` | Rounding: -0.761424 |
| 55 | 654 | Field stars | `R36/outputs/stars.json#/catalogue_rows` | Exact |
| 56 | 325 | Predicted radius-2 leak (ppm) | `R36/PLAN.json#/L/1` | Rounding: 325.360542 |
| 57 | 0.105 | Leak half-window (d) | `R36/PROTOCOL.md#L17` | Exact |
| 58 | 2 | Circular aperture radius (px) | `R32/outputs/attachment_review.json#/dv/page1/selected_aperture_pixels` | Exact |
| 59 | 58 | Leak / catalog depth (%) | `DERIVED.json#/leak_fraction_percent` | Rounding: 58.308341 |
| 60 | 558 | Catalog depth (ppm) | `R27/PLAN.json#/targets/12/depth_ppm` | Exact |
| 61 | two | The two depth estimators | `R36/REPORT.md#L24` | Exact count of the two compared estimators |
| 62 | 2 | LC2 processing label | `R36/REPORT.md#L26` | Exact label |
| 63 | 4 | Secondary aperture radius (px) | `R32/outputs/attachment_review.json#/dv/page2/aperture_radii_pixels/3` | Exact |
| 64 | 2538 | Radius-4 secondary depth (ppm) | `R36/outputs/real.json#/S/stack/depth/0` | Rounding: 2537.877918 |
| 65 | 350 | Radius-4 secondary error (ppm) | `R36/outputs/real.json#/S/stack/error/0` | Rounding: 350.417401 |
| 66 | 2235 | Predicted radius-4 secondary (ppm) | `R36/outputs/real.json#/S/prediction` | Rounding: 2234.811421 |
| 67 | 2 | Primary-dip aperture radius (px) | `R32/outputs/attachment_review.json#/dv/page1/selected_aperture_pixels` | Exact |
| 68 | 2 | Model-comparison aperture radius (px) | `R32/outputs/attachment_review.json#/dv/page1/selected_aperture_pixels` | Exact |
| 69 | 2 | LC2 processing label | `R37/REPORT.md#L17` | Exact label |
| 70 | 5 | CTOI box duration (h) | `R27/PLAN.json#/targets/12/duration_hours` | Rounding: 5.04 |
| 71 | 200 | EB injections; none as negative as real | `DERIVED.json#/lc2_eb_injection_count`; `DERIVED.json#/lc2_eb_injections_at_least_as_negative_as_real`; `R37/outputs/injections_summary.json#/trials/EB` | Exact: 0 of 200; 199 valid fits |
| 72 | 1.2 | LC2 EB amplitude / predicted leak | `R37/outputs/real.json#/bars/F/ratio` | Rounding: 1.193875 |
| 73 | 0.4 | LC2 injection-scatter error | `R37/outputs/real.json#/bars/F/sigma` | Rounding: 0.369996 |
| 74 | 1 | LC1 processing label | `R38/REPORT.md#L1` | Exact label |
| 75 | 3 | LC1 EB amplitude / predicted leak | `R38/outputs/real.json#/bars/F2/ratio` | Rounding: 3.045112 |
| 76 | 326 | Joint residual box depth (ppm) | `R38/outputs/real.json#/joints/2/combined/amplitudes/1` | Rounding: 326.191407 |
| 77 | 110 | Joint residual box error (ppm) | `R38/outputs/real.json#/joints/2/combined/sigma/1` | Rounding: 110.419838 |
| 78 | 3 | Joint residual box significance | `R38/outputs/real.json#/joints/2/combined/z/1` | Rounding: 2.954101 |
| 79 | 11 | Illustrative eclipse width (h) | `DERIVED.json#/catalogue_primary_width_hours` | Rounding: 11.128094 |
| 80 | 2 | LC2 processing label | `R38/REPORT.md#L5` | Exact label |
| 81 | 3 | Untested systematics timescale (d) | `R36/PLAN.json#/EB/period` | Rounding: 3.0305267; period scale, not a systematics measurement |
| 82 | 1 | Centroid target radii, start | `R32/outputs/attachment_review.json#/dv/page2/aperture_radii_pixels/0` | Exact |
| 83 | 3 | Centroid target radii, end | `R32/outputs/attachment_review.json#/dv/page2/aperture_radii_pixels/2` | Exact |
| 84 | 2 | Selected centroid aperture radius | `R32/outputs/attachment_review.json#/dv/page1/selected_aperture_pixels` | Exact |
| 85 | 4 | Other-source centroid radius | `R32/outputs/attachment_review.json#/dv/page2/aperture_radii_pixels/3` | Exact |
| 86 | 36 | Other-stars aperture light (%) | `DERIVED.json#/other_stars_light_percent` | Rounding: 36.080653 |
| 87 | eight | Later sectors including off-array cases | `DERIVED.json#/later_sector_count` | Exact count |
| 88 | 40 | Later sectors, first | `DERIVED.json#/later_sectors/0` | Exact |
| 89 | 82 | Later sectors, last | `DERIVED.json#/later_sectors/7` | Exact |
| 90 | two | Off-science-array sectors | `DERIVED.json#/off_science_sector_count` | Exact count |
| 91 | 40 | Expired-ephemeris example sector | `DERIVED.json#/later_sectors/0` | Exact |
| 92 | 0.3 | Propagated CTOI epoch error (d) | `DERIVED.json#/sector40_epoch_sigma_days` | Rounding: 0.284672; quadrature, cycle 239 |
| 93 | 10014830 | EB KIC identifier | `R29/RESULT.json#/targets/9/eb_within_R/0/kic` | Exact |

Untraced numbers: none. Original round-29 PROTOCOL.md and verify.py do not exist: its protocol is in PLAN.json; the provided pointer and saved-result verifier are explicitly packaging additions. Original full scientific verification has dependencies omitted from this small public bundle; see PACKAGING.md.
