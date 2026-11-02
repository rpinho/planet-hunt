Gate C2 **FAIL**; Bar T2 **UNSCORED**; Bar F2 **FAIL**. Real LC1 radius 2 after D descriptively favours the EB template: delta chi² = **+36.040808**. EB amplitude a/a0 = **3.045112**, 4.856876 injection-scatter sigma from 1. The requested attribution condition is **not met**.

**Protocol deviation: the literal pre-flux blind freeze was not achieved.** While inspecting local provenance, a shell command mistakenly displayed the full raw DIAmante file before PROTOCOL.md was committed. That display included LC1 flux. No LC1 fit, LC1 summary statistic, or flux-based parameter choice preceded the protocol. The deviation was disclosed immediately and recorded in the frozen protocol; the inherited random phase schedule predates the display. The protocol, new time mask and exact implementation were committed before deliberate numeric LC1 parsing. Numerical verification does not repair or certify blinding.

**Already-seen data disclosure.** The LC2 radius-2 results of r36 and r37 are seen; LC1 in-window flux had not been fitted in this lane before this run, but it is a different reduction of the same pixels. r37's verified LC2 radius-2 result had delta chi² -17.894316, with 0/200 EB-template injections reaching that negative value, and its calibration gate failed. The unchanged r37 fit was reproduced locally before constructing D: -17.894316088547043 on 15 events, identical within 1e-9 to its saved result. This lane uses the same models on LC1 and injects before a shared detrend; it is not an independent set of pixels.

**Verification.** Independent `verify.py`, importing no measurement, model or detrending code, returns **VERIFIED** for the numerical audit. It checks 413 profiles, 4472 events and 181120 fitting rows, plus 404 separately recomputed D series / 594319 detrending rows. Its centered-regression D reconstruction differs from saved trends by at most 2.66e-15 in normalized flux. It reconstructs source LC1, each multiplicative pre-D injection, acceptance, time masks, nuisance projection, joint algebra and the synthetic spline fit; regenerates the inherited time-only random phase stream; checks frozen bytes and source hashes; and recomputes all gate/bar decisions and joint products with 60-digit Decimal arithmetic. Floating-array checks use explicit tolerances in verify.py; this is a numerical audit, not an interval proof of physical transfer assumptions. `outputs/verification.json` explicitly records `prefreeze_blinding_requirement_met: false`.

**Gate C2.** Exactly 200 paired shifts, inherited unchanged from r36 through r37, were frozen before deliberate LC1 parsing. The random seed is 2026110237; the preceding R stream consumes 550 candidate draws and the E stream accepts its 200 phases after 263 candidates. The independent verifier regenerates both streams using timestamps alone. Each model is injected multiplicatively into original accepted LC1 radius-2 OOT flux, separately, before recalculating D and fitting both candidates. Errors are unchanged by injection and divided by the fitted trend by D. Strict positive delta prefers EB; strict negative delta prefers CTOI. No threshold calibration or retries.

| Injected model | Actual injections | Valid fits | Correct preference | Opposite / false preference | Unavailable or tie |
|---|---:|---:|---:|---:|---:|
| EB | 200/200 | 181/200 | 129/200 (64.5%) | 51/200 (25.5%) | 20/200 |
| CTOI | 200/200 | 181/200 | 164/200 (82.0%) | 17/200 (8.5%) | 19/200 |

The minimum actual-injection condition is met for both types. Required correct counts are at least 160/200 per type; required false counts are at most 10/200 per type. EB misses the correct-preference requirement and both types exceed their false-preference limit. All 400 D operations succeed; 19/200 trials per type have no eligible paired event after OOT/quality selection. They remain in the denominator. One additional EB trial is a valid-fit tie (both nonnegative amplitudes zero); no CTOI trial ties. At a zero preference threshold, correct and false outcomes are complementary among valid non-tied trials; the separate 5% false condition is retained exactly as requested. Actual injection counts use |fractional dimming|>1e-12, including measured out-of-eclipse template structure; they do not claim complete eclipse coverage in every trial.

EB trials: eligible events median 14, range 0–16; fitting rows median 533.5, range 0–726; nonzero injected cadences range 1463–1463. No missing row is interpolated into a fitting event.
CTOI trials: eligible events median 14, range 0–16; fitting rows median 533.5, range 0–726; nonzero injected cadences range 13–160. No missing row is interpolated into a fitting event.

**T2 and F2.** T2 is UNSCORED because C2 fails. The descriptive radius-2 chi² values are EB 727.467784757, CTOI 763.508592734, common null 964.302588211; delta +36.040807977, on 15 events / 778 fitting rows. Each event profiles an intercept and slope, and each candidate has one shared nonnegative amplitude. The EB sign preference alone does not meet the attribution rule.

F2 compares a/a0=3.045111504 with 1 using injection sigma=0.421075523. The scatter is the sample SD (ddof=1) from all 181/200 valid EB injections, without clipping, bias correction or selection by preference; mean=0.937337236, bias=-0.062662764. The discrepancy is 4.856876 sigma (>3), so F2 fails separately. The real EB mean-window-equivalent depth is 990.759 ppm; the single-box fit is 922.819 ppm. The EB coefficient multiplies a non-unit full template and is not itself a box depth.

**Joint descriptive fit.** After including the EB template, the radius-2 box coefficient is **326.191 ± 110.420 ppm**, z=**2.954101**. The joint EB coefficient is 0.00137769779798, box/EB coefficient correlation -0.807560, joint chi² 718.741069205. This uses the same event rows, weights and projected baselines. The signed unconstrained two-coefficient estimate and inverse-normal-matrix covariance are reported as frozen; the nonnegative joint optimum is identical in the combined radius-2 fit. The z is conditional on D and r37's local scatter weights, with no separate injection calibration for the joint diagnostic. It is not an extra attribution bar.

**Aperture and sector profiles, descriptive.** EB depth denotes amplitude times the fixed continuous template-window mean; CTOI is the single-model box depth. All fits use each aperture's own accepted LC1 data after its own D. Sector fits share amplitudes only within that sector; combined fits share amplitudes across sectors. No F2 bar is applied to radii 1/3 or sectors.

| Radius | Sector | Events | Rows | chi² EB | chi² box | delta chi² | EB equiv. ppm | Box ppm |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | combined | 15 | 778 | 749.0442 | 746.3414 | -2.7028 | 126.53 | 163.06 |
| 1 | 14 | 8 | 403 | 385.3613 | 387.1759 | +1.8146 | 197.09 | 158.03 |
| 1 | 15 | 7 | 375 | 362.4198 | 359.1589 | -3.2610 | 61.20 | 168.01 |
| 2 | combined | 15 | 778 | 727.4678 | 763.5086 | +36.0408 | 990.76 | 922.82 |
| 2 | 14 | 8 | 403 | 367.2145 | 379.7991 | +12.5847 | 1025.07 | 959.54 |
| 2 | 15 | 7 | 375 | 360.0219 | 383.4316 | +23.4098 | 962.80 | 890.73 |
| 3 | combined | 15 | 778 | 405.3483 | 912.1744 | +506.8260 | 4993.51 | 4323.07 |
| 3 | 14 | 8 | 403 | 219.7726 | 381.3711 | +161.5985 | 4659.28 | 4150.72 |
| 3 | 15 | 7 | 375 | 177.5708 | 528.6053 | +351.0344 | 5305.75 | 4492.61 |

| Radius | Sector | Joint EB equiv. ppm | Joint box ppm | Box sigma ppm | Box z | Coefficient correlation |
|---:|---|---:|---:|---:|---:|---:|
| 1 | combined | -9.56 | 170.89 | 103.78 | +1.6467 | -0.8073 |
| 1 | 14 | 208.29 | -13.54 | 154.01 | -0.0879 | -0.8253 |
| 1 | 15 | -180.76 | 314.97 | 140.78 | +2.2374 | -0.7903 |
| 2 | combined | 730.35 | 326.19 | 110.42 | +2.9541 | -0.8076 |
| 2 | 14 | 713.30 | 374.43 | 169.38 | +2.2106 | -0.8264 |
| 2 | 15 | 741.68 | 287.04 | 145.98 | +1.9663 | -0.7919 |
| 3 | combined | 4429.11 | 704.43 | 196.19 | +3.5906 | -0.8091 |
| 3 | 14 | 3892.99 | 915.90 | 296.57 | +3.0883 | -0.8338 |
| 3 | 15 | 4843.62 | 601.64 | 262.42 | +2.2927 | -0.7855 |

**Fixed models and D.** The exact form is frozen in PROTOCOL.md and implemented in detrend.py. EB P=3.0305267 d, T0 BJD 2454967.124438, full Kepler periodic interpolant and separately measured TESS sector offsets -0.038 d each. CTOI P=3.02553 d, T0 BJD 2458682.5504, duration 0.210 d (5.04 h). The complete interpolant is retained, including its secondary and out-of-eclipse structure. No target timing or shape search. a0=0.0006137423556910796 gives the specified 325.3605416475134 ppm primary-window leak prediction using the continuous template mean 0.5301256115543052, exactly as r37; this convention differs from r36's direct event-stack depth.

D is independent per sector, Gaussian local-linear regression with sigma=1.4 d, three Tukey reweight updates (four fits), evaluated on a fixed 0.25 d grid and linearly interpolated at existing cadences. Gaussian sigma exceeds three full primary widths (3 × 0.1530 × PE = 1.3910117553 d). Its training masks both EB eclipses using their full catalog widths plus 0.05 phase on either side, and the CTOI box plus half a duration on either side. Primary full width=0.4636705851 d (~11.13 h); secondary full width=0.41942489528 d. The smoother is recalculated for every injection, using shifted model masks and the OOT timestamps. Both flux and error are divided by its own trend. No smoothing parameter was tuned after LC1 fitting.

LC1 uses zero-based flux/error/flag columns 3/4/7, 10/11/14 and 17/18/21. Require flag=0 and finite positive flux/error, without clipping. Total timestamps: 2438; accepted radii 1/2/3: [2227, 2227, 2227]; accepted OOT radii 1/2/3: [1463, 1463, 1463]. Sector 14 accepted: [1123, 1123, 1123] / 1241; sector 15 accepted: [1104, 1104, 1104] / 1197. The timestamp-only OOT pool has 1602/2438 rows before quality acceptance. It intersects r36's original mask with full primary/secondary windows at the TESS offsets and CTOI windows, each padded by half the median cadence. This retimed intersection is an explicit frozen change from r37's OOT pool.

The unchanged r37 fitter retains ±0.60 d original paired-midpoint support, original two-box event eligibility and scatter mask, event out-flux median normalization, nonnegative single-model amplitudes, and weighted intercept/slope profiling on identical support for both models. Some broad EB wings remain in its scatter-estimation mask. D changes the input series; it does not change these fit conventions.

**Synthetic spline check.** Fixed before deliberate numeric LC1 parsing and saved before the first numeric-read receipt: f=1−a0*h_EB on all inherited timestamps, error=500 ppm, no added noise or target flux. The unmasked spline is a per-sector cubic least-squares spline with fixed internal knot spacing 0.5 d; flux/errors are divided by the fitted trend. The retained Montalto text specifies spline normalization for LC2 but not enough settings to reproduce its exact operator. This is therefore an illustrative spline check with one frozen setting, not a DIAmante LC2 reconstruction or a fitted explanation of r37.

| Synthetic processing | delta chi² | Change from unnormalized | Fitted EB a/a0 |
|---|---:|---:|---:|
| none | +21.979657875 | +0.000000000 | 1.000006466 |
| D | +21.981510120 | +0.001852245 | 1.000051237 |
| spline | +1.038918997 | -20.940738877 | 0.430032089 |

The D synthetic amplitude ratio differs from 1 by about 0.0051%; the illustrative spline retains about 43.00% of the input EB coefficient in the unchanged r37 fit. The pre-freeze raw display remains a chronology deviation for this check as well: it precedes deliberate parsing, not all exposure to LC1 values.

**Freeze and deviations.** Initial protocol commit `cb80016951fab9577b968ebd9b6d0e458ab7cfcc` contains PROTOCOL.md, PRIOR_ART.md and STATUS.md beginning FROZEN. Design/implementation commit `694bb003d63ecc65a26fba0566f9bfd30c54e79d` seals models, D, the inherited phase schedule and new OOT mask. Synthetic results were saved at 2026-11-02T14:47:35.316289+00:00; deliberate numeric LC1 parsing began after receipt 2026-11-02T14:47:39.657964+00:00. The inherited phase stream predates this lane. The raw pre-freeze display is the material protocol deviation; no model, width, timing, phase schedule, detrending parameter, gate threshold, decision rule or support criterion was changed after the real LC1 fit. No blocked dependency or external alert was needed.

Other limits: injections condition on an OOT source with gaps and only 181/200 usable paired fits; their cadence distribution differs from real in-window data. D is recomputed on each series but its masks and missing rows differ with phase, as specified. D estimation uncertainty is not added to formal errors; the empirical gate includes its nonlinear response. The local scatter mask and nonnegative boundary are inherited from r37. The joint z has only its stated conditional covariance. TESS offset-grid uncertainty, Kepler-to-TESS eclipse transfer, PRF/color/dilution uncertainties and DIAmante extraction transfer are not refitted or propagated. The upstream eigen-light-curve correction already exists in LC1; these injections test D and the fit after LC1, not the upstream extraction. C2 failure prevents a scored T2 inference regardless of its descriptive sign, and F2 fails independently.

**Prior-art record and credits.** PRIOR_ART.md records the local files inspected, reproduced baseline and all requested search categories that were not searched because the user prohibited network. No live scoreboard, public-repository, AI-lab/agent-result, arXiv or community query was made. Credit Montalto/DIAmante for the retained light curves, product documentation and methods; the Villanova Kepler EB Catalog team: Prsa et al. 2011, Slawson et al. 2011, Kirk et al. 2016; TESS/SPOC PRF for the inherited leak prediction; Kepler for the EB photometry/template; and TESS for the inherited sector timing data. No novelty, priority or importance judgment is made. No network access, message to another party, upload, or posting occurred.

**Reproduction.** Read-only inputs are hashed in DESIGN.json. The local scientific venv and all required libraries were already installed. From this directory:

```sh
./run.sh verify.py
./run.sh make_report.py
```

Historical full sequence, whose measurement stages refuse overwrite:

```sh
# PRIOR_ART.md first; inspect local prior work.
./run.sh baseline.py
# Write protocol + FROZEN status; initial local freeze commit.
./run.sh prepare.py
# Seal implementation/design in a second local commit.
./run.sh synthetic.py
./run.sh measure.py
./run.sh verify.py
./run.sh make_report.py
```

For a full rerun use a fresh local job directory and freeze, rather than removing these retained artifacts. run.sh fixes all numerical thread pools at 1 and actual nice at 10; one numerical worker runs at a time, within the two-core budget. Process snapshots before launches and at completion are in logs/. No worker is left running. Retained artifacts total about 39 MB, below 50 MB. Sufficient statistics, projected fitting rows, raw/injected/D series, smoother grids and robust weights are retained in outputs/. ARTIFACTS.sha256.json seals stable outputs and code; mutable verification timestamps/logs are excluded. Verification regenerates its own output. No source directory was modified.

The next step, if another experiment is requested, is a separate frozen design; no further selection, refit or adjusted gate is part of this report.
