Gate C **FAIL**. Bar T **UNSCORED**. The real radius-2 comparison descriptively prefers the CTOI 0.210 d box: delta chi² = **-17.894316**. Bar F **PASS**: a/a0 = **1.193875**, injection sigma **0.369996** (199/200 valid EB fits). The attribution rule is **not met**.

The descriptive result is: **"the dip does not follow the EB's TESS-era eclipses"**, which counts against the NEB reading under the requested wording. Because Gate C fails, this preference has no scored Bar T interpretation. No FPP or planet statement is made.

**Verification.** Independent `verify.py` imports no measurement code. `outputs/verification.json` is VERIFIED. It regenerates the inherited random-phase stream (550 R candidates then 263 E candidates), checks all 409 profiles / 4956 events / 200890 saved fitting rows, and reconstructs source flux, multiplicative injections, timestamp support, normalization, weights and nuisance projections. It recalculates sufficient-statistic arithmetic, gate decisions, amplitude sample scatter, bars and secondary stacks with 60-digit Decimal arithmetic, and independently reconstructs secondary extraction. It checks source hashes, frozen git bytes and commit order. This certifies numerical calculations and input/design identity; it does not validate the astrophysical template/PRF assumptions. Floating-array comparisons use explicit numerical tolerances; this is not an interval proof of the physical model.

**Already-seen data disclosure.** The target radius-2 flux was already inspected in r36. Its old two-box Bar E reported 354.57 ppm (about 355 ppm) at the Kepler-extrapolated EB timing and 469.70 ppm (about 470 ppm) at the CTOI timing, delta chi² = -26.054196 (about -26.05). This job is a follow-up with fixed external EB timing/shape, not a new blinded target data set. The current job froze before new numeric target-flux access.

**Baseline reproduction.** The unchanged local copy of r36 `measurement.py` reproduces its real Bar E delta chi² exactly to <1e-9: -26.054196387, with 15 events. Old fitted box depths: 354.570409 / 469.702316 ppm. This reproduction precedes interpretation of the new real-model comparison and is saved in `outputs/baseline.json`.

**Freeze.** Protocol commit `6fe3d8d26acb760524fd701132f36dc51ba709a0` contains PROTOCOL.md, PRIOR_ART.md and STATUS.md beginning FROZEN. Model/design/implementation commit `3cd90b660f394978f547943044f0fc56e9cbbd52` seals DESIGN.json/npz and measurement code before first numeric target-flux access at 2026-11-02T14:36:27.369873+00:00. The design was generated from inherited timestamps, the r36 Kepler template, r36 EB-centred TESS fits and r30 sector metadata. PROTOCOL.md contains the gate and bars verbatim. No model, phase, support, threshold, or decision rule changed after flux.

**Fixed models and common fit.** EB P=3.0305267 d, T0 BJD 2454967.124438. Its complete 1500-bin periodic Kepler-folded template is interpolated without width/shape fitting; the primary is about 11 h wide. The separate EB-centred S14 and S15 fits in r36 EB_TESS.json both supply offset -0.038 d (-54.72 min). The combined -0.036 d fit is not used. The EB-aperture fitted amplitude near 0.97 is not multiplied in again. CTOI P=3.02553 d, T0 BJD 2458682.5504, box duration 0.210 d (5.04 h). Each model has one free nonnegative amplitude, as in r36 Bar E.

The template's continuous mean dimming over ±0.105 d about its primary is 0.530125612, giving a0 = 0.000613742355691 from the specified 325.360542 ppm leak mean-window prediction. This template mean differs from r36's direct Kepler event-stack mean 0.444456586, because the folded template and event estimator have different normalization and stacking. The fixed conversion explicitly enforces the requested 325.36 ppm mean-window model depth; it is not fitted to target flux. The full interpolated template includes its measured out-of-eclipse structure. Its scaling prediction carries no fitted PRF/color uncertainty.

LC2 radii 1–3 use the original zero-flag, finite positive flux/error rule: 2227/2438 accepted per radius; 1520/2438 accepted OOT rows. S14 and S15 accept 1123/1241 and 1104/1197 rows. No clipping. OOT pool is r36's timestamp mask excluding full catalogued EB primary/secondary widths and CTOI boxes with half-cadence padding; it is deliberately not retimed after learning the EB offset.

Both candidate models use exactly r36 Bar E support: ±0.60 d about the midpoint of each paired original EB/CTOI cycle centre. The original two-box in-count and out-count eligibility tests and scatter-estimation mask remain unchanged. On each event, flux/errors are normalized by out-flux median; formal errors are inflated by local residual scatter, never reduced. The same weighted intercept/slope nuisance basis is profiled against both shapes on all common fitting rows. One amplitude is shared across events. Delta chi²=chi²_CTOI−chi²_EB=gain_EB−gain_CTOI. No timing/shape parameters are searched in target flux. The new real comparison uses the same 15 events as the reproduced baseline.

**Gate C.** Exactly the 200 paired E shifts from r36 timestamp-only supplement design are reused. Separately inject EB at a0 and CTOI at 558 ppm multiplicatively into original accepted radius-2 OOT flux, leaving supplied errors unchanged. Shift both trial ephemerides by the same amount, and fit both candidates on identical support. Threshold is exactly zero: strict positive prefers EB and strict negative prefers CTOI; no calibrated thresholds or retries. All requested trials stay in the denominator, including unavailable fits.

| Injected model | Actual injections | Valid paired fits | Correct preference | Opposite/false preference | No preference / unavailable |
|---|---:|---:|---:|---:|---:|
| EB | 200/200 | 199/200 | 169/200 (84.5%) | 30/200 (15.0%) | 1/200 |
| CTOI | 200/200 | 199/200 | 184/200 (92.0%) | 15/200 (7.5%) | 1/200 |

Both correct-preference rates exceed the required 160/200 (80%). Both false-preference rates exceed the allowed 10/200 (5%), so Gate C fails. At a fixed zero threshold, correct and opposite preferences are complements apart from ties/unavailable fits: meeting 80% alone cannot satisfy the separate 5% condition. All 200 trials contain actual nonzero injections of each type; one paired fit per type is unavailable after quality acceptance and counts against power. Full-template EB injections include out-of-eclipse structure; actual-injection counts use |fractional injected dimming|>1e-12 and do not imply every event has a complete observed eclipse.

EB trials: fitted events median 14, range 0–16; fitting rows median 557.5, range 0–752. No absent rows are interpolated or replaced.
CTOI trials: fitted events median 14, range 0–16; fitting rows median 557.5, range 0–752. No absent rows are interpolated or replaced.

**Bar T.** UNSCORED because Gate C fails. Descriptive radius-2 chi²: EB 766.999804, CTOI 749.105488; common null chi² 808.800410. Delta chi² -17.894316 favours CTOI. There are 15 events and 778 fitting rows; nuisance dimension is 30, plus one model amplitude. These are the frozen TESS-era EB timing and full shape, not the earlier Kepler-extrapolated box comparison.

**Bar F.** PASS: a/a0=1.193874591; sigma=0.369995768 from the sample standard deviation (ddof=1) of all 199/200 valid EB injections, regardless of their preferred model. Injection mean=1.150843983, bias from 1=0.150843983, without bias subtraction. |a/a0−1|/sigma=0.523991≤3. The real EB model's mean-window-equivalent depth is 388.440 ppm; the fitted CTOI box depth is 469.702 ppm. The shared EB coefficient a=0.000732731403932 multiplies a non-unit shape and is not itself a box depth. F is evaluated separately even though T is unscored; it cannot establish attribution.

**Descriptive aperture and sector comparisons.** The mean-window-equivalent EB depth below is a times the same fixed template-window mean. No amplitude bar is applied to radii 1 or 3. All sector fits have their own shared amplitude and per-event baselines. Negative delta favours CTOI.

| Radius | Sector | Events | Fitting rows | chi² EB | chi² CTOI | delta chi² | EB mean-window equiv. (ppm) | CTOI box (ppm) |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | combined | 15 | 778 | 754.0239 | 752.1805 | -1.8434 | 70.33 | 108.58 |
| 1 | 14 | 8 | 403 | 391.6119 | 391.6119 | +0.0000 | 0.00 | 0.00 |
| 1 | 15 | 7 | 375 | 360.5978 | 357.2027 | -3.3951 | 148.92 | 219.88 |
| 2 | combined | 15 | 778 | 766.9998 | 749.1055 | -17.8943 | 388.44 | 469.70 |
| 2 | 14 | 8 | 403 | 386.4254 | 377.1325 | -9.2929 | 277.48 | 388.93 |
| 2 | 15 | 7 | 375 | 377.9153 | 370.4970 | -7.4183 | 474.95 | 537.24 |
| 3 | combined | 15 | 778 | 918.5998 | 921.7260 | +3.1262 | 2301.08 | 2316.70 |
| 3 | 14 | 8 | 403 | 450.3275 | 425.5021 | -24.8254 | 1829.51 | 1935.36 |
| 3 | 15 | 7 | 375 | 444.4972 | 480.1193 | +35.6220 | 2714.08 | 2669.07 |

**Radius-2 TESS-era secondary, descriptive.** Centre each event at EB phase 0.5007 plus its sector offset -0.038 d. Use r36 local-linear ±0.105 d mean-window extraction and original flank windows; inverse-variance stack and inflate for event scatter, never deflate. There is no secondary gate or pass criterion here.

| Sector | Events | Depth (ppm) | Error (ppm) | Observed depth/error | Fixed r36 L prediction (ppm) | (measured−predicted)/error |
|---|---:|---:|---:|---:|---:|---:|
| combined | 16 | 81.321 | 73.470 | 1.1069 | 144.292 | -0.8571 |
| 14 | 8 | 43.320 | 134.464 | 0.3222 | 144.292 | -0.7509 |
| 15 | 8 | 119.012 | 91.567 | 1.2997 | 144.292 | -0.2761 |

The combined stack contains 16 events (8 in each sector), event chi²=19.234525 for 15 degrees of freedom, variance inflation=1.282302. The prediction in every table row is the requested r36 combined L radius-2 secondary prediction, not a newly sector-refitted prediction.

**Deviations and limits.** No post-flux changes or protocol deviations were required. The frozen implementation choices matter: retain the original r36 midpoint support and two-box eligibility/scatter mask while changing the physical EB template; use the complete interpolant rather than replacing its 11 h eclipse by a box; define a0 by its continuous primary-window mean; leave the inherited OOT mask unretimed. EB wings may lie in the old scatter mask, and missing OOT regions truncate injected eclipse coverage. Those effects are included in the measured gate performance rather than corrected after seeing flux. The timestamp-conditioned supplement selects support, not noise or recovery; outcomes are conditional on its real retained cadence/noise distribution. The nonnegative boundary, correlated residuals, and a small number of sector events limit chi² interpretation. Sector EB offsets are 0.002 d grid estimates whose uncertainty is not propagated. Kepler-to-TESS shape/contrast, PRF/dilution and DIAmante extraction transfer remain forward-model assumptions. Injection scatter is phase-dependent empirical scatter, with the reported nonzero bias; it is not a fitted prediction error budget. Gate failure prevents the sign preference from satisfying attribution regardless of F.

**Prior-art record and credits.** PRIOR_ART.md records exactly what was searched and found locally and what was not searched because network was prohibited. The baseline is the existing r36 profile, template and supplement design, using r30 retained sector metadata and r32 retained product context. No live leaderboard, GitHub/AI-lab repository, arXiv, or community search was performed. No novelty, priority or importance judgment is made. Credit Montalto/DIAmante for the light curve and existing validation/multiple-aperture methods; the Villanova Kepler EB Catalog team: Prsa et al. 2011, Slawson et al. 2011, Kirk et al. 2016; TESS/SPOC PRF and TESS data; Kepler for EB photometry, and the retained MAST/TESScut products for the prior EB-centred TESS fits. All inputs were already on this machine. No network operation, upload, message to another party, or posting occurred.

**Reproduction.** From this job directory, verify all retained results with:

```sh
./run.sh verify.py
./run.sh make_report.py
```

The wrapper uses the shared local scientific venv, reaches actual nice 10, and fixes OPENBLAS/OMP/MKL/VECLIB/NUMEXPR/BLIS threads to 1. All numerical stages used one process; no more than two cores were available to this job. No package installation was needed. Sources in r36/r30/r32 are read-only, and all files created by this job are here. The historical ordering was:

```sh
# PRIOR_ART.md first; PROTOCOL.md + FROZEN STATUS.md, then local freeze commit.
./run.sh prepare.py
# Commit DESIGN.json/npz and implementation; save design_commit.txt.
./run.sh measure.py
./run.sh verify.py
./run.sh make_report.py
# Final SCORED status and artifact manifest, then verify once more.
```

`prepare.py` and `measure.py` refuse to overwrite sealed runs. For a full repeat, use a fresh local job directory with a fresh freeze rather than deleting these artifacts. Read-only dependencies and SHA256s are recorded in DESIGN.json; `outputs/first_flux.json` seals implementation bytes and resource settings. `outputs/injections.json`, `outputs/real.json`, `outputs/profile_certificates.npz` and `outputs/baseline.json` retain the gate, bars, per-event equations and raw/projected fitting arrays. `verify.py` recomputes decisions from those saved outputs and independently checks them against source bytes. ARTIFACTS.sha256.json seals the final stable artifacts; verification timestamps/logs are excluded from that mutable-output loop. Retained size is below 50 MB. Process snapshots are in logs, including final worker inspection.

The requested measurement is complete. Any later change to model, OOT support, noise treatment or decision threshold requires a separate frozen follow-up; no such change is part of this result.

Measurements completed 2026-11-02T14:36:29.203770+00:00.
