Gate R failed on the completed 200-injection-per-case supplement: L 0/200 (0%), T 0/200 (0%). **Bar A is UNSCORED. No real four-aperture primary stack or R41 was computed.** The test therefore supplies no scored aperture-growth statement about the primary dip's source.

Bar S is SCORED: radius-4 secondary 2537.88 ± 350.42 ppm, z=7.2424, versus the fixed leak prediction 2234.81 ppm, a 0.8649 sigma residual. Bar E is DESCRIPTIVE: delta chi2=-26.054196 favours the CTOI timing template in radius 2. The EB's own TESS primary eclipses are earlier than the Kepler folded template by approximately 0.038 d in each sector. These are measurements of timing and the dip's light, with no FPP or planet conclusion.

**Verification.** A separate `verify.py`, importing no measurement code, independently recomputes every original and supplementary gate decision and all bar decisions from saved events and profile sufficient statistics with 60-digit Decimal arithmetic. It also regenerates the random-phase designs from timestamps, checks input SHA256s and git freezes, and verifies the EB timing profiles. `outputs/verification.json` reports VERIFIED. This certifies saved arithmetic, chronology, and input identity; extraction and the astrophysical PRF/color model are not independently validated by that certificate.

**Freeze and inputs.** Protocol commit `ef8289d57501670a5d0e493a1ce4c6a8eb124485` preceded all DIAmante flux access. Prediction/design/implementation commit `a77d71fe11cdfc836535e71eae40ca3828aea3dd` sealed PLAN.json before the first flux read at 2026-11-02T14:14:57.454902+00:00. `STATUS.md` began FROZEN. Supplement commit `eeaeb010838b3f5c3787e05195b6ade058420610` registered the time-support correction after the initial OOT injections and before true-event measurements; this is the material deviation described below. The original protocol and outputs are unchanged.

The DIAmante attachment is [TIC275488151L-mm20221020.txt](https://exofop.ipac.caltech.edu/tess/get_file.php?id=887179), retained from r32 and copied byte-for-byte into `sources/diamante.txt`. SHA256 `a762f968e4a5a5c7653b83062f48a0f5481cd3386bd0f76d30731ac82e5778f7` matches r32. No network request was needed or made in this job. The retained light curve covers BTJD 1683.3671387–1737.387918385, 2438 cadences. r32's attachment review, the DIAmante header and product guide were inspected before freezing. The public [DV](https://exofop.ipac.caltech.edu/tess/get_file.php?id=887180) already gives selected aperture radius 2, duration 0.21 d, apstat 4.33, and centroid Mahalanobis ranks (1,1,1,0). These are recorded prior findings.

Fixed ephemerides throughout: EB P=3.0305267 d, T0 BJD 2454967.124438; CTOI P=3.02553 d, T0 BJD 2458682.5504; secondary EB phase 0.5007. The scored windows and injections were not retimed after the TESS comparison.

**Predictions before flux.** The r30 TESS/SPOC UPDATED_2.0 PRF code and four-grid interpolation were reproduced using retained official surfaces. All 654 TIC/Gaia field stars have Tmag and all were included, with no brightness cutoff; positions with proper motions were propagated as in r30. Predicted stellar light was summed in pixel-centre circular apertures at the target, with no DIAmante background/pipeline transfer fit. S14/S15 were weighted 0.509023790/0.490976210 by 1241/1197 timestamp counts before quality/flux access. Target–EB separations are 3.95960/3.97100 pixels. All masks, star positions/fluxes, throughputs and sector vectors are in `PLAN.json` and `outputs/stars.json`.

At the fixed ephemeris the retained Kepler PDCSAP data give primary depth 444456.59 ± 3511.31 ppm over 418 events, and phase-0.5007 depth 197108.96 ± 758.45 ppm over 416 events, using this job's ±0.105 d windows. These are mean-window depths, not central eclipse depths. Multiplying them by the predicted EB flux fraction gives L and secondary vectors below. T is the target dilution vector normalized to 558 ppm at radius 2.

| Radius (px) | L primary (ppm) | T primary (ppm) | L secondary (ppm) | Median OOT stack error (ppm) | Real primary mean-window depth |
|---:|---:|---:|---:|---:|---|
| 1 | 75.57 | 752.14 | 33.52 | 77.56 | UNSCORED |
| 2 | 325.36 | 558.00 | 144.29 | 94.42 | UNSCORED |
| 3 | 1793.55 | 453.40 | 795.41 | 148.28 | UNSCORED |
| 4 | 5039.23 | 326.86 | 2234.81 | 337.30 | UNSCORED |

Predicted R41: L=66.678983; T=0.434572. The T vector decreases because additional catalogue stars dilute an on-target signal as the aperture grows. It is the computed target prediction, rather than an assumed unit ratio. The table's OOT uncertainties use the 200 support-conditioned, uninjected placebo stacks; their median event count is 14, range 1–16. The expected L radius-1 depth is only 0.974 times its median OOT error, so the ratio denominator has limited precision. No table entry is an undisclosed real A depth. Radius-2 E profile depths below use a different, explicitly defined shared box fit.

**Measurement.** Use normalized LC2_AP1–AP4 (1-based columns 6,13,20,27), ELC2 columns 7,14,21,28, and flags 8,15,22,29. The documented flags combine TESS quality bits with the DIAmante bit 4096. Requiring flag==0 rejects all nonzero flags: 2227/2438 accepted in each aperture, 211 rejected, with no further clipping. A common acceptance mask was used across apertures for R. There are 1664/2438 timestamp-only OOT cadences, 1520/2438 after quality/finite-positive filtering.

Each event uses an error-weighted local linear baseline in [-0.504,-0.252] and [0.252,0.504] d and the arithmetic mean of 1-flux/baseline in ±0.105 d. At least three in and three out on each side are required. Noise propagates in-window flux and fitted-baseline covariance, with OOT residual cross-aperture correlation and marginal noise at least the quoted scale. Inverse-variance event stacks inflate each marginal variance by max(1,event chi2/(N−1)); covariance is scaled accordingly. There is no error deflation. Synthetic linear baselines cancel and multiplicative boxes recover their known means to <1e-12 before flux (`outputs/synthetic.json`).

R41 is d4/d1 for positive d1. Its delta-method error uses the full covariance and gradient (−d4/d1²,1/d1). Success requires R41>2 and (R41−2)/sigma_R>=3. A nonpositive denominator or unavailable stack is a failure. No contrast statistic or alternate ratio definition replaced this frozen statistic.

**Gate R and the phase-support deviation.** The original uniform phase design did not reject empty OOT injection trains. Its 200 nominal draws per case contain only 168 actual nonempty injections, with only 82 extractable stacks; the other draws remained failures. This does not meet the user's minimum of 200 actual injections, despite its conservative 0/200 detection counts. These original artifacts are retained, and are not described as 200 successful injections.

`SUPPLEMENT_PROTOCOL.md` corrects eligibility using timestamps only: seed 2026110237, uniform candidate phases conditioned on at least one event meeting the frozen in/out cadence counts. It consumed 550 candidates to obtain 200 R phases and 263 for 200 E phases. Independent verification regenerates the identical stream without reading observed flux or quality to select phases. This is a post-OOT-flux protocol correction, not a fully pre-flux registered second experiment. It changes neither the estimator nor any physical prediction/window/threshold, and cannot promote a failed original decision.

Completed R supplement: L **0/200** successes versus required >=160/200; T **0/200** false detections versus allowed <=10/200. Both have 200/200 actual nonempty injections and 200/200 extractable stacks; positive-denominator ratios exist in 182/200 L and 200/200 T trials. All 200 trials count. Gate R FAIL. Bar A UNSCORED. The flat-depth wording is not triggered because real R41 was deliberately not measured.

**Bar E.** Shared nonnegative 558 ppm box injections at each ephemeris, with paired random phase shifts, calibrate delta chi2=chi2_CTOI−chi2_EB. A separate linear intercept/slope per paired event is projected out on identical support for both fits; the shared box amplitude is free. The local residual scatter is never smaller than its formal scale. Power is conditional on 558 ppm, not on a subsequently fitted amplitude. The EB/CTOI thresholds are respectively max(0,190th sorted delta under CTOI) and min(0,11th sorted delta under EB), with strict >/< comparisons. Calibration false rates are in-sample, not held-out validation.

| Calibration | Actual EB/CTOI injections | Valid EB/CTOI fits | EB correct | CTOI correct | EB false on CTOI | CTOI false on EB | Power criterion |
|---|---:|---:|---:|---:|---:|---:|---|
| Original 200 draws | 175/177 | 152/152 | 145/200 (72.5%) | 142/200 (71.0%) | 5/200 (2.5%) | 2/200 (1.0%) | FAIL |
| Timestamp-support supplement | 200/200 | 199/199 | 184/200 (92.0%) | 184/200 (92.0%) | 6/200 (3.0%) | 6/200 (3.0%) | PASS numerically |

Both runs' thresholds are zero. The supplement reaches the requested numerical power/false-preference bars, but E remains **DESCRIPTIVE** under the registered supplement's conservative rule requiring both calibrations to pass; it does not silently promote the original underpowered experiment after OOT data access. Real delta chi2=-26.054196 over 15 paired events; chi2 EB=775.159684, CTOI=749.105488, from a common null chi2=808.800410. Fitted radius-2 box depths are 354.57 ppm (EB) and 469.70 ppm (CTOI). Descriptive timing preference: CTOI. No source location follows from that timing comparison alone.

**Bar S.** The expected radius-4 secondary is 2234.81 ppm. Before reading its true cadences, the median radius-4 error from the original valid OOT placebo stacks was 326.02 ppm (82/200 valid stacks), yielding expected z=6.854748>=3. Thus S is SCORED. The independent support-conditioned placebo median is 337.30 ppm, also leaving expected z above 3; it does not replace the registered expected-z value.

The real phase-0.5007 stack is **2537.88 ± 350.42 ppm**, 15 events, observed z=7.242443; measured-minus-predicted residual 0.864873 sigma. Event chi2=90.263555, dof=14, variance inflation=6.447397. The detected secondary is compatible with the L prediction at this precision. The primary aperture-growth statement remains unscored.

**EB-centred TESS comparison.** Extract radius-1 pixel-centre apertures centred on TIC 275488133 (Tmag 13.838) within r30's retained 15×15 TESScut S14/S15 stamps. The stamps themselves are target-centred; the photometric apertures are EB-centred and fully on-stamp. Use quality==0, finite aperture pixels and finite outer-border median subtraction, with propagated flux/background errors. No download is needed. This satisfies the EB-centred aperture request using existing public TESScut bytes. These are diluted aperture depths; different sectors/Kepler passbands and background subtraction can change their amplitude.

| Sector | Accepted/total cadences | Primary events | Primary depth (%) | Secondary events | Secondary depth (%) | Kepler-template O−C (d) | Timing events |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 14 | 1237/1241 | 9 | 49.655 ± 0.221 | 8 | 22.093 ± 0.200 | -0.038 | 7 |
| 15 | 1190/1197 | 7 | 59.412 ± 0.307 | 8 | 26.911 ± 0.279 | -0.038 | 7 |

Combined primary: 52.992 ± 1.204% (16 events; variance inflation 45.095). Combined secondary: 23.732 ± 0.598% (16 events). A folded Kepler template, local linear nuisance baselines and shared nonnegative template amplitude give an all-sector best grid offset -0.036 d (-51.84 min), from 14 events. Each sector separately gives −0.038 d (−54.72 min). The combined shared-amplitude fit differs by one 0.002 d grid step from the per-sector fits; it is descriptive. The formal delta-chi2=1 width is unresolved at this grid spacing and is not reported as a zero timing uncertainty.

The EB eclipses remain present at the extrapolated cycle sequence, with a measured displacement relative to the fixed Kepler template. They do not sit exactly at zero template offset. The inherited catalogue flags include ETV context; its listed epoch error is 0.091901 d and period error 0.0000052 d, and catalogue error covariance is not propagated here. White-noise profile chi2, sector-dependent dilution, Kepler/TESS contrast differences and the 30-minute integrations limit a precise timing interpretation. This descriptive timing check does not refit any scored ephemeris.

**Limitations and deviations.** The material phase-support error and its post-OOT correction are disclosed above, including actual-injection denominators. The supplement makes the requested 200 actual injections available while retaining the conservative registered A/E decisions. The on-target prediction is diluted rather than flat; its T[2]=558 ppm normalization is exact. Pixel-centre circular masks, a finite PRF surface, the finite inherited catalogue field, no explicit diffuse light model, time-count rather than inverse-variance sector prediction weighting, difference-image pipeline transfer, and transfer of Kepler eclipse contrast to the TESS passband are uncalibrated forward-model assumptions. Formal prediction uncertainty does not include those effects. No thresholds or fixed ephemerides changed after flux. The secondary is evidence compatible with EB light in radius 4, but it cannot substitute for the deliberately withheld primary aperture-growth statistic.

No generic scoreboard, GitHub/AI-lab repository, arXiv or community network search was run: those domains are outside this job's allowlist. `PRIOR_ART.md` records the permitted local searches and findings, rather than asserting absence elsewhere. No novelty or priority judgment is made. The shared venv instruction was superseded by the user's permitted read-only r30 venv; its NumPy/SciPy and local FITS reader were reused with no installation or external service. All computations were single-process, at actual nice 10, BLAS/OMP threads 1, within the two-core budget. After measurement, the launcher wrapper was made portable to an ordinary nice-0 shell by calculating the increment needed to reach nice 10; it preserves an inherited nice 10. This changes no measurement code or saved results. No data, messages or uploads left this machine, and no outside party was contacted.

**Reproduction.** Run from this directory:

```sh
./run.sh verify.py
```

`run.sh` uses the permitted r30 read-only CPython environment and sets OPENBLAS_NUM_THREADS, OMP_NUM_THREADS, MKL_NUM_THREADS, VECLIB_MAXIMUM_THREADS and NUMEXPR_NUM_THREADS to 1. It ensures actual nice 10 (including when invoked from a normal terminal) and refuses an inherited nice greater than 10. Verification reads the retained source hashes in r30/r32 and the local saved outputs. The historical extraction order is below; scientific stages refuse to overwrite sealed outputs and should not be rerun over this completed record:

```sh
# First: PRIOR_ART.md; PROTOCOL.md + FROZEN STATUS.md; local freeze commit.
./run.sh prepare.py
# Seal PLAN.json, time-only design and code in prediction commit.
./run.sh inject.py
./run.sh supplement.py design
# Commit SUPPLEMENT_PROTOCOL.md and supplementary phases.
./run.sh supplement.py measure
./run.sh real.py
./run.sh eb_tess.py
./run.sh verify.py
./run.sh make_report.py
# Final STATUS.md, artifact seal, then saved-output verification again.
./run.sh seal.py
./run.sh verify.py
```

Before every measurement launch the process list was checked for older workers. No older instance was relaunched concurrently. Saved events, baseline/cadence counts, covariance, profile normal equations, source hashes, phase designs and exact local commits are retained. The final artifact SHA256 manifest excludes mutable verification timestamps and logs; its bytes are checked by verify.py when present. Total retained size is below the RULES 50 MB cap. Final process inspection is saved under logs. A future primary aperture-growth measurement would require a separately frozen ratio experiment whose injected denominator is measured precisely enough; this record does not change its statistic to obtain a score.

**Credit.** Montalto/DIAmante for the light curve, validation document and established multiple-aperture/centroid methods (the locally retained primary descriptions are Montalto et al. 2020, arXiv:2008.09832, and Montalto et al. 2022, arXiv:2210.14559). The Villanova Kepler EB Catalog team: Prsa et al. 2011, Slawson et al. 2011, Kirk et al. 2016, for catalogue ephemerides, eclipse widths, phase separation and flags. TESS/SPOC for PRFs and FFI data; MAST/TESScut for cutout products; Kepler for the EB light curves. The r29–r32 local records supply the retained catalogue/attachment context and baseline implementations reproduced here.

Report generated 2026-11-02T14:29:49.303641+00:00.
