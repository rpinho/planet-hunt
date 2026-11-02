# planets_r28_1102 — repaired periodic Kepler batch

**Verified saved-output result: 0/12 recovered archival leads, 0/12 mechanical refutations, 12/12 OPEN.** The same twelve periodic rows were 0 recovered, 0 refuted, 12 OPEN in r27. All five recall-qualified templates miss the >=0.9 coverage bar; the other seven fail recall. No threshold is relaxed. Stage 1 gate passes for 5/12 templates; all 1200 injections and all 600 frozen control slots are included. Overall matched-control recall is 954/1200; gating uses each template separately. Known native Kepler signals: 5/5 detected and localized in native unscaled flux. Only 1/5 passes a literal catalogue-central-depth versus fitted-box-depth comparison; this limits the known-planet validation.

A recovered CTOI is a lead, never a planet. DIAmante binary classifications and Boyajian’s Star context remain in force. No FPP, claim, note, posting, contact or ExoFOP write was made.

## Bars and gate

r27 batch bar, verbatim:

> report the recovered fraction and its denominator against 50 matched controls each. A CTOI is refuted if Kepler coverage >= 0.9 at its period with SNR_expected >= 15 and no event. It is recovered if the folded Kepler signal has z >= 7.1, depth within 3 sigma of the dilution-corrected CTOI depth, and exceeds the 99th percentile of its 50 controls. Anything else is open. Injection recall per target is reported beside each verdict; a refutation needs recall >= 0.9 at the frozen depth.

> Gate: injection recall >= 0.9 at the frozen depth on matched controls for every target, and >= 4 of 5 known planets recovered. A target whose recall is below 0.9 stays OPEN and is reported as such.

Principal recovered fraction: **0/12**. Among gate-passing templates: **0/5**. Validation uses two injections into each of 50 frozen matched controls per template: **1200/1200 trials, 600/600 controls**, no post-flux replacements. A failed injection gate vetoes both recovery and refutation. The known gate uses native detection/localization; physical catalogue central depths and fitted fixed-box depths are distinguished below.

## Per-CTOI results

Coverage is covered/predicted nominal-ephemeris transits; minimum is the least coverage over the entire registered period grid at that epoch. The 99th percentile is linear interpolation over all 50 null search maxima. Controls passing means a z/depth-compatible event, before the empirical null threshold.

| CTOI | KIC | r27 → r28 | Recall | Nominal coverage | Grid min | Expected SNR | Selected z | Control q99 | Control matches | Context |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---|
| 27233531.01 | 11768756 | OPEN → **OPEN** | 90/100 | 87/108 (0.806) | 0.757 | 474.53 | 4.573 | 12.730 | 0/50 | DIAmante binary=1 |
| 63069576.01 | 10274090 | OPEN → **OPEN** | 96/100 | 737/903 (0.816) | 0.808 | 482.26 | 3.654 | 12.683 | 0/50 | DIAmante binary=0 |
| 63120478.01 | 7118558 | OPEN → **OPEN** | 95/100 | 575/670 (0.858) | 0.838 | 561.18 | 4.355 | 7.679 | 0/50 | DIAmante binary=0 |
| 120355606.01 | 3934997 | OPEN → **OPEN** | 90/100 | 278/314 (0.885) | 0.869 | 99.40 | 2.698 | 8.657 | 2/50 | DIAmante binary=1 |
| 120963388.01 | 2695129 | OPEN → **OPEN** | 77/100 | 460/528 (0.871) | 0.834 | 973.45 | 4.006 | 26.900 | 3/50 | DIAmante binary=0 |
| 122673562.01 | 1718360 | OPEN → **OPEN** | 76/100 | 1065/1201 (0.887) | 0.882 | 9.44 | 8.948 | 17.215 | 0/50 | not DIAmante |
| 123125197.01 | 7254760 | OPEN → **OPEN** | 86/100 | 100/110 (0.909) | 0.909 | 38.60 | 2.315 | 7.494 | 0/50 | DIAmante binary=1 |
| 158984285.01 | 10528564 | OPEN → **OPEN** | 60/100 | 1449/1624 (0.892) | 0.876 | 74.54 | 5.293 | 34.486 | 0/50 | DIAmante binary=0 |
| 185336364.01 | 8462852 | OPEN → **OPEN** | 64/100 | 301/339 (0.888) | 0.866 | 2800.97 | 36.115 | 19.866 | 0/50 | not DIAmante; Boyajian’s Star, prior FP |
| 275488151.01 | 10014894 | OPEN → **OPEN** | 74/100 | 417/482 (0.865) | 0.849 | 430.19 | 3.943 | 25.432 | 3/50 | DIAmante binary=0 |
| 351189201.01 | 10873968 | OPEN → **OPEN** | 94/100 | 1432/2101 (0.682) | 0.679 | 374.00 | 2.777 | 10.225 | 0/50 | DIAmante binary=0 |
| 399793914.01 | 3424541 | OPEN → **OPEN** | 52/100 | 593/664 (0.893) | 0.877 | 25.27 | 3.504 | 19.940 | 2/50 | DIAmante binary=0 |

| CTOI | P (days) | Duration (hours) | Expected depth (aperture ppm) | Selected depth ± fit sigma (ppm) | Matches | Empirical rank |
|---|---:|---:|---:|---:|---:|---:|
| 27233531.01 | 13.61082 | 3.648 | 2003.531 | 19.519 ± 4.269 | 0 | 29/51 |
| 63069576.01 | 1.61694 | 2.496 | 4064.077 | 30.467 ± 8.338 | 0 | 31/51 |
| 63120478.01 | 2.17725 | 1.368 | 2344.868 | 17.997 ± 4.133 | 0 | 17/51 |
| 120355606.01 | 4.64364 | 2.712 | 738.772 | 20.122 ± 7.459 | 0 | 42/51 |
| 120963388.01 | 2.76544 | 5.808 | 1290.368 | 5.366 ± 1.340 | 0 | 37/51 |
| 122673562.01 | 1.2156 | 7.92 | 1199.293 | 981.084 ± 109.648 | 1 | 6/51 |
| 123125197.01 | 13.21099 | 5.928 | 1328.098 | 177.382 ± 76.608 | 0 | 50/51 |
| 158984285.01 | 0.89849 | 2.616 | 1466.265 | 109.114 ± 20.615 | 0 | 18/51 |
| 185336364.01 | 4.316 | 7.92 | 17828.594 | 231.337 ± 6.406 | 0 | 1/51 |
| 275488151.01 | 3.02553 | 5.04 | 872.975 | 7.912 ± 2.006 | 0 | 35/51 |
| 351189201.01 | 0.69463 | 1.32 | 692.172 | 5.274 ± 1.899 | 0 | 34/51 |
| 399793914.01 | 2.19755 | 4.224 | 153.608 | 21.595 ± 6.163 | 0 | 28/51 |

Missing epoch/period errors preserve r27’s conservative refutation veto. CTOI 122673562.01 lacks ephemeris errors; 185336364.01 lacks epoch and errors. Their all-phase searches can recover a lead if the recall gate passes, but cannot refute. The latter’s tabulated nominal coverage/expected SNR use the fixed Kepler midpoint as a descriptive placeholder. Missing depth errors contribute zero catalogue uncertainty; none is invented.

Why each OPEN stays open:

- 27233531.01: coverage below 0.9; no depth-compatible significant event.
- 63069576.01: coverage below 0.9; no depth-compatible significant event.
- 63120478.01: coverage below 0.9; no depth-compatible significant event.
- 120355606.01: coverage below 0.9; no depth-compatible significant event.
- 120963388.01: injection gate 77/100 <90/100; coverage below 0.9; no depth-compatible significant event.
- 122673562.01: injection gate 76/100 <90/100; missing ephemeris information vetoes refutation; coverage below 0.9; expected SNR below 15; matching event does not strictly exceed the control 99th percentile.
- 123125197.01: injection gate 86/100 <90/100; no depth-compatible significant event.
- 158984285.01: injection gate 60/100 <90/100; coverage below 0.9; no depth-compatible significant event.
- 185336364.01: injection gate 64/100 <90/100; missing ephemeris information vetoes refutation; coverage below 0.9; no depth-compatible significant event.
- 275488151.01: injection gate 74/100 <90/100; coverage below 0.9; no depth-compatible significant event.
- 351189201.01: coverage below 0.9; no depth-compatible significant event.
- 399793914.01: injection gate 52/100 <90/100; coverage below 0.9; no depth-compatible significant event.

## Diagnosis and repaired method

Synthetic Kepler-span ablation with all twelve frozen templates: r27’s selected fits pass its depth gate for **3/12**; refining only period/phase, while keeping its original SG filter, noise and depth gate, gives **12/12**. Fits at the injected truth recover the frozen depth. Old selected depths can be roughly 2–10% low despite very large z, while the nearest-period exact-phase error is much smaller. Removing only the epoch phase restriction does not repair these failures. This establishes grid quantization plus the strict depth gate as a sufficient failure mechanism in the reproduction; it does not uniquely apportion every historical raw-data failure among all implementation details.

A separate fixed-mask diagnostic on four preserved deep-injection trials gives injection-only transfer gain 0.99998923–1.00002262, while background box-depth offsets span -26.07 to 6.59 ppm at fit errors around 1.8 ppm. Those residual strict-depth failures persist even at the injected truth; they are not lost phase-grid peaks. The earlier normalization hypothesis did not explain these cases. Such failures stay in the recall denominator; no noise threshold is inflated to force a pass. See results/depth_confounds.json.

Paired actual-control ablation, the first frozen control for each of the first three templates, two identical raw injections per template: r27 finds a depth-compatible event in **0/6**, repaired search **6/6**. All injected periods lie within r27’s own interval; target flux is never used in this diagnosis. Exact old-template fits at the injected truth separate attenuation from search mismatch. Full evidence: results/diagnosis_synthetic.json, diagnosis_refinement.json, diagnosis_controls.json.

r27’s historical injections were on target flux; this repair’s registered recall population is matched controls, so historical and new rates have different data populations. The paired same-control experiment above isolates the algorithm comparison. The repaired search scans P ±5σP, or ±0.1% where uncertainty is absent, and every phase. Period steps accumulate <=0.1 duration drift across 1459.5 d; phase bins are <=0.025 duration. No period-grid cap and no finite shortlist. Exposure-integrated fixed-duration boxes use heteroscedastic weighted fits with a floating constant. Every independent grid maximum above z=7.1, plus the global maximum, gets raw-time period/phase refinement and its own transit-masked baseline. All templates, Amendment 1 durations, dilution ratios, targets and ordered control IDs are retained from r27.

Detrending is median-initialized quadratic SG with window max(2 days, 8 durations), negative-excursion masks, then periodic masks. Flux is divided by the masked continuum before fitting fractional depth; additive subtraction was rejected during repair because it mismatches multiplicative injections. Interpolation across masked transits keeps their flux out of the baseline. Adjacent-difference robust noise is inflated by duration-scale robust scatter, fixed from uninjected control flux for injections. The native Astropy BLS baseline was reproduced before repaired-flux validation; the first repair-code draft preceded that executable check, an ordering deviation recorded in DEVIATIONS.md. Compiled analytic fit derivatives agree with independent finite differences, and 100 parameter tests verify exact zero-template cadence pruning against the full fit. These are numerical checks, not a new detector or performance record. Full method: PROTOCOL.md and PLAN.json.

## Known-planet gate and its depth limitation

| Confirmed KOI / KIC | Matched CTOI | P (days) | Catalogue depth (ppm) | Native box depth ± sigma (ppm) | z | Located | Strict box-depth gate |
|---|---|---:|---:|---:|---:|---|---|
| K02017.01 / 8750043 | 399793914.01 | 2.2950351 | 172.1 | 195.69 ± 10.83 | 18.08 | yes | pass |
| K00439.01 / 12470954 | 63120478.01 | 1.9022086 | 2343.5 | 2030.65 ± 6.65 | 305.58 | yes | fail |
| K00961.01 / 8561063 | 158984285.01 | 1.2137704 | 1821.0 | 1656.04 ± 17.00 | 97.44 | yes | fail |
| K00874.01 / 7134976 | 120355606.01 | 4.6018118 | 755.9 | 664.13 ± 10.66 | 62.30 | yes | fail |
| K00952.05 / 9787239 | 351189201.01 | 0.74295791 | 843.3 | 292.02 ± 17.88 | 16.33 | yes | fail |

These five hosts were chosen from r27’s frozen CONFIRMED KOI catalogue before flux, by log-period/log-depth proximity with distinct hosts; none was replaced. Native flux was not deepened. The known gate tests blind detection and period/phase localization, not equality of catalogue physical central depth and fixed-box average depth. Several native known signals fail that stricter equality, particularly K00952.05. Therefore the 5/5 detection result does **not** establish completeness for realistic transit shapes under the CTOI depth gate. Only the frozen-depth integrated-box injections validate that gate. The 50-control empirical rank is bounded below by 1/51; interpolation of a 99th percentile does not establish a calibrated 1% false-alarm probability. This limitation remains material to interpreting mechanical refutations.

## Verification, registration and operational record

Freeze commit `949cd501952ce059f963cb671f7d1f671e344938`, plan SHA256 `2979f756cc01766e1c378bf282b2f664379a0a1e328dfe23f5b925ab71aadcd0`, method SHA256 `6e431fb3a6b1f4a9ed0816ec62b37d50759e1315667a1921b2b67056a5ce08fb`. STATUS.md first line FROZEN and the committed protocol preceded every r28 CTOI target search. Target-start receipts are saved per target. The sole score was written at 2026-11-02T10:59:11.086383+00:00; no target was searched or scored a second time. All final control and known outputs match the frozen implementation hash. Earlier repair runs are retained under repair_trials/, never scored.

Independent verify.py imports no measurement or scoring code. It validates the output seal, frozen source rows/hash, all **1200/1200** injection decisions, **600/600** control scores/percentiles, **5/5** known localizations and **12/12** verdicts. Independent accepted-exposure coverage reconstruction: True; 3961482 predicted transits checked over the period grids. It checks the final decision inequalities with exact rational arithmetic on saved decimal measurements. An independent pre-target gate certificate covers all 1200 trials. It certifies saved numerical decisions, not independent astrophysical truth, PDC completeness or a second photometric extraction.

Verifier-only amendment: the original registered reader failed because the percentile-position scalar q shadowed the exact-rational constructor q. Commit ad313cf corrects that name and adds an explicit hash bridge to the archived original reader. Measurement code, thresholds, original output seal, freeze and sole scores.json are unchanged; audit/verifier_amendment.json and the final certificate record both reader hashes and the unchanged score hash. No target was re-searched and no second score was written.

Deviations are recorded in DEVIATIONS.md: fresh pod/dependency installation; short-segment masking fallback; I/O locks; pre-freeze exact-fit optimization and preserved partial validation; segment noise instead of pooled-quarter noise; native known-depth distinction; initial small code transfers through the authorized proxy rather than the prescribed uploader. Final science is unchanged after freeze.

Only the named pod performed flux searches. Four scientific workers, each at nice 10 with four Numba threads, BLAS/OMP one thread. Data come directly from public stpubdata S3 with content-length and SHA256 receipts. Compact outputs were hash-pulled before the deadline; no FITS was downloaded to the Mini. results/shutdown.json records stopped job workers and the pod left up.

## Reproduction

For local saved-output verification (lightweight arithmetic/certificate check):

```sh
../.venv/bin/python audit/freeze_check.py
../.venv/bin/python verify.py . --arithmetic-only
```

The pod contains the completed run; do not rerun score.py there. The following records the validation and subsequent measurement/scoring commands on the named pod. It is historical provenance, not an instruction to score this workspace again. Dependencies are pinned in results/runtime.json; proxy, nohup and CPU limits are documented in PROTOCOL.md:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=4
# Stage 1, before the target freeze
nice -n 10 python pod_tasks/numerical_check.py
nice -n 10 python pod_tasks/validate.py baseline
nice -n 10 python pod_tasks/validate.py known
nice -n 10 python pod_tasks/validate.py controls
# After local freeze commit 949cd50
nice -n 10 python pod_tasks/measure_targets.py
nice -n 10 python pod_tasks/seal.py
nice -n 10 python score.py .  # executed exactly once at the saved score time
nice -n 10 python verify.py .
```

The registered verifier initially failed after the score; the final verify.py command was then repeated with the documented reader-only amendment, without repeating any measurement or scoring command. Saved-output reproduction uses the first two local commands above and the retained data; the corrected reader requires audit/verifier_amendment.json and the archived frozen reader. A fresh photometry replication needs a separately authorized workspace/receipt and an explicitly new score ledger; its result must not be represented as the original once-only score. Detailed proxy/upload/pull scripts are in podtools/.

## Prior art, credit and next steps

The baseline is established [BLS](https://arxiv.org/abs/astro-ph/0206099), [Astropy](https://docs.astropy.org/en/stable/api/astropy.timeseries.BoxLeastSquares.html), and established [transit masking](https://wotan.readthedocs.io/en/latest/Interface.html); see PRIOR_ART.md for official sources, recent-update public-code/issue/PR checks, AI-agent repositories, recent papers and community checks. [Dallas et al.](https://arxiv.org/abs/2610.32949) already uses Kepler ground truth for TESS; its different-pipeline counts are not comparable to this panel. No novelty or record claim is made.

Credit the original CTOI submitters montalto, roche and lynch for these periodic rows, and the broader frozen source panel’s chiarello, french, dragomir and TSTPC working group; Marco Montalto and the [DIAmante](https://archive.stsci.edu/hlsp/diamante) contributors; NASA Kepler and TESS teams, NASA Exoplanet Archive, MAST/STScI and public data contributors; the authors of BLS, Astropy, SciPy and Numba. [DIAmante](https://arxiv.org/abs/2210.14559) already classifies the three binary-flagged rows. [Boyajian et al.](https://arxiv.org/abs/1509.03622) establishes the irregular dipping-star context. These results audit submitted signals and do not claim their discoveries.

Low-recall rows remain unresolved by this search. Before stronger astrophysical interpretation, a separately registered study would need realistic transit-shape/PDC injections, target-specific sensitivity and a justified Kepler/TESS depth mapping. Keep any recovered archival signal as a lead only.
