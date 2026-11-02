# Kepler reverse-vetting of TESS community signals

**0/15 recovered leads, 1/15 refuted under the registered bars, 14/15 OPEN.** Batch recovery is **0/14**, or 0/12 with executable periodic templates. The special single CTOI 269033446.01 (KIC 9488989) is **REFUTED**. These are conditional archival-template outcomes; no planet, new candidate, FPP, or discovery claim is made.

Verification: **VERIFIED_VERDICTS**, all 15/15 verdicts independently recomputed from saved statistics. Coverage independently recomputed for 13/13 executable targets using a separate exposure-integral implementation. 1384 saved files passed the SHA256/size seal. This does not independently re-extract or validate the photometry.

## Registered bars (verbatim)

- Reverse-vetting bar for TIC 269033446 = KIC 9488989: freeze depth 3988 ppm and T14 4.7 h from the CTOI row; search Kepler Q1-Q17 PDCSAP with that template only. Report (a) every event with SNR >= 7.1 and depth within 3 sigma of 3988 ppm x (Kepler/TESS dilution ratio), and (b) the fraction of P in [1, 1460] d for which a transit would have fallen in Kepler data. Pass = >= 1 matching Kepler event (comb lead), or a coverage fraction >= 0.9 with no event (CTOI demoted as P < 1460 d false alarm). Control: the same template on 50 Tmag-matched Kepler targets on the same module.
- Batch bar for the other CTOIs: report the recovered fraction and its denominator against 50 matched controls each. A CTOI is refuted if Kepler coverage >= 0.9 at its period with SNR_expected >= 15 and no event. It is recovered if the folded Kepler signal has z >= 7.1, depth within 3 sigma of the dilution-corrected CTOI depth, and exceeds the 99th percentile of its 50 controls. Anything else is open.
- Injection recall per target is reported beside each verdict; a refutation needs recall >= 0.9 at the frozen depth.

## The single: TIC 269033446 / KIC 9488989

Frozen TESS depth 3988 ± 2974 ppm, duration 4.7 h and epoch BJD 2458713.905. Dilution ratio 1.05125405; expected aperture-equivalent Kepler depth 4192.40 ppm. Matching Kepler event groups: **0**. Uniform-period coverage: **94314/100000 = 94.31400%** on P ∈ [1,1460] d. Injection recall: **90/100** including gaps. Controls passing the SNR/depth gates: **5/50**; control maximum z 99th percentile 10.672. The single bar does not require exceeding this percentile.

The period-space fraction uses fixed TESS epoch and uniform midpoint quadrature in period, rather than a period posterior. A demotion under the requested bar does not exclude every period below 1460 d. Injection recall covers the input PDCSAP and this detrender; it does not measure upstream PDC preservation. The large submitted depth uncertainty makes the depth gate broad.

Matching-event list: **empty** for the registered search. The finite proposal limitation described below applies.

## Per-CTOI outcomes

Coverage below is nominal covered/predicted; minimum period-grid coverage is shown separately. The single uses the period-space denominator. z is the selected exact fit, which can be depth-inconsistent when no matching event exists. Controls show SNR/depth-passing stars out of 50, followed by their unrestricted maximum-z 99th percentile. Missing-template rows have no executed controls or injections.

| CTOI | KIC | Verdict | z | Controls SNR/depth pass/50; q99 | Recall | Nominal coverage | Min coverage | Expected SNR |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| 26493474.01 | 11137395 | OPEN: metadata | — | 0/50 executed | 0/100 executed; N/A | — | — | — |
| 26542657.01 | 12013550 | OPEN: metadata | — | 0/50 executed | 0/100 executed; N/A | — | — | — |
| 27233531.01 | 11768756 | OPEN | 4.295 | 0/50; 7.679 | 0/100 | 87/108 | 0.7757 | 446.35 |
| 63069576.01 | 10274090 | OPEN | 3.296 | 0/50; 4.654 | 0/100 | 737/903 | 0.8073 | 470.88 |
| 63120478.01 | 7118558 | OPEN | 3.679 | 0/50; 4.943 | 0/100 | 575/670 | 0.8376 | 633.34 |
| 120355606.01 | 3934997 | OPEN | 1.547 | 0/50; 3.503 | 32/100 | 278/314 | 0.8694 | 95.21 |
| 120963388.01 | 2695129 | OPEN | 2.952 | 0/50; 5.371 | 0/100 | 460/528 | 0.8393 | 908.77 |
| 122673562.01 | 1718360 | OPEN | 1.025 | 0/50; 3.504 | 100/100 | 1065/1201 | 0.8868 | 20.08 |
| 123125197.01 | 7254760 | OPEN | 0.533 | 0/50; 2.969 | 60/100 | 100/110 | 0.9091 | 15.24 |
| 158984285.01 | 10528564 | OPEN | 1.637 | 0/50; 5.043 | 1/100 | 1449/1624 | 0.8764 | 187.33 |
| 185336364.01 | 8462852 | OPEN | 29.671 | 0/50; 4.251 | 20/100 | 301/339 | 0.8879 | 2603.36 |
| 269033446.01 | 9488989 | REFUTED | 3.648 | 5/50; 10.672 | 90/100 | 94314/100000 | — | 25.40 |
| 275488151.01 | 10014894 | OPEN | 2.930 | 0/50; 5.607 | 0/100 | 417/482 | 0.8568 | 400.55 |
| 351189201.01 | 10873968 | OPEN | 2.444 | 0/50; 5.063 | 1/100 | 1432/2101 | 0.6794 | 360.89 |
| 399793914.01 | 3424541 | OPEN | 1.118 | 0/50; 4.030 | 67/100 | 593/664 | 0.8810 | 27.51 |

Frozen durations and fitted depths (aperture-equivalent ppm):

| CTOI | P (d) | Primary T14 (h) | Literal Montalto hours, descriptive only | Expected depth | Fitted depth ± σ | Depth gate |
|---|---:|---:|---:|---:|---:|---|
| 26493474.01 | 1.56000 | — | — | 160.47 | — | unavailable |
| 26542657.01 | — | — | — | — | — | unavailable |
| 27233531.01 | 13.61082 | 3.648 | 0.152 | 2003.53 | 19.17 ± 4.46 | FAIL |
| 63069576.01 | 1.61694 | 2.496 | 0.104 | 4064.08 | 28.05 ± 8.51 | FAIL |
| 63120478.01 | 2.17725 | 1.368 | 0.057 | 2344.87 | 13.64 ± 3.71 | FAIL |
| 120355606.01 | 4.64364 | 2.712 | 0.113 | 738.77 | 12.00 ± 7.76 | FAIL |
| 120963388.01 | 2.76544 | 5.808 | 0.242 | 1290.37 | 4.23 ± 1.43 | FAIL |
| 122673562.01 | 1.21560 | 7.920 | — | 1199.29 | 61.11 ± 59.64 | FAIL |
| 123125197.01 | 13.21099 | 5.928 | 0.247 | 1328.10 | 46.19 ± 86.72 | FAIL |
| 158984285.01 | 0.89849 | 2.616 | 0.109 | 1466.26 | 12.86 ± 7.86 | FAIL |
| 185336364.01 | 4.31600 | 7.920 | — | 17828.59 | 203.02 ± 6.84 | FAIL |
| 269033446.01 | — | 4.700 | — | 4192.40 | 723.96 ± 198.45 | PASS |
| 275488151.01 | 3.02553 | 5.040 | 0.210 | 872.98 | 6.32 ± 2.16 | FAIL |
| 351189201.01 | 0.69463 | 1.320 | 0.055 | 692.17 | 4.69 ± 1.92 | FAIL |
| 399793914.01 | 2.19755 | 4.224 | 0.176 | 153.61 | 6.27 ± 5.61 | FAIL |

Reasons and existing classifications:

- **26493474.01** — Missing frozen duration/depth/contamination. Missing Duration (hrs) Missing Period (days) Error No duration: requested template, 100 injections and 50 template-control scores unavailable; OPEN.
- **26542657.01** — Missing frozen duration/depth/contamination. Missing Duration (hrs) Missing Period (days) Error No duration: requested template, 100 injections and 50 template-control scores unavailable; OPEN. TIC contamination is null; no frozen dilution ratio.
- **27233531.01** — At least one recovery/refutation condition fails. Published DIAmante duration is 0.152 DAYS; CTOI gives same number in HOURS. Amendment 1 applied before flux: primary hours = 24 x CTOI value; literal hours descriptive only. DIAmante binflag=1
- **63069576.01** — At least one recovery/refutation condition fails. Published DIAmante duration is 0.104 DAYS; CTOI gives same number in HOURS. Amendment 1 applied before flux: primary hours = 24 x CTOI value; literal hours descriptive only. DIAmante binflag=0
- **63120478.01** — At least one recovery/refutation condition fails. Published DIAmante duration is 0.057 DAYS; CTOI gives same number in HOURS. Amendment 1 applied before flux: primary hours = 24 x CTOI value; literal hours descriptive only. DIAmante binflag=0
- **120355606.01** — At least one recovery/refutation condition fails. Published DIAmante duration is 0.113 DAYS; CTOI gives same number in HOURS. Amendment 1 applied before flux: primary hours = 24 x CTOI value; literal hours descriptive only. DIAmante binflag=1
- **120963388.01** — At least one recovery/refutation condition fails. Published DIAmante duration is 0.242 DAYS; CTOI gives same number in HOURS. Amendment 1 applied before flux: primary hours = 24 x CTOI value; literal hours descriptive only. DIAmante binflag=0
- **122673562.01** — At least one recovery/refutation condition fails. Missing Period (days) Error
- **123125197.01** — At least one recovery/refutation condition fails. Published DIAmante duration is 0.247 DAYS; CTOI gives same number in HOURS. Amendment 1 applied before flux: primary hours = 24 x CTOI value; literal hours descriptive only. DIAmante binflag=1
- **158984285.01** — At least one recovery/refutation condition fails. Published DIAmante duration is 0.109 DAYS; CTOI gives same number in HOURS. Amendment 1 applied before flux: primary hours = 24 x CTOI value; literal hours descriptive only. DIAmante binflag=0
- **185336364.01** — At least one recovery/refutation condition fails. Missing Transit Epoch (BJD) Missing Period (days) Error Boyajian star; already user FP; irregular dust dips known. No planetary novelty. CTOI submitter already marked FP.
- **269033446.01** — Single period coverage >=0.9, no event, recall >=0.9. Missing Period (days) Error
- **275488151.01** — At least one recovery/refutation condition fails. Published DIAmante duration is 0.210 DAYS; CTOI gives same number in HOURS. Amendment 1 applied before flux: primary hours = 24 x CTOI value; literal hours descriptive only. DIAmante binflag=0
- **351189201.01** — At least one recovery/refutation condition fails. Published DIAmante duration is 0.055 DAYS; CTOI gives same number in HOURS. Amendment 1 applied before flux: primary hours = 24 x CTOI value; literal hours descriptive only. DIAmante binflag=0
- **399793914.01** — At least one recovery/refutation condition fails. Published DIAmante duration is 0.176 DAYS; CTOI gives same number in HOURS. Amendment 1 applied before flux: primary hours = 24 x CTOI value; literal hours descriptive only. DIAmante binflag=0

## Denominators, sample and prior art

Pre-freeze snapshots: ExoFOP CTOI 5139 rows; TOI 8151; NASA cumulative KOI 9564; DR25 TCE 34032. Independent global TIC/KIC and Kepler long-cadence joins found **63 CTOI rows / 15 selected**, versus design 64/15. No selected host has a KOI, DR25 TCE, promotion flag, or current TOI. The missing original design join snapshot prevents assigning a cause to the one-row overlap difference. This run uses its audited 15, without changing the sample after flux.

Template-ready: 13/15. Requested injections 1500; executed **1300/1500**, 100 for each of 13 templates, and 0/100 for each missing-duration row. Requested control slots 750; frozen IDs 700; executed **650/750**, 50 per executable target. The missing-duration 26493474.01 has 50 unused frozen IDs; 26542657.01 has no frozen panel because its reference spring-module assignment is unavailable. These represent 650 distinct control stars; all executed IDs are distinct. All 15 remain in the principal denominator, including the two metadata-limited OPEN rows.

CTOI 26493474.01 and 26542657.01 lack duration; the latter also lacks contamination. Neither gets an inferred duration. CTOI 122673562.01 lacks ephemeris uncertainties; CTOI 185336364.01 lacks epoch. Their full-phase diagnostics cannot refute. Missing depth errors contribute zero catalog uncertainty, as frozen. Ten Montalto durations and errors were multiplied by 24 under Amendment 1 after numeric matching to DIAmante days. Literal-hours entries above are descriptive metadata; no second search uses them.

DIAmante already flags TIC 27233531, 120355606 and 123125197 as binaries (`binflag=1`). TIC 185336364 / KIC 8462852 is Boyajian’s Star, already CTOI-disposition FP with known irregular dips. A mechanical recovery would not override these established classifications. See [DIAmante source](https://archive.stsci.edu/hlsp/diamante), [Montalto catalog paper](https://arxiv.org/abs/2210.14559), and [Boyajian et al.](https://arxiv.org/abs/1509.03622).

The broad method is established: [Dallas et al. (2026)](https://arxiv.org/abs/2610.32949) uses Kepler ground truth for TESS. Their table-3 retention 665/748 and no-KOI surviving-event count 65/311381 describe their different TESS pipeline, not comparable scores here. The baseline reproduced here is official catalog federation. Full DR25/TGLC reproduction is outside this job. [Planet Hunters TESS](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/about/results), submitter sources and exact-identifier searches were checked before the freeze; absence of further indexed matches is not evidence of novelty. Full search provenance, inaccessible sources and method credits are in PRIOR_ART.md.

## Method, scope and limitations

Only final Kepler Q1–Q17 PDCSAP long cadence is used, quality exactly zero, finite positive flux/time, fixed BKJD [131.5,1591], no flux clipping. Normalize each quarter by its median. The registered detrender is a degree-2 Savitzky–Golay filter at the nearest odd 2-day window, split at >0.5-day gaps. Both data and exposure-integrated box template receive the same filter. Quarter noise is robust adjacent-difference scatter with frozen duration-scale inflation. Expected aperture depth is CTOI depth × (1 + TIC contamination ratio) × frozen median CROWDSAP; PDC crowding is not applied twice. The submitters’ upstream TESS correction conventions and bandpass differences are assumptions.

Periodic searches use the frozen period/phase grid and propagated epoch window, with a fixed shortlist of 16 strongest approximate and 16 strongest approximate depth-compatible proposals plus nominal ephemeris, fitted on actual times. Single searches use frozen local-maximum proposals and exact fits. This is a finite-shortlist search, not an exhaustive exact fit at every grid point. Frozen-depth signals can be missed or fail the depth gate even at high nominal expected SNR; measured injection recall is essential and is not adjusted after seeing results. Several periodic rows returned 0/100 injection successes despite expected SNR above 400. This is an observed sensitivity failure of the registered pipeline for those rows, so their null results do not refute the CTOIs. The precise contribution of phase-grid mismatch, shortlist selection and depth gating has not been separately established in this frozen run. No thresholds, grid, noise method, target or control IDs were tuned. Low recall prevents a refutation even if no event is found.

Coverage follows the registered accepted-exposure convention. Refutation uses both nominal and minimum period-grid coverage at nominal epoch, expected SNR at nominal ephemeris, and recall ≥0.9. This is not a proof over all continuously varying ephemerides. Injections are multiplicative into original input PDCSAP before this detrender, using fixed normalization/noise; gap injections count in the denominator. The implementation draws trial epochs uniformly over the fixed Kepler span, yielding approximately uniform folded phases for the short periodic rows, rather than exactly uniform orbital phase. It propagates the injection ephemeris acceptance window from that drawn epoch; this is not necessarily the same window width as back-propagation from the TESS epoch. Recall characterizes this frozen implementation. Recovered signals are leads only. Fifty controls give minimum empirical rank 1/51 ≈1.96%; interpolating a 99th percentile does not establish a calibrated 1% false-alarm probability.

Control matching summary:

| CTOI | Spring module | Max absolute ΔKp | Max absolute ΔTmag for single |
|---|---:|---:|---:|
| 26493474.01 | 17 | 0.2530 | — |
| 26542657.01 | 0 | — | — |
| 27233531.01 | 22 | 0.0360 | — |
| 63069576.01 | 17 | 0.0080 | — |
| 63120478.01 | 14 | 0.0290 | — |
| 120355606.01 | 4 | 0.0600 | — |
| 120963388.01 | 4 | 0.0740 | — |
| 122673562.01 | 10 | 0.0310 | — |
| 123125197.01 | 2 | 0.0860 | — |
| 158984285.01 | 12 | 0.0720 | — |
| 185336364.01 | 24 | 0.0570 | — |
| 269033446.01 | 24 | 0.4990 | 0.1251 |
| 275488151.01 | 18 | 0.0640 | — |
| 351189201.01 | 23 | 0.0430 | — |
| 399793914.01 | 4 | 0.6110 | — |

Full IDs, actual magnitude distances and primary-header module evidence are frozen in PLAN.json and catalogs/validated_headers.json. All recorded quarters, cadence counts, noise scales, fits, control maxima and injection truth/recovery records are retained under results/, coverage/ and receipts/. Public FITS are downloaded directly on the pod; fresh HTTP response lengths and all saved SHA256 receipts are checked. No TESS flux is analyzed.

## Freeze, deviations and verification

Registered commit **f9cd2d9**, ancestor 7ce2cbe, plan SHA256 `e016d9d16ed1aa8dc2a2e7722c2a04b469c27d544af4d7d904ec6c2fae7396e1`. FROZEN-v2 commit time 09:02:00 UTC predates first superseded flux at 09:02:50.886758 UTC. Implementation reference b8f07cb was committed 09:03:52 UTC, after that first old-run flux; implementation is not misrepresented as preceding it. The plan/protocol did precede flux. V2 measurements began after the sole-owner resume and passed synthetic validation first. Score time: **2026-11-02T09:34:53.899176+00:00**; one scores.json was written, with no rescoring.

Superseded v1 plan hash `895cfd3450ca25afd583640f1d72d5f650c91f85616c62be60148e94ee67d1c7` belongs only to quarantined outcomes. All old pod results were moved, never deleted, to results_superseded_literal. Old coverage, receipts, logs and FITS were likewise preserved separately. They were never scored or reused. audit/quarantine_manifest.json records file hashes, sizes and original numerical outcomes; DEVIATIONS.md preserves their statistics. An uncommitted quadratic-projection protocol draft was preserved and excluded. The working protocol was restored exactly to the registered text. Cache plan-hash guards and download length checks change provenance protection only.

Independent verify.py imports no measurement/scoring code. It checks the output seal, arithmetic for depths and norms, dilution, all 1300 injection success/localization decisions, 650 control maxima and percentiles, and all 15 verdicts. pod_tasks/verify_coverage.py independently reconstructed the coverage flags/counts from saved accepted times; its certificate is hash-bound to those arrays. audit/check_freeze.py checks registration ancestry and amendment isolation; audit/verify_sources.py verifies all nine frozen source hashes. Verification completed 2026-11-02T09:39:08.191969+00:00. This verifies the saved-output inference, not independent astrophysical truth or a second photometry implementation.

## Reproduction

Use the saved score and independent verifier without flux:

```sh
../.venv/bin/python audit/check_freeze.py
../.venv/bin/python audit/verify_sources.py
../.venv/bin/python verify.py .
```

For a fresh measurement reproduction, make a new results directory/workspace on the authorized pod rather than overwriting scored outputs. Restore PLAN.json, PROTOCOL.md, pod_tasks/ and metadata_headers.json from catalogs/metadata_headers.json to ../r27. Dependencies: Python, numpy, scipy, requests, astropy, numba; exact recorded versions are in results/runtime.json. The historical results/freeze_receipt.json and quarantine_manifest.json are required provenance inputs for the fresh working copy. Keep OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. The following commands run on the pod; use podtools/pod.sh for the interactive proxy. Do not rerun score.py in the completed workspace:

```sh
cd ../r27
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
nice -n 10 python pod_tasks/synthetic_check.py
nice -n 10 python pod_tasks/measure.py targets
nice -n 10 python pod_tasks/measure.py controls
nice -n 10 python pod_tasks/verify_coverage.py
nice -n 10 python pod_tasks/seal_outputs.py
python score.py .  # exactly once in a fresh unscored workspace
python verify.py .
```

Downloads use each frozen metadata S3 key at https://stpubdata.s3.amazonaws.com/. Pull compact products with podtools/pod_pull.py (SHA256 checked). No bulk flux is retained on the Mini. Raw measurement JSON, coverage arrays and receipts are retained in results/measurement_outputs.tar.gz; verify.py reads their exact bytes directly from that archive when individual files are absent. Frozen catalog bytes and auxiliary metadata are compressed in catalogs/retained_sources.tar.gz and retained_auxiliary.tar.gz; audit/restore_sources.py restores their original paths if needed. Worker shutdown and final retention evidence are in results/shutdown.json and audit/retention.json. The pod is left up.

## Credit and next steps

Credit the original CTOI submitters chiarello, french, montalto, roche, lynch, dragomir and the TSTPC working group; Marco Montalto and the DIAmante contributors; Planet Hunters TESS volunteers/authors; NASA Kepler and TESS, their mission teams and public-data contributors; NASA Exoplanet Archive; and MAST/STScI. This job audits their submitted signals and does not claim their discoveries.

Keep missing-metadata rows OPEN until authoritative durations/contamination are available. Treat any low-recall or ephemeris-limited OPEN row as unresolved by this test. Any future methodological improvement needs a new registered plan and separately labeled outputs; the present thresholds and outcomes remain fixed. Nothing was sent to ExoFOP, posted or communicated externally.
