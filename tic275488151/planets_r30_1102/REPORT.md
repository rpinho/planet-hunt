# TIC 275488151.01 versus Kepler EB KIC 10014830

**Frozen outcome: no attribution (0/1 CTOI). Bar A FAIL; injection gate FAIL; Bar B UNSCORED. Saved-output decisions independently VERIFIED.**

**Gate validity limitation:** the frozen injection exclusion was too narrow for this EB’s catalogued eclipses. Strict out-of-eclipse-only validation is INCOMPLETE; the retained 0/20 and 14/20 are the original failed registered trial counts. Bar B remains unscored.

The EB-ephemeris aperture stack is **171.50 ± 102.11 ppm, z=1.679567, 71 events**. The CTOI stack is **20.58 ± 71.37 ppm, z=0.288402, 72 events**. Thus z(EB)−z(CTOI)=1.391165; both z(EB)>=7.1 and the difference>=3 fail. This test does not establish that the CTOI is a planet, reject the EB blend hypothesis, or reproduce DIAmante's optimized photometry.

Gate: **target 0/20; EB 14/20**, against >=18/20 separately. All 40 frozen trials count. 28/40 reached an injected cube, 26/40 returned a valid position. Bar B was set UNSCORED before reading/fitting real eclipse locations; no real free-position fit was run. No thresholds, apertures, periods or epochs changed after flux. The frozen protocol additionally requires valid localization in all ten sectors; this extra conservative condition was never reached because the gate failed.

**Registration and data.**

Freeze commit `f991d095689455bd69b9461578970d1b8c0feae6` preceded every field-flux read. PROTOCOL.md contains the bars verbatim, full method and ten-sector list. Aperture SHA256 `b31f3ff290244fa4744158c9e0127f0d62508a0c09255af945d54d7be67d0566` was sealed before first flux; implementation hashes and chronology are in outputs/implementation_before_flux.json and outputs/first_flux_read.json. Metadata and all cutouts were retrieved with serialized public GETs; official TESScut delivery redirects to its MAST-managed staged-cutout storage are logged. Nothing was posted, uploaded or sent to ExoFOP; no FPP, note or novelty claim is made.

Target TIC 275488151 = KIC 10014894: TIC/Gaia RA 293.203758997702°, Dec 46.9094215256791°, Tmag 11.2155. EB TIC 275488133 = KIC 10014830: RA 293.176228612672°, Dec 46.9223884462716°, Tmag 13.838. Positions were propagated with catalogue proper motions from J2015.5 to each sector midpoint. Pixel separations below use the cutout WCS. Official SPOC UPDATED_2.0 PRFs, four-grid-point interpolation and the frozen catalogue-fraction ranking selected four pixels per sector without flux. Detector coordinates and minor grid-edge clamping are retained in aperture_manifest.json.

Quality==0, finite aperture values and a one-exposure expansion around momentum-dump bit 32 were required. Per-cadence finite-border median subtraction precedes aperture extraction. Local out windows [−0.504,−0.252] and [+0.252,+0.504] d define a linear baseline; the primary in window is ±0.126 d. Events require >=3 in and >=3 out on each side. Inverse-variance stacking inflates the formal error by event-to-event scatter. The EB all-event variance inflation is 10.086 (chi²=706.031, dof=70); these are mean-window depths, not central eclipse depths.

**Per-sector ephemeris results.**

| Sector | Accepted/total cadences | EB events | EB depth ± error (ppm) | z EB | CTOI events | CTOI depth ± error (ppm) | z CTOI | Target–EB (px) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 14 | 1229/1241 | 9 | 887.17 ± 242.28 | 3.662 | 9 | 865.28 ± 237.20 | 3.648 | 3.9596 |
| 15 | 1181/1197 | 7 | 363.69 ± 219.55 | 1.657 | 8 | 408.99 ± 154.14 | 2.653 | 3.9710 |
| 40 | 3912/3927 | 9 | -179.61 ± 353.26 | -0.508 | 9 | -84.70 ± 150.04 | -0.565 | 4.1297 |
| 41 | 3652/3668 | 5 | 110807.91 ± 575703.73 | 0.192 | 4 | 276282.48 ± 679795.55 | 0.406 | 3.9429 |
| 54 | 3565/3605 | 8 | 16.11 ± 77.19 | 0.209 | 8 | 33.94 ± 322.47 | 0.105 | 3.9162 |
| 55 | 3765/3781 | 8 | 238.55 ± 124.06 | 1.923 | 8 | -235.68 ± 106.45 | -2.214 | 3.9134 |
| 74 | 11264/11284 | 8 | -78.71 ± 512.47 | -0.154 | 9 | -250.80 ± 149.21 | -1.681 | 3.9690 |
| 75 | 11674/11694 | 9 | 372.84 ± 535.01 | 0.697 | 9 | -306.56 ± 228.31 | -1.343 | 4.0094 |
| 81 | 11029/11076 | 8 | 286.65 ± 97.65 | 2.935 | 8 | 257.37 ± 97.63 | 2.636 | 3.9034 |
| 82 | 10868/10898 | 0 | — | — | 0 | — | — | 3.9121 |

Downloaded sectors: 10/10; sectors with usable EB-event measurements: 9/10. Sector 82 has 10,868 accepted cadences but zero valid normalized events: its aperture values and border-subtracted continuum are zero in buffer rows. It stays in the sector and injection denominators. Counts and reasons are retained in outputs/photometry_s0082.json. Individual event inputs, continuum values and formal variances are retained for independent restacking.

**Science-array coverage limitation.** The header-derived target coordinate in sector 41 is FFI column 2101.789, beyond the illuminated column range 45–2092; its median raw four-pixel flux is only 0.0297 e−/s and median border-subtracted flux 0.00291 e−/s. Its enormous fractional depth/error reflects an effectively blank aperture. Sector 82 places the target at row 2051.058; the FFI header marks science rows 1–2048 and buffer rows 2049–2058. Its aperture flux is zero. TESScut metadata therefore includes geometric CCD/collateral coverage that is not useful illuminated stellar coverage. These sectors were retained without post-flux exclusion. The edge diagnosis is descriptive (outputs/edge_diagnostics.json) and changes no bar or trial count.

**Ephemeris drift.**

Both ephemerides stay fixed: EB P=3.0305267 d, BJD0=2454967.124438; CTOI P=3.02553 d, BJD0=2458682.5504. The period difference is −0.0049967 d per cycle (−7.195248 min, −0.1649% relative to EB). At cycle 1226, the CTOI epoch minus the EB primary is about +0.0002278 d. The cycle-aligned drift below continues from this pairing; the nearest-primary offset wraps modulo P_EB. Negative means the CTOI prediction arrives earlier. These are deterministic central-ephemeris predictions, without re-fitting or uncertainty broadening.

| Sector | Midpoint BJD | CTOI cycle from epoch | Cycle-aligned CTOI−EB (d) | Nearest EB offset (d) |
|---:|---:|---:|---:|---:|
| 14 | 2458696.28400 | 5 | -0.024756 | -0.024756 |
| 15 | 2458723.84650 | 14 | -0.069726 | -0.069726 |
| 40 | 2459404.29062 | 239 | -1.193983 | -1.193984 |
| 41 | 2459433.84267 | 248 | -1.238954 | -1.238954 |
| 54 | 2459783.60935 | 364 | -1.818571 | +1.211956 |
| 55 | 2459811.13710 | 373 | -1.863541 | +1.166985 |
| 74 | 2460326.33775 | 543 | -2.712980 | +0.317546 |
| 75 | 2460353.52291 | 552 | -2.757951 | +0.272576 |
| 81 | 2460520.14918 | 607 | -3.032769 | -0.002242 |
| 82 | 2460546.42344 | 616 | -3.077739 | -0.047213 |

**Gate and location.**

Trials use RNG seed 301102, random phases drawn once before recovery, real cadences outside the frozen EB and CTOI primary/secondary windows, and the same blind single-PRF fitter used for the proposed real location test. Two trials per star per sector yield 20 per star. Injection amplitude reproduces the observed aperture loss separately at each position. Nonpositive sector loss and absent event coverage count as failures; no phase redraws or substitutions occurred. Same-model injections do not validate actual PRF mismatch.

| Sector | Target recovered/planned | EB recovered/planned | Target positions returned | EB positions returned |
|---:|---:|---:|---:|---:|
| 14 | 0/2 | 2/2 | 1/2 | 2/2 |
| 15 | 0/2 | 2/2 | 1/2 | 2/2 |
| 40 | 0/2 | 0/2 | 0/2 | 0/2 |
| 41 | 0/2 | 2/2 | 2/2 | 2/2 |
| 54 | 0/2 | 2/2 | 2/2 | 2/2 |
| 55 | 0/2 | 2/2 | 2/2 | 2/2 |
| 74 | 0/2 | 0/2 | 0/2 | 0/2 |
| 75 | 0/2 | 2/2 | 2/2 | 2/2 |
| 81 | 0/2 | 2/2 | 2/2 | 2/2 |
| 82 | 0/2 | 0/2 | 0/2 | 0/2 |

**Material gate deviation:** the frozen exclusion window is too narrow to establish the requested out-of-eclipse-only injection backgrounds. Villanova gives full eclipse widths 0.1530 and 0.1384 in phase (0.46367 and 0.41942 d), versus the registered 0.252 d exclusions. Eclipse wings/variability may remain. The original counts above are retained as a failed registered experiment; a compliant clean-OOE validation is INCOMPLETE. No post-flux mask change or reassessment was made. Six trials per position already fail for nonpositive/unavailable observed signal, placing their maximum possible recall at 14/20 under the remaining frozen choices. The original failed gate and unscored Bar B remain unchanged.

Trials without positions: {"ValueError('no valid events')": 2, "ValueError('nonpositive signal/throughput')": 12}. Full 40-trial phases, amplitudes, coordinates, covariance matrices, distances, failures and outcomes are in outputs/gate.json and outputs/injection_design.json. For the failed target trials that did return positions, the saved fits often land several pixels away; these are injection recoveries, not real-source localizations. The real EB difference images remain retained but uninspected/unfitted. Therefore there is no measured mean source position, no measured distance to either star, and no Bar B significance to report.

**Descriptive eclipse depths and dilution.**

| TESS quantity | Mean-window depth ± error (ppm) | z | Events |
|---|---:|---:|---:|
| Primary | 171.50 ± 102.11 | 1.680 | 71 |
| Phase 0.5 | 83.41 ± 70.35 | 1.186 | 73 |
| Odd primary cycles | 178.29 ± 116.01 | 1.537 | 37 |
| Even primary cycles | 163.22 ± 176.52 | 0.925 | 34 |

Kepler PDCSAP was read only for KIC 10014830: 17/17 quarters, 56685/69302 accepted/total cadences. The same fixed-window estimator gives:

| Kepler EB quantity | Depth ± error | z | Events |
|---|---:|---:|---:|
| Primary | 42.712 ± 0.346% | 123.527 | 419 |
| Phase 0.5 | 18.892 ± 0.075% | 253.255 | 417 |
| Odd primary cycles | 42.565 ± 0.515% | 82.708 | 206 |
| Even primary cycles | 42.822 ± 0.469% | 91.225 | 213 |

Quarter primary-window depths vary substantially (38.7%–52.7%; all-event chi²/dof=822.1). The combined formal weighting and scatter inflation do not turn heterogeneous PDCSAP into a unique physical eclipse depth. The Villanova catalogue central-depth values and our mean-window depths need not agree. No on-target KIC 10014894 Kepler flux was read.

TIC Tmag predicts EB/target total flux ratio 0.08933; Kp 14.827 versus target 11.67 gives 0.05460. The following model uses the frozen border-subtracted aperture and all finite-Tmag catalogue stars. Kp values are an approximate ratio rescaling, without reconstructing all contaminant colors. Implied undiluted depths are descriptive ratios of noisy aperture depth to predicted EB fraction, not EB-depth measurements or bar criteria.

| Sector | EB aperture flux fraction (Tmag, %) | Kp approximation (%) | Implied undiluted depth (Tmag, %) | Kp approximation (%) |
|---:|---:|---:|---:|---:|
| 14 | 0.16238 | 0.09925 | 54.63 | 89.39 |
| 15 | 0.10394 | 0.06353 | 34.99 | 57.25 |
| 40 | 0.03614 | 0.02209 | -49.71 | -81.32 |
| 41 | 0.08986 | 0.05492 | 12331.37 | 20174.90 |
| 54 | 0.04349 | 0.02658 | 3.70 | 6.06 |
| 55 | 0.08006 | 0.04893 | 29.80 | 48.75 |
| 74 | 0.02087 | 0.01276 | -37.71 | -61.69 |
| 75 | 0.13137 | 0.08029 | 28.38 | 46.43 |
| 81 | 0.04739 | 0.02897 | 60.48 | 98.96 |
| 82 | 0.07690 | 0.04700 | — | — |

Negative or >100% ratios flag noise/model inconsistency, rather than physical negative or super-total eclipses. PRF wing uncertainties, differential TESS/Kepler colors, background-median response, catalogue incompleteness and Kepler PDC crowding corrections limit this comparison. No interpretation overrides the failed bars.

**Prior art and credit.**

PRIOR_ART.md was written before analysis code. The post-hoc r29 match was the starting lead; r29’s frozen 0/12 result remains unchanged. r28 left all twelve periodic CTOIs OPEN; its target result was not a mechanical Kepler refutation. A permitted live Villanova overview check confirms the EB ephemeris and morphology, and flags FB (flat-bottomed eclipses), TM (ETV evidence of a multiple), DV (eclipse-depth variations); the catalogue documentation/papers are retained in sources/. ExoFOP evidence is inherited from r29 REPORT.md (2026-11-02 PC, no TFOP disposition/notes/time series); it was not refreshed. Live TESS EB catalogue (Prsa et al. 2022), GitHub, arXiv and community searches are UNCHECKED because their sites fall outside this job’s explicit network allowlist. This limitation precludes a novelty assessment; none is claimed.

Credit Montalto (2022) and DIAmante for the CTOI and duration in DAYS; the Villanova Kepler EB Catalog team, Prsa et al. (2011), Slawson et al. (2011), Kirk et al. (2016), for the known EB and ephemeris; Prsa et al. (2022) for the TESS EB catalogue context (not live checked here); NASA TESS, SPOC PRF commissioning/model team, TESScut/MAST/STScI, Kepler and the public archive contributors; TIC/Gaia, NumPy and SciPy. See [DIAmante at MAST](https://archive.stsci.edu/hlsp/diamante), [Villanova EB](https://keplerebs.villanova.edu/overview/?k=10014830), [official PRF documentation](https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/start_s0004/00README.txt).

**Deviations, verification and reproduction.**

DEVIATIONS.md records the material narrow-OOE-mask limitation with the original gate outcome, pre-flux preparation failures, the official TESScut delivery-host handling, scheduler pause/resumption, inherited nice correction, offline environment and the unchanged Kepler-script dependency-order rerun. None changed a bar threshold or replaced a sector/trial. The custom FITS/TAN-SIP reader was validated with direct TIME bytes, celestial/full-FFI round trips (<0.002-pixel local approximation differences), and pre-flux synthetic PRF/box checks. It is not a full Astropy replacement. Scientific code hashes at first flux match the final seal. A report-only postprocessor corrects the initially mistaken zero-cadence description for sector 82 and puts the gate-validity limitation up front; original measurement and verifier files remain unchanged.

Independent verify.py certifies saved arithmetic, not OOE purity or full methodological compliance. It imports no measurement code. It recomputes per-sector and all-sector stacks and Bar A, all 40 gate decisions with exact rational half-pixel comparisons, the conditional Bar B decision and final 0/1 attribution. It also checks frozen protocol/aperture hashes and 160 sealed files. Decimal restacking uses 70-digit arithmetic. This verifies saved numerical decisions, not independent astrophysical truth or photometric extraction. outputs/verification.json is the certificate.

Run saved-output verification in this directory:

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 nice -n 0 .venv/bin/python verify.py
```

The launcher already has nice 10; `nice -n 0` retains actual nice 10. Dependencies: CPython 3.12.13, NumPy 2.5.4, SciPy 1.18.1, installed offline into the job-local venv; commands are in OFFLINE_ENV.md. No Astropy/network package installation was used. Historical extraction sequence (do not overwrite the original outputs; independent replication needs a fresh authorized directory):

```sh
# First write PRIOR_ART.md; fetch metadata; commit PROTOCOL.md and FROZEN STATUS.md.
.venv/bin/python metadata.py
.venv/bin/python more_metadata.py
.venv/bin/python prf_metadata.py
.venv/bin/python download.py
.venv/bin/python fetch_kepler.py
.venv/bin/python ffi_headers.py
.venv/bin/python prepare.py
.venv/bin/python validate_local.py
.venv/bin/python measure.py photometry
.venv/bin/python measure.py injections
.venv/bin/python measure.py location
.venv/bin/python kepler_compare.py
.venv/bin/python prior_finish.py
.venv/bin/python score.py
.venv/bin/python verify.py
.venv/bin/python make_report.py
.venv/bin/python finalize_report.py
```

All scientific commands use the environment/thread settings and actual nice10 shown above; no more than two CPU cores are used, and every network request is serialized at >=5-second starts with backoff. All measurement stages refuse to overwrite their final outputs. Full file/data hashes, GET URLs, version/quality/position/aperture metadata and logs are retained. Retained files stay below the 2 GB cap; no worker remains running at finalization. The only out-of-directory write was the explicitly requested scheduler-block alert.

A future test would need a correctly masked clean-OOE gate and a separately frozen method with demonstrable target-position localization sensitivity and adequate control of background/stellar variability. The present result remains no attribution, with Bar B unscored. No stronger scientific or public claim follows.

Report rendered 2026-11-02T12:08:12.098816+00:00.
