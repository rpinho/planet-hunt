# TIC 275488151.01: Kepler fixed-EB-ephemeris measurement

**SCORED — final all-quarter gate PASS; Bar K PASS; saved-output arithmetic independently VERIFIED.**

**Kepler excludes the TESS dip on KIC 10014894 at the ephemeris of KIC 10014830.** The primary depth is **0.008 ± 2.571 ppm**, with a **99.7% one-sided upper limit of 7.072 ppm**, below 291 ppm. With r29's already post-hoc timing match, the known EB **KIC 10014830 is the likely source (a nearby eclipsing binary)**. This is the whole claim; no FPP or planet statement is made.

A material sequencing deviation is disclosed: the inherited receipt omitted Q0. The initial Q1–Q17 real result was already seen when the archive completeness check found Q0. The same frozen 400 trials were completed incrementally with Q0, and the completed gate was independently verified before Q0 real-window flux was read. No original event was remeasured and no threshold, model or phase changed. The final result includes **18/18 available long-cadence products, Q0–Q17**. It must not be described as a fully blind, single-pass all-quarter execution. See the correction section below.

## Prior knowledge and registration

**r28 already searched this target's Kepler flux over P ± 5 sigma_P and all phases and found no depth-compatible event (best 7.9 ± 2.0 ppm); this job is not blind to that. It is a fixed-ephemeris measurement with its own power demonstration.** r28's 74/100 matched-control recall and 417/482 nominal coverage left that earlier result OPEN. Its frozen expected Kepler aperture depth was 872.975 ppm; the requested r32 test uses exactly 873 ppm and a 291 ppm limit.

[PRIOR_ART.md](PRIOR_ART.md) was the first artifact. [PROTOCOL.md](PROTOCOL.md) was committed locally as `4d44f169b675842e113db3adcc807d794527aa9e` at **2026-11-02T07:38:28-05:00**, before target-flux extraction; its SHA256 is `823f33953414104972294d2775c52edcfea7e92372169cb9f76719d9f257b729`. STATUS then began FROZEN. The seed-321102 design was saved at 2026-11-02T12:40:06.992267+00:00; initial gate flux access began at 2026-11-02T12:40:53.583078+00:00. The original protocol remains unchanged. All earlier r28–r31 files and the r30 environment were used read-only.

## Gate and bars

Bars, verbatim:

- Bar K. At the EB primary ephemeris, the stacked Kepler depth on KIC 10014894 has a 99.7% one-sided upper limit below 291 ppm (one third of the 873 ppm expected if the TESS signal were on target).
- If the gate and Bar K pass, report: "Kepler excludes the TESS dip on KIC 10014894 at the ephemeris of KIC 10014830." Together with the post-hoc r29 timing match, the report may say the known EB KIC 10014830 is the likely source (a nearby eclipsing binary). That is the whole claim. No FPP, no planet statement.

| Test | Result | Requirement | Decision |
|---|---:|---:|---|
| 873 ppm, 0.210 d box injections: recovered positive depth at ≥5σ | **200/200 (100%)** | ≥190/200 (95%) | PASS |
| Uninjected allowed phases: U99.7 <291 ppm | **200/200 (100%)** | ≥190/200 (95%) | PASS |
| Bar K: real primary upper limit | **7.072 ppm** | <291 ppm | PASS |

The 200 injection phases and 200 separate null phases were generated with PCG64 seed 321102 before flux access; every phase is ≥0.15 in circular phase from primary 0 and secondary 0.5007. There were no new draws, replacements, failed-phase exclusions or recovery-driven changes. The original EB primary/secondary ±0.126 d cadences were excluded using TIME before gate PDCSAP_FLUX access. The original gate and its Q0 completion both used precisely those 400 phases.

Final injection recovery significance: **198.0–260.1σ** (median 227.3σ). Recovered mean-window depth: **719.328–734.077 ppm**, median **727.788 ppm**. The mean-window depth is intentionally smaller than 873 ppm because the 0.210 d box is measured in a 0.252 d window; no transfer correction was applied. Null upper limits span **6.161–11.662 ppm**, median **8.083 ppm**. Valid events per trial span 366–423 for injections and 369–423 for nulls. All 400/400 trials were valid. These phase ensembles share observational data; no independent-binomial interpretation is asserted.

Final gate: [completed_gate.json](outputs/completed_gate.json). Initial gate and events remain [gate.json](outputs/gate.json), [gate_events.npz](outputs/gate_events.npz); Q0-only additions are [q0_gate_events.npz](outputs/q0_gate_events.npz). Independent pre-real certificates exist for both stages.

## Real fixed-phase measurements

| Phase | Depth ± inflated σ (ppm) | Signed significance | U99.7 (ppm) | Valid / quarter-intersecting predicted events |
|---|---:|---:|---:|---:|
| Primary, 0 | **0.008 ± 2.571** | 0.0032σ | **7.072** | **424/472** |
| Secondary, 0.5007 | 1.946 ± 2.522 | 0.7715σ | 8.875 | **417/473** |

Use **57,334/69,778** eligible/total LLC cadences. Every quality value must be zero; TIME, flux and error must be finite, with flux/error positive. No flux clipping or fitted ephemeris occurs. Per-quarter predicted denominators count centers whose ±0.126 d measurement window intersects that quarter's selected finite-time span. Across the full first-to-last finite archive TIME span, including inter-quarter gaps, the corresponding primary/secondary window-intersecting predictions are **485/486**. Thus 424/485 and 417/486 are the full-span usable fractions; neither is a separate threshold in Bar K. Missing/gap events contribute no invented measurement.

All events require ≥3 in points and ≥3 out points on each side. Event-to-event variance inflation is **2.729167** at primary and **2.532142** at secondary; errors are never deflated. The secondary is descriptive and cannot rescue a failed primary bar. Saved real outputs: [completed_real.json](outputs/completed_real.json), original [real_events.npz](outputs/real_events.npz), and Q0-only [q0_real_events.npz](outputs/q0_real_events.npz).

## Per-quarter photometry and aperture geometry

CROWDSAP and FLFRCSAP are archived LLC header values, not new fitted dilution estimates. Optimal pixels use bit 2 of each TPF APERTURE image. EB edge distance is the minimum distance to the union of closed half-pixel optimal-pixel squares, in Kepler pixels. All **18/18** quarterly masks and **4/4** seasons were checked. The target is inside the optimal aperture in **18/18**; the EB is outside in **18/18**.

| Q | Season | Eligible/total cadences | CROWDSAP | FLFRCSAP | Primary valid/pred. | Primary ppm ±σ | Secondary valid/pred. | Secondary ppm ±σ | Optimal pixels | EB edge distance (px) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 3 | 438/476 | 1.0000 | 0.9663 | 3/3 | -7.24 ± 43.84 | 3/4 | 33.22 ± 21.75 | 12 | 18.135 |
| 1 | 3 | 1420/1639 | 0.9999 | 0.9658 | 9/11 | 16.84 ± 17.41 | 10/11 | -1.95 ± 16.27 | 12 | 18.172 |
| 2 | 0 | 2829/4354 | 1.0000 | 0.9787 | 25/30 | 0.55 ± 12.63 | 26/29 | -8.13 ± 11.13 | 13 | 17.951 |
| 3 | 1 | 3786/4370 | 1.0000 | 0.9762 | 27/30 | 4.70 ± 9.84 | 26/29 | 6.14 ± 8.32 | 14 | 18.557 |
| 4 | 2 | 3342/4397 | 0.9999 | 0.9915 | 26/30 | 7.17 ± 8.98 | 28/30 | -12.75 ± 10.29 | 27 | 17.631 |
| 5 | 3 | 4121/4634 | 1.0000 | 0.9656 | 28/31 | -3.59 ± 11.32 | 29/31 | 4.91 ± 9.84 | 12 | 18.141 |
| 6 | 0 | 3916/4398 | 1.0000 | 0.9755 | 29/30 | -0.81 ± 9.55 | 28/29 | 1.00 ± 7.04 | 12 | 17.943 |
| 7 | 1 | 3757/4375 | 1.0000 | 0.9829 | 29/30 | -4.00 ± 8.05 | 25/29 | 16.64 ± 10.48 | 16 | 17.979 |
| 8 | 2 | 2693/3279 | 1.0000 | 0.9795 | 20/22 | 0.75 ± 12.41 | 21/23 | -16.39 ± 13.78 | 20 | 17.625 |
| 9 | 3 | 4231/4768 | 1.0000 | 0.9875 | 32/32 | 18.44 ± 9.10 | 29/33 | -13.18 ± 9.34 | 20 | 16.833 |
| 10 | 0 | 3746/4573 | 1.0000 | 0.9884 | 29/31 | -1.02 ± 8.73 | 29/31 | 4.27 ± 7.71 | 18 | 17.424 |
| 11 | 1 | 3645/4754 | 1.0000 | 0.9922 | 28/32 | 10.85 ± 8.74 | 26/32 | 10.14 ± 10.99 | 22 | 17.958 |
| 12 | 2 | 3053/4044 | 0.9999 | 0.9885 | 21/27 | -8.88 ± 12.66 | 21/28 | -5.46 ± 9.80 | 24 | 17.680 |
| 13 | 3 | 3853/4421 | 1.0000 | 0.9805 | 28/30 | -11.81 ± 9.89 | 27/30 | 8.26 ± 11.66 | 16 | 17.227 |
| 14 | 0 | 3923/4757 | 1.0000 | 0.9782 | 29/32 | -2.87 ± 9.09 | 27/32 | 12.50 ± 8.86 | 13 | 17.932 |
| 15 | 1 | 4128/4780 | 1.0000 | 0.9799 | 30/33 | 3.89 ± 10.07 | 28/32 | 19.92 ± 10.37 | 15 | 18.529 |
| 16 | 2 | 3238/4203 | 0.9999 | 0.9783 | 23/28 | -33.60 ± 12.15 | 24/29 | -11.85 ± 10.73 | 19 | 17.737 |
| 17 | 3 | 1215/1556 | 1.0000 | 0.9568 | 8/10 | 15.13 ± 22.05 | 10/11 | 2.83 ± 17.18 | 11 | 18.145 |


[quarter_table.csv](quarter_table.csv) retains full numerical values. Public TIC coordinates and proper motions were propagated from J2015.5 to each quarter midpoint and transformed with that TPF APERTURE TAN WCS. The EB lies about 20.5 pixels from the target; its nearest optimal-aperture edge is **16.833–18.557 pixels** away. This is an extrapolation of the small-stamp WCS, with a very large separation margin; it is not an EB centroid measured in this TPF. Mask membership does not imply zero PRF-wing or scattered-light contribution. Headers, masks and source hashes are retained in [apertures.json](outputs/apertures.json) and `data/tpf/*.aperture.json`; no TPF flux was analyzed. Download/decompression temporary files were removed, leaving the compact masks and headers.

## What the two Montalto attachments already say

The [public light curve](https://exofop.ipac.caltech.edu/tess/get_file.php?id=887179), `TIC275488151L-mm20221020.txt`, contains **2,438 rows and 39 columns**, over BTJD 1683.367139–1737.387918. Its header defines aperture **radii** 1, 2, 3 and 4 TESS pixels. For each it supplies raw and two normalized light curves, errors and flags, plus sky, eight flux-weighted RA/Dec centroid columns and parent FFI filenames. There are **2,227/2,438** zero-flag, nonzero-centroid rows per aperture; the other rows include zero centroid placeholders. The file carries no narrative on-target/contaminant verdict. These existing TESS centroids were inventoried, not used to form a new r32 TESS attribution test.

The [public DV PDF](https://exofop.ipac.caltech.edu/tess/get_file.php?id=887180), `TIC275488151O-mm20221020.pdf`, has two raster pages, both rendered and visually inspected. Its first page states **“Aperture=2 pix”**, apstat=**4.33**, TIC contamination ratio=**0.5645**, duration=0.210±0.012 d, and 17 transits. Its displayed fit P=3.025346 d and epoch 1682.547196 BTJD differ slightly from the CTOI summary; neither is substituted for the frozen EB ephemeris.

The second page shows the Gaia stellar field, the four aperture circles and aperture-dependent centroid-motion ellipses, with Mahalanobis ranks **(1,1,1,0)**. The [DIAmante 2020 method](https://arxiv.org/abs/2008.09832) defines rank 1 when the target is the nearest Gaia source to that aperture's centroid-motion distribution, and rank 0 when another source is selected. **The public DV therefore already has a non-target 4-pixel centroid association**, while the smaller three apertures select the target. It does not name that other source as KIC 10014830. [Montalto 2022](https://arxiv.org/abs/2210.14559), equation 2, defines apstat and selected values <5: this target's 4.33 passes that published cut. It is not interpreted as a new r32 sigma detection or an FPP. These are meaningful pre-existing aperture/contamination diagnostics, and must accompany any discussion of the result. Review metadata and source hashes: [attachment_review.json](outputs/attachment_review.json).

## Estimator, verification and limits

Use fixed P=3.0305267 d, T0=2454967.124438 BJD(TDB), and secondary T0+0.5007P. TIME is converted with the FITS BJDREFI/BJDREFF. Each quarter is processed separately. The in window is ±0.126 d; baseline intervals are [-0.504,-0.252] and [+0.252,+0.504] d. Ordinary least squares fits y=a+b dt only to the two baseline windows. Fractional depth is (a+b mean(dt_in)−mean(y_in))/a. Event variance follows the local residual variance with the frozen formal PDCSAP-error floor and baseline-fit prediction uncertainty. The inverse-variance stack is inflated by max(1,Q/(N−1)). The one-sided normal limit is U=max(0,D)+2.7477813854449926σ; it is a 99.7% one-sided Gaussian convention, rather than exactly 3σ.

[verify.py](verify.py) imports no measurement code. It reconstructs every event depth and variance from saved sufficient statistics, restacks with **70-digit Decimal arithmetic**, and applies **exact rational squared threshold comparisons to the reconstructed decimal numbers**. It checks all 400 initial and 400 completed trial decisions, both primary bar decisions, all per-quarter stacks, 400 regenerated phases, 18 LLC hashes/archive completeness, unchanged original-event prefixes, two local registration commits and chronology, and all 18 aperture transforms with a separately implemented scalar WCS calculation. It also checks the serialized GET log and final artifact seal. The final combined dataset has **164,445 gate-event records** and **841 real-event records**. This verifies saved arithmetic and provenance, not a second extraction from raw flux or an interval proof of astrophysical noise calibration. Certificate: [verification.json](outputs/verification.json).

The baseline-normalization error uses the frozen first-order fixed-continuum convention inherited from r30. Cadence boxes ignore exposure integration; real transit shape and bandpass evolution are not inferred. A Kepler-era, strictly periodic on-target signal of the frozen scale is the tested hypothesis; time-dependent signals or ephemeris drift outside this fixed test are not excluded. Random-phase injections validate this estimator on this target and these windows, not universal transit completeness. Event scatter inflates diagonal errors but is not a complete correlated-noise model. No period/phase/duration scan, r28 search rerun, new fitted centroid, or additional astrophysical classification is hidden in this test.

## Completeness correction and other deviations

[COMPLETENESS_CORRECTION.md](COMPLETENESS_CORRECTION.md) was committed as `f54502e721166542ece5036b737a055bd956df3e` before any Q0 flux extraction. The initial Q1–Q17 gate had already passed 200/200 in both arms and been independently verified before its real measurements. Its primary was **0.069±2.575 ppm**, U99.7=**7.145 ppm**; its secondary was **1.624±2.534 ppm**. All original outputs and hashes remain unchanged. A later direct LLC prefix listing found Q0, which the inherited r28 receipt omitted; the first implementation had wrongly asserted Q1–Q17 were complete.

Q0 was then added to the existing 400 trials without touching any original event: **1,200 additional gate-event records**. The completed gate again passed 200/200 in both arms and was independently verified and saved before Q0 EB-window flux access. Q0 contributed 3 valid primary and 3 secondary real events. Final primary is **0.008±2.571 ppm**, U99.7=**7.072 ppm**. This is an input-completeness repair, not a second set of random trials or a failed-gate rescue. Nevertheless, the final all-quarter gate completion followed viewing the original 17-quarter real result, so the strict all-quarter gate-before-all-real execution order was not fully met. The numeric exclusion is robust to the correction, but the registration deviation is material and is not concealed.

Other choices—OLS baseline, formal-error floor, nonnegative upper-limit convention, narrow EB-window sealing for gate rows, cadence-center box injection and per-quarter boundaries—were frozen before flux and were not tuned. One archive directory URL returned 403 on three serialized attempts; the public stpubdata listing resolved it. The job was not blocked and no external alert was necessary. No other data/model/threshold changes occurred.

## Reproduction, operations, prior-art scope and credit

To verify saved outputs without network or remeasurement:

```sh
cd ../planets_r32_1102
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ps -axo pid,ppid,ni,etime,command
# This session already inherits nice 10. From an ordinary nice-0 shell, prefix nice -n 10.
../planets_r30_1102/.venv/bin/python verify.py
```

[reproduce.sh](reproduce.sh) creates a fresh directory inside this job, uses read-only links to the retained LLC inputs, makes new local registrations and replays the exact initial-plus-Q0 historical calculation, then compares the final stacks. It performs no network request. It is provided for review and was not run as an extra scoring pass. Its shell syntax was checked. Actual scientific order: design.py → validate_synthetic.py → measure.py gate → verify.py → measure.py real → verify.py → registered completeness correction → complete_q0.py gate → verify.py → complete_q0.py real → verify.py. `measure.py` deliberately preserves its original 17-quarter assertion; use the reproduction wrapper for historical replay rather than rerunning it on all 18 files in this directory.

Only public GETs were made: **56/59** succeeded, **3/59** failed with 403; **69,665,319 bytes** were downloaded, with minimum request-start separation **5.001 s**. Sources were only authorized MAST/Kepler stpubdata, ExoFOP-TESS public pages/attachments, and arXiv. No redirects, credentials, contacts, posts, uploads or ExoFOP submissions occurred. No local data were uploaded. All writes stayed in this job. The launcher supplied nice 10; all numerical code used one process with numerical thread counts fixed to 1, within the 2-core budget. Retained artifacts are below 50 MB; compressed TPF download and decompression temporaries were session-created and removed. Final process inspection is retained in logs/process_final.txt; no analysis worker remains.

Prior-art scope: exact-identifier arXiv searches returned no results, and broader newest-first DIAmante/Kepler-EB pages plus the primary references were checked. Generic GitHub/open-PR, AI-lab and community searches required by generic RULES were outside the user's narrower network authorization and were not made. No novelty or comprehensive literature claim follows. Existing Villanova values were inherited from r29–r31; no forbidden new Villanova request was made. The two Montalto attachments and their diagnostics are prior art. A future review should assess the disclosed sequencing correction and the assumptions of the fixed-ephemeris exclusion before reusing the wording; no additional scientific test is required to finish this registered measurement.

Credit: **Montalto/DIAmante** (Montalto et al. 2020; Montalto 2022) for the original detection, photometry and diagnostics; the **Villanova Kepler EB Catalog team: Prsa et al. 2011, Slawson et al. 2011, Kirk et al. 2016**, for the known EB and ephemeris; **Kepler**, **TESS**, NASA, MAST/STScI and archive contributors. r28 supplied the earlier search result and incomplete source receipt; r29 supplied the post-hoc timing context; r30 supplied the local-linear baseline method, FITS reader lineage and read-only environment; r31 supplied documentary context. NumPy, Python Decimal/Fraction and PyMuPDF supported computation and review.
