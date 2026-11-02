# TIC 275488151.01 versus leaked light from KIC 10014830

**POWER_FAIL — saved-output decisions independently VERIFIED. No attribution (0/1 CTOI). Bar 1 UNSCORED; Bar 2 UNSCORED. No real EB- or CTOI-ephemeris difference image was formed.**

The registered fitter has ample power to recover the injected EB component, but insufficient power to exclude a target signal at the required level. This outcome does not establish or reject the leak hypothesis, classify a planet, or compute an FPP.

Disclosure: **r30's aperture depths and noise levels were already seen; its per-pixel difference images were not.** No r30 difference-image file was opened in this job. The supplied post-hoc leak interpretation motivated this test; it is not treated as a new result.

## Registration and power gate

[PROTOCOL.md](PROTOCOL.md) was committed locally as `99ffe523506fb6b83dd7f2ed17979e3b8cb8fe0f` at **2026-11-02 12:19:24 UTC**, before any r31 pixel difference image or field FLUX read. STATUS.md then began FROZEN. The entire seed-311102 phase design was saved before OOT flux access at **12:22:05.743 UTC**. The completed gate was saved at **12:22:06.390 UTC**. Implementation/protocol/input hashes, actual nice 10 and all numerical thread limits are in [outputs/before_oot_flux.json](outputs/before_oot_flux.json).

Exactly 200 matched-background trials per arm pooled all eight usable sectors. All **400/400** fits were valid; each sector contributed to **200/200** phases in both arms. Every trial used only original cadences outside both ephemerides' broad primary/secondary exclusions. All predicted injected centers were also away from both sets of excluded windows. Phase rejection used times and ephemerides only, before reading FLUX. No phases were replaced for gaps, recovered signal, or fit quality.

| Registered requirement | Result | Required | Decision |
|---|---:|---:|---|
| EB injection: positive EB component at ≥5σ | **200/200 (100%)** | ≥180/200 (90%) | PASS |
| EB injection: target 95% upper limit <279 ppm | **0/200 (0%)** | ≥180/200 (90%) | FAIL |
| Target injection: false EB detection at ≥5σ | **3/200 (1.5%)** | ≤10/200 (5%) | PASS |
| Target injection: positive target detection at ≥3σ | **0/200 (0%)** | ≥160/200 (80%) | FAIL |

The joint EB recovery-and-target-limit condition holds in **0/200** trials. Failure of either arm requires stopping; both arms failed a required target-sensitivity condition. **POWER_FAIL was written before any real difference image or fit.**

### Injection estimates, not measurements at the real ephemerides

Here r=b/b0; a is target fractional depth, reported in ppm. Quantiles summarize the 200 registered recoveries in each arm; they are not real-signal estimates or posterior intervals.

| Quantity | EB-only arm: median [min, max] | Target-only arm: median [min, max] |
|---|---:|---:|
| Injected a (ppm) | 0 | 558 |
| Injected r | 1 | 0 |
| Recovered a (ppm) | 206.32 [−1462.51, 1091.56] | 764.32 [−904.51, 1649.56] |
| Formal σ(a), ppm | 638.09 [607.01, 733.35] | 638.09 [607.01, 733.35] |
| Target 95% upper limit, ppm | 1263.02 [998.45, 2185.39] | 1821.02 [998.45, 2743.39] |
| Recovered r | 0.99632 [0.91768, 1.09492] | −0.00368 [−0.08232, 0.09492] |
| Formal σ(r) | 0.01712 [0.01632, 0.01941] | 0.01712 [0.01632, 0.01941] |
| EB significance r/σ(r) | 57.79 [49.17, 65.41] | −0.21 [−4.81, 5.44] |
| Target significance a/σ(a) | 0.325 [−2.362, 1.641] | 1.176 [−1.448, 2.481] |

Under the frozen one-sided rule U95=max(0,a)+1.6448536269514722 σ(a), the **smallest registered σ(a), 607.01 ppm, already implies U95≥998.45 ppm**. Therefore these trials cannot meet 279 ppm even at a≤0. For scale, U95<279 ppm at a=0 requires σ(a)<169.62 ppm. This diagnosis uses only the registered OOT trials and does not score real eclipse pixels. Residual chi-square covariance inflation was 1 in 400/400 fits; no errors were deflated.

Trials share the same finite OOT data and the two arms share phase vectors. Thus the 200 recoveries per arm are an empirical phase ensemble, not 200 independent observational datasets; binomial confidence intervals are not asserted.

## Bars, verbatim, and their outcome

- Bar 1, leak. On the real data, b > 0 at >= 5 sigma and b/b0 within 3 sigma of 1.
- Bar 2, target. a's 95% upper limit is below 279 ppm, and a is below 3 sigma.
- Attribution only if the power gate passes and Bars 1 and 2 pass. Otherwise report the outcome per bar. A pass means the CTOI signal is attributed to KIC 10014830 (a nearby eclipsing binary); it is not a statement about any planet and no FPP is computed.

**Bar 1: UNSCORED. Bar 2: UNSCORED. Attribution: no.** The 200/200 successful EB injections are not a real Bar 1 detection.

The requested real fit at the CTOI ephemeris, EB/CTOI fits restricted to sectors 14–15, real per-sector a and b, and real residual images were **not produced**, as required by the failed gate's stop rule. These are conditional withheld outputs, not zero-valued measurements. No real-image visualization or real-location inspection was performed.

## Sector geometry, cadences and fixed flux predictions

All **10/10** inherited TESScut files match r30's saved SHA256 hashes. No cutout was downloaded. Geometry uses inherited proper-motion-propagated TIC/Gaia positions and parent FFI TAN-SIP headers, with 1-based illuminated detector coordinates. Require each star's distance from the science-area boundaries (columns 45–2092, rows 1–2048) to exceed 2 px. **8/10 sectors retained; 2/10 dropped before flux.** Every retained cutout has **225/225 illuminated pixels**, and all pixels enter the two-component fit.

- Sector 41 dropped: target (2101.789,1062.755), EB (2099.760,1066.136), both outside science columns.
- Sector 82 dropped: target (962.575,2051.058), outside science rows; EB (963.656,2047.298), only 0.702 px from the upper science edge.

In the eight retained sectors, quality/finite-time selection accepts **47,619/47,805** cadences; the broad mask retains **16,915/47,619** for the gate. The table's minimum edge margin is across both stars. Ftarget and b0 are OOT-calibrated total stellar flux and predicted total EB eclipse loss in e−/s, not aperture fluxes or observed eclipse depths. Valid phases count once per sector and are shared by the two arms.

| Sector | Minimum edge margin (px) | Accepted/total cadences | OOT/accepted | Ftarget (e−/s) | b0 (e−/s) | Events per trial (range) | Valid phases | Real a | Real b/b0 |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 14 | 280.437 | 1229/1241 | 578/1229 | 5238.36 | 199.87 | 7–8 | 200/200 | UNSCORED | UNSCORED |
| 15 | 675.803 | 1181/1197 | 532/1181 | 5246.13 | 200.17 | 7–8 | 200/200 | UNSCORED | UNSCORED |
| 40 | 357.476 | 3912/3927 | 1146/3912 | 5190.98 | 198.06 | 6–8 | 200/200 | UNSCORED | UNSCORED |
| 54 | 189.588 | 3565/3605 | 1078/3565 | 5505.87 | 210.08 | 6–7 | 200/200 | UNSCORED | UNSCORED |
| 55 | 179.075 | 3765/3781 | 1025/3765 | 5372.15 | 204.98 | 6–8 | 200/200 | UNSCORED | UNSCORED |
| 74 | 819.473 | 11264/11284 | 3312/11264 | 5266.12 | 200.93 | 8–9 | 200/200 | UNSCORED | UNSCORED |
| 75 | 653.704 | 11674/11694 | 3921/11674 | 5313.04 | 202.72 | 8–9 | 200/200 | UNSCORED | UNSCORED |
| 81 | 25.527 | 11029/11076 | 5323/11029 | 5428.80 | 207.14 | 7–9 | 200/200 | UNSCORED | UNSCORED |

Exact detector positions, all candidate sectors, 52 inherited input hashes, and the fixed Kepler/Tmag prediction are in [outputs/metadata.json](outputs/metadata.json). Event counts, per-sector sufficient statistics and invalid-trial fields are in `outputs/gate_sector_sNNNN.json`; there were no invalid trials or nonfinite-pixel omissions.

## Frozen method and its limits

The EB ephemeris is P=3.0305267 d, T0=2454967.124438 BJD(TDB); CTOI P=3.02553 d, T0=2458682.5504. For **both** ephemerides, baseline/gate exclusions have primary half-width 0.1530/2+0.05=0.1265 phase and secondary half-width 0.1384/2+0.05=0.1192, centered at Villanova's measured secondary phase 0.5007. Quality must be zero; momentum-dump bit 32 is expanded by one cadence. In the gate, FLUX is read from only the retained OOT rows, rather than loading real eclipse rows and masking them afterward.

The primary measurement/injection box is ±0.126 d, matching the already-seen Kepler depth estimator. Local OOT baseline samples lie within ±P/2 and outside an additional ±(0.126+0.05P) guard. Each event requires ≥3 in points and ≥3 baseline points on each side. Every pixel's local temporal constant is subtracted. Event weighting and overlap covariance are retained in the final cadence contrast coefficients; per-pixel noise is the original OOT sample variance.

The two fixed-position official SPOC PRFs give D=in−out=−a Ftarget PRFtarget−r b0 PRFEB+gamma_sector. Positive a and r represent eclipse losses. A per-pixel temporal constant cancels in the contrast; fitting a separate free offset to every difference-image pixel would make the stellar components unidentifiable. A spatially uniform difference-background nuisance is profiled per sector. WLS pools every pixel and sector; covariance is the inverse profiled normal matrix, inflated only if reduced residual chi-square exceeds one. No fitted centroid, real-data aperture selection, flux clipping, period search or tuning occurred.

Tmag gives EB/target total flux ratio **0.08933054837332963**, from 13.838 versus 11.2155. OOT WLS of r30's complete finite-Tmag catalogue scene plus uniform sky fixes each sector's target flux scale. b0=0.42712322107786777×0.08933054837332963×Ftarget. The Kepler factor is **42.7123±0.3458%**, from **419 primary events in 17/17 quarters**, **56,685/69,302** accepted/total Kepler cadences, inherited from r30. The verifier independently restacks those saved event measurements; no new Kepler photometric extraction or on-target KIC 10014894 Kepler flux read occurs.

Additive PRF eclipses were injected into the cadence model before the same linear contrast. Distributed application of that linear operator avoids hundreds of cube copies; **16/16 explicit injected-cube equivalence checks** pass, covering both arms in every sector. Header-time synthetic tests recover a known simultaneous two-star eclipse and cancel a constant background in all **8/8** sectors before any OOT flux read.

Limitations remain: injections and recovery share the same PRFs; Tmag/color/calibration and Kepler-depth uncertainties are held fixed; the Kepler mean-window depth is not the catalogue's central depth; the wider real eclipse wings are excluded rather than modeled; the EB has catalogue depth-variation/multiple-system timing flags; constant local baselines do not remove arbitrary target variability or instrumental trends; diagonal pixel/cadence noise does not model all correlations. The empirical gate evaluates those OOT artifacts for this frozen fitter and fails. It is not legitimate to tune a better detrending method after this failure and score the same real images under this registration.

## Prior art and DIAmante aperture

[PRIOR_ART.md](PRIOR_ART.md), [prior/summary.json](prior/summary.json) and [network.jsonl](network.jsonl) retain the exact permitted search scope, responses, times and hashes. Checks were on 2026-11-02; no novelty assessment or completeness claim is made.

- Refreshed [ExoFOP target](https://exofop.ipac.caltech.edu/tess/target.php?id=275488151&json) and [EB TIC](https://exofop.ipac.caltech.edu/tess/target.php?id=275488133&json). Target is still CTOI PC with empty TFOP disposition and the original Montalto (2022) note; EB TIC has no CTOI. ExoFOP confirms both TIC/KIC cross-identifications. Target overview lists two Montalto attachments, a light curve and DV PDF, despite an empty structured time_series array. No explicit blend attribution appears in the returned overview fields; this does not rule out statements in uninspected attachments. Obsolete download_target.php returned migration notices; the documented target.php endpoints supplied the actual refreshed JSON.
- [Villanova KIC 10014830](https://keplerebs.villanova.edu/overview/?k=10014830) confirms the pre-existing EB, ephemeris, widths, secondary separation and flags. The known EB and its catalogue are prior art.
- [Prsa et al. 2022 TESS EB HLSP](https://archive.stsci.edu/hlsp/tess-ebs): neither TIC occurs among **4584** rows of the linked [v1.0 catalogue snapshot](https://archive.stsci.edu/hlsps/tess-ebs/hlsp_tess-ebs_tess_lcf-ffi_s0001-s0026_tess_v1.0_cat.csv). This is a version-limited archival nonmatch, not absence from a current live EB catalogue. The CSV has TIC rather than KIC identifiers; ExoFOP supplies the cross-identifications.
- arXiv all-field, newest-first searches for **TIC 275488151, TIC 275488133, KIC 10014894, KIC 10014830** each return no results. All-date searches include the last 12 months but do not cover every paper's full text or machine-readable catalogue attachment. Exact requests are saved.
- ADS searches for the same four quoted identifiers returned **0-byte bodies**. ADS results remain **UNCHECKED**, not zero matches. No alternate forbidden API or site was used. RULES generic GitHub/recent-code, AI-lab-repository, and community searches were outside the explicit allowlist and were not performed.

The live [DIAmante HLSP](https://archive.stsci.edu/hlsp/diamante) identifies the established difference image analysis pipeline and Montalto 2020/2022 references. Its [v2 catalogue](https://archive.stsci.edu/hlsps/diamante/hlsp_diamante_tess_lightcurve_catalog_tess_v2_cat.csv) has a row for TIC 275488151: **aperture=2**, **apstat=4.330**, depth=558 ppm, duration=0.210 **days**, 17 transits, periodLin=3.02553 d, t0Lin=1682.5504 BTJD. The [catalogue column documentation](https://archive.stsci.edu/hlsps/diamante/hlsp_diamante_tess_lightcurve_all_tess_v2_readme-catalog.txt) labels aperture in pixels and apstat as the adjacent-aperture depth difference normalized by the combined errors. The [product documentation](https://archive.stsci.edu/hlsps/diamante/hlsp_diamante_tess_lightcurve_all_tess_v2_readme-products.pdf) describes four concentric photometric apertures centered on the target. These sources document the selected 2-pixel aperture parameter; no exact cadence-dependent pixel footprint or mask is reconstructed here. apstat is pre-existing catalogue context, not an r31 attribution test. DIAmante detrending/photometry is not reproduced by this frozen constant-baseline PRF gate.

## Verification and reproduction

[verify.py](verify.py) imports no measurement code. It reconstructs **400/400** pooled WLS fits from the saved per-sector normal equations with **70-digit Decimal arithmetic** and **exact rational threshold inequalities**, recomputes signed significance/upper-limit tests and all gate counts, and returns both real bars UNSCORED and attribution false. It also independently checks the geometric edge rule, all OOT times and injected-center exclusions, event count minimums, Tmag/b0 calibration arithmetic, explicit injection equivalence, inherited Kepler event stack, input/protocol/implementation hashes and freeze chronology. The final artifact seal is checked on a final verifier run. This certifies saved arithmetic and recorded masks, not an independent raw-pixel extraction or astrophysical truth.

Run the retained independent verifier without changing measurements:

```sh
cd ../planets_r31_1102
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=../planets_r30_1102/.venv/bin/python
ps -axo pid,ppid,ni,etime,command
nice -n 0 "$PY" verify.py
```

The launcher already sets nice 10; nice 0 preserves it. From a shell with nice 0, use `nice -n 10` instead. CPython 3.12.13, NumPy 2.5.4, SciPy 1.18.1 are reused read-only from the r30 job-local venv; no packages or new cutouts were downloaded.

The original scientific sequence was `metadata.py`, local protocol commit, `design.py`, `validate_synthetic.py`, `measure.py gate`, then `verify.py`; `measure.py real` was never executed. Frozen outputs refuse replacement. For exact measurement reproduction, [reproduce.sh](reproduce.sh) creates a fresh subdirectory within this job, copies only code/protocol metadata, records a new local freeze, reruns that sequence on the same read-only r30 inputs and checks the original gate counts. It never invokes the real stage. It is provided for review and **was not run after POWER_FAIL**.

[outputs/gate.json](outputs/gate.json) contains all 400 recoveries; [outputs/design.json](outputs/design.json) contains phases and time-mask hashes; [outputs/verification.json](outputs/verification.json) is the independent certificate. [outputs/summary.json](outputs/summary.json) supplies descriptive quantiles and per-sector counts. Retained new artifacts are under 50 MB; cutouts, original PRFs and Kepler data remain only in the read-only predecessor directory.

## Deviations, operations and next step

There were no post-freeze changes to the protocol, thresholds, phases, model or measurement implementation; no recovery reruns, early stopping or substitutions occurred. The per-pixel temporal-constant interpretation, signed eclipse model, OOT flux calibration, narrow measurement window, one-sided upper-limit convention and extra per-sector spatial background were explicitly frozen, rather than silently inferred after looking at results. The local FITS reader was extended before measurement to read selected OOT rows; the predecessor was untouched.

Operational issues: the first public-source refresh briefly inherited nice 20 after adding nice 10 to an already nice-10 launcher, then was corrected to 10; all scientific design, validation and measurement ran at nice 10. The old ExoFOP endpoint and empty ADS responses are the source limitations recorded above. The optional DIAmante column-documentation follow-up completed after the gate; it affected only the report and did not change the frozen method. No external decision or resource was needed, so no blocked-job alert was required.

Network audit: **21** serialized public GETs, **5.030 s** minimum between request starts, no redirects or retries, only the five allowed source domains; about **3.42 MB** downloaded, all documentation/catalogue/search material. No contacts, credentials, uploads, posts or ExoFOP submissions. All writes stayed in this directory. At most two compute cores were available; actual scientific processing used one process with one numerical thread. STATUS checkpoints exceeded the requested hourly cadence. No worker remains after finalization; the final ps receipt is retained.

A future, separately registered test would need substantially better validated target sensitivity, for example an OOT-only treatment of temporal variability plus an appropriate correlated-noise model, and would need to pass the same power requirements before real-image scoring. This job stops at the required negative power result.

Credit: Marco Montalto and DIAmante (Montalto et al. 2020; Montalto 2022) for the original detection, catalogue and pipeline; the Villanova Kepler EB Catalog team, **Prsa et al. 2011, Slawson et al. 2011, Kirk et al. 2016**, for EB identification, ephemerides and morphology; **Prsa et al. 2022** for TESS EB catalogue context; NASA **TESS**, **TESScut/MAST/STScI**, the **SPOC PRF** model team, **Kepler** and archive contributors; TIC/Gaia; NumPy, SciPy and PyMuPDF (documentation text extraction). r30 supplied the input geometry, normalized PRF implementation, hashed cutouts, local environment and saved Kepler event measurements. No planet claim or FPP is made.
