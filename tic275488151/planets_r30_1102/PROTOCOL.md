# Frozen protocol — planets_r30_1102

Written and locally committed before any field TESS flux/pixels are opened. Metadata only have been read. No bar threshold changes after freeze. No period/epoch/duration search. All deviations retain the original outcome. This local commit is explicitly authorized by the task despite the general RULES no-commit rule.

## Bars (verbatim)

- Bar A, ephemeris. On the target's aperture (frozen before flux: the pixels whose PRF-weighted target flux fraction is highest, up to 4 px), the stacked primary-eclipse depth over all sectors has z >= 7.1 on the EB ephemeris, and z(EB) - z(CTOI) >= 3. Report per-sector z on both ephemerides and the CTOI ephemeris's predicted drift against the EB's.
- Bar B, location. The per-sector in-minus-out difference image at the EB ephemeris, fitted with a single PRF of free position and amplitude, gives an inverse-variance mean position within 0.5 px of KIC 10014830 and at least 3 sigma from TIC 275488151. Positions from the TIC/Gaia catalogue; give the separations in pixels per sector.
- Gate, before Bar B is read: in at least 18 of 20 injections of an eclipse of the observed size into real cutouts (out-of-eclipse cadences only, random phases), the same fitter returns a position within 0.5 px of the injected star, for injections at the target position and, separately, at the EB position. If the gate fails, Bar B is unscored.
- Attribution only if Bar A and Bar B both pass and the gate passes. Otherwise report the outcome per bar. A pass means the TESS signal is attributed to KIC 10014830 for this CTOI; it is not a statement about any planet.

## Inputs, catalogue and sectors

Ephemerides, fixed BJD(TDB): EB P=3.0305267 d, epoch=2454967.124438; CTOI P=3.02553 d, epoch=2458682.5504. Uncertainties are descriptive only (EB period 0.0000052 d; CTOI period 0.00119 d, epoch 0.0122 d). Duration 0.210 d. Every predicted-primary selection has half-width 0.6*0.210=0.126 d. Do not optimize phases.

Sector list from sources/sectors.json (TESScut metadata only):

| Sector | Camera | CCD |
|---:|---:|---:|
|14|2|4|
|15|2|3|
|40|2|4|
|41|2|3|
|54|3|3|
|55|3|4|
|74|4|1|
|75|4|2|
|81|3|2|
|82|3|4|

Retrieve one 15x15 FFI TESScut per sector, centered at target TIC coordinate (293.203758997702,46.9094215256791 deg). No sector replacement. Missing sectors remain in denominators and are explicitly unscored. Header/time metadata do not authorize reading flux until this protocol is committed.

TIC v8/Gaia positions and Tmag from sources/tic_field.json (654 stars within 0.07 deg). Target TIC 275488151 (KIC 10014894), EB TIC 275488133 (KIC 10014830). Propagate catalogue proper motions from Gaia epoch J2015.5 to sector midpoint where finite; otherwise retain catalogue position, flag missing PM. Use Tmag for all relative fluxes, with Kp descriptive crosscheck. Primary coordinate/wcs transformations: FITS TAN WCS, checked round trips. Save catalogue and header-derived pixel positions before reading field FLUX arrays.

## PRF and target aperture (freeze before flux)

Use the official SPOC UPDATED_2.0 PRF FITS (start_s0004, valid sector >=4). Bilinearly interpolate the four CCD grid points bracketing the cutout, with physical WCS/header coordinates accounting for the documented collateral-column convention. Each is a 117x117 surface of 9 samples per physical pixel over 13x13 pixels. Interpolate this surface at the requested pixel-to-star offsets; it is already a pixel response, do not integrate it again. Normalize the sum over its support to unity at each subpixel phase. No PRF shape fit or tuning to the eclipse. Document PRF range limits and boundary handling.

For each sector compute predicted target and all other catalogue stellar PRF contributions across the 15x15 cutout, scaled by 10**(-0.4*(Tmag-Tmag_target)); omit stars missing Tmag and record count. For each pixel define target flux fraction = target contribution / total catalogue stellar contribution. Rank all cutout pixels by that fraction, tie-break by larger target PRF then row then column. Take exactly the highest four finite pixels with positive target contribution. No signal-to-noise/flux-based aperture changes, no manual exclusions. Save all fractions, selected indices, positions, WCS and PRF receipts to aperture_manifest.json and hash it before opening FLUX. This literal fraction-ranking rule may select faint wing pixels; report their throughput and retain the outcome. Background estimation cannot change the aperture.

## Cadences, local photometry and z

Use finite TIME and finite FLUX on all aperture pixels; QUALITY==0, rejecting every flagged cadence including momentum dumps. Also exclude cadences within one exposure of a nonzero-quality momentum-dump bit (TESS bit 32) where QUALITY provides it. No flux clipping, no adaptive detrending. Background in aperture: subtract per-cadence median of outer cutout border (outermost row/column); preserve the same linear operation for PRF throughput and injection amplitudes.

For each ephemeris and each event: in cadences are |t-primary|<=0.126 d; out cadences are 0.252<=|t-primary|<=0.504 d. Fit a local constant plus time slope to out flux, linearly predict the in continuum, measure mean continuum-minus-in loss. Require >=3 in points and >=3 out points on EACH side. Use exact cadence times, no sigma clipping. Event normalization is fitted local continuum at primary; reject nonpositive continuum and retain counts. Point errors use out residual standard deviation, with local linear-fit covariance propagated to the in mean and the point scatter/n_in term. The per-event fractional depth and error are retained. The fixed box/window is a mean-window depth, not physical central eclipse depth.

Combine events per sector and across sectors by inverse formal variance. Error of the combined depth is the formal inverse-weight error multiplied by sqrt(max(1, weighted residual chi2/(Nevents-1))). Thus event-to-event scatter can only inflate errors. z=depth/error, signed. No trials correction: ephemerides are fixed. Each ephemeris uses its own valid event set with all counts reported. Bar A uses all events across all sectors, no sector selection. Save event measurements so verify.py recomputes stack and z independently. Report CTOI-minus-EB nearest-primary timing offsets at each sector midpoint and accumulated cycle-aligned drift since the CTOI epoch. Describe primary odd/even depths and phase-0.5 depth with the same window/out method; these never affect bars.

## Difference imaging and position

Use raw cutout FLUX with the same accepted cadence masks, without the aperture border subtraction. For each valid EB event apply the same out-only constant+slope interpolation per pixel; image is mean(in-minus-predicted-out), negative for an eclipse. Combine event images by inverse aperture-event formal variance, and estimate per-pixel image variance from event-to-event scatter of the weighted mean, floored at propagated local pixel noise. Save images and noise arrays before localization. No examination or free-position fit of the real EB image until the injection gate is saved and scored.

The single-PRF fitter uses a free (x,y), free signed amplitude and constant background. Search starts at every pixel of the interior 13x13, profiles linear amplitude/background, takes the smallest weighted chi2, then bounded nonlinear least squares refines position and amplitude over the whole cutout. Do not initialize solely on either named star. Position covariance comes from the fit Jacobian, multiplied by max(1,reduced chi2), with a 0.1-pixel WCS/PRF floor on each coordinate. Failed, singular or boundary solutions are recorded as failures; no removal after seeing position. Transform each sector's offset from its catalogue EB to a common tangent plane in units of 21 arcsec per pixel using WCS derivatives. Inverse-covariance combine all valid sector fits, retain failed-sector counts; Bar B requires usable positions for all ten sectors. Within-EB distance is Euclidean norm of mean offset; target rejection is the Mahalanobis distance to the target offset, including the saved covariance. Save per-sector detector and sky offsets, covariance, separations from both stars and all decisions.

## Gate, performed before real location fits

40 trials: exactly 20 injections at target and 20 separately at EB. Deterministic RNG seed 301102. Assign trial j to sector j modulo 10 (two trials per sector for each star). Every trial starts with real cutouts restricted to out-of-eclipse cadences: exclude all points within 0.126 d of EB primary AND EB secondary, and CTOI primary AND CTOI secondary. Draw a uniform phase over [0,P_EB), inject the fixed 0.252-d-wide negative box at that phase. Draw all 40 phases before any recovery; do not redraw for gaps or failure. Same event/window construction, same free-position fitter, same noise procedure as the real image. An injection is an additive PRF eclipse, applied to actual cadence pixels before image extraction.

Observed size for a sector is its EB-ephemeris aperture depth times its measured aperture continuum. For each injected position scale total PRF amplitude so that its loss in the frozen, border-subtracted target aperture equals this observed loss. Nonpositive observed depth, negligible/nonpositive effective throughput, no usable injected events, or fit failure count as failures. Save the actual amplitude and effective throughput. This tests apparent loss in this target aperture, rather than assuming identical intrinsic stellar eclipses. Count success only if fitted position is within 0.5 detector pixels of the injected catalogue location. Both sets must have >=18 successes out of 20; all 40 count, no substitutions. Do not read/fix Bar B if gate fails; report UNSCORED, optionally retain difference arrays without real-location fit. Injection uses the same official PRF as recovery, so it does not independently validate PRF shape mismatch; report this limitation.

## Kepler descriptive comparison

Download all available DR25 Kepler long-cadence PDCSAP for ONLY KIC 10014830 via MAST. Do not read KIC 10014894 Kepler flux. QUALITY==0, finite positive PDCSAP, same EB window and local-linear out baseline, per-quarter/event weighting. Record quarter and accepted-cadence denominators, primary/secondary/odd/even mean-window depths. Estimate EB fraction of frozen target aperture with TIC Tmag+PRF; give Kp-based comparison (14.827 vs 11.67) and implied intrinsic eclipse depth from aperture depth / EB fraction. These are descriptive approximations: TESS/Kepler bandpasses, PRF wings, catalogue incompleteness, border subtraction and PDC crowding correction differ. Never use them to tune the gate or bars.

## Operations, verification and retained outputs

2 cores maximum, all compute nice 10, OPENBLAS/OMP/MKL/VECLIB/NUMEXPR threads=1. No parallel downloads; only public HTTPS GETs to MAST/TESScut and keplerebs.villanova.edu; >=5 seconds between request starts, exponential backoff 5/10/20 s, four attempts maximum. RULES.md contains no numeric rate; these conservative limits are adopted. No remote Python/package index calls. Reuse pinned locally available NumPy/SciPy through a job-local venv with recorded interpreter and package paths. A minimal FITS/TAN reader avoids a prohibited package download if Astropy is unavailable. Record dependency versions and reader validation.

Retain under 2 GB (explicit task overrides generic 50 MB rule); stop downloads before exceeding it. All normal writes inside this job directory. User explicitly authorizes the named alert file outside it only if blocked; record that exception. No deletions outside; use reversible Trash inside for pre-existing files, rm only temporary files created this session. No posting, contacts, uploads or ExoFOP operations. Prior-art restrictions are in PRIOR_ART.md, with no novelty claim.

STATUS.md first line FROZEN until final SCORED; update at least hourly. Record start/finish UTC, PID/process checks before launches/relaunches and before finishing. At end write REPORT.md with all bars, gate numerators/denominators, per-sector tables, drift, descriptive depths, deviations, exact commands and credits. Independent verify.py imports no measurement code and uses saved event measurements, injection fits and covariance outputs to recompute all applicable bar decisions; verify freeze/aperture/output hashes. It checks numerical decisions, not independent astrophysical truth. Thresholds stay exactly as above even if they produce a null or an unscored result.
