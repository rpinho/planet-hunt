# planets_r14_1030 frozen protocol

Frozen before any new measurement. Explicit task authorization permits this local protocol commit despite RULES.md's general no-commit rule. Nothing is pushed or uploaded. All writes stay here, except the explicitly requested blocker alert if required. No FPP. No candidate claim. Nothing leaves this machine.

Started 2026-10-31 01:28 UTC (2026-10-30 21:28 ET); approximately three hours, deadline 04:28 UTC. At most two single-thread workers at nice 10; wait through SIGSTOP. Use ../.venv/bin/python with r13's read-only dependency path. Footprint <5 GB; stop if free space <30 GB. Retained certificate/data target <50 MB. Do not modify any r13 file, process or shared lock. Prefer read-only cached curves. Any network GET must use the holdout limiter with <=180,000 B/s outside 23:00–06:00 ET and logged verification; no independent limiter competing with r13. If safely sharing the existing limiter requires modifying an input lock, acquisition is blocked and cached work continues.

## Decision (verbatim)

Decision, which resolves the block in planets_r13_1030/BLOCKED.md: the frozen pre-hold-out ephemeris is the right input and the note value is not a valid comparator. The comparators here are the r13 frozen values, TIC 101653309 4.97 and TIC 262843927 3.42. All r13 Track A bars are otherwise unchanged.

## Track C (verbatim)

For each signal below, find the ephemeris and duration that existed before its hold-out sectors entered any fit, and cite the file and field. Then state, from code or data, which sectors went into the published note ephemeris and whether any hold-out sector was among them.

- TIC 144439305 0.940 d (rnaas_check.py globals; expected to match the note)
- TIC 101653309 0.703 d (r13 frozen result, reuse)
- TIC 262843927 24.11 d (r13 frozen result, reuse)
- TIC 262843927 14.52 d (find the input: planets_audit_1030/evidence/legacy/AAS82326 tic262843927_1452 files and the round 8 worktree; if no input predates the hold-out, say so and do not pick one)
- TOI-4405 3.601 d (planets_holdout_1030/results/old_input_comparison.json, source toi4405_lcgate.json cand)

Use the same sectors, products, estimator, detrending and known-signal masks as planets_holdout_1030 and planets_r13_1030/legacy_photometry.py. Report unmasked, mask A and mask B depth, scatter z, and null z from 2,000 random ephemerides, plus exceedance p. Also report the same statistics at the published note ephemeris and published duration, reproduced within 0.1 z of planets_holdout_1030 (5.43, 4.89, 4.42, 3.59, 4.35). If a note ephemeris does not reproduce, stop that row and report it.

Bar: a hold-out supports the note only at its pre-hold-out ephemeris. A signal stays supported at mask-A null z of at least 3. Below 3, report it plainly as lost hold-out support; this is a statement about the hold-out, not about the planet.

Give a table of published text, corrected number and whether the note's wording ("at an ephemeris fixed beforehand", "the ephemeris fixed from the detection data") is true for that signal.

## Track A (verbatim)

Use planets_r13_1030/PROTOCOL.md Track A and planets_r13_1030/controls.py as they stand, at the frozen round 8 ephemerides and durations.

1. For each of the 17 survivors: same hold-out sectors, mask A and mask B depth and z, null z from 2,000 random ephemerides, edge SNR-squared share in the unmasked data.
2. Bar per signal: it stays if null z under mask A is at least 2.5 and the edge share is at most 50%. Otherwise it comes out of the count. The two controls must reproduce r13's 4.97 and 3.42 within 0.1, or stop.
3. Report k of 17 that stay.
4. Control: the inverted-light-curve control under mask A for all 4,021 signals with data if cached or fetchable inside the budget, otherwise the frozen 600-signal subsample in planets_r13_1030/results/inverted_sample.json. Report the expected count of false survivors with its uncertainty, and whether the frozen subsample was used.

Independent verification: a separate verify.py, importing no measurement code, recomputes every count, z and bar decision from saved outputs.

REPORT.md: each track's bar, numbers with denominators, verification status, reproduction commands, deviations. A null is fine. STATUS.md hourly. If blocked, write alerts/planets_r14_1030-needs-claude.md and continue with whatever isn't blocked.

## Frozen implementation conventions

The inherited r13 protocol is saved verbatim in PROTOCOL_r13.md. Its Track A controls' obsolete 4.89/4.42 comparator is overridden ONLY by the explicit decision above. All other bars and implementation conventions are unchanged: real.his P,t0,dur, sectors his=false, HOST masked biweight detrend (0.5 d for P<6 otherwise 1 d; gap 0.5, trim 0.25), red-noise beta per sector, weighted local box depth, null z=depth/std(null,ddof=1), 2,000 valid nulls with at most 40,000 attempts, null P uniform [0.8,1.2]P excluding within 2%, epoch uniform over trial period, remove real and known transit windows at 0.8 duration plus half cadence. Mask before detrending. Edge share is sum event (depth/error)^2 weighted by event in-transit inverse-variance edge fraction / sum event (depth/error)^2; use mask A for the Track A edge bar.

Seeds remain r13's 20261030+1000*idx+j, controls idx=1/2 and j=0/1/2 for unmasked/A/B. Other survivors idx=10+index in sorted survivor IDs. Track C note seeds use holdout's original signal index and j; pre-holdout uses the same seed as that signal's note for a controlled comparison. Additional unmasked nulls use j=0. Reuse existing identical arrays read-only where possible, with SHA256 provenance.

Inverted fallback is the already frozen 600 IDs; no replacements for missing cached data. Reproduce the original full inverted search and gates after mask A, not a sign flip at the real ephemeris. Report available/attempted denominator, k and Jeffreys equal-tail 95% Beta(k+.5,n-k+.5) interval scaled by 4,021. Record missing and failed rows; incomplete population is not a measured zero. All deviations retain original outcomes. No numerical thresholds may change after this commit.
