This is a mechanical copy, scrub and provenance bundle. NOTE.md has a three-line supplied header followed by the byte-identical version-13 note. PROVENANCE.json records the original note digest, source digests of copied files, oversized JSON exclusions and original files that do not exist. The copied reports retain their scientific limitations and registration deviations. No novelty or worth assessment was performed for this bundle.

Rounds 27 and 28 supply frozen target metadata, dilution and earlier Kepler search context; round 29 supplies the nearby-EB match and its rejected registered period tolerance; round 30 supplies Kepler EB depth and sector geometry; round 31 supplies the pixel-fit gate and magnitude metadata; round 32 supplies the completed target Kepler measurement, aperture checks and attachment inventory; rounds 36 through 38 supply leakage predictions, EB timing and light-curve model comparisons. Rounds 33 through 35 concern other targets and contribute no numbers to this note.

Local path prefixes are rewritten to relative paths and local account names to neutral descriptions. Scientific numerical values are not edited. PROVENANCE.json distinguishes copies changed by this scrub. Historical source hashes are retained as historical records, so a scrubbed file may not match its old source seal. BUNDLE.sha256.json is the integrity seal for the public bytes. Nothing in the bundle is a log, cache, environment, FITS file, binary array or downloaded source dataset. JSON/CSV result files are individually below the decimal 2 MB limit.

The original round-29 job has neither PROTOCOL.md nor verify.py, including in its tracked file list. PLAN.json contains its original registered protocol. The bundled PROTOCOL.md is an explicit pointer to that plan; verify.py is a new, read-only saved-result arithmetic check. The original match.py and posthoc_epoch.py are copied for inspection, require omitted catalog inputs, and write results when run.

The original injection JSON files in rounds 37 and 38 exceed the per-file limit. Each injections_summary.json is an exact key/value projection of its original, dropping only the per-trial events lists; trial amplitudes, normal equations, chi-squared values, validity, counts and calibration summaries are retained unchanged. The source SHA256 and projection are recorded in PROVENANCE.json. The oversized round-36 supplementary event JSON is omitted, with its hash and size recorded. DERIVED.json contains transparent arithmetic on the copied results, recomputed by verify_bundle.py. No new fit or measurement is performed.

The original full verify.py scripts check more than the small result tables: omitted event/cadence certificates, FITS/PRF and DIAmante source bytes, implementation source files, original source seals, input receipts and original registration commits. Those original git histories are not on this branch. The scripts are retained for inspection and for use in a restored original job checkout; path rewriting alone does not supply those dependencies. Rounds 27/28 additionally use their original catalog, receipt and search-code trees; round 30 uses its aperture manifest and implementation/source seals; rounds 31/32 require the source archives and event arrays; rounds 36–38 require their source data, event/cadence arrays, inherited jobs and registration histories. Consult each copied report and the file reads in its verifier for the exact dependency names. The original scientific verifiers were not rerun or claimed to pass in this reduced bundle. Their saved verification JSON files are historical certificates. The portable bundle checks do run on this branch.

To repeat the portable checks from this directory:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 verify_bundle.py
PYTHONDONTWRITEBYTECODE=1 python3 planets_r29_1102/verify.py
cd ..
egress-gate scan tic275488151
```

For full historical audits, use the original job checkout and restore omitted inputs rather than running a partial verifier over the sealed public files: original verifiers may rewrite their saved verification JSON. This bundle does not claim to provide a complete raw-data reproduction environment.
