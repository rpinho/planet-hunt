# Deviations and superseded freezes

2026-11-02T09:02:00.845216+00:00 — Before any Kepler flux was read, FROZEN-v2 superseded freeze commit `7ce2cbe` (retained as an ancestor). The original freeze used literal hours for the ten Montalto CTOIs and failed to apply design/AMENDMENT_1.md. No original measured or scored outcome exists: no flux had been read. Amendment 1 is quoted verbatim in PROTOCOL.md. Primary durations and their errors are now converted from DIAmante days to hours = 24 x CTOI values, after exact numeric matching to the DIAmante v2 source table. The original literal-hours values remain beside them as descriptive sensitivity metadata only; they never decide verdicts. The two missing-duration CTOIs stay OPEN. All other frozen settings and all bars remain unchanged.

## Authorized resume and quarantine — 2026-11-02T09:12:12.980316+00:00

The user confirmed exclusive ownership and authorized FROZEN-v2. Pod process inspection found no remaining measurement workers. Every file under pod results was moved to `../r27/results_superseded_literal`; coverage, receipts, logs and data were likewise moved to sibling directories ending `_superseded_literal`. No files were deleted. Full paths, sizes, SHA256s and original saved statistics are in `audit/quarantine_manifest.json`.

Superseded plan SHA256: `895cfd3450ca25afd583640f1d72d5f650c91f85616c62be60148e94ee67d1c7`. Registered plan SHA256: `e016d9d16ed1aa8dc2a2e7722c2a04b469c27d544af4d7d904ec6c2fae7396e1`. These are superseded outcomes, never scored or reused, including the unchanged single template. No old score exists. The following numbers are preserved original measurements, not verdicts:

| Superseded file | Original selected z | Matching events | Saved injection successes / trials |
|---|---:|---:|---:|
| 120355606.01_measurement.json | 3.7443202687401795 | 0 | — |
| 120963388.01_injections.partial.json | — | — | 10/30 |
| 120963388.01_measurement.json | 4.763011424419191 | 0 | — |
| 122673562.01_injections.partial.json | — | — | 13/60 |
| 122673562.01_measurement.json | 14.120384098230769 | 0 | — |
| 269033446.01_injections.partial.json | — | — | 90/100 |
| 269033446.01_measurement.json | 4.4609678637606 | 0 | — |
| 269033446.01_target.json | 4.4609678637606 | 0 | 90/100 |
| 399793914.01_injections.partial.json | — | — | 41/60 |
| 399793914.01_measurement.json | 5.344487580277919 | 0 | — |

Partial injection JSON lists have no embedded plan hash; their v1 provenance follows from the enclosing run and matching measurement files, rather than an invented embedded hash. Full original coverage values and fitted parameters remain in the manifest and quarantined files.

The uncommitted working PROTOCOL had been replaced by an unfrozen quadratic-block protocol inconsistent with PLAN.json. It was preserved in `audit/resume_unregistered_working_files/`, then restored byte-for-byte from b8f07cb (same registered v2 protocol). No quadratic-block analysis is run. The registered Savitzky–Golay algorithm and finite search shortlist remain unchanged. Cache guards now reject foreign plan hashes; fresh downloads check HTTP Content-Length. These are provenance checks only.

FROZEN-v2 commit f9cd2d9 is timestamped 09:02:00 UTC, before first old-run flux at 09:02:50.886758 UTC. Implementation commit b8f07cb is timestamped 09:03:52 UTC; do not describe that implementation commit as preceding first flux. The plan and protocol did precede flux. All old FITS were moved aside, so no old downloaded FITS is reused without an independent MAST hash/size record. V2 starts with a fresh data directory.

Verifier strengthening before scoring: added independent coverage recomputation from saved accepted cadence times, an output completeness/hash seal, and assertions of fitted-depth/norm arithmetic. No measurement gates or score thresholds changed. The coverage verifier imports no measurement code. Numerical photometry itself is not independently re-extracted. The measurement pipeline remains the registered finite proposal shortlist, including its known search limitation.

Frozen-ID denominator clarification: direct inspection of registered PLAN.json shows 700 frozen control IDs, all distinct, rather than the 750 asserted in some concurrent drafts. CTOI 26493474.01 has 50 frozen but unexecuted IDs; CTOI 26542657.01 has no frozen controls and module_spring=0, in addition to missing duration/contamination. The executable 13 targets each have 50 controls (650 executed slots). Both missing-template rows stay OPEN. Requested denominator remains 750; no post-flux controls are added and PLAN.json remains unchanged.

Pre-score implementation clarification from inspection of b8f07cb: periodic injection epochs are drawn uniformly on the full fixed Kepler span, then periods drawn in the frozen interval. Their folded phases are approximately rather than exactly uniform on one orbit. The injection acceptance mask is propagated from the drawn Kepler epoch (epoch_override), so its uncertainty width can differ from that for a TESS epoch propagated back to Kepler. The registered implementation is retained, with its original measured recall; no rerun, threshold or verdict adjustment is made. This limits interpretation of recall as an exact replication of the target ephemeris window. Report this beside the finite search shortlist and upstream PDC limitation.

Download scheduling only: after the target stage completed, an I/O-only helper prefetches the already frozen periodic control-product S3 keys with eight bounded HTTP threads and persistent connections. Complete files are published atomically with separate temporary paths; no FITS decoding or measurement runs in this helper. Fresh response length and SHA256 are logged. The four measurement workers and scientific pipeline are unchanged. All old FITS remain excluded. The final seal waits for both measurement and prefetch completion. Verification can also read the exact saved output bytes from a compact tar archive, to meet local retention limits; its gates are unchanged.

Superseded result retention independently checked: all 12 files moved out of pod results are byte-identical to entries already in the hash-pulled local results/v1_checkpoint.tar.gz. See audit/superseded_retention_verification.json. They will survive pod disk termination and remain never scored.

Final outcome (one scoring pass, 2026-11-02T09:34:53.899176+00:00): 0 recovered leads, 1 bar-level refutation, 14 OPEN /15. The single is demoted under the registered period-space bar (94314/100000 coverage, 90/100 recall, no matching event; 5/50 controls pass SNR/depth). Several periodic rows have zero injection recall despite high expected SNR; these are registered-pipeline sensitivity failures and stay OPEN. No threshold, detrender, target/control ID or scored outcome was adjusted. Independent coverage verification and all 15 saved-output verdict checks pass. A second local verification is a certificate check, not a second score.

Local retention: exact frozen catalog bytes and auxiliary metadata are compressed, with per-file SHA256 manifests; redundant originals were moved reversibly to the macOS Trash, never irreversibly removed. audit/restore_sources.py restores original paths. All raw numerical outputs are retained in the hash-checked measurement archive; verify.py can read them directly. No bulk FITS were pulled to the Mini. Pod workers exited, with results/shutdown.json confirming no remaining job processes and the pod left up.
