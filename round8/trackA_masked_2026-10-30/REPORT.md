# Packaging verification

Source: `planets_r14_1030`, dated 2026-10-30, local protocol commit `a530985`.
Packaging date: 2026-11-02. This is a local evidence copy for AAS82326 v7.

Verification: all 17 per-signal results, the cohort, and status are included.
All copied JSON numbers are identical to their source values. The complete
Track A section is verbatim after its three-line header. All 51 referenced
arrays match their saved SHA256 values and contain only numerical/boolean data.
The independent table script matches the README, all status rows, all report
rows, and the reported 14/17 kept and 17/17 measured. Mask-A z is checked from
saved depth and null standard deviation, and from the saved null arrays.
Copied Python files parse successfully; the measurement was not rerun.

Reproduce the table from this directory:

```sh
python3 table_trackA.py
```

From the repository root, scan the complete bundle:

```sh
egress-gate scan round8/trackA_masked_2026-10-30
```

External literature was not queried for this mechanical task; local sources
are recorded in PRIOR_ART.md. The README documents privacy edits and the
external inputs required for a full scientific rerun. Next step: readers can
run the table script and inspect the result JSONs; no publication or transfer
was performed by this task.
