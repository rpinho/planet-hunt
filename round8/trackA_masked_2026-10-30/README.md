# Track A masked recalculation

Track A remeasures the 17 original round 8B hold-out survivors at their frozen
round-8 ephemerides and durations. Mask A removes the first 1.0 and last 0.5 days
of QLP orbits in sectors 94–104 before detrending. The bar is `z_A >= 2.5`;
the frozen rule also requires unmasked mask-A edge SNR-squared share <=50%.
`Kept` below is the saved `bar_stays` decision, checked against both conditions.
Period is `signal.P` in days; `z_A` is `cases.A.z_null`, without rounding.

| TIC | Period (days) | z_A | Kept |
|---:|---:|---:|:---:|
| 47393760 | 1.8317291077466848 | 3.0877176731454288 | yes |
| 101653309 | 0.7034399053390048 | 4.972589815885718 | yes |
| 189827896 | 0.997345728044974 | 3.3902338202623 | yes |
| 201752790 | 2.512650570870438 | 3.8044189136900273 | yes |
| 262843927 | 24.108771213642704 | 3.4247204801545372 | yes |
| 317300777 | 2.8210766304506065 | 2.8957891520680397 | yes |
| 348760393 | 12.388635480281128 | 2.9591833799614 | yes |
| 410418450 | 6.8202239655237875 | 1.9640825055924842 | no |
| 423840582 | 16.16750076129867 | 4.903815322610392 | yes |
| 447558713 | 5.557604899603225 | 8.088628559188844 | yes |
| 254643343 | 0.6997152308529201 | 1.6344654161932741 | no |
| 261117566 | 0.6231819495303007 | 2.5859539054300655 | yes |
| 267546678 | 2.2325561399863254 | 1.7748906202586365 | no |
| 363119909 | 0.5537644772326115 | 10.309001875966066 | yes |
| 421279377 | 2.057467986290972 | 10.028997634763499 | yes |
| 421279377 | 7.226765599613811 | 2.7692698645265215 | yes |
| 198787224 | 0.6467256752614395 | 18.80925083878296 | yes |

14/17 stay; 17/17 measured.

TIC 254643343 has z_A = 1.6344654161932741 < 2.5, the masked-recalculation
failure cited in AAS82326 v7. [TRACK_A.md](TRACK_A.md) contains the verbatim
Track A section of the source report, preceded by its provenance header.

Recompute this table from the JSONs with the standalone standard-library script:

```sh
python3 table_trackA.py
```

All 51 arrays referenced by the 17 result JSONs are included, including the
reused control arrays; their recorded SHA256 values are unchanged. The bundle
is below 25 MB. No logs, runtime contents, vendor directory, network module,
or raw input curves are included.

The historical `track_A.py` and its local helpers are provided for inspection.
A full measurement rerun additionally needs the original acquisition plan,
control JSONs, raw products, and external hold-out method/data archive, plus
NumPy, SciPy, pandas, Astropy, and wotan. Those inputs are outside this bundle;
`python3 track_A.py` is therefore not a standalone reproduction command.
Relative archive roots in `legacy_photometry.py` must be supplied by a rerunner.
Only the table command above is self-contained.

Privacy edits replace local roots with relative paths and redact source
usernames in cohort attribution. The external pipeline directory is called
`upstream_pipeline` here. `disk_guard.py` retains the original disk-space helper;
its two imports replace the network-module imports. Scientific numbers and
measurement logic are unchanged. [DEVIATIONS.md](DEVIATIONS.md) is the scrubbed
historical record; [REPORT.md](REPORT.md) records packaging checks.
