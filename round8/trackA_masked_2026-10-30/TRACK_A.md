Source job: planets_r14_1030
Date: 2026-10-30
Local protocol commit: a530985

## Track A: 17 survivors

Frozen bar: mask-A null z ≥ 2.5 **and** unmasked mask-A edge SNR-squared share ≤50%. The cohort remains the original 17 round 8B survivors; this count applies the prescribed hold-out and edge checks. Control gate: 4.9725898159 vs 4.97 and 3.4247204802 vs 3.42; both within 0.1. The printed/refined ephemeris is never substituted for the frozen round-8 input.

|Signal|Hold-out sectors|A depth ± error ppm|A scatter z|A null z|B depth ± error ppm|B scatter z|B null z|Edge share A|Decision|
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
|C4_47393760_1.8317|[45]|451.17 ± 124.82|3.614|3.088|451.17 ± 124.82|3.614|3.106|0.00%|stays|
|E1_101653309_0.7034|[94]|598.57 ± 121.68|4.919|4.973|611.59 ± 127.85|4.784|4.810|1.19%|stays|
|E1_189827896_0.9972|[103]|2625.10 ± 691.25|3.798|3.390|2698.37 ± 742.48|3.634|3.393|17.51%|stays|
|E1_201752790_2.5127|[96, 97]|1145.51 ± 287.68|3.982|3.804|1057.83 ± 297.22|3.559|3.318|33.56%|stays|
|E1_262843927_24.1087|[95, 96, 97]|793.61 ± 216.08|3.673|3.425|792.72 ± 214.09|3.703|3.351|0.00%|stays|
|E1_317300777_2.8211|[94]|587.30 ± 213.89|2.746|2.896|505.36 ± 238.10|2.122|2.240|19.05%|stays|
|E1_348760393_12.3887|[94, 95]|908.30 ± 249.80|3.636|2.959|901.55 ± 245.73|3.669|2.790|0.00%|stays|
|E1_410418450_6.8202|[94]|2065.44 ± 778.53|2.653|1.964|2076.34 ± 950.80|2.184|1.793|17.04%|out|
|E1_423840582_16.1678|[103]|1390.51 ± 240.73|5.776|4.904|1390.51 ± 260.52|5.337|4.940|0.00%|stays|
|E1_447558713_5.5577|[101, 103]|3002.25 ± 345.61|8.687|8.089|3242.44 ± 382.53|8.476|8.642|0.88%|stays|
|E2_254643343_0.6990|[103]|115.86 ± 69.52|1.666|1.634|97.17 ± 62.80|1.547|1.445|20.59%|out|
|E2_261117566_0.6232|[102]|1034.63 ± 340.61|3.038|2.586|985.18 ± 372.39|2.646|2.200|11.40%|stays|
|E2_267546678_2.2325|[95]|7169.65 ± 3120.66|2.297|1.775|5291.97 ± 3494.12|1.515|1.336|76.63%|out|
|E2_363119909_0.5538|[103, 104]|2541.94 ± 317.04|8.018|10.309|2508.47 ± 349.16|7.184|8.951|12.71%|stays|
|E2_421279377_2.0575|[103]|16202.85 ± 1950.90|8.305|10.029|13540.10 ± 1872.10|7.233|7.888|14.53%|stays|
|E2_421279377_7.2268|[103]|12433.73 ± 3541.37|3.511|2.769|8341.76 ± 3951.60|2.111|1.771|0.00%|stays|
|E3_198787224_0.6467|[74]|13203.43 ± 817.78|16.145|18.809|13203.43 ± 817.78|16.145|18.852|0.00%|stays|

Result: **14 of 17 stay**. 17/17 measured; 14 meet both bars among measured. Every A/B null has 2,000 trials. [A_status.json](results/A_status.json) preserves missing/error rows rather than treating them as measured failures.

### Inverted control and uncertainty

The full original procedure is rerun after mask A: same product/sector selections, original known-signal masks, rotation assessment, detrending, inverted flux, discovery-data BLS, hold-out measurement, full-data BLS, odd/even, half-phase secondary, consistency and season gates. Fixed/window modes use original round-8 conventions. This is not a fixed-real-ephemeris sign flip. The population count remains 4,021, and the reference unmasked control is 3/4,021. On the ten fully cached controls whose masks do not change any cadence, the reconstructed procedure exactly reproduces the original ephemerides, hold-out z, full-data SNR and decisions (one undefined hold-out z remains undefined).

Frozen fallback sample used: **yes**, the exact 600 IDs in r13/results/inverted_sample.json, independently regenerated with seed 20261030. Full-pool acquisition requires approximately 20,427 uncached products and 14.3 GB; that raw footprint exceeds the 5 GB acquisition budget for this workflow. [Budget audit](results/full_population_budget.json).

Measured **0/600** inverted passes; target denominator **600**. Scaled MLE expected false count: **0.000 of 4,021**; Jeffreys posterior-mean expected count: **3.345**; equal-tail 95% Jeffreys interval scaled by 4,021: **[0.0032893777279054985, 16.792029687631786]**.

Finite hold-out z exists for **543/600** controls; **32/600** have zero new cadences after the prescribed processing. All 600 remain in the frozen denominator and an undefined hold-out z cannot pass. The same frozen sample had **0/600** passes in the original unmasked outputs. Thus this sample does not demonstrate a mask-induced reduction in the control rate. The requested full-population reference remains **3/4,021** unmasked.

A final product-binding audit found three controls with 11 historical-sector URLs pointing to TESS-SPOC although the frozen original sector metadata required QLP. Those three full procedures were rerun on the matching QLP products from the contemporaneous MAST index; all three still fail. The original outcomes and certificates are archived separately. The final 0/600 result includes the corrections; all 600 controls and all 17 real signals pass the independent product/source/parameter audit. The URL correction overlay is automatic in acquisition and replay scripts. See [deviations](DEVIATIONS.md) and [input audit](results/input_verification.json).

The Jeffreys interval quantifies finite-subsample sampling uncertainty; it does not model changes in target distribution or detection power caused by the mask.

