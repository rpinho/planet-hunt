# planets_r29_1102: are the r28 periodic CTOIs blends of known Kepler EBs?

**Frozen bar: 0 of 12 attributed (0 of 5 high-SNR rows). No target is itself in the Kepler EB catalogue.** Plan frozen in git (fc74e00) before the catalogue was opened. One deviation: the Villanova online export (rev 2019-08-08) replaced VizieR J/AJ/151/68, which was unreachable (DEVIATIONS.md).

Rule: an EB within 120 arcsec with P_EB/(h P_CTOI) within 0.1% of 1 for h in {1/2, 1, 2}, and chance < 0.01. Null: 0.0943 catalogue EBs within 120 arcsec on average over 700 r27 control positions; per-target chance 0 to 1.9e-4. EBs within 120 arcsec of any target: 3 (around 63120478, 275488151, 351189201), none at a matching period inside 0.1%.

## Post-hoc lead (found after the frozen result; descriptive, not a bar outcome)

TIC 275488151.01 (KIC 10014894, P 3.02553 ± 0.00119 d, 558 ppm, Montalto 2022) has Villanova EB KIC 10014830 at 82.7 arcsec, Kp 14.83 (target 11.67), P 3.0305267 d, morph 0.61. The periods differ by 0.165% (4.2 sigma of the CTOI error), outside the frozen tolerance. But the CTOI epoch, BJD 2458682.5504, lands 0.0002 d from the EB's predicted primary eclipse at cycle 1226 (phase offset 0.0001). A wider descriptive net (EBs within 300 arcsec, period within 1%, all 12 targets) finds only this pair. Rough chance of a phase hit that close for any of the 12: about 1e-5. r28 searched the target's Kepler flux over P ± 5 sigma (includes 3.0305) and found nothing on target, which fits a source 83 arcsec away and outside the Kepler aperture.

Reading: the TESS signal on TIC 275488151 is most likely the eclipse of KIC 10014830 leaking into the TESS aperture (21 arcsec pixels; 83 arcsec is about 4 pixels), with a short-baseline TESS period that is 0.005 d off. Not yet tested in TESS pixels. ExoFOP shows the CTOI as PC with no TFOP disposition, notes or time series (prior/exofop_275488151.json, 11-02). One web search found nothing on the three IDs; that is a search limit, not proof of absence.

Next: a frozen TESS test (centroid shift toward KIC 10014830 and which ephemeris the dips follow) before anything is said publicly.

Credit: Villanova Kepler EB Catalog (Prsa et al. 2011, Slawson et al. 2011, Kirk et al. 2016); Montalto (DIAmante) for the CTOI; ExoFOP-TESS.

Reproduce: `../.venv/bin/python -I match.py` then `posthoc_epoch.py`.

## Correction (11-02, after review)

"Rough chance of a phase hit that close for any of the 12: about 1e-5" is wrong. It treats the 0.0002 d offset as the matching window, but the CTOI epoch error is 0.0122 d, so a random EB phase lands inside that error about 1% of the time (Claude Science review, frame dc7acbc1). The stronger evidence is the catalogued EB at 82 arcsec with a period within 0.165%, roughly 1e-4 to 1e-3 by chance after trials. The frozen result (0/12) is unchanged.
