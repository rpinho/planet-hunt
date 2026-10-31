Reruns for the TIC 144439305 correction (2026-10-30).

FPP (TRICERATOPS): selected 200-second cadence at its true exposure, original apertures and cleared-neighbour list, five runs per configuration. MOLUSC runs use five companion samples at N=1,000,000. `verified_padding_probes.csv` shows that N changes the result: TRICERATOPS zero-pads the surviving MOLUSC companions to N, so the draw count sets the companion prior. Values in the note are means and sample standard deviations of each `tri_*.json`.

Hold-out: `tic144_holdout.json` has the QLP comparison with the public orbit-edge mask (first 1.0 d and last 0.5 d of each orbit, Sectors 94-104), z calibrated against 2,000 random ephemerides.
