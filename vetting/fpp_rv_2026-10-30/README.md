Reruns for the TOI-4342 correction (2026-10-30).

FPP (TRICERATOPS): the original nine-sector 2-minute fold, 2-minute bins, 120 s modeled exposure, five runs per configuration. MOLUSC configurations use five companion samples at N=1,000,000. `verified_padding_probes.csv` shows that N changes those results: TRICERATOPS zero-pads the surviving MOLUSC companions to N, so the draw count sets the companion prior. The values in the note are the means and sample standard deviations of each `tri_*.json`.

Radial velocities: `rv_fit.py` fits a circular orbit at P=18.908226 d, T0=59904.5666 RJD to the 55 public ESPRESSO LBL velocities of Parc et al. 2026 (DACE), with b and c at published ephemerides and the published activity kernel at fixed hyperparameters. `rv.json` holds the results. The limits are conditional: GP hyperparameters are not marginalized.

1.96 d: `short4342.json` holds the depth with and without 0.5 d around orbit boundaries, the orbit-locking controls, and the free-search period check, which fails for the literal period.

`rnaas_fig.py` is the figure script with the panel (a) label corrected to eight sectors.
