# Deviations and explicit operational choices

No post-freeze protocol or measurement changes; no recovery reruns or threshold/phase/model tuning. Gate failed, so real EB/CTOI images, 14–15 fits, real a/b and residuals were withheld exactly as required. Failure of power is the scientific stopping result, not an external blocker.

Frozen operational definitions: positive-depth sign for negative in-minus-out images; per-pixel temporal nuisance constants eliminated by event contrasts; one spatial background nuisance per sector; OOT complete-scene Tmag calibration; inherited ±0.126-d Kepler mean-depth window with broad Villanova exclusions; secondary center 0.5007; conservative one-sided U95=max(0,a)+1.6448536269514722 sigma. These were explicit before any image/flux read. See PROTOCOL.md.

The first network refresh briefly inherited nice 20 (nice 10 added to a nice-10 launcher), then was reniced to 10. Scientific stages ran at actual nice 10, one process and one BLAS/OMP thread. No source directory was modified; bytecode writes were disabled. No cutouts or dependencies were downloaded.

ExoFOP's obsolete endpoint returned notices; the documented target.php endpoint was used for the successful refresh. ADS returned empty search bodies and remains UNCHECKED. Generic GitHub/AI-lab/community prior-art checks are prohibited by the narrower user network list. DIAmante column documentation was read after the gate for the report only; catalogue/product documents already existed before it.

The independent verifier was implemented after the gate and imports no measurement modules. It separately recomputes all 400 WLS fits and decisions from saved sufficient statistics, rather than rerunning extraction or changing the estimator. Verification code adds an inherited Kepler event restack; it does not alter measurement outputs. The optional reproduce.sh was written for review, syntax-checked only, and not executed after POWER_FAIL.
