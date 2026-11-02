#!/usr/bin/env python3
"""Packaging-added saved-result audit: the original job had no verify.py.

Original calculation scripts match.py and posthoc_epoch.py are retained for
inspection. They require omitted source catalogs and overwrite result files.
This read-only check does not reproduce catalog matching or source positions.
"""
from pathlib import Path
import json
import math

root = Path(__file__).resolve().parent
plan = json.loads((root / 'PLAN.json').read_text())
result = json.loads((root / 'RESULT.json').read_text())
for target in result['targets']:
    expected_lambda = result['mean_eb_within_R_null'] * target['period_match_fraction']
    assert math.isclose(target['lambda'], expected_lambda, rel_tol=1e-12)
    assert math.isclose(target['chance'], -math.expm1(-expected_lambda), abs_tol=1e-14)
    assert target['blend_attributed'] == bool(target['matches'] and target['chance'] < 0.01)
assert sum(t['blend_attributed'] for t in result['targets']) == result['attributed_12']
target = next(t for t in result['targets'] if t['ctoi'] == '275488151.01')
eb = next(t for t in target['eb_within_R'] if t['kic'] == 10014830)
assert eb['sep_arcsec'] < plan['radius_arcsec']
assert abs(eb['P_eb'] / target['P_ctoi'] - 1) > plan['period_ratio_tolerance']
assert not target['blend_attributed']
print('PASS: saved chance arithmetic and registered rejection of the post-hoc pair. Catalog matching not rerun.')
