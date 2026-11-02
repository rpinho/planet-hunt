#!/usr/bin/env python3
"""Read-only checks of this public bundle; no photometry or network access."""
from pathlib import Path
import hashlib
import json
import math
import re
import statistics

ROOT = Path(__file__).resolve().parent
NUMBER = re.compile(r'\d+(?:\.\d+)?|\beight\b|\btwo\b')


def read(name):
    return json.loads((ROOT / name).read_text())


def job(n, name):
    return read(f'planets_r{n}_1102/{name}')


def derived():
    """Arithmetic on copied results only, with all inputs identified below."""
    target = next(t for t in job(27, 'PLAN.json')['targets']
                  if t['ctoi'] == '275488151.01')
    plan = job(36, 'PLAN.json')
    eb = plan['EB']
    ctoi = plan['CTOI']
    meta = job(30, 'aperture_manifest.json')
    tess = job(36, 'outputs/EB_TESS.json')
    review = job(32, 'outputs/attachment_review.json')
    start = review['light_curve']['btjd_start'] + 2457000
    end = review['light_curve']['btjd_end'] + 2457000
    anchor_cycle = round((ctoi['epoch_BJD'] - eb['epoch_BJD']) / eb['period'])
    anchor = eb['epoch_BJD'] + anchor_cycle * eb['period']
    offset = tess['sectors'][0]['timing']['best_offset_days']

    def drift(t):
        # Continuous accumulated phase drift, evaluated at boundary cadences.
        k = (t - ctoi['epoch_BJD']) / ctoi['period']
        return (ctoi['epoch_BJD'] - anchor
                + k * (ctoi['period'] - eb['period']) - offset) * 24

    mids = {}
    for s in meta['sectors']:
        if s['sector'] in [14, 15]:
            k = round((s['mid_bjd'] - ctoi['epoch_BJD']) / ctoi['period'])
            mids[str(s['sector'])] = (ctoi['epoch_BJD'] - anchor
                + k * (ctoi['period'] - eb['period']) - offset) * 24
    s40 = next(s for s in meta['sectors'] if s['sector'] == 40)
    k40 = round((s40['mid_bjd'] - ctoi['epoch_BJD']) / ctoi['period'])
    inj = job(37, 'outputs/injections_summary.json')['trials']['EB']
    real = job(37, 'outputs/real.json')['profiles']['2']['combined']['delta_chi2']
    valid = [t['delta_chi2'] for t in inj if t['valid']]
    apertures = job(32, 'outputs/apertures.json')['quarters']
    return {
        'tmag_difference': job(31, 'outputs/metadata.json')['eb_tmag']
                           - job(31, 'outputs/metadata.json')['target_tmag'],
        'kepler_primary_percent': job(30, 'outputs/kepler.json')['combined']['primary']['depth'] * 100,
        'catalogue_primary_width_hours': 0.1530 * eb['period'] * 24,
        'ctoi_duration_days': target['duration_hours'] / 24,
        'kepler_start_year': int(apertures[0]['table_header']['DATE-OBS'][:4]),
        'kepler_end_year': int(apertures[-1]['table_header']['DATE-END'][:4]),
        'planned_period_tolerance_percent': job(29, 'PLAN.json')['period_ratio_tolerance'] * 100,
        'period_difference_percent': abs(ctoi['period'] / eb['period'] - 1) * 100,
        'period_difference_sigma': abs(ctoi['period'] - eb['period']) / target['period_error_days'],
        'tess_to_kepler_primary_depth_ratio': tess['combined']['primary']['depth'][0]
                                              / plan['kepler']['primary']['depth'][0],
        'boundary_timing_hours': {'start': drift(start), 'end': drift(end)},
        'midpoint_cycle_timing_hours': mids,
        'leak_fraction_percent': plan['L'][1] / plan['T'][1] * 100,
        'other_stars_light_percent': target['tic_contratio'] / (1 + target['tic_contratio']) * 100,
        'later_sectors': [s['sector'] for s in meta['sectors'] if s['sector'] >= 40],
        'later_sector_count': sum(s['sector'] >= 40 for s in meta['sectors']),
        'off_science_sector_count': len(job(31, 'outputs/metadata.json')['dropped_sectors']),
        'sector40_ctoi_cycle': k40,
        'sector40_epoch_sigma_days': math.hypot(target['epoch_error_days'],
                                               k40 * target['period_error_days']),
        'lc2_eb_injection_count': len(inj),
        'lc2_eb_valid_fit_count': len(valid),
        'lc2_eb_injections_at_least_as_negative_as_real': sum(d <= real for d in valid),
        'lc2_eb_minimum_delta_chi2': min(valid),
        'lc2_real_delta_chi2': real,
    }


def resolve(reference):
    name, sep, key = reference.partition('#')
    path = ROOT / name
    assert path.is_file(), reference
    if not sep:
        return path.read_text()
    if key.startswith('L'):
        return path.read_text().splitlines()[int(key[1:]) - 1]
    value = read(name)
    for part in key.strip('/').split('/'):
        value = value[int(part)] if isinstance(value, list) else value[part]
    return value


def main():
    provenance = read('PROVENANCE.json')
    note = (ROOT / 'NOTE.md').read_bytes()
    header, body = note.split(b'\n\n', 1)
    assert len(header.splitlines()) == 3
    assert hashlib.sha256(body).hexdigest() == provenance['source_note_sha256']
    assert len(body) == provenance['source_note_bytes']
    rows = read('EVIDENCE.json')
    tokens = list(NUMBER.finditer(note.decode()))
    assert len(rows) == len(tokens)
    for row, token in zip(rows, tokens):
        assert row['number'] == token.group()
        assert row['offset'] == token.start()
        for reference in row['references']:
            resolve(reference)
        if row.get('numeric_reference'):
            value = resolve(row['numeric_reference'])
            value = float(value) * row.get('scale', 1)
            if row.get('absolute'):
                value = abs(value)
            expected = {'eight': 8, 'two': 2}.get(row['number'], row['number'])
            expected = float(expected)
            tolerance = row.get('tolerance', 1e-12)
            assert abs(value - expected) <= tolerance, (row['context'], value, expected)
    actual = derived()
    assert actual == read('DERIVED.json')
    assert actual['lc2_eb_injections_at_least_as_negative_as_real'] == 0
    assert actual['lc2_eb_injection_count'] == 200
    seal = read('BUNDLE.sha256.json')
    for name, digest in seal.items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    files = [p for p in ROOT.rglob('*') if p.is_file()]
    assert set(seal) == {p.relative_to(ROOT).as_posix() for p in files
                         if p.name != 'BUNDLE.sha256.json'}
    for path in files:
        assert path.suffix not in {'.fits', '.npz', '.npy', '.log'}
        assert path.stat().st_size < 2000000, path
        assert not {'logs', 'cache', '.venv', '__pycache__'} & set(path.parts)
        text = path.read_text()
        assert not re.search(r'/(?:Users|home|root|workspace|private|tmp)/', text), path
    untraced = [r for r in rows if r['match'].startswith('UNTRACED')]
    print(f'PASS: source note bytes; {len(rows)} numeric occurrences; {len(seal)} sealed files; '
          f'derived arithmetic; size/type/path checks. Untraced numbers: {len(untraced)}.')
    print('Scope: public saved-result bundle. Original full scientific verifiers require omitted inputs/history.')


if __name__ == '__main__':
    main()
