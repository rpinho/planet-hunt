"""Reproduce the round-8 pixel vet and TRICERATOPS fold for TIC 262843927's second signal, using cached sectors only.
Run with .venv/bin/python round8/vet_1452.py from the repository root.
Then run (PYTHONDONTWRITEBYTECODE=1, nice -n 15):
.trivenv/bin/python round6/run_tri.py --tic 262843927
  --fold runs/round8_crowd/tri/tic262843927_1452_fold.npz
  --lc-dir runs/round8_crowd/tri/E1_262843927_24.1087_apertures.npz
  --runs 3 --n 1000000 --out runs/round8_crowd/tri
  --name triceratops_tic262843927_1452
  --local-cutouts cache/tesscut/52.57321_-63.50180
The 24.1-day reference had no MOLUSC file, contrast curve or cleared stars.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import json
import re
import numpy as np
import deep_vet as D
from planet_hunt import pixelvet as pv
import holdout as H
from holdout import hs, ms, r5

OUT = Path('runs/round8_crowd')
NAME = 'tic262843927_1452'
old = json.loads((OUT / 'deepvet/E1_262843927_24.1087.json').read_text())
context = json.loads((OUT / 'context_262843927_14.517.json').read_text())
f = context['real']['full']
spec = dict(old['spec'], name=NAME, period=f['P'], t0=f['t0'],
            duration_h=f['dur'] * 24, depth_ppm=f['depth_ppm'])
cdir = pv.cut_dir(spec['ra'], spec['dec'])
cuts = sorted((int(re.search(r'-s(\d{4})-', p.name)[1]), p)
              for p in cdir.glob('*_11x11_astrocut.fits'))
assert len(cuts) == 8 and all(pv._readable(p) for _, p in cuts)
spec['sectors'] = [s for s, _ in cuts]
(OUT / f'{NAME}_targets.json').write_text(json.dumps({'targets': [spec]}, indent=2))
# No live sector discovery or download: use exactly the original cached cutouts.
pv.available_sectors = lambda *args: spec['sectors']
pv.fetch_cutouts = lambda *args, **kwargs: cuts
D.OUT = OUT / 'deepvet'
if not (D.OUT / f'{NAME}.json').exists():
    res = D.run_target(spec)
    (D.OUT / f'{NAME}.json').write_text(json.dumps(D._jsonable(res), indent=1, allow_nan=False))

# Same fold preparation as round8/trifold.py, with the second signal's saved ephemeris.
plans = {int(k): v for k, v in json.loads((OUT / 'plans.json').read_text()).items()}
secs = H.load_sectors(spec['tic'], plans[spec['tic']], {'rule': context['row']['rule'], 'src_sectors': []})
assert sorted(s['meta']['sector'] for s in secs) == sorted(s['sector'] for s in context['sectors'])
P, t0, dur = spec['period'], spec['t0'], spec['duration_h'] / 24
kl = hs.known_list(r5.known_ephemerides(spec['tic'], *r5.load_catalogs()), spec['rstar'], spec['mstar'])
km = lambda t: hs.mask_known(t, kl)
rot = hs.rotation(secs, km)
lc = hs.detrend(secs, hs.W_SHORT if P < hs.SPLIT_P else hs.W_LONG,
                lambda t: km(t) | ms.in_tr(t, P, t0, dur, 1.6), rot)
ph = ms.phase(lc['t'], P, t0)
near = np.abs(ph) < max(3 * dur, 0.25 * P if P < 1 else 0.25)
np.savez(OUT / 'tri' / f'{NAME}_fold.npz', phase_d=ph[near], flux=lc['flat'][near],
         err=lc['dy'][near], period=P, depth=spec['depth_ppm'] * 1e-6,
         sectors=np.array(sorted(set(lc['sector'].astype(int)))))
print('Saved pixel vet and fold; reuse E1_262843927_24.1087_apertures.npz unchanged.', flush=True)
