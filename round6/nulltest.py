"""Noise null: same box-depth comb statistic on the data with the real 18.9-d windows removed,
max over phase, for 2000 trial periods in 12-30 d. How often does noise reach the observed value?"""
from common import *
import json
z = np.load(OUT / 'data' / 'single_event.npz'); g, D, S = z['grid'], z['D'], z['S']
E = json.load(open(OUT / 'data' / 'ephem.json'))
good = np.isfinite(D) & (S > 0) & ~in_tr(g, E['d18']['P'], E['d18']['t0'], E['d18']['dur'], widen=4)
g, D, S = g[good], D[good], S[good]; w = 1 / S ** 2
step = 10 / 1440
rng = np.random.default_rng(3)
Ps = rng.uniform(12, 30, 3000)
mx = []
for P in Ps:
    ph = ((g - 2900.0) / P) % 1.0; nb = int(np.ceil(P / step)); idx = np.floor(ph * nb).astype(int)
    sw = np.bincount(idx, w, minlength=nb); sd = np.bincount(idx, w * D, minlength=nb)
    nep = np.bincount(idx, minlength=nb)
    with np.errstate(invalid='ignore', divide='ignore'):
        s = sd / np.sqrt(sw)
    s[nep < 5] = np.nan  # require >=5 single-event samples (like >=5 windows)
    mx.append(np.nanmax(s))
mx = np.array(mx)
obs = json.load(open(OUT / 'data' / 'alias.json'))['peaks'][0][1]
print('observed comb SNR at 18.908 d:', obs)
print('noise: max-over-phase per period: median %.2f, 99%% %.2f, 99.9%% %.2f, max %.2f' % (np.median(mx), *np.percentile(mx, [99, 99.9]), mx.max()))
print('fraction of trial periods with noise >= observed:', np.mean(mx >= obs))
json.dump(dict(observed=obs, median=float(np.median(mx)), p99=float(np.percentile(mx, 99)), p999=float(np.percentile(mx, 99.9)), max=float(mx.max()), n=len(mx), frac_ge_obs=float(np.mean(mx >= obs))), open(OUT / 'data' / 'nulltest.json', 'w'), indent=1)
