import json, csv, math
import numpy as np
# ICRS->Galactic rotation (Hipparcos/ICRS, as in astropy); icrs = T.T @ gal
T = np.array([[-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
              [ 0.4941094278755837, -0.4448296299600112,  0.7469822444972189],
              [-0.8676661490190047, -0.1980763734312015,  0.4559837761750669]])
def vec(a, d):
    a, d = np.radians(a), np.radians(d)
    return np.stack([np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)], -1)
def sep_arcsec(v1, v2):  # v1 (n,3), v2 (m,3) -> (n,m)
    return np.degrees(np.arccos(np.clip(v1 @ v2.T, -1, 1))) * 3600
plan = json.load(open('PLAN.json')); co = json.load(open('sources/tic_coords.json'))
R = plan['radius_arcsec']; tol = plan['period_ratio_tolerance']; H = plan['harmonics']
eb = [r for r in csv.reader(l for l in open('sources/villanova_ebs.csv') if not l.startswith('#')) if r and r[0].strip()]
kic = np.array([int(r[0]) for r in eb]); per = np.array([float(r[1]) for r in eb])
gl = np.array([float(r[6]) for r in eb]); gb = np.array([float(r[7]) for r in eb])
pos = gl != -1
E = vec(gl[pos], gb[pos]) @ T  # rows: (T.T @ g).T = g.T @ T
ekic, eper = kic[pos], per[pos]
def pmatch(P, p):  # p array
    return np.any([np.abs(p/(h*P) - 1) <= tol for h in H], axis=0)
# null: EB counts within R around control positions, excluding the control's own KIC
C = plan['null_positions']
CS = vec(np.array([c['ra'] for c in C]), np.array([c['dec'] for c in C]))
S = sep_arcsec(CS, E)
idx_c, idx_e = np.nonzero(S <= R)
own = ekic[idx_e] == np.array([C[i]['kic'] for i in idx_c])
mean_count = np.sum(~own) / len(C)
out = {'n_eb': int(len(kic)), 'n_eb_no_position': int((~pos).sum()), 'n_null_positions': len(C),
       'mean_eb_within_R_null': mean_count, 'targets': []}
for t in plan['targets']:
    c = co[str(t['tic'])]; P = t['period_days']
    s = sep_arcsec(vec(np.array([c['ra']]), np.array([c['dec']])), E)[0]
    near = s <= R
    selfin = bool(np.any(kic == t['kic']))
    frac = float(np.mean(pmatch(P, per)))
    lam = mean_count * frac; chance = 1 - math.exp(-lam)
    m = near & pmatch(P, eper) & (ekic != t['kic'])
    matches = [{'kic': int(k), 'sep_arcsec': float(x), 'P_eb': float(p), 'ratio': float(p/P)} for k, x, p in zip(ekic[m], s[m], eper[m])]
    neigh = [{'kic': int(k), 'sep_arcsec': round(float(x),1), 'P_eb': float(p)} for k, x, p in zip(ekic[near], s[near], eper[near]) if k != t['kic']]
    attributed = bool(matches) and chance < 0.01
    out['targets'].append({'ctoi': t['ctoi'], 'kic': t['kic'], 'P_ctoi': P, 'target_in_eb_catalogue': selfin,
        'self_P_eb': [float(p) for p in per[kic == t['kic']]], 'eb_within_R': neigh, 'matches': matches,
        'period_match_fraction': frac, 'lambda': lam, 'chance': chance, 'blend_attributed': attributed})
hi = set(plan['high_snr_five'])
out['attributed_12'] = sum(x['blend_attributed'] for x in out['targets'])
out['attributed_high5'] = sum(x['blend_attributed'] for x in out['targets'] if x['ctoi'] in hi)
out['self_in_catalogue'] = [x['ctoi'] for x in out['targets'] if x['target_in_eb_catalogue']]
json.dump(out, open('RESULT.json', 'w'), indent=1)
print('EBs', out['n_eb'], 'nopos', out['n_eb_no_position'], 'mean EB within R at null positions', round(mean_count, 4))
for x in out['targets']:
    print(x['ctoi'], 'P', x['P_ctoi'], 'self', x['target_in_eb_catalogue'], x['self_P_eb'], 'nearEB', x['eb_within_R'], 'match', x['matches'], 'chance %.2e' % x['chance'], 'ATTR' if x['blend_attributed'] else '')
print('attributed', out['attributed_12'], '/12; high5', out['attributed_high5'], '/5; self', out['self_in_catalogue'])
