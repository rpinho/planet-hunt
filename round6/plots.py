from fitlib import *
import json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
store = load_store()
E = json.load(open(OUT / 'data' / 'ephem.json')); rho = E['_rho']['rho']
mc = json.load(open(OUT / 'data' / 'mcmc.json'))
PD = OUT / 'plots'
def others_removed(kind, pl):
    lc = concat(store, kind); f = lc['f'].copy()
    for o in ['b', 'c', 'd18', 'd2']:
        if o == pl: continue
        oe = E[o]; f -= transit_model(lc['t'], oe['t0'], oe['P'], oe['k'], oe['b'], rho) - 1
    lc['f'] = f; return lc
def binned(x, y, w):
    edges = np.arange(x.min(), x.max() + w, w); idx = np.digitize(x, edges)
    xb, yb, eb = [], [], []
    for i in np.unique(idx):
        s = idx == i
        if s.sum() >= 3: xb.append(x[s].mean()); yb.append(y[s].mean()); eb.append(y[s].std() / np.sqrt(s.sum()))
    return np.array(xb), np.array(yb), np.array(eb)
# ---------------- per-window plot for d18
pt = json.load(open(OUT / 'data' / 'per_transit.json'))['d18']
ep = E['d18']; P, t0 = ep['P'], ep['t0']
wins = sorted(set(r['epoch'] for r in pt))
fig, axes = plt.subplots(2, 5, figsize=(16, 6.5), sharey=True)
lc2 = others_removed('2min_loose', 'd18'); lcq = others_removed('qlp', 'd18'); lc20 = others_removed('20s', 'd18')
for ax, n in zip(axes.ravel(), wins):
    tc = t0 + n * P
    for lc, col, lab, bw in [(lc2, 'k', '2-min PDCSAP', 20), (lcq, 'tab:orange', 'QLP FFI', 20)]:
        s = np.abs(lc['t'] - tc) < 0.5
        if s.sum() < 10: continue
        xb, yb, eb = binned((lc['t'][s] - tc) * 24, lc['f'][s], bw / 60)
        ax.errorbar(xb, (yb - 1) * 1e6 - (0 if col == 'k' else 2500), eb * 1e6, fmt='o', ms=3, color=col, label=lab)
    tt = np.linspace(tc - 0.5, tc + 0.5, 500)
    m = transit_model(tt, tc, P, ep['k'], ep['b'], rho)
    ax.plot((tt - tc) * 24, (m - 1) * 1e6, 'r-', lw=1)
    ax.plot((tt - tc) * 24, (m - 1) * 1e6 - 2500, 'r-', lw=1)
    r = [x for x in pt if x['epoch'] == n and x['kind'] in ('2min_loose',)]
    sec = r[0]['sector'] if r else '?'
    a = f"amp {r[0]['amp_lin']:.2f}$\\pm${r[0]['amp_lin_err']:.2f}" if r else ''
    flag = ' (SPOC-flagged cadences)' if n == -66 else (' (partial)' if n == 67 else '')
    ax.set_title(f'S{sec} E{n}{flag}\n{a}', fontsize=9)
    ax.axvline(0, color='0.8', lw=0.5); ax.set_xlim(-12, 12)
axes[0, 0].set_ylim(-5000, 2000); axes[0, 0].legend(fontsize=7, loc='lower left')
for ax in axes[1]: ax.set_xlabel('hours from predicted mid-transit')
for ax in axes[:, 0]: ax.set_ylabel('flux - 1 (ppm); QLP offset -2500')
fig.suptitle('TOI-4342 18.908-d candidate: every TESS window (linear ephemeris, 20-min bins)')
fig.tight_layout(); fig.savefig(PD / 'd18_every_window.png', dpi=130); plt.close(fig)
# --------------- S101 detail
fig, axes = plt.subplots(4, 1, figsize=(10, 9), sharex=True)
tc = t0 + 64 * P
for ax, (kind, lab) in zip(axes, [('2min', '2-min PDCSAP'), ('2min_sap', '2-min SAP (local detrend)'), ('20s', '20-s PDCSAP'), ('qlp', 'QLP FFI 200-s')]):
    d = store[(kind, 101)]
    s = np.abs(d['t'] - tc) < 2.5
    xb, yb, eb = binned(d['t'][s] - tc, d['fd'][s], 20 / 1440)
    ax.errorbar(xb * 24, (yb - 1) * 1e6, eb * 1e6, fmt='.', color='k', ms=3)
    ax.axvspan(-mc['d18']['T14_h'][1] / 2, mc['d18']['T14_h'][1] / 2, color='r', alpha=0.15)
    ax.set_ylabel('ppm'); ax.set_title(lab, fontsize=9); ax.set_ylim(-3500, 2500)
axes[-1].set_xlabel('hours from predicted mid-transit (S101, E64)')
fig.suptitle('The "S101 null": the transit sits 0.2 d after a 2.4-d gap in the 2-min data; round 5 trimmed 0.25 d after gaps')
fig.tight_layout(); fig.savefig(PD / 'd18_s101_detail.png', dpi=130); plt.close(fig)
# --------------- alias / period scan
z = np.load(OUT / 'data' / 'period_scan.npz'); al = json.load(open(OUT / 'data' / 'alias.json'))
fig, ax = plt.subplots(1, 2, figsize=(13, 4))
ax[0].plot(z['P'], z['snr'], 'k-', lw=0.4); ax[0].set_xscale('log'); ax[0].set_xlabel('period (d)'); ax[0].set_ylabel('box-depth SNR (best phase)')
for m_, lab in [(0.5, 'P/2'), (1, 'P'), (1.5, '3P/2'), (2, '2P'), (3, '3P')]:
    ax[0].axvline(P * m_, color='r', alpha=0.3); ax[0].text(P * m_, ax[0].get_ylim()[1] * 0.97, lab, color='r', fontsize=8)
sel = (z['P'] > 18.0) & (z['P'] < 19.8)
ax[1].plot(z['P'][sel], z['snr'][sel], 'k-', lw=0.5); ax[1].set_xlabel('period (d)'); ax[1].set_title('zoom: cycle-count aliases near 18.9 d', fontsize=9)
fig.tight_layout(); fig.savefig(PD / 'd18_period_scan.png', dpi=130); plt.close(fig)
# --------------- folds d18 and d2
fig, axes = plt.subplots(1, 2, figsize=(13, 4.2))
for ax, pl, bw, lim in [(axes[0], 'd18', 10, 6), (axes[1], 'd2', 10, 4)]:
    e_ = E[pl]; lc = others_removed('2min', pl)
    ph = ((lc['t'] - e_['t0'] + 0.5 * e_['P']) % e_['P'] - 0.5 * e_['P']) * 24
    s = np.abs(ph) < lim
    ax.plot(ph[s], (lc['f'][s] - 1) * 1e6, ',', color='0.75')
    xb, yb, eb = binned(ph[s], lc['f'][s], bw / 60)
    ax.errorbar(xb, (yb - 1) * 1e6, eb * 1e6, fmt='o', color='k', ms=3)
    m = mc[pl]
    tt = np.linspace(-lim, lim, 500) / 24
    mod = transit_model(tt, 0, m['P'][1], m['k'][1], m['b'][1], m['rho_circ'][1])
    ax.plot(tt * 24, (mod - 1) * 1e6, 'r-')
    ax.set_ylim(-2500 if pl == 'd18' else -900, 1200 if pl == 'd18' else 600); ax.set_xlabel('hours from mid-transit'); ax.set_ylabel('ppm')
    ax.set_title(f"{pl}: P={m['P'][1]:.5f} d, depth {m['depth_ppm'][1]:.0f} ppm, 9 sectors 2-min", fontsize=9)
fig.tight_layout(); fig.savefig(PD / 'folds_d18_d2.png', dpi=130); plt.close(fig)
# --------------- O-C
ptall = json.load(open(OUT / 'data' / 'per_transit.json'))
fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
for ax, pl in zip(axes, ['b', 'c', 'd18']):
    sel = [r for r in ptall[pl] if r['kind'] == '2min' and r['coverage'] > 0.5]
    tt = np.array([r['post_tc'] for r in sel]); n = np.array([r['epoch'] for r in sel]); e_ = np.array([r['post_err_min'] for r in sel])
    A = np.vstack([np.ones_like(n), n]).T / (e_[:, None] / 1440)
    c_, *_ = np.linalg.lstsq(A, tt / (e_ / 1440), rcond=None)
    oc = (tt - (c_[0] + c_[1] * n)) * 1440
    ax.errorbar(tt, oc, e_, fmt='o', ms=3, color='k')
    for name, col in [('b+c', 'tab:blue'), ('b+c+d18', 'tab:red')]:
        f_ = OUT / 'dyn' / f'tt_{name}.npz'
        if f_.exists():
            z2 = np.load(f_)
            if pl in z2.files:
                T = z2[pl]; nn = np.arange(len(T)); B = np.vstack([np.ones_like(nn), nn]).T
                cc, *_ = np.linalg.lstsq(B, T, rcond=None)
                ax.plot(T, (T - B @ cc) * 1440, '-', color=col, lw=0.8, label=f'N-body {name} (circular, M-R masses)')
    ax.axhline(0, color='0.7'); ax.set_ylabel(f'{pl} O-C (min)'); ax.set_ylim(-60, 60)
axes[0].legend(fontsize=7); axes[-1].set_xlabel('BTJD')
fig.tight_layout(); fig.savefig(PD / 'ttv.png', dpi=130); plt.close(fig)
# --------------- pixel maps
loc = json.load(open(OUT / 'data' / 'tpf_template_loc.json'))
grid = np.arange(-110, 110.1, 2.0)
nb = json.load(open(OUT / 'data' / 'neighbours_local.json'))
fig, axes = plt.subplots(1, 4, figsize=(17, 4.3))
for ax, pl in zip(axes, ['b', 'c', 'd18', 'd2']):
    m = np.load(OUT / 'data' / f'locmap_{pl}.npy')
    im = ax.imshow(np.sqrt(m), origin='lower', extent=[grid[0], grid[-1], grid[0], grid[-1]], cmap='viridis_r', vmax=15)
    ax.contour(grid, grid, m, levels=[2.3, 11.8], colors=['w', 'w'], linestyles=['-', '--'])
    for x in nb:
        xi, eta = x['sep'] * np.sin(np.radians(x['pa'])), x['sep'] * np.cos(np.radians(x['pa']))
        mim = (x['need1150'] if pl != 'd2' else x['need284']) < 1
        ax.plot(xi, eta, 'o', mfc='none', mec='r' if mim else '0.6', ms=4 + max(0, 18 - x['tmag']) * 1.5)
    ax.plot(0, 0, '*', color='gold', ms=12)
    ax.set_title(f"{pl}: $\\Delta\\chi^2$(target) = {loc[pl]['dchi2_target']:.1f}", fontsize=10)
    ax.set_xlabel('East offset (arcsec)'); ax.invert_xaxis()
axes[0].set_ylabel('North offset (arcsec)')
fig.colorbar(im, ax=axes, shrink=0.8, label='sqrt(Δχ²)')
fig.suptitle('Difference-image localisation (empirical-PRF template, 9 sectors of 2-min TPFs). Red circles: Gaia stars bright enough to mimic the depth.', fontsize=10)
fig.savefig(PD / 'pixel_localisation.png', dpi=130, bbox_inches='tight'); plt.close(fig)
print('plots done')
