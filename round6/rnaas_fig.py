from fitlib import *
import json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8.5})
store = load_store(); E = json.load(open(OUT / 'data' / 'ephem.json')); rho = E['_rho']['rho']
mp = json.load(open(OUT / 'data' / 'mcmc_rhoprior.json'))['d18']
pt = json.load(open(OUT / 'data' / 'per_transit.json'))['d18']
lc = concat(store, '2min'); f = lc['f'].copy()
for o in ['b', 'c', 'd2']:
    oe = E[o]; f -= transit_model(lc['t'], oe['t0'], oe['P'], oe['k'], oe['b'], rho) - 1
P, t0 = mp['P'][1], mp['t0'][1]
ph = ((lc['t'] - t0 + 0.5 * P) % P - 0.5 * P) * 24
fig, ax = plt.subplots(1, 3, figsize=(11, 3.3), gridspec_kw=dict(width_ratios=[1.2, 1.2, 1]))
s = np.abs(ph) < 6
ax[0].plot(ph[s], (f[s] - 1) * 1e3, '.', color='0.8', ms=1, rasterized=True)
edges = np.arange(-6, 6.01, 0.25); idx = np.digitize(ph[s], edges)
xb = [ph[s][idx == i].mean() for i in np.unique(idx)]; yb = [(f[s][idx == i].mean() - 1) * 1e3 for i in np.unique(idx)]
eb = [f[s][idx == i].std() / np.sqrt((idx == i).sum()) * 1e3 for i in np.unique(idx)]
ax[0].errorbar(xb, yb, eb, fmt='o', color='k', ms=3)
tt = np.linspace(-6, 6, 600) / 24
ax[0].plot(tt * 24, (transit_model(tt, 0, P, mp['k'][1], mp['b'][1], mp['rho_circ'][1]) - 1) * 1e3, 'r-', lw=1.2)
ax[0].set_ylim(-2.2, 1.2); ax[0].set_xlabel('Hours from mid-transit'); ax[0].set_ylabel('Relative flux - 1 (ppt)')
ax[0].text(0.03, 0.05, '(a) 9 sectors, 2-min, 15-min bins', transform=ax[0].transAxes)
dep = mp['depth_ppm'][1] / 1e3
for kind, mk, col, lab in [('2min_loose', 'o', 'k', 'SPOC 2-min'), ('qlp', 's', 'tab:orange', 'QLP FFI')]:
    rows = [r for r in pt if r['kind'] == kind]
    x = np.array([r['tc_lin'] for r in rows]) + (0 if kind == '2min_loose' else 25)
    ax[1].errorbar(x, [r['amp_lin'] * dep for r in rows], [r['amp_lin_err'] * dep for r in rows], fmt=mk, color=col, ms=4, label=lab, mfc='none' if kind == 'qlp' else col)
ax[1].axhline(dep, color='r', lw=0.8); ax[1].axhline(0, color='0.6', lw=0.6)
ax[1].set_xlabel('BTJD (QLP offset +25 d)'); ax[1].set_ylabel('Depth at predicted time (ppt)'); ax[1].legend(loc='lower right', fontsize=7)
ax[1].text(0.03, 0.92, '(b) every predicted window', transform=ax[1].transAxes)
ax[1].set_ylim(-0.6, 2.4)
grid = np.arange(-110, 110.1, 2.0); m = np.load(OUT / 'data' / 'locmap_d18.npy')
im = ax[2].imshow(np.sqrt(m), origin='lower', extent=[grid[0], grid[-1], grid[0], grid[-1]], cmap='viridis_r', vmax=12)
ax[2].contour(grid, grid, m, levels=[11.8], colors='w', linestyles='--')
for x in json.load(open(OUT / 'data' / 'neighbours_local.json')):
    xi, eta = x['sep'] * np.sin(np.radians(x['pa'])), x['sep'] * np.cos(np.radians(x['pa']))
    if x['need1150'] < 1: ax[2].plot(xi, eta, 'o', mfc='none', mec='r', ms=4)
ax[2].plot(0, 0, '*', color='gold', ms=10); ax[2].invert_xaxis()
ax[2].set_xlim(80, -80); ax[2].set_ylim(-80, 80)
ax[2].set_xlabel('East offset (arcsec)'); ax[2].set_ylabel('North offset (arcsec)')
ax[2].text(0.03, 0.92, '(c) difference-image $\\sqrt{\\Delta\\chi^2}$', transform=ax[2].transAxes, color='w')
fig.tight_layout()
O = OUT.parent / 'round6_toi4342_rnaas'
fig.savefig(O / 'fig1.pdf', dpi=200); fig.savefig(O / 'fig1.png', dpi=200)
print('ok')
