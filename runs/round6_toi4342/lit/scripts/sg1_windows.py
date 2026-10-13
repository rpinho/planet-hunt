import sys, glob, os, json
import numpy as np
d = sys.argv[1]
OFF = 2457000.0
eph = {
 'A_18.9083': (1675.535, 18.9083, 1.89),
 'A_half_9.45415': (1675.535, 18.9083/2, 1.89),
 'B_1.95769': (1658.952, 1.95769, 1.26),
 'b_5.5382592': (2458654.53479-OFF, 5.5382592, 2.182),
 'c_10.688662': (2458659.3486-OFF, 10.688662, 2.854),
}
def read_time(path):
    if path.endswith('.dat'):
        a = np.loadtxt(path, comments='#')
        return a[:,0], a[:,1], a[:,2]
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        i = hdr.index('BJD_TDB')
        t = []
        for line in fh:
            p = line.rstrip('\n').split('\t')
            try: t.append(float(p[i]))
            except: pass
    t = np.array(t); return t, None, None
out = []
files = sorted(glob.glob(os.path.join(d, '*.dat'))) + sorted(glob.glob(os.path.join(d, '*measurements.tbl'))) + sorted(glob.glob(os.path.join(d,'*measurements_5px.tbl')))
for f in files:
    t, fl, er = read_time(f)
    if t.size == 0: continue
    if t[0] > 2400000: tb = t - OFF
    else: tb = t
    rec = {'file': os.path.basename(f), 'n': int(t.size), 'bjd_first': float(t.min()), 'bjd_last': float(t.max()),
           'btjd_first': round(float(tb.min()),5), 'btjd_last': round(float(tb.max()),5), 'span_h': round(float((tb.max()-tb.min())*24),3), 'events': {}}
    for k,(t0,P,dur) in eph.items():
        mid = 0.5*(tb.min()+tb.max())
        n = round((mid - t0)/P)
        evs = []
        for nn in (n-1, n, n+1):
            tc = t0 + nn*P
            h = dur/48.
            # fraction of transit [tc-h, tc+h] inside window
            lo, hi = max(tc-h, tb.min()), min(tc+h, tb.max())
            frac = max(0., hi-lo)/(2*h)
            dist_h = 0.0 if tb.min() <= tc <= tb.max() else min(abs(tc-tb.min()), abs(tc-tb.max()))*24
            if frac > 0 or dist_h < 12:
                evs.append({'epoch': int(nn), 'tc_btjd': round(tc,4), 'frac_covered': round(frac,3), 'tc_outside_window_by_h': round(dist_h,2)})
        if evs: rec['events'][k] = evs
    if fl is not None:
        rec['rms_ppt'] = round(float(np.std(fl)*1e3),3)
        # rough in-transit test for A and B nominal events (linear baseline, mask b/c transits)
        for k in ('A_18.9083','B_1.95769'):
            for ev in rec['events'].get(k, []):
                if ev['frac_covered'] <= 0: continue
                t0,P,dur = eph[k]; tc = ev['tc_btjd']; h = dur/48.
                intr = np.abs(tb-tc) < h
                mask_known = np.zeros_like(intr)
                for kk in ('b_5.5382592','c_10.688662'):
                    t0k,Pk,dk = eph[kk]
                    nk = np.round((tb-t0k)/Pk); mask_known |= np.abs(tb-(t0k+nk*Pk)) < (dk/48.+0.03)
                oot = (~intr) & (~mask_known)
                if intr.sum() < 5 or oot.sum() < 10: continue
                c = np.polyfit(tb[oot], fl[oot], 1)
                res = fl - np.polyval(c, tb)
                dep = -np.mean(res[intr & ~mask_known]) if (intr & ~mask_known).sum()>=5 else float('nan')
                err = np.std(res[oot])/np.sqrt(max(1,(intr & ~mask_known).sum()))
                ev['rough_depth_ppt_linear_baseline'] = round(float(dep*1e3),3)
                ev['rough_err_ppt_white'] = round(float(err*1e3),3)
                ev['n_in'] = int((intr & ~mask_known).sum())
    out.append(rec)
json.dump(out, open(sys.argv[2],'w'), indent=1)
for r in out:
    print(r['file'], r['n'], r['btjd_first'], r['btjd_last'], r['span_h'], r.get('rms_ppt'))
    for k,v in r['events'].items(): print('   ', k, v)
