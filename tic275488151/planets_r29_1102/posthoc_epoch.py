# POST-HOC, descriptive only (written after RESULT.json). Wider net: EBs within 300", P within 1% at harmonics 1/2,1,2; epoch offset of the CTOI vs EB ephemeris.
import json, csv, math, numpy as np
exec(open('match.py').read().split("eb = [r")[0])  # T, vec, sep_arcsec, plan, co
eb = [r for r in csv.reader(l for l in open('sources/villanova_ebs.csv') if not l.startswith('#')) if r and r[0].strip()]
kic=np.array([int(r[0]) for r in eb]); per=np.array([float(r[1]) for r in eb]); t0=np.array([float(r[3]) for r in eb])+2400000
gl=np.array([float(r[6]) for r in eb]); gb=np.array([float(r[7]) for r in eb]); km=np.array([float(r[8]) for r in eb]); pos=gl!=-1
E=vec(gl[pos],gb[pos])@T
src={str(t['tic']):t for t in json.load(open('../planets_r28_1102/R27_PLAN.json'))['targets']}
rows=[]
for t in plan['targets']:
    c=co[str(t['tic'])]; P=t['period_days']; s=sep_arcsec(vec(np.array([c['ra']]),np.array([c['dec']])),E)[0]
    tc=src[str(t['tic'])]['epoch_bjd']
    for i in np.nonzero(s<=300)[0]:
        j=np.nonzero(pos)[0][i]
        for h in (0.5,1,2):
            r=per[j]/(h*P)
            if abs(r-1)<=0.01:
                Pe=per[j]; n=(tc-t0[j])/Pe; d=(n-round(n))*Pe
                rows.append(dict(ctoi=t['ctoi'],eb=int(kic[j]),sep=round(float(s[i]),1),kp_eb=float(km[j]),kp_t=src[str(t['tic'])]['kepmag'],P_ctoi=P,P_eb=float(Pe),h=h,ratio_minus_1=round(float(r-1),5),epoch_offset_d=round(float(d),4),phase_offset=round(float(d/Pe),4)))
for r in rows: print(r)
json.dump(rows,open('POSTHOC_epoch.json','w'),indent=1)
