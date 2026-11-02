"""Independent saved-output verifier. Imports NO r31/r30 measurement modules.
Decimal arithmetic reconstructs every normal-equation fit and exact threshold decisions.
Only numpy is used for saved OOT arrays and phase audits; no real images are accessed.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from decimal import Decimal as D,getcontext
from fractions import Fraction as Q
import json,hashlib,subprocess,datetime
import numpy as np
getcontext().prec=70
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
def read(path,decimal=False):return json.loads(Path(path).read_text(),parse_float=D) if decimal else json.loads(Path(path).read_text())
def dec(v):return v if isinstance(v,D) else D(str(v))
def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def close(x,y,name):
    x,y=dec(x),dec(y)
    assert abs(x-y)<=D('2e-8')*max(abs(x),abs(y),D('1e-14')),(name,str(x),str(y))

def recompute_fit(stats,saved):
    h=[[sum(dec(s['H'][i][j]) for s in stats) for j in range(2)] for i in range(2)]
    g=[sum(dec(s['g'][i]) for s in stats) for i in range(2)]
    assert abs(h[0][1]-h[1][0])<D('1e-8')*max(abs(h[0][1]),D(1))
    determinant=h[0][0]*h[1][1]-h[0][1]*h[1][0];assert determinant>0
    inv=[[h[1][1]/determinant,-h[0][1]/determinant],[-h[1][0]/determinant,h[0][0]/determinant]]
    beta=[sum(inv[i][j]*g[j] for j in range(2)) for i in range(2)]
    dof=sum(int(s['n_pixels'])-1 for s in stats)-2
    yy=sum(dec(s['yy']) for s in stats);chi=yy-sum(beta[i]*g[i] for i in range(2));assert chi>=-D('1e-7')*max(yy,D(1))
    chi=max(D(0),chi);infl=max(D(1),chi/dec(dof));cov=[[v*infl for v in row] for row in inv]
    a,r=beta;sa=cov[0][0].sqrt();sr=cov[1][1].sqrt();ul=max(D(0),a)+D('1.6448536269514722')*sa
    values={'a':a,'r':r,'sigma_a':sa,'sigma_r':sr,'chi2':chi,'inflation':infl,'z_a':a/sa,'z_r':r/sr,'upper95_a':ul,'yy':yy}
    assert saved['dof']==dof and saved['n_sectors']==len(stats)
    for key,val in values.items():close(val,saved[key],key)
    for i in range(2):
        close(g[i],saved['normal_g'][i],'g')
        for j in range(2):close(h[i][j],saved['normal_H'][i][j],'H');close(cov[i][j],saved['covariance'][i][j],'covariance')
    # Recompute the threshold-bearing fit as exact rationals of saved decimal inputs.
    # Even the upper limit is decided without evaluating a rounded square root:
    # max(0,a)+Z*sqrt(V)<L iff max(0,a)<L and Z^2*V<(L-max(0,a))^2.
    qh=[[sum(Q(dec(s['H'][i][j])) for s in stats) for j in range(2)] for i in range(2)]
    qg=[sum(Q(dec(s['g'][i])) for s in stats) for i in range(2)]
    qdet=qh[0][0]*qh[1][1]-qh[0][1]*qh[1][0]
    qi=[[qh[1][1]/qdet,-qh[0][1]/qdet],[-qh[1][0]/qdet,qh[0][0]/qdet]]
    qa,qr=[sum(qi[i][j]*qg[j] for j in range(2)) for i in range(2)]
    qchi=max(Q(0),sum(Q(dec(s['yy'])) for s in stats)-qa*qg[0]-qr*qg[1])
    qinfl=max(Q(1),qchi/dof);va=qi[0][0]*qinfl;vr=qi[1][1]*qinfl
    lower=max(Q(0),qa);limit=Q('0.000279');z95=Q('1.6448536269514722')
    decisions={'leak_detect':qr>0 and qr*qr>=25*vr,
               'target_limit':lower<limit and z95*z95*va<(limit-lower)**2,
               'target_detect':qa>0 and qa*qa>=9*va,
               'false_eb':qr>=0 and qr*qr>=25*vr}
    bars={'bar1':decisions['leak_detect'] and (qr-1)**2<=9*vr,
          'bar2':decisions['target_limit'] and (qa<0 or qa*qa<9*va)}
    assert decisions['target_limit']==(ul<D('0.000279'))
    return decisions,bars,values

meta=read('outputs/metadata.json');before=read('outputs/before_oot_flux.json');gate=read('outputs/gate.json',True);design=read('outputs/design.json')
commit=Path('freeze_commit.txt').read_text().strip();assert commit==before['freeze_commit']==gate['freeze_commit']
frozen=subprocess.check_output(['git','show',commit+':PROTOCOL.md']);assert hashlib.sha256(frozen).hexdigest()==digest('PROTOCOL.md')
frozen_meta=subprocess.check_output(['git','show',commit+':outputs/metadata.json']);assert hashlib.sha256(frozen_meta).hexdigest()==digest('outputs/metadata.json')
for p,h in before['sha256'].items():assert digest(p)==h,('implementation changed',p)
assert before['nice']==10 and all(v=='1' for v in before['threads'].values())
assert design['seed']==311102 and design['trials_per_arm']==200 and gate['denominator_each']==200
freeze_time=datetime.datetime.fromisoformat(subprocess.check_output(['git','show','-s','--format=%cI',commit],text=True).strip())
assert freeze_time<=datetime.datetime.fromisoformat(design['utc'])<datetime.datetime.fromisoformat(before['utc'])<datetime.datetime.fromisoformat(gate['utc'])
for r in meta['receipts']:assert digest(Path(meta['source_root'])/r['relative_path'])==r['sha256'],r['relative_path']
# Reproduce the inherited, already-seen Kepler prediction from saved event measurements.
kepler=read(Path(meta['source_root'])/'outputs/kepler.json',True)
events=[e for q in kepler['quarters'] for e in q['measurements']['primary']['events']]
weights=[1/dec(e['variance']) for e in events];ws=sum(weights)
depth=sum(w*dec(e['depth']) for w,e in zip(weights,events))/ws
kchi=sum(w*(dec(e['depth'])-depth)**2 for w,e in zip(weights,events))
ksigma=(max(D(1),kchi/dec(len(events)-1))/ws).sqrt()
assert len(events)==meta['kepler_primary']['n']==419 and len(kepler['quarters'])==17
close(depth,meta['kepler_primary']['depth'],'inherited Kepler depth')
close(ksigma,meta['kepler_primary']['sigma'],'inherited Kepler error')
geom=[]
for row in meta['geometry']:
    margins={}
    for star,xy in row['physical_xy'].items():
        x,y=map(dec,xy);m=min(x-45,2092-x,y-1,2048-y);margins[star]=m;close(m,row['edge_margin_px'][star],'edge margin')
    keep=all(m>2 for m in margins.values());assert keep==row['usable'];geom.append((row['sector'],keep))
assert [s for s,k in geom if k]==meta['usable_sectors']==[14,15,40,54,55,74,75,81]
assert [s for s,k in geom if not k]==meta['dropped_sectors']==[41,82]
assert all(row['illuminated_pixels']==225 for row in meta['geometry'] if row['usable'])

# Explicit independent modulo phase mask; saved arrays contain OOT TIME, never real images.
def outside(t):
    t=np.asarray(t);good=np.ones(t.shape,bool)
    for per,epoch in [(3.0305267,2454967.124438),(3.02553,2458682.5504)]:
        ph=(t-epoch)/per
        good &= abs((ph+.5)%1-.5)>.1265
        good &= abs((ph-.5007+.5)%1-.5)>.1192
    return good
sectors={};phase_audits=[]
for d in design['sectors']:
    sec=d['sector'];z=read(f'outputs/gate_sector_s{sec:04d}.json',True);sectors[sec]=z
    assert z.get('error') is None and len(z['trials'])==200 and all(t['valid'] for t in z['trials'])
    with np.load(f'outputs/oot_summary_s{sec:04d}.npz') as n:
        times=n['t'];assert len(times)==d['oot'] and np.all(outside(times));assert hashlib.sha256(times.astype('<f8').tobytes()).hexdigest()==d['oot_times_sha256']
        v=n['variance'];scene=n['scene'];mu=n['mean'];A=np.column_stack([scene,np.ones(225)]);H=A.T@((1/v)[:,None]*A);rhs=A.T@(mu/v);cal=np.linalg.solve(H,rhs)
        close(cal[0],z['calibration']['F_target_e_per_s'],'Fcal');close(cal[1],z['calibration']['scene_background'],'bgcal')
        F=dec(z['calibration']['F_target_e_per_s']);R=dec(meta['tmag_flux_ratio']);dk=dec(meta['kepler_primary']['depth']);close(F*R,z['calibration']['F_EB_e_per_s'],'Tmag EB flux');close(F*R*dk,z['calibration']['b0_e_per_s'],'b0')
        assert np.all(v>0) and len(v)==225
    for trial,ph in enumerate(d['phase']):
        assert ph is not None
        lo,hi=d['min_bjd'],d['max_bjd'];P=3.0305267;T=2454967.124438
        ks=np.arange(int(np.floor((lo-T)/P))-1,int(np.ceil((hi-T)/P))+2)
        centers=T+(ks+ph)*P;centers=centers[(centers>=lo)&(centers<=hi)];assert len(centers)>0 and np.all(outside(centers))
        tr=z['trials'][trial];assert tr['trial']==trial and dec(tr['phase'])==dec(ph);con=tr['contrast'];assert con['n_events']==len(con['events'])>0
        assert abs(dec(con['contrast_sum']))<D('1e-10') and abs(dec(con['injection_response'])-1)<D('1e-10')
        for ev in con['events']:assert ev['n_in']>=3 and ev['n_out_left']>=3 and ev['n_out_right']>=3
        # Injection shifts are exactly the requested linear components after profiling background.
        L=tr['leak'];U=tr['target'];H=np.array(L['H'],float)
        assert np.allclose(H,np.array(U['H'],float),rtol=1e-12,atol=1e-8)
        delta_g=np.array(L['g'],float)-np.array(U['g'],float);assert np.allclose(delta_g,H@np.array([-558e-6,1]),rtol=1e-8,atol=1e-8)
    phase_audits.append({'sector':sec,'oot_cadences':d['oot'],'valid_trials':200,'all_centers_away_from_both_windows':True})
for check in gate['explicit_injection_checks']:assert dec(check['max_abs_error_e_per_s'])<=D('1e-8')+D('1e-8')*dec(check['scale_e_per_s'])
assert len(gate['explicit_injection_checks'])==16

counts={};bounds={};verified=0
for arm in ['leak','target']:
    counts[arm]={k:0 for k in ['leak_detect','target_limit','target_detect','false_eb','valid']};values=[]
    for j,row in enumerate(gate['results'][arm]):
        assert row['trial']==j and row['fit'] is not None
        stats=[sectors[sec]['trials'][j][arm] for sec in meta['usable_sectors']]
        decisions,bars,v=recompute_fit(stats,row['fit']);assert decisions==row['decisions']
        for k,x in decisions.items():counts[arm][k]+=int(x)
        counts[arm]['valid']+=1;values.append(v);verified+=1
    bounds[arm]={k:{'min':str(min(v[k] for v in values)),'max':str(max(v[k] for v in values))} for k in ['sigma_a','upper95_a','z_a','z_r']}
assert counts==gate['counts']
requirements={'leak_detection':counts['leak']['leak_detect']>=180,'leak_target_upper_limit':counts['leak']['target_limit']>=180,'target_false_eb':counts['target']['false_eb']<=10,'target_detection':counts['target']['target_detect']>=160}
assert requirements==gate['requirements'];passed=all(requirements.values());assert passed==gate['gate_pass']
if passed:
    real=read('outputs/real.json',True)
    for name,data in real['fits'].items():
        stat=[];stat14=[]
        for row in data['sectors']:
            recompute_fit([row['stats']],row['fit']);stat.append(row['stats'])
            if row['sector'] in [14,15]:stat14.append(row['stats'])
        _,bars,_=recompute_fit(stat,data['pooled']);recompute_fit(stat14,data['sectors14_15'])
        if name=='EB':assert real['bar1']==bars['bar1'] and real['bar2']==bars['bar2'];bar_result=bars
    attribution=passed and bar_result['bar1'] and bar_result['bar2'];assert attribution==real['attribution'];status='SCORED'
else:
    assert not Path('outputs/real.json').exists() and not Path('outputs/before_real_flux.json').exists()
    assert not list(Path('outputs').glob('real_*.npz')) and gate['real_difference_images_formed'] is False
    bar_result={'bar1':'UNSCORED','bar2':'UNSCORED'};attribution=False;status='POWER_FAIL'
assert Path('STATUS.md').read_text().splitlines()[0]==status
# Optional final artifact seal is outside its own hash graph.
sealpath=Path('outputs/final_seal.json')
if sealpath.exists():
    for name,expected in read(sealpath)['sha256'].items():assert digest(name)==expected,('final seal mismatch',name)
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'verified':True,'precision_decimal_digits':70,'decision_arithmetic':'exact rational inequalities, including a squared upper-limit comparison','fits_recomputed':verified,'trial_denominator_each':200,'counts':counts,'requirements':requirements,'gate_pass':passed,'bars':bar_result,'attribution':attribution,'status':status,'freeze_commit':commit,'phase_audits':phase_audits,'bounds':bounds,'input_hashes_checked':len(meta['receipts']),'final_seal_checked':sealpath.exists(),'limits':'Saved arithmetic, masks, geometry thresholds, injection equivalence and hashes; not independent raw-pixel noise extraction, PRF truth, or an astrophysical validation.'}
Path('outputs/verification.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:result[k] for k in ['verified','fits_recomputed','counts','requirements','bars','attribution','status']},indent=2))
