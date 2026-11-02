"""Independent saved-output verifier. Imports no measurement or extraction code.

Rebuilds timestamp support and checks raw fit inputs and nuisance projections;
uses Decimal for sufficient-statistic arithmetic, scatter, gate and bar decisions.
No network. This is a numerical certificate, not an astrophysical model proof.
"""
from pathlib import Path
import json,hashlib,datetime,subprocess,os
from decimal import Decimal,localcontext
import numpy as np

def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dec(x): return Decimal(str(x))
def close(a,b,rtol=2e-9,atol=2e-10):
    assert np.allclose(np.asarray(a,float),np.asarray(b,float),rtol=rtol,atol=atol),(a,b)
def git(*args): return subprocess.check_output(['git',*args],text=True).strip()
def regression(xx,yy,ww):
    # Closed-form centered weighted regression, separate from measurement's solve.
    xm=np.dot(ww,xx)/ww.sum();ym=np.dot(ww,yy)/ww.sum()
    slope=np.dot(ww*(xx-xm),yy-ym)/np.dot(ww,(xx-xm)**2)
    return ym-slope*xm,slope

meta=read('DESIGN.json');d=np.load('DESIGN.npz');first=read('outputs/first_flux.json')
inj=read('outputs/injections.json');real=read('outputs/real.json');raw=np.load('outputs/profile_certificates.npz')
PE=3.0305267;TE=2454967.124438-2457000;PC=3.02553;TC=2458682.5504-2457000;H=.105
for p,h in meta['input_hashes'].items(): assert sha(p)==h,p
for p,h in first['hashes'].items(): assert sha(p)==h,p
for p,h in read('outputs/completed.json')['output_hashes'].items(): assert sha(p)==h,p
assert sha('DESIGN.npz')==meta['design_sha256']
assert sha('legacy_r36.py')==meta['input_hashes']['../planets_r36_1102/measurement.py']
freeze=Path('freeze_commit.txt').read_text().strip();seal=Path('design_commit.txt').read_text().strip()
assert first['freeze_commit']==meta['freeze_commit']==freeze and first['design_commit']==seal
assert git('rev-parse',seal+'^')==freeze
assert hashlib.sha256(subprocess.check_output(['git','show',freeze+':PROTOCOL.md'])).hexdigest()==sha('PROTOCOL.md')
assert subprocess.check_output(['git','show',freeze+':STATUS.md']).splitlines()[0]==b'FROZEN'
for p,h in first['hashes'].items():
    assert hashlib.sha256(subprocess.check_output(['git','show',seal+':'+p])).hexdigest()==h,p
assert datetime.datetime.fromisoformat(git('show','-s','--format=%cI',seal))<=datetime.datetime.fromisoformat(first['utc'])
assert datetime.datetime.fromisoformat(meta['utc'])<=datetime.datetime.fromisoformat(first['utc'])
assert first['nice']==10 and set(first['threads'].values())=={'1'}

t=d['time'];sector=d['sector'];clean=d['clean'];shifts=d['shifts']
oldtime=np.load('../planets_r36_1102/outputs/time_design.npz')
assert np.array_equal(t,oldtime['time']) and np.array_equal(clean,oldtime['clean'])
def dist(times,p,ep): return (times-ep+p/2)%p-p/2
pad=np.median(np.diff(t))/2
clean2=((abs(dist(t,PE,TE))>.1530*PE/2+pad)&
        (abs(dist(t,PE,TE+.5007*PE))>.1384*PE/2+pad)&(abs(dist(t,PC,TC))>H+pad))
assert np.array_equal(clean,clean2)
man=read('../planets_r30_1102/aperture_manifest.json');ss=np.zeros(len(t),int)
for s in man['sectors']:
    if s['sector'] in [14,15]:
        z=s['header'];ss[(t>=z['TSTART'])&(t<=z['TSTOP'])]=s['sector']
assert np.array_equal(sector,ss)
eb=read('../planets_r36_1102/outputs/EB_TESS.json')
for s in eb['sectors']: assert meta['offsets'][str(s['sector'])]==s['timing']['best_offset_days']
oldph=read('../planets_r36_1102/outputs/supplement_phases.json')
assert np.array_equal(shifts,oldph['phases']['E'])
tt=t[clean];rng=np.random.default_rng(2026110237);attempts={}
for name in ['R','E']:
    accepted=[];n=0
    while len(accepted)<200:
        shift=float(rng.uniform(0,PE));n+=1;ok=False
        if name=='R':
            ep=TE+shift
            centres=ep+np.arange(np.ceil((tt.min()-ep)/PE),np.floor((tt.max()-ep)/PE)+1)*PE
            for c in centres:
                x=tt-c
                if ((abs(x)<=H).sum()>=3 and ((x>=-.504)&(x<=-.252)).sum()>=3
                        and ((x>=.252)&(x<=.504)).sum()>=3): ok=True;break
        else:
            ref=TE+round((TC-TE)/PE)*PE
            for k in range(-3,23):
                a=ref+k*PE+shift;b=TC+k*PC+shift;mid=(a+b)/2
                v=tt[abs(tt-mid)<=.60];x=v-mid
                be=abs(v-a)<=H;bc=abs(v-b)<=H;out=~(be|bc)
                if len(v)>=12 and be.sum()>=3 and bc.sum()>=3 and (out&(x<0)).sum()>=3 and (out&(x>0)).sum()>=3:
                    ok=True;break
        if ok: accepted.append(shift)
    attempts[name]=n
    assert np.array_equal(accepted,oldph['phases'][name])
assert attempts==oldph['attempts']
kt=np.load('../planets_r36_1102/outputs/kepler_template.npz')
assert np.array_equal(d['phase'],kt['phase']) and np.array_equal(d['flux'],kt['flux'])
# Exact integral of the piecewise-linear interpolant, with saved knot checks.
xx=d['integral_x'];hh=d['integral_h']
close(hh,1-np.interp((xx/PE)%1,kt['phase'],kt['flux'],period=1),atol=1e-15)
with localcontext() as ctx:
    ctx.prec=60
    integral=sum((dec(xx[i+1])-dec(xx[i]))*(dec(hh[i+1])+dec(hh[i]))/2 for i in range(len(xx)-1))/(2*dec(H))
    close(integral,meta['template_mean_window'],atol=1e-15)
    close(dec(meta['predicted_primary'])/integral,meta['a0'],atol=1e-15)

lc=np.loadtxt('../planets_r36_1102/sources/diamante.txt',usecols=range(29))
assert np.array_equal(lc[:,0],t)
flux=lc[:,[5,12,19]];error=lc[:,[6,13,20]];flags=lc[:,[7,14,21]]
good=(flags==0)&np.isfinite(flux)&np.isfinite(error)&(flux>0)&(error>0)
ac=read('outputs/acceptance.json');assert ac['accepted']==good.sum(axis=0).tolist()==[2227]*3
assert ac['oot_accepted']==(good&clean[:,None]).sum(axis=0).tolist()
for s in [14,15]: assert ac['per_sector'][str(s)]==good[sector==s].sum(axis=0).tolist()

def shape(times,sectors,shift):
    displacement=np.array([meta['offsets'][str(int(s))] for s in sectors])
    phase=((times-(TE+shift+displacement))/PE)%1
    return np.column_stack([1-np.interp(phase,kt['phase'],kt['flux'],period=1),
                            (abs(dist(times,PC,TC+shift))<=H).astype(float)])

checked_events=0;checked_rows=0;checked_profiles=0
def check_profile(p,key,radius,mask,shift=0,injected=None):
    global checked_events,checked_rows,checked_profiles
    arr=raw[key];assert arr.shape==(p['n_rows'],9)
    # Rebuild every qualifying event using only accepted timestamps.
    ref=TE+round((TC-TE)/PE)*PE;expected=[]
    for k in range(-3,23):
        a=ref+k*PE+shift;b=TC+k*PC+shift;mid=(a+b)/2
        ids=np.flatnonzero(mask&(abs(t-mid)<=.60));v=t[ids];x=v-mid
        be=abs(v-a)<=H;bc=abs(v-b)<=H;out=~(be|bc)
        if len(ids)<12 or be.sum()<3 or bc.sum()<3 or (out&(x<0)).sum()<3 or (out&(x>0)).sum()<3: continue
        expected.append((a,b,ids,be,bc,out))
    assert len(expected)==p['n_events']==len(p['events'])
    sums=[[Decimal(0)]*3 for _ in range(2)]
    for ev,(a,b,ids,be,bc,out) in zip(p['events'],expected):
        z=arr[ev['start']:ev['end']];x,y,w,h1,h2,yr,r1,r2,_=z.T
        assert ev['certificate_key']==key and np.array_equal(z[:,8].astype(int),ids)
        close([ev['eb_center'],ev['ctoi_center']],[a,b],atol=1e-10)
        close(x,t[ids]-(a+b)/2,atol=1e-10)
        assert ev['in_eb']==int(be.sum()) and ev['in_ctoi']==int(bc.sum()) and ev['n_out']==int(out.sum())
        h=shape(t[ids],sector[ids],shift)
        close(z[:,3:5],h,rtol=2e-8,atol=3e-11)
        orig=flux[ids,radius].copy();ee=error[ids,radius]
        if injected is not None:
            j,amp=injected;orig*=1-amp*h[:,j]
        norm=float(np.median(orig[out]));beta=regression(x[out],orig[out],1/ee[out]**2)
        resid=orig[out]-beta[0]-beta[1]*x[out]
        scale=max(1,float(np.sum((resid/ee[out])**2)/(out.sum()-2)))
        close(norm,ev['norm'],atol=1e-12);close(beta,ev['baseline'],atol=1e-11)
        close(scale,ev['scale'],rtol=1e-8)
        close(y,orig/norm,atol=1e-12)
        close(w,(norm/ee)**2/scale,rtol=1e-8,atol=1e-5)
        for source,res in [(y,yr),(h1,r1),(h2,r2)]:
            b0,b1=regression(x,source,w)
            close(res,source-b0-b1*x,rtol=1e-8,atol=2e-12)
        local=[]
        with localcontext() as ctx:
            ctx.prec=60
            wd=[dec(v) for v in w];yd=[dec(v) for v in yr]
            yy=sum(ww*vv*vv for ww,vv in zip(wd,yd))
            for j,r in enumerate([r1,r2]):
                rd=[dec(v) for v in r]
                hh=sum(ww*vv*vv for ww,vv in zip(wd,rd))
                hy=sum(ww*vv*zz for ww,vv,zz in zip(wd,rd,yd))
                terms=[hh,hy,yy];local.append(terms)
                for k,value in enumerate(terms): sums[j][k]+=value
        close(local,ev['equations'],rtol=2e-10,atol=1e-9)
        checked_events+=1;checked_rows+=len(ids)
    answers=[]
    with localcontext() as ctx:
        ctx.prec=60
        for eq,fit in zip(sums,p['fits']):
            hh,hy,yy=eq
            amp=max(Decimal(0),-hy/hh) if hh>0 else Decimal(0)
            gain=hy*hy/hh if hh>0 and hy<0 else Decimal(0)
            chi=yy-gain
            close(eq,[fit['hh'],fit['hy'],fit['yy']],rtol=2e-10,atol=2e-8)
            close([amp,gain,chi],[fit['amplitude'],fit['gain'],fit['chi2']],rtol=2e-10,atol=2e-8)
            answers.append(dict(amplitude=amp,gain=gain,chi2=chi))
        delta=answers[0]['gain']-answers[1]['gain']
    close(delta,p['delta_chi2'],atol=2e-8)
    assert p['valid']==bool(expected and min(q[0] for q in sums)>0)
    checked_profiles+=1
    return delta,answers

counts=dict(EB_correct=0,CTOI_correct=0,EB_false=0,CTOI_false=0)
actual={};valid={};ratios=[]
for name in ['EB','CTOI']:
    assert len(inj['trials'][name])==200;actual[name]=0;valid[name]=0
    j=0 if name=='EB' else 1;amp=meta['a0'] if name=='EB' else .000558
    mask=good[:,1]&clean
    for k,p in enumerate(inj['trials'][name]):
        assert p['trial']==k and p['shift']==shifts[k]
        n=int((abs(amp*shape(t[mask],sector[mask],p['shift'])[:,j])>1e-12).sum())
        assert n==p['injected_cadences'];actual[name]+=n>0
        delta,answer=check_profile(p,f'inj_{name}_{k:03d}',1,mask,p['shift'],(j,amp))
        if p['valid']:
            valid[name]+=1
            if delta>0: counts['EB_correct' if name=='EB' else 'EB_false']+=1
            if delta<0: counts['CTOI_false' if name=='EB' else 'CTOI_correct']+=1
            if name=='EB': ratios.append(answer[0]['amplitude']/dec(meta['a0']))
assert counts==inj['gate']['counts'] and actual==inj['gate']['actual_injections'] and valid==inj['gate']['valid_fits']
gate=(min(actual.values())>=200 and counts['EB_correct']>=160 and counts['CTOI_correct']>=160 and counts['EB_false']<=10 and counts['CTOI_false']<=10)
assert gate==inj['gate']['passed'] and inj['gate']['threshold']==0 and inj['gate']['denominator']==200
with localcontext() as ctx:
    ctx.prec=60
    mean=sum(ratios)/len(ratios)
    sigma=(sum((r-mean)**2 for r in ratios)/(len(ratios)-1)).sqrt()
close([mean,mean-1,sigma],[inj['amplitude_scatter'][k] for k in ['mean','bias','sigma']],atol=1e-12)
assert len(ratios)==inj['amplitude_scatter']['n']
computed_real={}
for radius,profiles in real['profiles'].items():
    computed_real[radius]={}
    for sec,p in profiles.items():
        mask=good[:,int(radius)-1].copy()
        if sec!='combined': mask&=sector==int(sec)
        computed_real[radius][sec]=check_profile(p,f'real_r{radius}_{sec}',int(radius)-1,mask)
delta,answer=computed_real['2']['combined'];ratio=answer[0]['amplitude']/dec(meta['a0'])
pref='EB' if delta>0 else 'CTOI' if delta<0 else 'tie'
fpass=abs(ratio-1)<=3*sigma
bars=real['bars']
assert bars['T']['status']==('SCORED' if gate else 'UNSCORED') and bars['T']['preference']==pref
close(delta,bars['T']['delta_chi2'],atol=2e-8)
assert bars['F']['status']==('PASS' if fpass else 'FAIL')
close([ratio,sigma,abs(ratio-1)/sigma],[bars['F'][k] for k in ['ratio','sigma','residual_sigma']],atol=1e-12)
assert bars['attribution']==bool(gate and pref=='EB' and fpass)

# Independently reconstruct the secondary event estimator and covariance.
all_events=[]
for sec in [14,15]:
    sec_events=[];mask=good[:,1]&(sector==sec);st=t[mask];ep=TE+meta['offsets'][str(sec)]+.5007*PE
    for c in ep+np.arange(np.ceil((st.min()-ep)/PE),np.floor((st.max()-ep)/PE)+1)*PE:
        dt=t-c;inside=np.flatnonzero(mask&(abs(dt)<=H))
        left=mask&(dt>=-.504)&(dt<=-.252);right=mask&(dt>=.252)&(dt<=.504)
        outside=np.flatnonzero(left|right)
        if len(inside)<3 or left.sum()<3 or right.sum()<3: continue
        xi=dt[inside];xo=dt[outside];fo=flux[outside,1];fi=flux[inside,1];eo=error[outside,1];ei=error[inside,1]
        w=1/eo**2;s0=w.sum();s1=np.dot(w,xo);s2=np.dot(w,xo*xo);det=s0*s2-s1*s1
        inv=np.array([[s2,-s1],[-s1,s0]])/det
        mapping=inv@np.vstack([w,w*xo]);beta=mapping@fo;base=beta[0]+beta[1]*xi
        dep=float(np.mean(1-fi/base))
        din=-1/len(inside)/base
        db=np.array([np.mean(fi/base**2),np.mean(fi*xi/base**2)])
        dout=db@mapping
        noise=max(1,float(np.sum(((fo-beta[0]-beta[1]*xo)/eo)**2)/(len(outside)-2)))
        variance=float((np.sum((din*ei)**2)+np.sum((dout*eo)**2))*noise)
        sec_events.append(dict(center=float(c),depth=[dep],covariance=[[variance]],n_in=len(inside),n_left=int(left.sum()),n_right=int(right.sum())))
    saved=real['secondary'][str(sec)]['events'];assert len(saved)==len(sec_events)
    for a,b in zip(saved,sec_events):
        close(a['center'],b['center'],atol=1e-10);close(a['depth'],b['depth'],atol=1e-12)
        close(a['covariance'],b['covariance'],rtol=1e-8,atol=1e-18)
        assert [a[k] for k in ['n_in','n_left','n_right']]==[b[k] for k in ['n_in','n_left','n_right']]
    all_events+=saved
assert real['secondary']['combined']['events']==all_events
secondary={}
for sec,z in real['secondary'].items():
    ev=z['events'];st=z['stack'];assert st['n_events']==len(ev)
    with localcontext() as ctx:
        ctx.prec=60
        w=[1/dec(e['covariance'][0][0]) for e in ev];v=[dec(e['depth'][0]) for e in ev]
        sw=sum(w);mu=sum(a*b for a,b in zip(w,v))/sw
        chi=sum(a*(b-mu)**2 for a,b in zip(w,v));infl=max(Decimal(1),chi/(len(ev)-1))
        var=infl/sw;err=var.sqrt()
    close([mu,var,1/sw,chi,infl,err],[st['depth'][0],st['covariance'][0][0],st['formal_covariance'][0][0],st['chi2'][0],st['inflation'][0],st['error'][0]],atol=1e-15)
    secondary[sec]=dict(depth_ppm=float(mu*1000000),error_ppm=float(err*1000000),events=len(ev),
        prediction_residual_sigma=float((mu-dec(meta['predicted_secondary']))/err))

# Previously seen two-box baseline: independent sufficient-statistic check.
base=read('outputs/baseline.json');prev=read('../planets_r36_1102/outputs/real.json')['E']
old=base['reproduced'];gains=[]
with localcontext() as ctx:
    ctx.prec=60
    for f in old['fits']:
        hh=dec(f['hh']);hy=dec(f['hy']);yy=dec(f['yy']);gain=hy*hy/hh if hh>0 and hy<0 else Decimal(0)
        close([max(Decimal(0),-hy/hh),gain,yy-gain],[f['depth'],f['gain'],f['chi2']],atol=1e-9)
        gains.append(gain)
close(gains[0]-gains[1],prev['delta_chi2'],atol=1e-9)
assert old['events']==prev['events']
for s in read('outputs/synthetic.json'):
    assert abs(s['fitted']-s['normalization_adjusted_expected'])<1e-12
    close(s['absolute_error'],s['fitted']-s['normalization_adjusted_expected'],atol=1e-18)

out=dict(status='VERIFIED',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),pid=os.getpid(),
    independent_of_measurement_imports=True,decimal_precision=60,profiles_checked=checked_profiles,
    events_checked=checked_events,rows_checked=checked_rows,phase_design_attempts=attempts,
    gate=dict(passed=gate,counts=counts,actual_injections=actual,valid_fits=valid,denominator=200),
    bars=dict(T_status=bars['T']['status'],preference=pref,delta_chi2=float(delta),
              F_pass=bool(fpass),amplitude_ratio=float(ratio),injection_sigma=float(sigma),attribution=bars['attribution']),
    secondary=secondary,scope='Saved arithmetic, raw fit inputs, nuisance projection, timestamp design, acceptance, secondary extraction, hashes and freeze chronology; no astrophysical forward-model validation.')
if Path('ARTIFACTS.sha256.json').exists():
    for p,h in read('ARTIFACTS.sha256.json').items(): assert sha(p)==h,p
    out['artifact_seal_checked']=True
Path('outputs/verification.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps(out,indent=2))
