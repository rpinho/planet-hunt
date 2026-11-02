"""Independent saved-output audit. Imports no measurement, model or detrend code.

Decimal arithmetic for every gate/bar and joint decision; independent centered
regression for D and nuisance fits; source/injection and frozen-byte checks.
Numerical tolerances are stated, not an interval proof of physical assumptions.
"""
from pathlib import Path
from decimal import Decimal,localcontext
import json,hashlib,datetime,subprocess,os,ast
import numpy as np
from scipy.interpolate import BSpline

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dec(x):return Decimal(str(x))
def close(a,b,rtol=2e-8,atol=2e-10,label=''):
    assert np.allclose(np.asarray(a,float),np.asarray(b,float),rtol=rtol,atol=atol),label or (a,b)
def git(*args):return subprocess.check_output(['git',*args],text=True).strip()
def regress(x,y,w):
    sw=w.sum();xm=(w*x).sum()/sw;ym=(w*y).sum()/sw
    slope=np.dot(w*(x-xm),y-ym)/np.dot(w,(x-xm)**2)
    return np.array([ym-slope*xm,slope])

m=read('DESIGN.json');d=np.load('DESIGN.npz');receipt=read('outputs/first_numeric_lc1.json')
inj=read('outputs/injections.json');real=read('outputs/real.json')
raw=np.load('outputs/profile_certificates.npz');dr=np.load('outputs/detrend_certificates.npz')
syn=read('outputs/synthetic.json');sr=np.load('outputs/synthetic_certificates.npz')
PE=3.0305267;TE=2454967.124438-2457000;PC=3.02553;TC=2458682.5504-2457000;H=.105
assert receipt['prefreeze_raw_lc1_display'] is True and receipt['blinded_freeze'] is False
assert m['prefreeze_raw_lc1_display'] is True
assert receipt['nice']==10 and set(receipt['threads'].values())=={'1'}
for p,h in m['input_hashes'].items():assert sha(p)==h,p
for p,h in receipt['hashes'].items():assert sha(p)==h,p
for p,h in read('outputs/completed.json')['output_hashes'].items():assert sha(p)==h,p
assert sha('outputs/synthetic.json')==receipt['synthetic_sha256']
assert sha('DESIGN.npz')==m['design_sha256']
freeze=Path('freeze_commit.txt').read_text().strip();seal=Path('design_commit.txt').read_text().strip()
assert receipt['freeze_commit']==m['freeze_commit']==freeze and receipt['design_commit']==seal
assert git('rev-parse',seal+'^')==freeze
assert subprocess.check_output(['git','show',freeze+':PROTOCOL.md'])==Path('PROTOCOL.md').read_bytes()
assert subprocess.check_output(['git','show',freeze+':STATUS.md']).splitlines()[0]==b'FROZEN'
for p,h in receipt['hashes'].items():
    assert hashlib.sha256(subprocess.check_output(['git','show',seal+':'+p])).hexdigest()==h,p
assert datetime.datetime.fromisoformat(git('show','-s','--format=%cI',seal))<=datetime.datetime.fromisoformat(receipt['utc'])
assert datetime.datetime.fromisoformat(syn['utc'])<=datetime.datetime.fromisoformat(receipt['utc'])
assert sha('r37_measurement.py')==sha('../planets_r37_1102/measurement.py')
assert sha('legacy_r36.py')==sha('../planets_r37_1102/legacy_r36.py')
# Assert the verifier has no measurement/detrending module imports.
for node in ast.walk(ast.parse(Path(__file__).read_text())):
    if isinstance(node,ast.ImportFrom):assert node.module not in ['measurement','r37_measurement','legacy_r36','detrend','measure','synthetic','prepare']
    if isinstance(node,ast.Import):
        for alias in node.names:assert alias.name not in ['measurement','r37_measurement','legacy_r36','detrend','measure','synthetic','prepare']

t=d['time'];sectors=d['sector'];shifts=d['shifts'];indices=np.arange(len(t))
old=np.load('../planets_r37_1102/DESIGN.npz');oldmeta=read('../planets_r37_1102/DESIGN.json')
for k in ['time','sector','shifts','phase','flux','integral_x','integral_h']:assert np.array_equal(d[k],old[k]),k
assert np.array_equal(d['inherited_clean'],old['clean'])
assert np.array_equal(shifts,read('../planets_r36_1102/outputs/supplement_phases.json')['phases']['E'])
assert len(shifts)==200 and m['denominator']==200 and m['threshold']==0
kt=np.load('../planets_r36_1102/outputs/kepler_template.npz')
assert np.array_equal(kt['phase'],d['phase']) and np.array_equal(kt['flux'],d['flux'])
for k in ['offsets','template_mean_window','predicted_primary','a0']:assert m[k]==oldmeta[k]
eb=read('../planets_r36_1102/outputs/EB_TESS.json')
for s in eb['sectors']:assert m['offsets'][str(s['sector'])]==s['timing']['best_offset_days']

def dist(tt,p,ep):return np.remainder(tt-ep+p/2,p)-p/2

def masks(tt,ss,shift=0,margin=True,pad=0):
    offset=np.array([m['offsets'][str(int(z))] for z in ss])
    extra=.05*PE if margin else pad
    return ((abs(dist(tt,PE,TE+offset+shift))<=PE*.153/2+extra)|
        (abs(dist(tt,PE,TE+offset+shift+PE*.5007))<=PE*.1384/2+extra)|
        (abs(dist(tt,PC,TC+shift))<=(.210 if margin else H+pad)))

pad=np.median(np.diff(t))/2
inherited=((abs(dist(t,PE,TE))>PE*.153/2+pad)&
    (abs(dist(t,PE,TE+PE*.5007))>PE*.1384/2+pad)&(abs(dist(t,PC,TC))>H+pad))
assert np.array_equal(inherited,d['inherited_clean'])
assert np.array_equal(d['clean'],inherited&~masks(t,sectors,margin=False,pad=pad))
assert m['oot_time_count']==int(d['clean'].sum()) and m['total_time_count']==len(t)
# Regenerate the inherited schedule from timestamps alone, including the R
# stream consumed before the E stream used by this job.
oldph=read('../planets_r36_1102/outputs/supplement_phases.json')
rng=np.random.default_rng(2026110237);phase_attempts={};pool=t[inherited]
for name in ['R','E']:
    chosen=[];attempt=0
    while len(chosen)<200:
        sh=float(rng.uniform(0,PE));attempt+=1;ok=False
        if name=='R':
            ep=TE+sh
            centers=ep+np.arange(np.ceil((pool.min()-ep)/PE),np.floor((pool.max()-ep)/PE)+1)*PE
            for center in centers:
                dt=pool-center
                if ((abs(dt)<=H).sum()>=3 and ((dt>=-.504)&(dt<=-.252)).sum()>=3 and ((dt>=.252)&(dt<=.504)).sum()>=3):
                    ok=True;break
        else:
            ref=TE+round((TC-TE)/PE)*PE
            for k in range(-3,23):
                a=ref+k*PE+sh;b=TC+k*PC+sh;mid=(a+b)/2
                v=pool[abs(pool-mid)<=.60];dt=v-mid
                be=abs(v-a)<=H;bc=abs(v-b)<=H;out=~(be|bc)
                if len(v)>=12 and be.sum()>=3 and bc.sum()>=3 and (out&(dt<0)).sum()>=3 and (out&(dt>0)).sum()>=3:
                    ok=True;break
        if ok:chosen.append(sh)
    phase_attempts[name]=attempt
    assert np.array_equal(chosen,oldph['phases'][name])
assert phase_attempts==oldph['attempts']
assert datetime.datetime.fromisoformat(inj['utc'])<=datetime.datetime.fromisoformat(real['utc'])
with localcontext() as ctx:
    ctx.prec=60
    xx=d['integral_x'];hh=d['integral_h']
    integ=sum((dec(xx[j+1])-dec(xx[j]))*(dec(hh[j+1])+dec(hh[j]))/2 for j in range(len(xx)-1))/(2*dec(H))
    close(integ,m['template_mean_window'],atol=1e-15)
    close(dec(m['predicted_primary'])/integ,m['a0'],atol=1e-15)

lc=np.loadtxt('../planets_r36_1102/sources/diamante.txt',usecols=[0,3,4,7,10,11,14,17,18,21])
assert np.array_equal(lc[:,0],t)
f=lc[:,[1,4,7]];error=lc[:,[2,5,8]];q=lc[:,[3,6,9]]
good=(q==0)&np.isfinite(f)&np.isfinite(error)&(f>0)&(error>0)
accept=read('outputs/acceptance.json')
assert accept['accepted']==good.sum(axis=0).tolist()
assert accept['oot_accepted']==(good&d['clean'][:,None]).sum(axis=0).tolist()
for sec in [14,15]:assert accept['per_sector'][str(sec)]==good[sectors==sec].sum(axis=0).tolist()

def shape(tt,ss,shift=0):
    offsets=np.array([m['offsets'][str(int(z))] for z in ss])
    ph=((tt-TE-offsets-shift)/PE)%1
    eb=1-np.interp(ph,kt['phase'],kt['flux'],period=1)
    box=(abs(dist(tt,PC,TC+shift))<=H).astype(float)
    return np.column_stack([eb,box])

nD=0;Drows=0;max_D_error=0.
def check_D(bundle,key,ids,fin,ein,shift=0):
    global nD,Drows,max_D_error
    z=bundle[key];assert z.shape==(len(ids),8)
    assert np.array_equal(z[:,0].astype(int),ids)
    close(z[:,1],fin,atol=2e-15);close(z[:,2],ein,atol=2e-15)
    times=t[ids];ss=sectors[ids]
    train=~masks(times,ss,shift)
    assert np.array_equal(z[:,6].astype(bool),train)
    trend=np.zeros(len(ids));robust=np.zeros(len(ids))
    sigma=max(1.4,3*.153*PE);assert sigma==m['smoother_sigma_days']
    for sec in np.unique(ss):
        loc=np.flatnonzero(ss==sec);it=loc[train[loc]];assert len(it)>=3
        grid=np.unique(np.r_[np.arange(times[loc].min(),times[loc].max(),.25),times[loc].max()])
        yy=fin[it];rw=np.ones(len(it))
        # Independent centered regression, with one vector per grid node.
        for stage in range(4):
            x=times[it][None,:]-grid[:,None]
            w=np.exp(-x*x/(2*sigma*sigma))*rw[None,:]/ein[it][None,:]**2
            s0=w.sum(axis=1);xm=(w*x).sum(axis=1)/s0;ym=(w*yy[None,:]).sum(axis=1)/s0
            dx=x-xm[:,None]
            slope=(w*dx*(yy[None,:]-ym[:,None])).sum(axis=1)/(w*dx*dx).sum(axis=1)
            intercept=ym-slope*xm
            if stage<3:
                resid=yy-np.interp(times[it],grid,intercept)
                median=np.median(resid);scale=1.4826*np.median(abs(resid-median))
                if scale>0:
                    u=(resid-median)/(4.685*scale)
                    rw=np.zeros(len(it));ok=abs(u)<1;rw[ok]=(1-u[ok]**2)**2
                else:rw=np.ones(len(it))
        saved=bundle[f'{key}_s{int(sec)}_grid']
        close(saved[:,0],grid,atol=1e-12)
        close(saved[:,1],intercept,atol=5e-12)
        close(saved[:,2],slope,atol=5e-12)
        trend[loc]=np.interp(times[loc],grid,intercept);robust[it]=rw
    err=float(np.max(abs(z[:,3]-trend)));max_D_error=max(max_D_error,err)
    close(z[:,3],trend,rtol=1e-11,atol=5e-12,label='D trend '+key)
    close(z[:,7],robust,rtol=2e-7,atol=2e-9,label='D robust weights '+key)
    close(z[:,4],fin/trend,atol=5e-12);close(z[:,5],ein/trend,atol=5e-15)
    nD+=1;Drows+=len(ids)
    return z[:,4],z[:,5]

nprofiles=0;nevents=0;nrows=0

def check_profile(p,rows,ids,fin,ein,shift=0):
    global nprofiles,nevents,nrows
    assert rows.shape==(p['n_rows'],9)
    ref=TE+round((TC-TE)/PE)*PE;expected=[]
    for k in range(-3,23):
        a=ref+k*PE+shift;b=TC+k*PC+shift;mid=(a+b)/2
        loc=np.flatnonzero(abs(t[ids]-mid)<=.60);ii=ids[loc];v=t[ii];x=v-mid
        be=abs(v-a)<=H;bc=abs(v-b)<=H;out=~(be|bc)
        if len(ii)<12 or be.sum()<3 or bc.sum()<3 or (out&(x<0)).sum()<3 or (out&(x>0)).sum()<3:continue
        if np.median(fin[loc][out])<=0:continue
        expected.append((a,b,loc,ii,be,bc,out))
    assert len(expected)==p['n_events']==len(p['events'])
    sums=[[Decimal(0)]*3 for _ in range(2)]
    for ev,(a,b,loc,ii,be,bc,out) in zip(p['events'],expected):
        z=rows[ev['start']:ev['end']];x=z[:,0];y=z[:,1];w=z[:,2]
        assert np.array_equal(z[:,8].astype(int),ii)
        close([ev['eb_center'],ev['ctoi_center']],[a,b],atol=1e-10)
        close(x,t[ii]-(a+b)/2,atol=1e-10)
        assert [ev['n'],ev['in_eb'],ev['in_ctoi'],ev['n_out']]==[len(ii),int(be.sum()),int(bc.sum()),int(out.sum())]
        h=shape(t[ii],sectors[ii],shift);close(z[:,3:5],h,atol=3e-11)
        norm=np.median(fin[loc][out]);beta=regress(x[out],fin[loc][out],1/ein[loc][out]**2)
        resid=fin[loc][out]-beta[0]-beta[1]*x[out]
        scale=max(1,float(np.sum((resid/ein[loc][out])**2)/(out.sum()-2)))
        close(norm,ev['norm'],atol=1e-12);close(beta,ev['baseline'],atol=1e-11)
        close(scale,ev['scale'],rtol=2e-8)
        close(y,fin[loc]/norm,atol=2e-12)
        close(w,(norm/ein[loc])**2/scale,rtol=2e-8,atol=1e-5)
        for source,res in [(y,z[:,5]),(h[:,0],z[:,6]),(h[:,1],z[:,7])]:
            beta=regress(x,source,w)
            close(res,source-beta[0]-beta[1]*x,atol=5e-12)
        with localcontext() as ctx:
            ctx.prec=60
            wd=list(map(dec,w));yd=list(map(dec,z[:,5]));yy=sum(ww*v*v for ww,v in zip(wd,yd));local=[]
            for j in range(2):
                hd=list(map(dec,z[:,6+j]));hh=sum(ww*v*v for ww,v in zip(wd,hd));hy=sum(ww*v*y0 for ww,v,y0 in zip(wd,hd,yd))
                terms=[hh,hy,yy];local.append(terms)
                for k,v in enumerate(terms):sums[j][k]+=v
        close(local,ev['equations'],rtol=2e-10,atol=1e-8)
        nevents+=1;nrows+=len(ii)
    answers=[]
    with localcontext() as ctx:
        ctx.prec=60
        for terms,ans in zip(sums,p['fits']):
            hh,hy,yy=terms;amp=max(Decimal(0),-hy/hh) if hh>0 else Decimal(0)
            gain=hy*hy/hh if hh>0 and hy<0 else Decimal(0);chi=yy-gain
            close(terms,[ans['hh'],ans['hy'],ans['yy']],rtol=2e-10,atol=5e-8)
            close([amp,gain,chi],[ans['amplitude'],ans['gain'],ans['chi2']],rtol=2e-10,atol=5e-8)
            answers.append(dict(amplitude=amp,gain=gain,chi2=chi))
        delta=answers[0]['gain']-answers[1]['gain']
    close(delta,p['delta_chi2'],atol=5e-8)
    assert p['valid']==bool(expected and min(a[0] for a in sums)>0)
    nprofiles+=1
    return delta,answers

def check_joint(j,rows):
    # Explicit independent Decimal 2x2 inverse and sufficient products.
    with localcontext() as ctx:
        ctx.prec=60
        w=list(map(dec,rows[:,2]));y=list(map(dec,rows[:,5]));h0=list(map(dec,rows[:,6]));h1=list(map(dec,rows[:,7]))
        g00=sum(ww*h*h for ww,h in zip(w,h0));g11=sum(ww*h*h for ww,h in zip(w,h1));g01=sum(ww*a*b for ww,a,b in zip(w,h0,h1))
        b0=-sum(ww*h*yy for ww,h,yy in zip(w,h0,y));b1=-sum(ww*h*yy for ww,h,yy in zip(w,h1,y));yy=sum(ww*v*v for ww,v in zip(w,y))
        det=g00*g11-g01*g01
        assert j['valid']==(det>0 and len(rows)>0)
        if not j['valid']:return
        c00=g11/det;c11=g00/det;c01=-g01/det
        a0=c00*b0+c01*b1;a1=c01*b0+c11*b1;s0=c00.sqrt();s1=c11.sqrt()
        chi=yy-a0*b0-a1*b1
        close(j['G'],[[g00,g01],[g01,g11]],rtol=2e-10,atol=1e-7)
        close(j['b'],[b0,b1],rtol=2e-10,atol=1e-8);close(j['yy'],yy,atol=1e-8)
        close(j['amplitudes'],[a0,a1],atol=1e-14)
        close(j['covariance'],[[c00,c01],[c01,c11]],atol=1e-18)
        close(j['sigma'],[s0,s1],atol=1e-14);close(j['z'],[a0/s0,a1/s1],atol=1e-10)
        close(j['correlation'],c01/(s0*s1),atol=1e-12);close(j['chi2'],chi,atol=5e-8)
        candidates=[(Decimal(0),Decimal(0)),(max(Decimal(0),b0/g00),Decimal(0)),(Decimal(0),max(Decimal(0),b1/g11))]
        if a0>=0 and a1>=0:candidates.append((a0,a1))
        def objective(v):
            a,b=v;return yy-2*(a*b0+b*b1)+a*a*g00+2*a*b*g01+b*b*g11
        winner=min(candidates,key=objective)
        close(j['nonnegative_amplitudes'],winner,atol=1e-14);close(j['nonnegative_chi2'],objective(winner),atol=5e-8)

counts=dict(EB_correct=0,CTOI_correct=0,EB_false=0,CTOI_false=0);actual={};valid={};ratios=[]
mask=good[:,1]&d['clean'];ids=indices[mask]
for name in ['EB','CTOI']:
    assert len(inj['trials'][name])==200;actual[name]=0;valid[name]=0
    typ=0 if name=='EB' else 1;amp=m['a0'] if name=='EB' else .000558
    for k,p in enumerate(inj['trials'][name]):
        assert p['trial']==k and p['shift']==shifts[k]
        sh=p['shift'];h=shape(t[ids],sectors[ids],sh)
        ni=int((abs(amp*h[:,typ])>1e-12).sum());assert ni==p['injected_cadences'];actual[name]+=ni>0
        key=f'inj_{name}_{k:03d}'
        # No D failures in this saved run; assert rather than silently skip.
        assert p['D_valid'] is True
        ff,ee=check_D(dr,key,ids,f[ids,1]*(1-amp*h[:,typ]),error[ids,1],sh)
        delta,ans=check_profile(p,raw[key],ids,ff,ee,sh)
        if p['valid']:
            valid[name]+=1
            if delta>0:counts['EB_correct' if name=='EB' else 'EB_false']+=1
            if delta<0:counts['CTOI_false' if name=='EB' else 'CTOI_correct']+=1
            if name=='EB':ratios.append(ans[0]['amplitude']/dec(m['a0']))
    print('verified injection type',name,flush=True)
gate=(min(actual.values())>=200 and counts['EB_correct']>=160 and counts['CTOI_correct']>=160 and counts['EB_false']<=10 and counts['CTOI_false']<=10)
assert counts==inj['gate']['counts'] and actual==inj['gate']['actual_injections'] and valid==inj['gate']['valid_fits']
assert gate==inj['gate']['passed'] and inj['gate']['threshold']==0 and inj['gate']['denominator']==200
with localcontext() as ctx:
    ctx.prec=60
    mean=sum(ratios)/len(ratios);sigma=(sum((x-mean)**2 for x in ratios)/(len(ratios)-1)).sqrt()
close([mean,mean-1,sigma],[inj['amplitude_scatter'][k] for k in ['mean','bias','sigma']],atol=1e-12)
assert len(ratios)==inj['amplitude_scatter']['n']
results={}
for radius,profiles in real['profiles'].items():
    ri=int(radius)-1;ids=indices[good[:,ri]]
    ff,ee=check_D(dr,'real_r'+radius,ids,f[ids,ri],error[ids,ri]);results[radius]={}
    for sec,p in profiles.items():
        sub=np.ones(len(ids),bool) if sec=='combined' else sectors[ids]==int(sec)
        key=f'real_r{radius}_{sec}'
        results[radius][sec]=check_profile(p,raw[key],ids[sub],ff[sub],ee[sub])
        check_joint(real['joints'][radius][sec],raw[key])
delta,ans=results['2']['combined'];ratio=ans[0]['amplitude']/dec(m['a0'])
fpass=abs(ratio-1)<=3*sigma;pref='EB' if delta>0 else 'CTOI' if delta<0 else 'tie'
bars=real['bars']
assert bars['T2']['status']==('SCORED' if gate else 'UNSCORED') and bars['T2']['preference']==pref
close(bars['T2']['delta_chi2'],delta,atol=5e-8)
assert bars['F2']['status']==('PASS' if fpass else 'FAIL')
close([ratio,sigma,abs(ratio-1)/sigma],[bars['F2'][k] for k in ['ratio','sigma','residual_sigma']],atol=1e-12)
assert bars['attribution']==bool(gate and delta>0 and fpass)

# Synthetic check: reconstruct no-noise EB series and both normalization curves.
sf=1-m['a0']*shape(t,sectors)[:,0];se=np.full(len(t),.0005)
ff,ee=check_D(sr,'D',indices,sf,se)
sp=np.zeros(len(t))
for sec in np.unique(sectors):
    sel=sectors==sec;tt=t[sel];knots=np.arange(tt.min()+.5,tt.max(),.5)
    close(knots,syn['spline']['knots'][str(int(sec))],atol=1e-12)
    fullknots=np.r_[[tt.min()]*4,knots,[tt.max()]*4]
    A=BSpline.design_matrix(tt,fullknots,3).toarray()
    coeff=np.linalg.lstsq(A/se[sel,None],sf[sel]/se[sel],rcond=None)[0]
    sp[sel]=A@coeff
for name,flux0,err0,trend in [('none',sf,se,np.ones(len(t))),('D',ff,ee,sr['D'][:,3]),('spline',sf/sp,se/sp,sp)]:
    saved=sr['series_'+name]
    close(saved[:,1],sf,atol=1e-14);close(saved[:,2],se,atol=1e-15)
    close(saved[:,3],trend,atol=1e-11);close(saved[:,4],flux0,atol=1e-11);close(saved[:,5],err0,atol=1e-14)
    p=syn['models'][name]['profile'];ds,aa=check_profile(p,sr['fit_'+name],indices,saved[:,4],saved[:,5])
    close(syn['models'][name]['amplitude_ratio'],aa[0]['amplitude']/dec(m['a0']),atol=1e-12)
    close(syn['models'][name]['delta_change'],ds-dec(syn['models']['none']['profile']['delta_chi2']),atol=5e-8)

base=read('outputs/baseline.json');br=np.load('outputs/baseline_certificate.npz')
x=np.loadtxt('../planets_r36_1102/sources/diamante.txt',usecols=[0,12,13,14])
g=(x[:,3]==0)&np.isfinite(x[:,1:3]).all(axis=1)&(x[:,1:3]>0).all(axis=1)
db,_=check_profile(base['reproduced'],br['baseline'],indices[g],x[g,1],x[g,2])
close(db,read('../planets_r37_1102/outputs/real.json')['profiles']['2']['combined']['delta_chi2'],atol=1e-9)
close(base['absolute_difference'],db-dec(base['previous_delta']),atol=1e-9)
size=sum(p.stat().st_size for p in Path('.').rglob('*') if p.is_file())
assert size<50_000_000,size
if Path('ARTIFACTS.sha256.json').exists():
    for p,h in read('ARTIFACTS.sha256.json').items():assert sha(p)==h,p
out=dict(status='VERIFIED',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),pid=os.getpid(),
    verifier_sha256=sha('verify.py'),
    artifact_seal_checked=Path('ARTIFACTS.sha256.json').exists(),
    independent_of_measurement_imports=True,decimal_precision=60,prefreeze_blinding_requirement_met=False,
    deliberate_parse_after_freeze=True,profiles_checked=nprofiles,events_checked=nevents,fit_rows_checked=nrows,
    D_series_checked=nD,D_rows_checked=Drows,max_independent_D_trend_absolute_error=max_D_error,
    inherited_phase_generation_attempts=phase_attempts,
    retained_bytes=size,gate=dict(passed=gate,counts=counts,actual_injections=actual,valid_fits=valid,denominator=200),
    bars=dict(T2_status=bars['T2']['status'],preference=pref,delta_chi2=float(delta),F2_pass=bool(fpass),
        amplitude_ratio=float(ratio),injection_sigma=float(sigma),attribution=bars['attribution']),
    scope='Saved numerical arithmetic, source and pre-D injection reconstruction, independent D and nuisance projection, joint algebra, synthetic spline check, hashes and post-display freeze chronology. Not a blinded experiment or a validation of physical forward assumptions.')
Path('outputs/verification.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps(out,indent=2))
