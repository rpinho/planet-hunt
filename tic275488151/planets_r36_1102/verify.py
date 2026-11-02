"""Independent saved-output verifier; imports no measurement/preparation code.

Decimal event restacking and decision arithmetic, independent profile normal-equation
evaluation, time-only phase regeneration, git freeze and SHA256 audit. No network.
"""
from pathlib import Path
from decimal import Decimal,localcontext
import json,hashlib,subprocess,datetime,statistics,math
import numpy as np

D=Decimal
ROOT=Path(__file__).resolve().parent
assert Path.cwd()==ROOT
checks=[]
def record(s):checks.append(s)
def load(p):return json.loads(Path(p).read_text(),parse_float=D)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def close(a,b,where,rtol=D('1e-8'),atol=D('1e-15')):
    a=D(a);b=D(b)
    assert abs(a-b)<=atol+rtol*max(abs(a),abs(b)),(where,a,b)

def restack(events):
    if not events:return None
    n=len(events);nr=len(events[0]['depth'])
    w=[[1/r['covariance'][j][j] for j in range(nr)] for r in events]
    sw=[sum(v[j] for v in w) for j in range(nr)]
    weights=[[v[j]/sw[j] for j in range(nr)] for v in w]
    mu=[sum(a[j]*ev['depth'][j] for a,ev in zip(weights,events)) for j in range(nr)]
    c0=[[sum(a[j]*a[k]*ev['covariance'][j][k] for a,ev in zip(weights,events)) for k in range(nr)] for j in range(nr)]
    chi=[sum(ww[j]*(ev['depth'][j]-mu[j])**2 for ww,ev in zip(w,events)) for j in range(nr)]
    inf=[max(D(1),v/max(1,n-1)) for v in chi]
    cov=[[c0[j][k]*(inf[j]*inf[k]).sqrt() for k in range(nr)] for j in range(nr)]
    return dict(n_events=n,depth=mu,covariance=cov,formal_covariance=c0,chi2=chi,inflation=inf,
                error=[cov[j][j].sqrt() for j in range(nr)])

def compare_stack(actual,saved,label):
    assert (actual is None)==(saved is None),label
    if actual is None:return
    assert actual['n_events']==saved['n_events'],label
    for name in ['depth','chi2','inflation','error']:
        for j,(a,b) in enumerate(zip(actual[name],saved[name])):close(a,b,f'{label}/{name}/{j}')
    for name in ['covariance','formal_covariance']:
        for j,row in enumerate(actual[name]):
            for k,a in enumerate(row):close(a,saved[name][j][k],f'{label}/{name}/{j}/{k}')
    assert all(x>=1 for x in actual['inflation'])

def ratio(s):
    if s is None or s['depth'][0]<=0:return dict(valid=False,passed=False)
    a,b=s['depth'][0],s['depth'][3];c=s['covariance'];r=b/a
    g=-b/a**2;h=1/a
    var=g*g*c[0][0]+h*h*c[3][3]+2*g*h*c[0][3]
    se=max(D(0),var).sqrt();z=(r-2)/se if se>0 else D('-Infinity')
    return dict(valid=True,ratio=r,error=se,z_above_2=z,denominator_z=a/c[0][0].sqrt(),
                passed=r>2 and z>=3)

def profile_result(rec,label):
    vals=[]
    for k,fit in enumerate(rec['fits']):
        hh,hy,yy=fit['hh'],fit['hy'],fit['yy']
        dep=max(D(0),-hy/hh) if hh>0 else D(0)
        gain=hy**2/hh if hh>0 and hy<0 else D(0)
        for name,v in [('depth',dep),('gain',gain),('chi2',yy-gain)]:close(v,fit[name],f'{label}/{k}/{name}')
        vals.append(gain)
    delta=vals[0]-vals[1];close(delta,rec['delta_chi2'],label+'/delta',atol=D('1e-10'))
    assert rec['n_events']==len(rec['events'])
    return delta

with localcontext() as ctx:
    ctx.prec=60
    plan=load('PLAN.json');first=load('outputs/first_flux.json')
    for p,h in first['hashes'].items():assert sha(p)==h,p
    assert first['source_sha256']==sha('sources/diamante.txt')=='a762f968e4a5a5c7653b83062f48a0f5481cd3386bd0f76d30731ac82e5778f7'
    for p,h in load('sources/input_hashes.json').items():assert sha(p)==h,p
    assert first['nice']==10 and all(v=='1' for v in first['threads'].values())
    for file,key in [('freeze_commit.txt','freeze_commit'),('prediction_commit.txt','prediction_commit')]:
        assert Path(file).read_text().strip()==first[key]
    for commit,paths in [(first['freeze_commit'],['PROTOCOL.md']),
                          (first['prediction_commit'],['PLAN.json','measurement.py','prepare.py','inject.py','real.py','outputs/phases.json','outputs/time_design.npz'])]:
        for path in paths:
            blob=subprocess.check_output(['git','show',f'{commit}:{path}'])
            assert hashlib.sha256(blob).hexdigest()==sha(path),path
        when=int(subprocess.check_output(['git','show','-s','--format=%ct',commit]))
        assert when<=datetime.datetime.fromisoformat(first['utc']).timestamp()
    assert plan['utc']<first['utc'] and plan['DIAmante_flux_read'] is False
    record('Pre-flux git protocol/prediction/code seals, chronology, input SHA256s, nice/thread settings')

    # Independently regenerate the ephemeris masks, seeded phases, and supplement eligibility.
    td=np.load('outputs/time_design.npz');tt=td['time'];rawt=np.loadtxt('sources/diamante.txt',usecols=[0])
    assert np.array_equal(tt,rawt)
    pe=3.0305267;pc=3.02553;te=2454967.124438-2457000;tc=2458682.5504-2457000
    def dist(t,p,t0):return (t-t0+p/2)%p-p/2
    pad=np.median(np.diff(tt))/2
    clean=(abs(dist(tt,pe,te))>.1530*pe/2+pad)&(abs(dist(tt,pe,te+.5007*pe))>.1384*pe/2+pad)&(abs(dist(tt,pc,tc))>.105+pad)
    assert np.array_equal(clean,td['clean'])
    # Read observed values only to verify eligibility, never to form an A stack.
    raw=np.loadtxt('sources/diamante.txt',usecols=range(29))
    flux=raw[:,[5,12,19,26]];error=raw[:,[6,13,20,27]];quality=raw[:,[7,14,21,28]]
    good=(quality==0)&np.isfinite(flux)&np.isfinite(error)&(flux>0)&(error>0)
    common_time=tt[clean&good.all(axis=1)];radius2_time=tt[clean&good[:,1]]
    ac=load('outputs/acceptance.json')
    assert len(tt)==ac['total'] and int(clean.sum())==ac['time_oot']
    assert good.sum(axis=0).tolist()==ac['accepted'] and len(common_time)==ac['common_oot']
    phases=load('outputs/phases.json');rng=np.random.default_rng(2026110236)
    for k in ['R','E']:assert np.array_equal(rng.uniform(0,pe,200),np.array(phases[k],float))
    st=tt[clean]
    def supported(x,name):
        if name=='R':
            epoch=te+x
            ns=np.arange(np.ceil((st.min()-epoch)/pe),np.floor((st.max()-epoch)/pe)+1)
            for center in epoch+ns*pe:
                dt=st-center
                if sum(abs(dt)<=.105)>=3 and sum((dt>=-.504)&(dt<=-.252))>=3 and sum((dt>=.252)&(dt<=.504))>=3:return True
            return False
        ref=te+round((tc-te)/pe)*pe
        for k in range(-3,23):
            a=ref+k*pe+x;b=tc+k*pc+x;mid=(a+b)/2;use=abs(st-mid)<=.6
            if use.sum()<12:continue
            time=st[use];dt=time-mid;ha=abs(time-a)<=.105;hb=abs(time-b)<=.105
            out=~(ha|hb)
            if min(ha.sum(),hb.sum())>=3 and sum(out&(dt<0))>=3 and sum(out&(dt>0))>=3:return True
        return False
    sp=load('outputs/supplement_phases.json');rng=np.random.default_rng(2026110237)
    for name in ['R','E']:
        accepted=[];tried=0
        while len(accepted)<200:
            x=float(rng.uniform(0,pe));tried+=1
            if supported(x,name):accepted.append(x)
        assert np.array_equal(accepted,np.array(sp['phases'][name],float))
        assert tried==sp['attempts'][name]
    sc=Path('supplement_commit.txt').read_text().strip()
    for p in ['SUPPLEMENT_PROTOCOL.md','supplement.py','outputs/supplement_phases.json']:
        assert hashlib.sha256(subprocess.check_output(['git','show',f'{sc}:{p}'])).hexdigest()==sha(p)
    record('Original and supplementary phase streams independently regenerated from timestamps only')

    # Reproduce Kepler mean-window baseline and both vectors from stored sector fractions.
    ke=load('outputs/kepler_events.json');ks=load('outputs/kepler_summary.json')
    kval={}
    for name in ['primary','secondary']:
        kval[name]=restack(ke[name]);compare_stack(kval[name],ks['combined'][name],'Kepler/'+name)
        compare_stack(kval[name],plan['kepler'][name],'PLAN/Kepler/'+name)
    tf=[D(0)]*4;leak=[D(0)]*4;sec=[D(0)]*4
    counts=[s['time_count'] for s in plan['sectors']];weights=[D(n)/sum(counts) for n in counts]
    for weight,s in zip(weights,plan['sectors']):
        for j in range(4):
            ll=s['eb_flux_fraction'][j]*kval['primary']['depth'][0]
            ss=s['eb_flux_fraction'][j]*kval['secondary']['depth'][0]
            close(ll,s['L'][j],f'predicted L sector {s["sector"]} r{j+1}')
            close(ss,s['secondary'][j],f'predicted S sector {s["sector"]} r{j+1}')
            leak[j]+=weight*ll;sec[j]+=weight*ss;tf[j]+=weight*s['target_flux_fraction'][j]
    for j in range(4):
        close(leak[j],plan['L'][j],f'L{j+1}');close(sec[j],plan['secondary'][j],f'S{j+1}')
        close(D('.000558')*tf[j]/tf[1],plan['T'][j],f'T{j+1}')
    close(leak[3]/leak[0],plan['L_R41'],'L ratio');close(tf[3]/tf[0],plan['T_R41'],'T ratio')
    record('Kepler 418/416 event stacks and all aperture prediction arithmetic')

    rresults={};noise_errors=[]
    for suffix,eventfile in [('', 'outputs/R_events.json'),('_supplement','outputs/R_supplement_events.json')]:
        saved=load(f'outputs/R{suffix}.json');raw=load(eventfile);counts={};actual={}
        for name in ['L','T','noise']:
            passed=0;nonempty=0
            assert len(raw[name])==len(saved['trials'][name])==200
            for k,(rr,ss) in enumerate(zip(raw[name],saved['trials'][name])):
                assert rr['trial']==ss['trial']==k
                source=phases['R'][k] if not suffix else sp['phases']['R'][k]
                assert rr['phase']==ss['phase']==source
                assert int(sum(abs(dist(common_time,pe,te+float(source)))<=.105))==ss['injected_cadences']
                a=restack(rr['events']);compare_stack(a,ss['stack'],f'R{suffix}/{name}/{k}')
                nonempty+=ss['injected_cadences']>0
                if name=='noise':
                    if not suffix and a:noise_errors.append(a['error'][3])
                    continue
                rt=ratio(a);assert rt['valid']==ss['ratio']['valid']
                assert rt['passed']==ss['ratio']['passed']
                if rt['valid']:
                    for key in ['ratio','error','z_above_2','denominator_z']:close(rt[key],ss['ratio'][key],f'R{suffix}/{name}/{k}/{key}')
                passed+=rt['passed']
            if name!='noise':counts[name]=passed;actual[name]=nonempty
        assert counts==saved['counts']
        passed=counts['L']>=160 and counts['T']<=10
        assert passed==saved['passed']
        if suffix:assert actual==saved['actual_injections']=={'L':200,'T':200}
        rresults[suffix or 'original']=dict(counts=counts,actual_injections=actual,passed=passed)
    original=load('outputs/R.json')
    noise=statistics.median(noise_errors);expected=plan['secondary'][3]/noise
    close(noise,original['secondary_noise_median'],'secondary noise')
    close(expected,original['expected_secondary_z'],'secondary expected z')
    assert original['secondary_noise_valid']==len(noise_errors)
    assert original['secondary_powered']==(expected>=3)
    record('All 800 L/T trial decisions independently recomputed from saved events; 400 placebo stacks restacked')

    eresults={}
    for file,key in [('outputs/E_injections.json','original'),('outputs/E_supplement.json','supplement')]:
        saved=load(file);deltas={};actual={}
        for name in ['EB','CTOI']:
            assert len(saved['trials'][name])==200
            for k,trial in enumerate(saved['trials'][name]):
                phase=phases['E'][k] if key=='original' else sp['phases']['E'][k]
                assert trial['trial']==k and trial['phase']==phase
                period,epoch=(pe,te) if name=='EB' else (pc,tc)
                assert int(sum(abs(dist(radius2_time,period,epoch+float(phase)))<=.105))==trial['injected_cadences']
            deltas[name]=[profile_result(r,f'E/{key}/{name}/{k}') for k,r in enumerate(saved['trials'][name])]
            actual[name]=sum(r['injected_cadences']>0 for r in saved['trials'][name])
        et=max(D(0),sorted(deltas['CTOI'])[189]);ct=min(D(0),sorted(deltas['EB'])[10])
        close(et,saved['threshold_EB'],f'E/{key}/threshold EB');close(ct,saved['threshold_CTOI'],f'E/{key}/threshold CTOI')
        counts=dict(EB_correct=sum(x>et for x in deltas['EB']),CTOI_correct=sum(x<ct for x in deltas['CTOI']),
                    EB_false=sum(x>et for x in deltas['CTOI']),CTOI_false=sum(x<ct for x in deltas['EB']))
        assert counts==saved['counts']
        power=counts['EB_correct']>=160 and counts['CTOI_correct']>=160 and counts['EB_false']<=10 and counts['CTOI_false']<=10
        assert power==saved['powered']
        if key=='supplement':assert actual==saved['actual_injections']=={'EB':200,'CTOI':200}
        eresults[key]=dict(counts=counts,actual_injections=actual,powered=power)
    record('All 800 ephemeris profile statistics, empirical thresholds and power/false-preference decisions')

    real=load('outputs/real.json')
    assert original['utc']<real['utc'] and load('outputs/E_injections.json')['utc']<real['utc']
    assert load('outputs/R_supplement.json')['utc']<real['utc'] and load('outputs/E_supplement.json')['utc']<real['utc']
    gate=rresults['original']['passed'] and rresults['_supplement']['passed']
    if not gate:
        assert real['A']['status']=='UNSCORED'
        assert 'stack' not in real['A'] and not Path('outputs/A_real_events.json').exists()
    else:
        a=restack(load('outputs/A_real_events.json'));compare_stack(a,real['A']['stack'],'real A')
        rt=ratio(a);assert real['A']['status']==('PASS' if rt['passed'] else 'FAIL')
        assert real['A']['T_consistent_2sigma']==(rt['valid'] and abs(rt['ratio']-plan['T_R41'])<=2*rt['error'])
    delta=profile_result(real['E'],'real E')
    ec=load('outputs/E_injections.json')
    preference='EB' if delta>ec['threshold_EB'] else 'CTOI' if delta<ec['threshold_CTOI'] else 'none'
    assert real['E']['preference']==preference
    power=eresults['original']['powered'] and eresults['supplement']['powered']
    assert real['E']['status']==('SCORED' if power else 'DESCRIPTIVE')
    ss=restack(load('outputs/S_real_events.json'));compare_stack(ss,real['S']['stack'],'real S')
    assert real['S']['status']==('SCORED' if expected>=3 else 'DESCRIPTIVE')
    close(ss['depth'][0]/ss['error'][0],real['S']['observed_z'],'S z')
    close((ss['depth'][0]-plan['secondary'][3])/ss['error'][0],real['S']['prediction_residual_z'],'S vs L')
    record('A withheld after failed gate, real E descriptive preference, S power/detection/prediction agreement')

    eb=load('outputs/EB_TESS.json');comb={name:[] for name in ['primary','secondary']}
    def verify_timing(tr,label):
        chis=[]
        for k,r in enumerate(tr['trials']):
            gain=r['hy']**2/r['hh'] if r['hy']<0 and r['hh']>0 else D(0)
            chi=r['yy']-gain;close(chi,r['chi2'],label+'/chi2',atol=D('1e-9'));chis.append(chi)
            close(max(D(0),-r['hy']/r['hh']) if r['hh'] else 0,r['amplitude'],label+'/amplitude')
        best=min(range(len(chis)),key=lambda i:chis[i]);lo=hi=best
        while lo>0 and chis[lo-1]<=chis[best]+1:lo-=1
        while hi<200 and chis[hi+1]<=chis[best]+1:hi+=1
        close(tr['trials'][best]['offset_days'],tr['best_offset_days'],label+'/offset')
        close(chis[best],tr['best_chi2'],label+'/best chi2',atol=D('1e-9'))
        close(chis[100]-chis[best],tr['delta_chi2_zero'],label+'/zero delta',atol=D('1e-9'))
        assert tr['grid_delta_chi2_1_interval']==[tr['trials'][lo]['offset_days'],tr['trials'][hi]['offset_days']]
    for s in eb['sectors']:
        ev=load(f'outputs/EB_TESS_s{s["sector"]:04d}_events.json')
        for name in ['primary','secondary']:
            comb[name]+=ev[name];compare_stack(restack(ev[name]),s[name],f'EB TESS {s["sector"]} {name}')
        assert sha(s['source'])==s['source_sha256']
        verify_timing(s['timing'],f'EB TESS {s["sector"]}')
    for name in ['primary','secondary']:compare_stack(restack(comb[name]),eb['combined'][name],'EB TESS combined '+name)
    verify_timing(eb['combined']['timing'],'EB TESS combined')
    record('EB-centred TESScut event depths, source hashes and all 603 timing-profile normal equations')
    if Path('ARTIFACTS.sha256.json').exists():
        manifest=load('ARTIFACTS.sha256.json')
        for path,h in manifest['files'].items():assert sha(path)==h,path
        record(f'Final artifact SHA256 seal: {len(manifest["files"])} files')
    out=dict(verified=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=checks,
        Gate_R=rresults,Bar_E_calibration=eresults,Bar_A=real['A']['status'],Bar_E=real['E']['status'],Bar_S=real['S']['status'],
        precision_digits=60,scope='Independent saved-event and profile arithmetic plus freezes and timestamps; not independent extraction or astrophysical truth.')
    Path('outputs/verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
