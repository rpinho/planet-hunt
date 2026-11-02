"""Independent saved-output verifier. No measurement or FITS-reader imports.
70-digit Decimal restacking; exact rational threshold inequalities on those decimals.
This verifies saved sufficient-statistic arithmetic, not independent photometric extraction.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,hashlib,datetime,subprocess,math,urllib.parse,xml.etree.ElementTree as ET
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
import numpy as np
Z=D('2.7477813854449926');LIMIT=D('0.000291')
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def close(a,b,atol=D('1e-12'),rtol=D('1e-10')):
    assert abs(a-D(str(b)))<=atol+rtol*abs(a),(str(a),b)
def ulpass(d,v):
    gap=LIMIT-max(D(0),d)
    return gap>0 and F(gap)**2>F(Z)**2*F(v)
def snrpass(d,v):return d>0 and F(d)**2>=25*F(v)

def independent_events(path):
    z=np.load(path);a=z['events'];names=list(z['columns']);assert names[-2:]==['depth','variance']
    groups={}
    for row in a:
        x={k:D(str(float(v))) for k,v in zip(names,row)}
        ni=x['n_in'];nl=x['n_left'];nr=x['n_right'];no=nl+nr
        assert min(ni,nl,nr)>=3 and x['a']>0
        sx=x['sum_dt_out'];sxx=x['sum_dt2_out'];det=no*sxx-sx*sx;assert det>0
        xm=x['mean_dt_in'];d=(x['a']+x['b']*xm-x['mean_flux_in'])/x['a']
        h=(sxx-2*xm*sx+xm*xm*no)/det
        v=max(x['residual_s2'],x['err2_mean'])*(1/ni+h)/(x['a']**2)
        assert v>0;close(d,x['depth']);close(v,x['variance'],D('1e-25'),D('1e-10'))
        trial=int(x['trial']);quarter=int(x['quarter']);groups.setdefault(trial,[]).append((d,v,quarter))
    return groups,len(a)

def independent_stack(ev):
    assert len(ev)>=2
    w=[1/v for d,v,q in ev];sw=sum(w);d=sum(wi*e[0] for wi,e in zip(w,ev))/sw
    chi=sum(wi*(e[0]-d)**2 for wi,e in zip(w,ev));infl=max(D(1),chi/(len(ev)-1));v=infl/sw
    return {'n':len(ev),'d':d,'v':v,'chi':chi,'infl':infl,'upper':max(D(0),d)+Z*v.sqrt(),'z':d/v.sqrt()}
def check_stack(s,r):
    assert s['n']==r['n'] and r['valid']
    for a,b in [('d','depth'),('upper','upper'),('z','z'),('chi','chi2'),('infl','inflation_variance')]:close(s[a],r[b])
    close(s['v'].sqrt(),r['sigma'])
def compact(s):return {k:float(v) if isinstance(v,D) else v for k,v in s.items()}

def geometry_and_network():
    # Separate scalar tangent-plane transform and explicit 2x2 solve.
    ap=load('outputs/apertures.json');nq=18 if Path('outputs/completed_gate.json').exists() else 17;assert len(ap['quarters'])==nq
    distances=[];seasons=set()
    net=[json.loads(l) for l in Path('network.jsonl').read_text().splitlines()]
    for q in ap['quarters']:
        h=q['aperture_header'];mask=q['mask'];pixels=[(i,j) for j,row in enumerate(mask) for i,v in enumerate(row) if int(v)&2]
        assert len(pixels)==q['optimal_pixels'];seasons.add(q['season'])
        assert any(r.get('sha256')==q['compressed_sha256'] and r['url']==q['source_url'] for r in net)
        for star,c in ap['coordinates'].items():
            ra=float(c['ra']);dec=float(c['dec']);years=(q['position_epoch_bjd']-(2451545+15.5*365.25))/365.25
            ra+=float(c['pm_ra'])*years/(3600000*math.cos(math.radians(dec)));dec+=float(c['pm_dec'])*years/3600000
            rr,dd,r0,d0=map(math.radians,[ra,dec,h['CRVAL1'],h['CRVAL2']]);dr=rr-r0
            den=math.sin(d0)*math.sin(dd)+math.cos(d0)*math.cos(dd)*math.cos(dr)
            u=math.degrees(math.cos(dd)*math.sin(dr)/den);v=math.degrees((math.cos(d0)*math.sin(dd)-math.sin(d0)*math.cos(dd)*math.cos(dr))/den)
            a=h['CDELT1']*h['PC1_1'];b=h['CDELT1']*h['PC1_2'];c2=h['CDELT2']*h['PC2_1'];d=h['CDELT2']*h['PC2_2'];det=a*d-b*c2
            x=(d*u-b*v)/det+h['CRPIX1']-1;y=(-c2*u+a*v)/det+h['CRPIX2']-1
            r=q['positions'][star];assert max(abs(x-r['cutout_xy0'][0]),abs(y-r['cutout_xy0'][1]))<1e-7
            dmin=min(math.hypot(max(abs(x-i)-.5,0),max(abs(y-j)-.5,0)) for i,j in pixels)
            assert abs(dmin-r['minimum_distance_to_optimal_pixel_union'])<1e-7
            inside=any(abs(x-i)<=.5 and abs(y-j)<=.5 for i,j in pixels);assert inside==r['inside']
            assert inside if star=='target' else not inside
            if star=='eb':distances.append(dmin)
    assert seasons==set(range(4))
    starts=[datetime.datetime.fromisoformat(r['utc']).timestamp() for r in net]
    assert min(b-a for a,b in zip(starts,starts[1:]))>=4.99
    for r in net:
        u=urllib.parse.urlparse(r['url'])
        assert u.scheme=='https' and u.hostname in {'arxiv.org','exofop.ipac.caltech.edu','archive.stsci.edu','stpubdata.s3.amazonaws.com','mast.stsci.edu'}
        if r.get('status'):
            assert r['status'] in (200,206) and r['url']==r['final_url']
            p=Path(r['path'])
            if p.exists():assert sha(p)==r['sha256']
    return {'quarters':nq,'seasons':4,'target_inside':nq,'eb_outside':nq,'eb_distance_min_px':min(distances),'eb_distance_max_px':max(distances),'public_GET_attempts':len(net),'minimum_start_spacing_s':min(b-a for a,b in zip(starts,starts[1:]))}

with localcontext() as ctx:
    ctx.prec=70
    design=load('outputs/design.json');gate=load('outputs/gate.json');before=load('outputs/before_gate_flux.json')
    commit=Path('freeze_commit.txt').read_text().strip();assert design['commit']==commit
    frozen=subprocess.check_output(['git','show',commit+':PROTOCOL.md']);assert hashlib.sha256(frozen).hexdigest()==sha('PROTOCOL.md')==design['protocol_sha256']==before['protocol_sha256']
    assert before['implementation_sha256']==sha('measure.py') and before['reader_sha256']==sha('fits_local.py')
    assert subprocess.check_output(['git','show',commit+':STATUS.md']).startswith(b'FROZEN\n')
    commit_time=subprocess.check_output(['git','show','-s','--format=%cI',commit],text=True).strip()
    dt=lambda s:datetime.datetime.fromisoformat(s)
    assert dt(commit_time)<dt(design['utc'])<dt(before['utc'])<dt(gate['utc'])
    assert before['nice']==10 and design['nice']==10 and all(v=='1' for v in design['threads'].values())
    # Independently regenerate every allowed phase from the registered seed.
    rng=np.random.Generator(np.random.PCG64(321102));regenerated=[]
    while len(regenerated)<400:
        p=float(rng.uniform())
        if min(abs((p-c+.5)%1-.5) for c in (0,.5007))>=.15:regenerated.append(p)
    assert regenerated==design['injection_phases']+design['null_phases']
    inherited=load('../planets_r28_1102/receipts/10014894.json')
    expect={Path(r['key']).name:r['sha256'] for r in inherited['files']}
    assert len(gate['manifest'])==17
    for r in gate['manifest']:
        assert sha(r['file'])==r['sha256']==expect[Path(r['file']).name]
    assert gate['events_sha256']==sha('outputs/gate_events.npz')
    groups,n_events=independent_events('outputs/gate_events.npz');assert set(groups)==set(range(400))
    counts={'injection':0,'null':0};recovered=[]
    assert len(gate['trials'])==400
    for i,t in enumerate(gate['trials']):
        assert t['trial']==i and t['phase']==regenerated[i]
        arm='injection' if i<200 else 'null';assert t['arm']==arm
        s=independent_stack(groups[i]);check_stack(s,t['stack'])
        passed=bool(snrpass(s['d'],s['v']) if arm=='injection' else ulpass(s['d'],s['v']))
        assert passed==t['passed'];counts[arm]+=passed
        assert sum(q['valid'] for q in t['quarters'])==s['n']
        for q in t['quarters']:
            assert q['predicted']==q['valid']+q['too_few']+q['invalid_fit']
            assert q['valid']==sum(e[2]==q['quarter'] for e in groups[i])
        recovered.append({'trial':i,'passed':passed,**compact(s)})
    passed=all(v>=190 for v in counts.values());assert counts==gate['counts'] and passed==gate['pass']
    assert gate['denominators']=={'injection':200,'null':200}
    cert={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'verified':True,'precision_digits':70,'thresholds':'exact rational squared inequalities from 70-digit reconstructed stack','gate_counts':counts,'gate_denominators':gate['denominators'],'gate_pass':passed,'gate_events':n_events,'bar_K':'UNSCORED','claim_permitted':False,'freeze_commit':commit,'phases_regenerated':400,'input_LLC_hashes':17,'gate_sha256':sha('outputs/gate.json'),'verifier_sha256':sha(__file__)}
    if Path('outputs/real.json').exists():
        assert passed
        real=load('outputs/real.json');br=load('outputs/before_real_flux.json')
        assert dt(gate['utc'])<dt(br['utc'])<dt(real['utc']) and br['gate_sha256']==sha('outputs/gate.json') and br['implementation_sha256']==sha('measure.py')
        assert real['events_sha256']==sha('outputs/real_events.npz')
        g2,n2=independent_events('outputs/real_events.npz');assert set(g2)=={0,1};cert['real']={}
        for j,name in enumerate(('primary','secondary')):
            r=real['results'][name];s=independent_stack(g2[j]);check_stack(s,r['stack']);cert['real'][name]=compact(s)
            for q in r['quarters']:
                eq=[e for e in g2[j] if e[2]==q['quarter']];assert len(eq)==q['valid']
                if len(eq)>=2:check_stack(independent_stack(eq),q['stack'])
                assert q['predicted']==q['valid']+q['too_few']+q['invalid_fit']
        p=independent_stack(g2[0]);bar=bool(ulpass(p['d'],p['v']));assert bar==real['bar_K']
        cert.update(bar_K=bar,claim_permitted=passed and bar,real_events=n2,real_sha256=sha('outputs/real.json'))
    else:
        assert not Path('outputs/before_real_flux.json').exists()
    if Path('outputs/completed_gate.json').exists():
        cg=load('outputs/completed_gate.json');bq=load('outputs/before_q0_gate_flux.json')
        cc=Path('completeness_commit.txt').read_text().strip();assert cc==bq['completion_commit']
        cct=subprocess.check_output(['git','show','-s','--format=%cI',cc],text=True).strip()
        assert dt(cct)<dt(bq['utc'])<dt(cg['utc']) and dt(real['utc'])<dt(cct)
        assert subprocess.check_output(['git','show',cc+':COMPLETENESS_CORRECTION.md'])==Path('COMPLETENESS_CORRECTION.md').read_bytes()
        assert bq['implementation_sha256']==sha('complete_q0.py') and bq['original_implementation_sha256']==sha('measure.py') and bq['nice']==10
        assert cg['initial_gate_sha256']==sha('outputs/gate.json') and cg['initial_real_was_seen']
        assert cg['q0_events_sha256']==sha('outputs/q0_gate_events.npz')
        gq,nq=independent_events('outputs/q0_gate_events.npz')
        assert all(e[2]==0 for v in gq.values() for e in v)
        assert len(cg['manifest'])==18 and {r['quarter'] for r in cg['manifest']}==set(range(18))
        assert cg['manifest'][:17]==gate['manifest']
        assert sha(cg['manifest'][-1]['file'])==cg['manifest'][-1]['sha256']
        ns={'s':'http://s3.amazonaws.com/doc/2006-03-01/'}
        listed={Path(k.text).name for k in ET.parse('prior/llc_list.xml').findall('s:Contents/s:Key',ns) if k.text.endswith('_llc.fits')}
        assert listed=={Path(r['file']).name for r in cg['manifest']}
        ccnt={'injection':0,'null':0}
        for i,t in enumerate(cg['trials']):
            orig=gate['trials'][i];assert t['phase']==orig['phase'] and t['quarters'][:17]==orig['quarters']
            ev=groups[i]+gq.get(i,[]);s=independent_stack(ev);check_stack(s,t['stack'])
            ok=bool(snrpass(s['d'],s['v']) if t['arm']=='injection' else ulpass(s['d'],s['v']))
            assert ok==t['passed'];ccnt[t['arm']]+=ok
            assert sum(q['valid'] for q in t['quarters'])==len(ev)
            assert t['quarters'][-1]['valid']==len(gq.get(i,[]))
        cp=all(v>=190 for v in ccnt.values());assert ccnt==cg['counts'] and cp==cg['pass']
        cert.update(initial_gate_counts=cert['gate_counts'],initial_bar_K=cert['bar_K'],initial_real=cert.get('real'),gate_counts=ccnt,gate_pass=cp,gate_events=n_events+nq,bar_K='UNSCORED',claim_permitted=False,archive_LLC_products=18,completion_commit=cc,sequencing_deviation='Q1-Q17 real result was seen before missing Q0 completion; no original event was remeasured',completed_gate_sha256=sha('outputs/completed_gate.json'))
        if Path('outputs/completed_real.json').exists():
            assert cp
            cr=load('outputs/completed_real.json');bqr=load('outputs/before_q0_real_flux.json')
            assert dt(cg['utc'])<dt(bqr['utc'])<dt(cr['utc']) and bqr['completed_gate_sha256']==sha('outputs/completed_gate.json')
            assert bqr['implementation_sha256']==sha('complete_q0.py')
            assert cr['initial_real_sha256']==sha('outputs/real.json') and cr['q0_events_sha256']==sha('outputs/q0_real_events.npz')
            qr,nqr=independent_events('outputs/q0_real_events.npz');cert['real']={}
            for j,name in enumerate(('primary','secondary')):
                ev=g2[j]+qr.get(j,[]);s=independent_stack(ev);check_stack(s,cr['results'][name]['stack']);cert['real'][name]=compact(s)
                assert cr['results'][name]['quarters'][:17]==real['results'][name]['quarters']
                for q in cr['results'][name]['quarters']:
                    eq=[e for e in ev if e[2]==q['quarter']];assert len(eq)==q['valid']
                    if len(eq)>=2:check_stack(independent_stack(eq),q['stack'])
                    assert q['predicted']==q['valid']+q['too_few']+q['invalid_fit']
            s=independent_stack(g2[0]+qr.get(0,[]));bar=bool(ulpass(s['d'],s['v']));assert bar==cr['bar_K']
            cert.update(bar_K=bar,claim_permitted=cp and bar,real_events=n2+nqr,completed_real_sha256=sha('outputs/completed_real.json'))
        else:assert not Path('outputs/before_q0_real_flux.json').exists()
    if Path('outputs/apertures.json').exists():cert['geometry_network']=geometry_and_network()
    if Path('outputs/summary.json').exists():
        summ=load('outputs/summary.json')
        assert summ['gate_counts']==cert['gate_counts'] and summ['gate_pass']==cert['gate_pass'] and summ['bar_K']==cert['bar_K']
        for name in ('primary','secondary'):
            close(D(str(cert['real'][name]['d'])),summ[name]['depth'])
            close(D(str(cert['real'][name]['upper'])),summ[name]['upper'])
        assert summ['primary_valid_predicted']==[sum(q['valid'] for q in cr['results']['primary']['quarters']),sum(q['predicted'] for q in cr['results']['primary']['quarters'])]
        for row in summ['quarters']:
            for name in ('primary','secondary'):
                q=next(q for q in cr['results'][name]['quarters'] if q['quarter']==row['quarter'])
                assert row[name+'_valid']==q['valid'] and row[name+'_predicted']==q['predicted']
                close(D(str(q['stack']['depth']))*1000000,row[name+'_depth_ppm'])
                close(D(str(q['stack']['sigma']))*1000000,row[name+'_sigma_ppm'])
                assert row['CROWDSAP']==q['crowdsap'] and row['FLFRCSAP']==q['flfrcsap']
        cert['report_summary_checked']=True
    if Path('ARTIFACTS.sha256.json').exists():
        seal=load('ARTIFACTS.sha256.json')
        for p,h in seal.items():assert sha(p)==h,p
        cert['sealed_artifacts_checked']=len(seal)
    Path('outputs/verification.json').write_text(json.dumps(cert,indent=2,allow_nan=False))
    print(json.dumps(cert,indent=2))
