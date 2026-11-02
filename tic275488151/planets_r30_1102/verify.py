#!/usr/bin/env python3
"""Independent saved-output verifier. Imports no measurement/FITS/PRF/scoring code.
Decimal arithmetic at 70 digits; exact decimal inputs and explicit inequalities.
This certifies decisions on saved measurements, not extraction or astrophysical truth.
"""
from pathlib import Path
import json,hashlib,subprocess,decimal,datetime
from decimal import Decimal as D
from fractions import Fraction as Q
decimal.getcontext().prec=70
ROOT=Path(__file__).resolve().parent

def read(name):return json.loads((ROOT/name).read_text(),parse_float=D)
def ds(x):return D(str(x))
def check_close(a,b,tol='1e-8'):
    assert abs(ds(a)-ds(b))<=D(tol)*max(D(1),abs(ds(a)),abs(ds(b))),(a,b)
def stack(ev):
    if not ev:return {'n':0,'depth':None,'sigma':None,'z':None}
    w=[D(1)/e['variance'] for e in ev];sw=sum(w);dep=sum(a*e['depth'] for a,e in zip(w,ev))/sw
    chi=sum(a*(e['depth']-dep)**2 for a,e in zip(w,ev));infl=max(D(1),chi/D(max(1,len(ev)-1)));sig=(infl/sw).sqrt()
    return {'n':len(ev),'depth':dep,'sigma':sig,'z':dep/sig}
def inv(c):
    a,b=c[0];cc,d=c[1];det=a*d-b*cc;assert det>0
    return [[d/det,-b/det],[-cc/det,a/det]]
def matvec(a,v):return [sum(a[i][j]*v[j] for j in range(2)) for i in range(2)]
def norm(v):return sum(x*x for x in v).sqrt()
def verify():
    commit=(ROOT/'freeze_commit.txt').read_text().strip()
    frozen=subprocess.check_output(['git','show',commit+':PROTOCOL.md'],cwd=ROOT)
    assert frozen==(ROOT/'PROTOCOL.md').read_bytes(),'frozen protocol changed'
    manifest=read('aperture_manifest.json');h=hashlib.sha256((ROOT/'aperture_manifest.json').read_bytes()).hexdigest()
    assert h==(ROOT/'aperture_manifest.sha256').read_text().strip()
    assert manifest['protocol_sha256']==hashlib.sha256(frozen).hexdigest()
    first=read('outputs/first_flux_read.json');assert first['aperture_sha256']==h
    assert first['utc']>=manifest['created_utc']
    implementation=read('outputs/implementation_before_flux.json');assert implementation['utc']<=first['utc']
    for p,expected in implementation['sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==expected,p
    assert len(manifest['sectors'])==10
    assert [r['sector'] for r in manifest['sectors']]==[14,15,40,41,54,55,74,75,81,82]
    for r in manifest['sectors']:
        assert len(r['aperture_flat'])==4 and len(set(r['aperture_flat']))==4
        tp=[v for row in r['target_prf'] for v in row];fp=[v for row in r['fractions'] for v in row]
        ii=[i for i in range(225) if tp[i]>0];ii.sort(key=lambda i:(-fp[i],-tp[i],i//15,i%15));assert ii[:4]==r['aperture_flat']
    photo=read('outputs/photometry.json');calc={};event_count=0
    for r in photo['sectors']:
        for name,m in r['measurements'].items():
            ev=m['events'];event_count+=len(ev)
            for e in ev:assert e['n_in']>=3 and e['n_out_left']>=3 and e['n_out_right']>=3 and e['variance']>0
            s=stack(ev)
            for k in ('depth','sigma','z'):
                if s[k] is not None:check_close(s[k],m['stack'][k])
            assert s['n']==m['stack']['n']
    for name in ('EB','CTOI','secondary','odd','even'):
        ev=[e for r in photo['sectors'] for e in r['measurements'][name]['events']];calc[name]=stack(ev)
        for k in ('depth','sigma','z'):
            if calc[name][k] is not None:check_close(calc[name][k],photo['combined'][name][k])
    z=calc['EB']['z'];zc=calc['CTOI']['z'];a=z is not None and zc is not None and z>=D('7.1') and z-zc>=3
    assert a==photo['bar_A']
    gate=read('outputs/gate.json');design=read('outputs/injection_design.json');assert len(design['trials'])==40 and len(gate['results'])==40
    counts={};trial_checks=0
    for star in ('target','eb'):
        trials=[r for r in gate['results'] if r['star']==star];assert sorted(r['trial'] for r in trials)==list(range(20));n=0
        for r in trials:
            trial_checks+=1;success=False
            frozen_trial=next(d for d in design['trials'] if d['star']==star and d['trial']==r['trial'])
            assert r['sector']==frozen_trial['sector'] and r['phase_days']==frozen_trial['phase_days']
            if r['fit'] is not None:
                rec=next(s for s in manifest['sectors'] if s['sector']==r['sector']);xy=rec[star+'_xy'];fit=r['fit']['xy']
                # Exact rational comparison to the half-pixel circle.
                sq=sum((Q(str(fit[i]))-Q(str(xy[i])))**2 for i in range(2));success=sq<=Q(1,4)
                check_close(D(sq.numerator).sqrt()/D(sq.denominator).sqrt(),r['distance_pixels'])
            assert success==r['success'];n+=int(success)
        counts[star]=n;assert n==gate['counts'][star]
    gatepass=all(v>=18 for v in counts.values());assert gatepass==gate['gate_pass']
    loc=read('outputs/location.json')
    if not gatepass:
        assert loc['bar_B']=='UNSCORED' and loc['sectors']==[];b='UNSCORED'
    else:
        valid=[r for r in loc['sectors'] if r['fit'] is not None];prec=[[D(0),D(0)],[D(0),D(0)]];rhs=[D(0),D(0)];trhs=[D(0),D(0)]
        for r in valid:
            c=inv(r['covariance_tangent']);v=matvec(c,r['offset_eb']);tv=matvec(c,r['target_offset'])
            for i in range(2):
                rhs[i]+=v[i];trhs[i]+=tv[i]
                for j in range(2):prec[i][j]+=c[i][j]
        cov=inv(prec);mean=matvec(cov,rhs);target=matvec(cov,trhs);dist=norm(mean);delta=[mean[i]-target[i] for i in range(2)];pd=matvec(prec,delta);sig=sum(delta[i]*pd[i] for i in range(2)).sqrt()
        b=len(valid)==10 and dist<=D('.5') and sig>=3
        check_close(dist,loc['distance_eb_pixels']);check_close(sig,loc['target_sigma']);assert b==loc['bar_B']
    attributed=a and gatepass and b is True
    score=read('outputs/score.json');assert score['bar_A']==a and score['bar_B']==b and score['gate_pass']==gatepass and score['attributed']==attributed
    seal=read('outputs/seal.json')
    for p,expected in seal['sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==expected,p
    report={'verified':True,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'freeze_commit':commit,'aperture_sha256':h,'bar_A':a,'z_EB':str(z),'z_CTOI':str(zc),'delta_z':str(z-zc) if z is not None and zc is not None else None,'gate_counts':counts,'gate_denominator_each':20,'gate_pass':gatepass,'bar_B':b,'attributed':attributed,'event_rows_checked_including_descriptive':event_count,'injection_trials_checked':trial_checks,'sealed_files_checked':len(seal['sha256']),'limits':'Saved-output arithmetic and hashes only; does not independently reproduce FITS extraction or PRF model.'}
    (ROOT/'outputs/verification.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
if __name__=='__main__':verify()
