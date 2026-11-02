"""Independent saved-output audit. Imports NO measurement/search/scoring modules."""
import pathlib,json,hashlib,math,datetime,sys,csv
from fractions import Fraction
def q(x):return Fraction(str(x))
R=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else pathlib.Path('.')
def read(name):return json.loads((R/name).read_text())
def sha(name):return hashlib.sha256((R/name).read_bytes()).hexdigest()
def close(a,b):assert math.isclose(a,b,rel_tol=2e-9,abs_tol=2e-9),(a,b)
def depth_ok(m,s):
 cat=q(s.get('depth_error_ppm') or 0)*q(s.get('kepler_tess_dilution_ratio',1.))
 delta=q(m['depth'])-q(s['expected_kepler_aperture_depth_ppm'])
 return delta*delta<=9*(q(m['sigma'])**2+cat*cat)

def check_fit(m):
 assert all(math.isfinite(m[k]) for k in ['z','depth','sigma','norm','period_days','epoch_bkjd'])
 if m['norm']>0:close(m['sigma'],1/math.sqrt(m['norm']));close(m['z'],m['depth']/m['sigma'])
def check_search(m,s):
 for f in [m['best'],m['selected']]+m['matches']:check_fit(f)
 for f in m['matches']:assert f['z']>=7.1 and depth_ok(f,s)
 if m['matches']:assert m['selected']==max(m['matches'],key=lambda f:f['z'])
 else:assert m['selected']==m['best']
 if 'local_maxima' in m:
  assert m['best']==max(m['local_maxima'],key=lambda x:x['z'])
  assert m['matches']==[x for x in m['local_maxima'] if x['z']>=7.1 and depth_ok(x,s)]
 d=s['duration_hours']/24;p=s['period_days'];half=5*s['period_error_days'] if s.get('period_error_days') else .001*p
 expected_n=max(2,math.ceil(2*half/(d*p/(10*1459.5)))+1)
 expected_phases=max(32,math.ceil((p+half)/(d/40)))
 assert m['period_grid_n']==expected_n and m['phase_grid_n']==expected_phases
 for item in [m['best'],m['selected']]+m['matches']:assert p-half-1e-12<=item['period_days']<=p+half+1e-12
 assert m['phase_step_days_max']<=d/5
 assert m['max_drift_days']<=d/5+1e-10 and not m['grid_capped']
 close(m['period_step_days'],2*half/(m['period_grid_n']-1))
 assert m['trial_count']==m['phase_grid_n']*m['period_grid_n']
plan=read('PLAN.json');old=read('R27_PLAN.json');freeze=read('freeze.json');seal=read('results/output_manifest.json');score=read('results/scores.json')
amendment=read('audit/verifier_amendment.json')
assert sha('results/scores.json')==amendment['sole_score_sha256']
assert sha('verify.py')==amendment['corrected_verifier_sha256']
assert sha('audit/verify_frozen_949cd50.py')==amendment['registered_verifier_sha256']==freeze['implementation_sha256']['verify.py']
assert amendment['freeze_commit']==freeze['commit'] and amendment['plan_sha256']==freeze['plan_sha256']
def registered_source(name):return 'audit/verify_frozen_949cd50.py' if name=='verify.py' else name
assert sha('PLAN.json')==freeze['plan_sha256']==score['plan_sha256']
assert sha('R27_PLAN.json')==plan['r27_plan_sha256']
assert sha('metadata_headers.json')==old['sources_sha256']['catalogs/metadata_headers.json']
assert plan['targets']==[s for s in old['targets'] if s['period_days'] and s['duration_hours']]
assert len(plan['targets'])==12
method=hashlib.sha256((R/'pod_tasks/search.py').read_bytes()+(R/'pod_tasks/io_utils.py').read_bytes()).hexdigest();assert method==freeze['method_sha256']==score['method_sha256']
for name,want in seal['sha256'].items():assert sha(registered_source(name))==want,name
for name,want in freeze['implementation_sha256'].items():assert sha(registered_source(name))==want,name
assert freeze['utc']<seal['utc']<score['scored_utc']
assert sha('results/gate.json')==freeze['gate_sha256']
assert sha('results/stage1_verification.json')==freeze['stage1_verification_sha256']
gate_certificate=read('results/gate.json')
assert sha('catalogs/koi_frozen.csv')==old['sources_sha256']['catalogs/koi.csv']
kois={r['kepoi_name']:r for r in csv.DictReader(open(R/'catalogs/koi_frozen.csv'))}
assert sha('catalogs/diamante_frozen.csv')==old['sources_sha256']['catalogs/diamante_v2.csv']
diamante={int(r['ticID']):r for r in csv.DictReader(open(R/'catalogs/diamante_frozen.csv'))}
chosen_hosts=set()
for k in plan['known_planets']:
 template=next(s for s in plan['targets'] if s['ctoi']==k['matched_to'])
 candidates=[r for r in kois.values() if r['koi_disposition']=='CONFIRMED' and r['koi_period'] and r['koi_depth'] and r['koi_duration'] and int(r['kepid']) not in chosen_hosts and .5<float(r['koi_period'])<20 and 100<float(r['koi_depth'])<6000]
 nearest=min(candidates,key=lambda r:math.log(float(r['koi_period'])/template['period_days'])**2+math.log(float(r['koi_depth'])/template['expected_kepler_pdcsap_depth_ppm'])**2)
 assert nearest['kepoi_name']==k['koi'] and int(nearest['kepid'])==k['kic']
 chosen_hosts.add(k['kic'])
known_success=0
for k in plan['known_planets']:
 source=kois[k['koi']];assert source['koi_disposition']=='CONFIRMED' and int(source['kepid'])==k['kic']
 for dest,src in [('period_days','koi_period'),('epoch_bkjd','koi_time0bk'),('depth_ppm','koi_depth'),('duration_hours','koi_duration')]:close(k[dest],float(source[src]))
 r=read(f"known/{k['kic']}.json");assert r['method_sha256']==method and r['native_flux_unscaled']
 for key in ['period_days','epoch_bkjd','depth_ppm','duration_hours','kic','koi']:assert r[key]==k[key]
 f=r['measurement']['best'];check_fit(f)
 half=.001*k['period_days'];d=k['duration_hours']/24;n=max(2,math.ceil(2*half/(d*k['period_days']/(10*1459.5)))+1);dp=2*half/(n-1)
 phase=abs((f['epoch_bkjd']-k['epoch_bkjd']+f['period_days']/2)%f['period_days']-f['period_days']/2)
 success=f['z']>=7.1 and abs(f['period_days']-k['period_days'])<=2*dp and phase<=max(.02043365,d)
 assert success==r['success'];close(phase,r['phase_distance_days']);known_success+=success
kg=len(plan['known_planets'])==5 and known_success>=4
assert score['known_recovered']==known_success and score['known_total']==5
counts={'RECOVERED':0,'REFUTED':0,'OPEN':0};gatecount=0;injection_count=0;control_count=0;coverage_events=0
# Coverage is independently rebuilt from accepted exposure times on the pod.
# --arithmetic-only is a local certificate check, never represented as re-extraction.
full_coverage='--arithmetic-only' not in sys.argv
if full_coverage:import numpy as np
for s in plan['targets']:
 sid=s['ctoi'];r=read(f'targets/{sid}.json');m=r['measurement'];cv=r['coverage'];row=next(x for x in score['rows'] if x['ctoi']==sid)
 assert r['plan_sha256']==freeze['plan_sha256'] and r['method_sha256']==method and r['started_utc']>freeze['utc']
 ratio=(1+s['tic_contratio'])*s['kepler_crowdsap_median_q1_q17'];close(ratio,s['kepler_tess_dilution_ratio']);close(s['depth_ppm']*ratio,s['expected_kepler_aperture_depth_ppm'])
 check_search(m,s);check_fit(r['nominal_fit']);close(r['expected_snr'],s['expected_kepler_aperture_depth_ppm']*math.sqrt(r['nominal_fit']['norm']))
 close(cv['nominal_fraction'],cv['nominal_covered']/cv['nominal_total']);close(cv['grid_min_fraction'],min(a/b if b else 0 for a,b in zip(cv['grid_covered'],cv['grid_total'])));close(cv['best_fraction'],cv['best_covered']/cv['best_total'])
 if full_coverage:
  arr=np.load(R/f'coverage/{sid}.npz');t=arr['times'];periods=arr['periods'];d=s['duration_hours']/24;lc=.02043365;e=s['epoch_bjd']-2454833 if s['epoch_bjd'] is not None else (131.5+1591)/2
  half=5*s['period_error_days'] if s['period_error_days'] else .001*s['period_days']
  assert np.allclose(periods,np.linspace(s['period_days']-half,s['period_days']+half,m['period_grid_n']),rtol=0.,atol=1e-12)
  starts=t-lc/2;ends=t+lc/2;prefix_start=np.r_[0.,np.cumsum(starts)];prefix_end=np.r_[0.,np.cumsum(ends)]
  def independent_counts(p,epoch):
   ns=np.arange(math.ceil((131.5-epoch)/p),math.floor((1591-epoch)/p)+1);centers=epoch+ns*p
   # Integral of the sum of accepted exposure indicator functions, including edges.
   def integral(x):
    a=np.searchsorted(starts,x,side='left');b=np.searchsorted(ends,x,side='left')
    return a*x-prefix_start[a]-(b*x-prefix_end[b])
   covered=integral(centers+d/2)-integral(centers-d/2);hit=covered>=d/2
   for ix in np.flatnonzero(np.abs(covered-d/2)<5e-7):
    center=float(centers[ix]);left=center-d/2;right=center+d/2
    a=np.searchsorted(t,center-(d+lc)/2);b=np.searchsorted(t,center+(d+lc)/2)
    overlaps=[max(0.,min(float(x)+lc/2,right)-max(float(x)-lc/2,left)) for x in t[a:b]]
    hit[ix]=math.fsum(overlaps)>=d/2
   return int(hit.sum()),len(ns)
  a,b=independent_counts(s['period_days'],e);assert(a,b)==(cv['nominal_covered'],cv['nominal_total'])
  for ip,p0 in enumerate(periods):
   a,b=independent_counts(float(p0),e);assert(a,b)==(cv['grid_covered'][ip],cv['grid_total'][ip]),(sid,ip,a,b)
   coverage_events+=b
  a,b=independent_counts(m['selected']['period_days'],m['selected']['epoch_bkjd']);assert(a,b)==(cv['best_covered'],cv['best_total'])
  assert len(t)==r['points'] and len(periods)==m['period_grid_n']
 rr=[read(f"validation/{sid}_{c['kic']}.json") for c in s['controls']];assert len(rr)==50
 injections=[];zs=[];passing=0
 for c,ctrl in zip(s['controls'],rr):
  assert ctrl['kic']==c['kic'] and ctrl['ctoi']==sid and ctrl['available'] and ctrl['method_sha256']==method
  check_search(ctrl['control'],s);zs.append(ctrl['control']['best']['z']);passing+=bool(ctrl['control']['matches']);assert len(ctrl['injections'])==2
  for tr in ctrl['injections']:
   f=tr['found'];check_fit(f);p=f['period_days'];phase=abs((f['epoch_bkjd']-tr['injected_epoch_bkjd']+p/2)%p-p/2);close(phase,tr['phase_distance_days'])
   located=phase<=max(.02043365,s['duration_hours']/24) and abs(p-tr['injected_period_days'])<=tr['period_tolerance_days']+1e-10
   dg=depth_ok(f,s);success=f['z']>=7.1 and dg and located
   assert located==tr['located'] and dg==tr['depth_consistent'] and success==tr['success']
   half=5*s['period_error_days'] if s['period_error_days'] else .001*s['period_days']
   assert s['period_days']-half<=tr['injected_period_days']<=s['period_days']+half
   assert abs(tr['injected_epoch_bkjd']-(131.5+1591)/2)<=tr['injected_period_days']/2
   close(tr['injected_depth_ppm'],s['expected_kepler_aperture_depth_ppm']);injections.append(success)
 assert len(injections)==100;recalled=sum(injections);gate=kg and recalled>=90;gatecount+=gate;injection_count+=len(injections);control_count+=len(rr)
 cert=next(x for x in gate_certificate['per_target'] if x['ctoi']==sid);assert cert['recalled']==recalled and cert['trials']==100 and cert['controls_available']==50 and cert['gate_pass']==(recalled>=90)
 zs.sort();percentile_position=.99*49;low=math.floor(percentile_position);q99=zs[low]+(percentile_position-low)*(zs[low+1]-zs[low]);close(q99,row['control_q99'])
 f=m['selected'];dg=depth_ok(f,s);match=bool(m['matches']);ephem=all(s[k] is not None for k in ['epoch_bjd','epoch_error_days','period_error_days'])
 rational_q99=q(zs[48])+Fraction(51,100)*(q(zs[49])-q(zs[48]))
 recovered=gate and q(f['z'])>=q('7.1') and dg and match and q(f['z'])>rational_q99
 coverage_pass=cv['nominal_covered']*10>=9*cv['nominal_total'] and all(a*10>=9*b for a,b in zip(cv['grid_covered'],cv['grid_total']))
 refuted=gate and ephem and coverage_pass and q(r['expected_snr'])>=15 and not match
 verdict='RECOVERED' if recovered else 'REFUTED' if refuted else 'OPEN';assert row['verdict']==verdict,(sid,row['verdict'],verdict);counts[verdict]+=1
 for key,value in [('z',f['z']),('depth_ppm',f['depth']),('depth_sigma_ppm',f['sigma']),('expected_depth_ppm',s['expected_kepler_aperture_depth_ppm']),('nominal_coverage',cv['nominal_fraction']),('grid_min_coverage',cv['grid_min_fraction']),('expected_snr',r['expected_snr'])]:close(row[key],value)
 assert row['nominal_covered']==cv['nominal_covered'] and row['nominal_total']==cv['nominal_total'] and row['matching_events']==len(m['matches']) and row['ephemeris_refutation_eligible']==ephem and row['known_gate_pass']==kg
 assert row['injection_successes']==recalled and row['injection_trials']==100 and row['gate_pass']==gate and row['control_count']==50 and row['control_passing_count']==passing
 assert row['empirical_rank_numerator']==1+sum(v>=f['z'] for v in zs) and row['empirical_rank_denominator']==51
 if s['submitter']=='montalto':
  source=diamante[s['tic']];assert row['diamante_binary_flag']==int(source['binflag']);close(s['duration_hours'],24*float(source['duration']))
 else:assert row['diamante_binary_flag'] is None
 assert row['boyajian_star']==(s['kic']==8462852)
 close(row['recall'],recalled/100)
assert score['recovered']==counts['RECOVERED'] and score['refuted']==counts['REFUTED'] and score['open']==counts['OPEN'] and score['total']==12 and score['gate_passing']==gatecount
out={'status':'VERIFIED_VERDICTS','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'plan_sha256':freeze['plan_sha256'],'score_sha256':sha('results/scores.json'),'verifier_amendment_sha256':sha('audit/verifier_amendment.json'),'effective_verifier_sha256':sha('verify.py'),'registered_verifier_sha256':amendment['registered_verifier_sha256'],'verdicts_verified':12,'injections_verified':injection_count,'controls_verified':control_count,'known_planets_verified':5,'coverage_recomputed':full_coverage,'exact_rational_verdict_inequalities':True,'predicted_transits_checked_across_grid':coverage_events,'counts':counts}
destination='results/verification.json' if full_coverage else 'results/verification_arithmetic.json'
(R/destination).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
