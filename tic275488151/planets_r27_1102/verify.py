"""Independent verdict verifier. Standard library only; no measurement/scoring imports.

Checks every gate against saved sufficient statistics, injection trial localization,
control order statistics, catalog dilution arithmetic, hashes and full denominators.
It does not claim to independently re-extract Kepler photometry.
"""
import json,pathlib,hashlib,math,sys,datetime,tarfile
R=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else pathlib.Path(__file__).resolve().parent
ARCHIVED=None
def raw(name):
 global ARCHIVED
 p=R/name
 if p.exists():return p.read_bytes()
 if ARCHIVED is None:
  with tarfile.open(R/'results/measurement_outputs.tar.gz','r:gz') as a:
   ARCHIVED={m.name:a.extractfile(m).read() for m in a.getmembers() if m.isfile()}
 return ARCHIVED[str(name)]
def readjson(name):return json.loads(raw(name))
P=readjson('PLAN.json');S=readjson('results/scores.json');H=hashlib.sha256(raw('PLAN.json')).hexdigest();assert S['plan_sha256']==H
manifest=readjson('results/output_manifest.json')
assert manifest['status']=='SEALED_BEFORE_SCORE' and manifest['plan_sha256']==H
assert (manifest['available_targets'],manifest['controls'],manifest['injections'])==(13,650,1300)
for name,record in manifest['files'].items():
 b=raw(name);assert len(b)==record['bytes'] and hashlib.sha256(b).hexdigest()==record['sha256'],name
coverage_certificate=readjson('results/coverage_verification.json')
assert coverage_certificate['status']=='PASS' and coverage_certificate['targets_checked']==13 and coverage_certificate['plan_sha256']==H
assert readjson('results/freeze_receipt.json')['plan_sha256']==H

def depth_ok(row,target):
 difference=abs(row['depth']-target['depth_ppm']*target['kepler_tess_dilution_ratio'])
 uncertainty=math.sqrt(row['sigma']**2+((target['depth_error_ppm'] or 0)*target['kepler_tess_dilution_ratio'])**2)
 return difference<=3*uncertainty

def close(a,b):return abs(a-b)<=1e-8*max(1.,abs(a),abs(b))

counts={'RECOVERED_LEAD':0,'REFUTED':0,'OPEN':0};checked=[]
for target,score in zip(P['targets'],S['rows']):
 assert target['ctoi']==score['ctoi'];c=target['ctoi'];T=readjson(f'results/{c}_target.json');assert T['plan_sha256']==H
 if not T['available']:
  assert score['verdict']=='OPEN' and len(T['injections'])==0;counts['OPEN']+=1;checked.append(c);continue
 assert close(target['kepler_tess_dilution_ratio'],target['kepler_crowdsap_median_q1_q17']*(1+target['tic_contratio']))
 successes=0
 assert len(T['injections'])==100
 for i in T['injections']:
  period=i['injected_period_days'];e=i['injected_epoch_bkjd'];found=i['found_epoch_bkjd']
  if period is None:phase=abs(found-e)
  else:
   true_mid=e+round(((131.5+1591.)/2-e)/period)*period
   found_period=i['found_period_days']
   phase=abs((true_mid-found+found_period/2)%found_period-found_period/2)
  located=phase<=max(target['duration_hours']/24,0.02043365) and (period is None or abs(i['found_period_days']-period)<=i['period_tolerance_days']+1e-10)
  assert close(phase,i['phase_distance_days']) and located==i['located']
  ok=i['z']>=7.1 and depth_ok(i,target) and located;assert ok==i['success'];successes+=ok
 assert successes==score['injection_successes'];recall=successes/100
 assert score['injection_total']==100 and close(score['recall'],recall)
 controls=[readjson(f"results/{c}_control_{x['kic']}.json") for x in target['controls']];assert len(controls)==50 and all(x['available'] and x['plan_sha256']==H for x in controls)
 values=sorted(x['measurement']['best']['z'] for x in controls);percentile=values[48]+.51*(values[49]-values[48]);assert close(percentile,score['control_p99_z'])
 falsepasses=sum(any(m['z']>=7.1 and depth_ok(m,target) for m in x['measurement']['matches']) for x in controls);assert falsepasses==score['controls_passing']
 m=T['measurement']['selected'];matches=[x for x in T['measurement']['matches'] if x['z']>=7.1 and depth_ok(x,target)];assert len(matches)==score['matching_events'];assert close(m['z'],score['z']);assert depth_ok(m,target)==score['depth_consistent']
 assert score['controls_total']==50 and close(score['control_rank_p'],(1+sum(v>=m['z'] for v in values))/51)
 assert close(score['depth_ppm'],m['depth']) and close(score['depth_sigma_ppm'],m['sigma'])
 assert close(score['expected_depth_ppm'],target['expected_kepler_aperture_depth_ppm'])
 assert close(score['expected_snr'],T['expected_snr']) and score['coverage']==T['coverage']
 for model in [T['measurement']['best'],m]+T['measurement']['matches']+T['measurement'].get('proposals',[]):
  assert close(model['z'],model['depth']/model['sigma'])
  assert close(model['norm'],1/model['sigma']**2)
 for control in controls:
  assert control['measurement']['best']['z']>=max((x['z'] for x in control['measurement']['matches']),default=-math.inf)-1e-9
 if target['period_days'] is not None:
  epoch=target['epoch_bjd']-2454833 if target['epoch_bjd'] is not None else (131.5+1591.)/2
  nominal=[x for x in T['measurement']['proposals'] if x['period_days']==target['period_days'] and x['epoch_bkjd']==epoch]
  assert nominal and close(T['expected_snr'],target['expected_kepler_aperture_depth_ppm']*math.sqrt(nominal[0]['norm']))
 coverage=T['coverage'];assert close(coverage['fraction'],coverage['covered']/coverage['total'])
 if target['period_days'] is None:
  recovered=len(matches)>0;refuted=not recovered and coverage['fraction']>=.9 and recall>=.9
 else:
  recovered=m['z']>=7.1 and depth_ok(m,target) and m['z']>percentile
  known=all(target[k] is not None for k in ['epoch_bjd','epoch_error_days','period_error_days']) and not T['measurement']['grid_capped'];assert known==score['ephemeris_complete']
  refuted=not recovered and not matches and coverage['nominal_fraction']>=.9 and coverage['grid_min_fraction']>=.9 and T['expected_snr']>=15 and recall>=.9 and known
 verdict='RECOVERED_LEAD' if recovered else 'REFUTED' if refuted else 'OPEN';assert verdict==score['verdict'],(c,verdict,score['verdict']);counts[verdict]+=1;checked.append(c)
assert len(checked)==len(P['targets'])==S['total']==15
assert counts=={'RECOVERED_LEAD':S['recovered'],'REFUTED':S['refuted'],'OPEN':S['open']}
result={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'VERIFIED_VERDICTS','scope':'Independent gate arithmetic and hashes from saved measurements; coverage independently recomputed by verify_coverage.py; not independent flux extraction','plan_sha256':H,'targets_checked':len(checked),'counts':counts,'measurement_code_imported':False,'sealed_files_checked':len(manifest['files']),'independent_coverage_targets':13}
(R/'results/verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
