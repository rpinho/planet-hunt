import bootstrap
from pathlib import Path
import json,time,traceback,sys
import numpy as np
from measure import one,save,mark_pid
from read_products import ready,read

def controls():
 rows=[]
 for key,ref in [('tic101',4.97),('tic262_24',3.42)]:
  r=json.loads(Path('results',key+'_round8.json').read_text());z=r['cases']['A']['z_null'];rows.append(dict(key=key,z=z,reference=ref,passes=abs(z-ref)<=.1))
 save('results/control_gate.json',dict(rows=rows,passes=all(x['passes'] for x in rows)))
 return all(x['passes'] for x in rows)
def run(poll=False):
 assert controls(),'Control gate failed; stop'
 plan=json.loads(Path('results/acquisition_plan.json').read_text());summary=[]
 for row in plan['A']:
  sid=row['id'];dest=Path('results','A_'+sid+'.json')
  try:
   if sid=='E1_101653309_0.7034' or sid=='E1_262843927_24.1087':
    key='tic101' if '101653309' in sid else 'tic262_24';r=json.loads(Path('results',key+'_round8.json').read_text());r.update(name='A_'+sid,source_control=key+'_round8.json');save(dest,r)
   elif dest.exists() and json.loads(dest.read_text()).get('complete'):r=json.loads(dest.read_text())
   elif not all(ready(p) for p in row['products']):
    summary.append(dict(id=sid,status='WAITING_FOR_PRODUCTS',missing=[p['destination'] for p in row['products'] if not ready(p)]));continue
   else:
    old=row['signal'];e=old['real']['his'];s=dict(key='round8',tic=old['tic'],P=e['P'],t0=e['t0'],duration=24*e['dur'],mask=[])
    raw=[read(pr) for pr in row['products']];r=one('A_'+sid,s,raw,row['idx'],'tic101',all_null=False);r.update(input_id=sid,input_sha256=row['input_sha256']);save(dest,r)
   summary.append(dict(id=sid,status='MEASURED',stays=r['bar_stays'],z_A=r['cases']['A']['z_null'],z_B=r['cases']['B']['z_null'],edge_share=r['cases']['unmasked']['edge_share']['A']))
  except Exception as ex:traceback.print_exc();summary.append(dict(id=sid,status='ERROR',error=repr(ex)))
 n=sum(x['status']=='MEASURED' for x in summary);k=sum(x.get('stays',False) for x in summary);save('results/A_status.json',dict(n_target=17,n_measured=n,k=(k if n==17 else None),measured_passes=k,rows=summary,complete=n==17));print('A status',n,k,flush=True)
 return n==17
if __name__=='__main__':
 mark_pid();run()
