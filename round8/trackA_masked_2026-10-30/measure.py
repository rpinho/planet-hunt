"""Continue the frozen estimator. No discovery/ephemeris refit on hold-out data."""
import bootstrap
from pathlib import Path
import json,hashlib,datetime,os,sys,traceback
import numpy as np
import legacy_photometry as p
from disk_guard import space
ROOT=Path.cwd();H=p.H

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,value):Path(path).write_text(json.dumps(value,indent=2))
def mark_pid():
 with Path('runtime/pids.jsonl').open('a') as f:f.write(json.dumps(dict(pid=os.getpid(),script=sys.argv[0],utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))+'\n')
def decorate(raw):
 for d in raw:
  start,end,orb,source=p.bounds(d);d.update(start=start,end=end,orb=orb+1000*d['sector']);ap=d['product']=='QLP' and 94<=d['sector']<=104
  for lab,b,e in [('A',1.,.5),('B',1.5,1.)]:d['edge'+lab]=ap&((d['t']<start+b)|(d['t']>end-e))
 return raw

def edge_stats(a,s):
 ep=np.rint((a['t']-s['t0'])/s['P']).astype(int);inn,_=p.masks(a,s)
 if inn.dtype!=bool:inn=inn>.25
 rows=[]
 for sec in np.unique(a['sec']):
  for epoch in np.unique(ep[inn&(a['sec']==sec)]):
   m=(ep==epoch)&(a['sec']==sec);v=p.sufficient(a,s,sel=m)
   if v is None:continue
   dep,err=p.from_suff(v,s);im=m&inn
   rows.append(dict(sector=int(sec),epoch=int(epoch),depth_ppm=float(dep*1e6),err_ppm=float(err*1e6),**{'edge'+l+'_fraction':float(np.average(a['edge'+l][im],weights=a['w'][im])) for l in ['A','B']}))
 sn=np.array([(r['depth_ppm']/r['err_ppm'])**2 for r in rows]);return rows,{l:float(np.dot(sn,[r['edge'+l+'_fraction'] for r in rows])/sn.sum()) for l in ['A','B']}

def one(name,s,raw,idx,method_key,all_null=True):
 outpath=Path('results',name+'.json')
 if outpath.exists() and json.loads(outpath.read_text()).get('complete'):return json.loads(outpath.read_text())
 space();raw=decorate(raw)
 if method_key=='tic144':p.R144=p.module('rnaas144',dict(P0=s['P'],T00=s['t0'],DUR=s['duration']/24))
 result=dict(name=name,signal=s,idx=idx,method_key=method_key,sectors=sorted(d['sector'] for d in raw),raw_manifest=[{k:d[k] for k in ['sector','source','sha256','product']} for d in raw],cases={},complete=False)
 for j,lab in enumerate(['unmasked','A','B']):
  space();parts=[]
  for d in raw:
   mm=np.ones(len(d['t']),bool) if lab=='unmasked' else ~d['edge'+lab]
   ans=p.flatten(p.subset(d,mm),s,method_key);ans=p.subset(ans,~p.known(ans['t'],s,ans['cad']));parts.append(ans)
  a=p.cat(parts);r=p.stat(a,s);r['n_cadences']=len(a['t']);arrays={k:a[k] for k in ['t','r','w','cad','sec','edgeA','edgeB']};r['sector_beta']={str(d['sector']):d['beta'] for d in parts}
  if all_null or lab!='unmasked':
   seed=20261030+idx*1000+j;ns,eph,meta=p.null(a,s,seed,n=2000)
   if len(ns)!=2000:raise RuntimeError('Exactly 2000 valid nulls required')
   exc=int((ns>=r['depth_ppm']).sum());r.update(null=meta,z_null=float(r['depth_ppm']/meta['std_ppm']),null_seed=seed,exceedances=exc,p=(1+exc)/2001);arrays.update(null_ppm=ns,null_ephemerides=eph)
  if lab=='unmasked':r['per_transit'],r['edge_share']=edge_stats(a,s)
  dest=Path('results',name+'_'+lab+'.npz');np.savez_compressed(dest,**arrays);r.update(arrays=str(dest),sha256=sha(dest));result['cases'][lab]=r;save(outpath,result)
  print(name,lab,'depth',r['depth_ppm'],'z_scatter',r['z_scatter'],'z_null',r.get('z_null'),'edges',r.get('edge_share'),flush=True)
 result['complete']=True;result['bar_stays']=bool(result['cases']['A']['z_null']>=2.5 and result['cases']['unmasked']['edge_share']['A']<=.5);result['holdout_supported']=bool(result['cases']['A']['z_null']>=3);save(outpath,result);return result

def run_C():
 inputs=json.loads(Path('results/C_inputs.json').read_text());rows=[]
 for row in inputs['rows']:
  key=row['key'];idx=row['idx'];s=dict(row['signal']);raw=p.raw_for(key)
  try:
   assert sorted(d['sector'] for d in raw)==row['holdout_sectors']
   note=one('C_'+key+'_note',s,raw,idx,key)
   reproduces=abs(note['cases']['A']['z_null']-row['reference_z'])<=.1
   status=dict(key=key,note_reproduces=reproduces,note_z=note['cases']['A']['z_null'],reference_z=row['reference_z'])
   if not reproduces:status.update(stopped=True,reason='note comparator fails');rows.append(status);save('results/C_status.json',rows);continue
   if key=='tic144':
    glob=row['ancillary']['globals_input'];sg=dict(s,P=glob['P'],t0=glob['t0'],duration=glob['duration_h']);one('C_'+key+'_globals_diagnostic',sg,raw,idx,key)
   sp=dict(s,**row['pre']);pre=one('C_'+key+'_pre',sp,raw,idx,key);status.update(stopped=False,pre_z=pre['cases']['A']['z_null'],supported=pre['holdout_supported']);rows.append(status);save('results/C_status.json',rows)
  except Exception as ex:
   traceback.print_exc();rows.append(dict(key=key,status='ERROR',error=repr(ex)));save('results/C_status.json',rows)

if __name__=='__main__':
 mark_pid()
 if sys.argv[1]=='C':run_C()
