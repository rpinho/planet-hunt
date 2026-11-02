import bootstrap
"""Frozen-ephemeris hold-outs. Original functions are imported from local AST extracts.
All outputs carry original source hashes and explicit statistic definitions.
"""
from pathlib import Path
import numpy as np
import json,hashlib,ast,sys,types,time,traceback
from astropy.io import fits
from astropy.time import Time
from disk_guard import space
ROOT=Path(__file__).resolve().parent
H=Path('../planets_holdout_1030')
AUD=Path('../planets_audit_1030')
W=Path('../worktrees')
SS={s['key']:s for s in json.loads((H/'signals.json').read_text())}
INV=json.loads((H/'results/inventory.json').read_text())

def phase(t,P,t0):return (t-t0+P/2)%P-P/2

def module(name,extra=None):
 from wotan import flatten as wf
 env={'np':np,'fits':fits,'re':__import__('re'),'flatten':wf};env.update(extra or {})
 exec(compile((H/'vendor'/(name+'_methods.py')).read_text(),name,'exec'),env)
 return types.SimpleNamespace(**env)
R144=module('rnaas144',dict(P0=SS['tic144']['P'],T00=SS['tic144']['t0'],DUR=SS['tic144']['duration']/24))
HOST=module('host');NOISE=module('noise');DEP=module('depth',{'ms':types.SimpleNamespace(phase=phase)})
GATE=module('lcgate');R6=module('round6')
src=(H/'evidence/upstream_pipeline/pipeline/common.py').read_text();tree=ast.parse(src)
for node in tree.body:
 if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='QLP_ORBIT_STARTS' for t in node.targets):START=np.array(ast.literal_eval(node.value))
edgeenv={'np':np,'QLP_ORBIT_STARTS':START,'QLP_ORBIT_MASK_DAYS':1.,'QLP_ORBIT_MASK_END_DAYS':.5,'QLP_ORBIT_MASK_SECTORS':range(94,105)}
exec(ast.get_source_segment(src,next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='mask_qlp_orbit_edges')),edgeenv)

# Scheduled physical boundaries: UTC -> BTJD_TDB; QLP grouping supersedes mislabeled S102/103.
import pandas as pd
ORB=pd.read_csv(H/'evidence/TESS_orbit_times.csv').dropna(subset=['Start of Orbit','End of Orbit'])
ORB['Sector']=pd.to_numeric(ORB['Sector'],errors='coerce')
ORB['start']=Time(ORB['Start of Orbit'].to_list(),format='iso',scale='utc').tdb.jd-2457000
ORB['end']=Time(ORB['End of Orbit'].to_list(),format='iso',scale='utc').tdb.jd-2457000
ORB=ORB.groupby(['Sector','Orbit'],as_index=False).agg(start=('start','min'),end=('end','max'))

def read(p):
 cached=H/'evidence/curve_cache'/(Path(p).name+'.npz')
 if cached.exists():
  meta=json.loads(cached.with_suffix('.json').read_text());assert hashlib.sha256(cached.read_bytes()).hexdigest()==meta['compact_sha256']
  with np.load(cached) as z:meta.update({k:z[k] for k in z.files})
  return meta
 with fits.open(p,memmap=False) as h:
  tab=h[1].data;names=tab.names;t=np.array(tab['TIME'],float)
  q=np.array(tab['QUALITY'],int) if 'QUALITY' in names else np.zeros(len(t),int)
  col=next(x for x in ['PDCSAP_FLUX','KSPSAP_FLUX','DET_FLUX'] if x in names)
  f=np.array(tab[col],float);e=np.array(tab[col+'_ERR'],float) if col+'_ERR' in names else np.full(len(t),np.nan)
  sec=int(h[0].header['SECTOR']);cad=float(h[1].header.get('TIMEDEL',np.nan))
  if not np.isfinite(cad):cad=np.nanmedian(np.diff(t))
  if np.nanmedian(t)>2e6:t-=2457000
  good=np.isfinite(t)&np.isfinite(f)&(q==0)
  med=np.median(f[good]);rt=t.copy();rq=q.copy()
  oid=np.array(tab['ORBITID'])[good] if 'ORBITID' in names else None
 return dict(t=t[good],f=f[good]/med,e=e[good]/med,cad=np.full(good.sum(),cad),sec=np.full(good.sum(),sec,int),sector=sec,product=('QLP' if 'qlp' in str(p) else 'TESS-SPOC' if 'hlsp' in str(p) else 'SPOC'),source=str(p),sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest(),empty_error=not np.any(np.isfinite(e[good])&(e[good]>0)),raw_t=rt,raw_q=rq,orbitid=oid)

def bounds(d):
 t=d['t'];sec=d['sector']
 if 94<=sec<=104:
  k=np.searchsorted(START,t+.02,side='right')-1
  first=np.full(len(START),np.inf);last=np.full(len(START),-np.inf)
  np.minimum.at(first,k[k>=0],t[k>=0]);np.maximum.at(last,k[k>=0],t[k>=0])
  kk=np.maximum(k,0);return first[kk],last[kk],k,'QLP_ORBIT_STARTS (group); first/last quality-selected cadences'
 o=ORB[ORB.Sector==sec]
 a=o.start.to_numpy();b=o.end.to_numpy()
 if not len(a):raise ValueError('No scheduled orbit '+str(sec))
 # Scheduled start identifies orbit, including cadences falling in small BJD/UTC gaps.
 k=np.clip(np.searchsorted(a,t,side='right')-1,0,len(a)-1)
 return a[k],b[k],k,'TESS official table (UTC->TDB), no QLP mask outside S94–104'

def subset(d,m):
 return {k:(v[m] if isinstance(v,np.ndarray) and v.shape==d['t'].shape and k not in ['raw_t','raw_q'] else v) for k,v in d.items()}

def raw_for(key):
 s=SS[key];ps=[Path(p) for p in INV[str(s['tic'])] if str(p).endswith('.fits')]
 if key=='tic144':ps=[p for p in ps if 'qlp' in p.name and any(f's{x:04d}-' in p.name for x in [96,103,104])]
 elif key=='tic101':ps=[p for p in ps if 'qlp' in p.name and 's0094-' in p.name]
 elif key.startswith('tic262'):ps=[p for p in ps if 'qlp' in p.name and any(f's{x:04d}-' in p.name for x in [95,96,97])]
 elif key=='toi4405':
  secs=[1,3,4,5,6,7,8,9,10,11,13];ps=[p for p in ps if 'tess-spoc' in p.name and any(f'-s{x:04d}_' in p.name for x in secs)]
 elif key=='toi4342_18':ps=[p for p in ps if 'qlp' in p.name]
 elif key=='toi4342_2':ps=[p for p in ps if p.name.startswith('tess') and p.name.endswith('-s_lc.fits')]
 out=[];seen=set()
 for p in sorted(ps):
  d=read(p)
  if d['sector'] in seen:continue
  seen.add(d['sector']);out.append(d)
 return sorted(out,key=lambda d:d['sector'])

def known(t,s,cad):
 m=np.zeros(len(t),bool)
 for p,t0,dh in s['mask']:m|=np.abs(phase(t,p,t0))<.8*dh/24+cad/2
 return m

def flatten(d,s,key):
 t=d['t'];f=d['f'];D=s['duration']/24;P=s['P'];t0=s['t0'];km=known(t,s,d['cad'])
 if key=='tic144':
  tt,r,sig=R144.detrend(t,f);ix=np.searchsorted(t,tt);ans=subset(d,np.isin(np.arange(len(t)),ix));ans.update(r=r,e=np.full(len(tt),sig),w=np.full(len(tt),1/sig**2));ans['beta']=1.;return ans
 mask=(np.abs(phase(t,P,t0))<.8*D+d['cad']/2)|km
 if key in ['tic101','tic262_14','tic262_24']:
  inp=dict(time=t,flux=f,cx=np.zeros(len(t)),cy=np.zeros(len(t)),meta={'sector':d['sector']})
  kmfn=lambda tt:known(tt,s,np.median(d['cad']))|(np.abs(phase(tt,P,t0))<.8*D)
  lc=HOST.detrend([inp],.5 if P<6 else 1.,kmfn,rot=None)
  if not lc:raise RuntimeError('Empty flatten')
  ii=np.searchsorted(t,lc['t']);ans=subset(d,np.isin(np.arange(len(t)),ii));be=NOISE.red_noise_beta(lc,D,np.ones(len(lc['t']),bool));ee=lc['dy']*be;ans.update(r=lc['flat']-1,e=ee,w=1/ee**2,beta=float(np.median(be)));return ans
 if key=='toi4405':
  mask=(np.abs(phase(t,P,t0))<.75*D+d['cad']/2)|km
  fd=GATE.flatten({'t':t,'f':f},mask,max(.5,3.5*D));tt=fd['t'];ii=np.searchsorted(t,tt);ans=subset(d,np.isin(np.arange(len(t)),ii));r=fd['r'];oot=~mask[ii];be=GATE.beta_at(tt[oot],r[oot],D,np.median(d['cad']));ee=np.full(len(tt),fd['sig']*be);ans.update(r=r,e=ee,w=1/ee**2,beta=be);return ans
 dd=R6.detrend({'t':t,'f':f},mask,window=.75);r=dd['fd']-1;mad=1.4826*np.nanmedian(np.abs(r-np.nanmedian(r)));keep=np.isfinite(r)&(r<4*mad)&((r>-8*mad)|mask);ans=subset(d,keep);rr=r[keep];ee=d['e'][keep].copy();ee[~np.isfinite(ee)|(ee<=0)]=mad
 # Same sector scatter convention as audit: formal errors at least adjacent-difference robust scatter.
 rms=1.4826*np.median(np.abs(np.diff(rr)-np.median(np.diff(rr))))/np.sqrt(2);ee*=max(1.,rms/np.median(ee));ans.update(r=rr,e=ee,w=1/ee**2,beta=float(max(1.,rms/np.median(d['e'][np.isfinite(d['e'])&(d['e']>0)]))) if np.any(np.isfinite(d['e'])&(d['e']>0)) else 1.);return ans

def cat(parts):
 out={k:np.concatenate([d[k] for d in parts]) for k in ['t','r','e','w','cad','sec','start','end','edgeA','edgeB','orb']}
 ii=np.argsort(out['t']);return {k:v[ii] for k,v in out.items()}

def masks(a,s,P=None,t0=None):
 P=P or s['P'];t0=s['t0'] if t0 is None else t0;D=s['duration']/24;x=phase(a['t'],P,t0)
 if s['key']=='toi4405':
  o=np.clip((np.minimum(x+a['cad']/2,D/2)-np.maximum(x-a['cad']/2,-D/2))/a['cad'],0,1)
  use=np.abs(x)<np.maximum(4*D,D+2*np.median(a['cad']));return o,use
 inn=np.abs(x)<D/2;oot=(np.abs(x)>D)&(np.abs(x)<(.5*P-D if s['key']=='tic144' else 3*D));return inn,oot

def sufficient(a,s,P=None,t0=None,sel=None):
 if sel is None:sel=np.ones(len(a['t']),bool)
 x,oo=masks(a,s,P,t0);w=a['w']
 if s['key']=='toi4405':
  m=sel&oo
  if ((x>.25)&m).sum()<3 or m.sum()<20:return None
  return np.array([w[m].sum(),(w[m]*x[m]).sum(),(w[m]*x[m]**2).sum(),(w[m]*a['r'][m]).sum(),(w[m]*a['r'][m]*x[m]).sum()])
 inn=sel&x;oot=sel&oo
 if inn.sum()<3 or oot.sum()<10:return None
 if s['key']=='tic144':return np.array([inn.sum(),a['r'][inn].sum(),oot.sum(),a['r'][oot].sum(),np.sum(a['r'][oot]**2)])
 return np.array([w[inn].sum(),np.sum(w[inn]*a['r'][inn]),w[oot].sum(),np.sum(w[oot]*a['r'][oot]),0.])

def from_suff(v,s):
 if v is None:return np.nan,np.nan
 if s['key']=='toi4405':
  sw,sx,sxx,sy,sxy=v;det=sw*sxx-sx*sx
  if det<=0:return np.nan,np.nan
  return -(sw*sxy-sx*sy)/det,np.sqrt(sw/det)
 wi,yi,wo,yo,_=v
 if wi<=0 or wo<=0:return np.nan,np.nan
 dep=yo/wo-yi/wi
 if s['key']=='tic144':
  sig=np.sqrt(max(0.,v[4]/wo-(yo/wo)**2));return dep,sig*np.sqrt(1/wi+1/wo)
 return dep,np.sqrt(1/wi+1/wo)

def stat(a,s,sel=None,P=None,t0=None):
 d,e=from_suff(sufficient(a,s,P,t0,sel),s)
 if s['key']=='tic144':
  sel=np.ones(len(a['t']),bool) if sel is None else sel
  if sel.sum()>=10:d,e,_,_=R144.box_depth(a['t'],a['r'],P=P or s['P'],T0=s['t0'] if t0 is None else t0,dur=s['duration']/24,sel=sel)
 return dict(depth_ppm=float(d*1e6),err_ppm=float(e*1e6),z_scatter=float(d/e))

def null(a,s,seed,n=2000):
 rng=np.random.default_rng(seed);clean=~known(a['t'],s,a['cad']);clean&=np.abs(phase(a['t'],s['P'],s['t0']))>(.8*s['duration']/24+a['cad']/2)
 aa={k:v[clean] for k,v in a.items()};ds=[];prs=[];attempts=0
 while len(ds)<n and attempts<n*20:
  attempts+=1;P=s['P']*rng.uniform(.8,1.2)
  if abs(P/s['P']-1)<.02:continue
  t0=float(aa['t'].min()+rng.uniform(0,P));d,e=from_suff(sufficient(aa,s,P,t0),s)
  if np.isfinite(d):ds.append(d*1e6);prs.append([P,t0])
 ds=np.array(ds)
 if len(ds)<1000:raise RuntimeError('Insufficient valid nulls '+str(len(ds)))
 return ds,np.array(prs),dict(n=len(ds),attempts=attempts,mean_ppm=float(ds.mean()),std_ppm=float(ds.std(ddof=1)),q025_975_ppm=np.quantile(ds,[.025,.975]).tolist(),clean_cadences=int(clean.sum()),removed_real_transit_cadences=int((~clean).sum()))

def transit_boot(a,s,seed):
 ep=np.rint((a['t']-s['t0'])/s['P']).astype(int);x,oo=masks(a,s);ints=x>.25 if s['key']=='toi4405' else x
 rows=[];suffs=[]
 for sec in np.unique(a['sec']):
  for epoch in np.unique(ep[ints&(a['sec']==sec)]):
   m=(ep==epoch)&(a['sec']==sec);v=sufficient(a,s,sel=m)
   if v is None:continue
   d,e=from_suff(v,s);suffs.append(v)
   forlab={};im=m&ints;wc=a['w'][im]
   for lab in ['A','B']:
    forlab['edge'+lab+'_fraction']=float(np.average(a['edge'+lab][im],weights=wc))
   rows.append(dict(sector=int(sec),epoch=int(epoch),tc=float(s['t0']+epoch*s['P']),depth_ppm=float(d*1e6),err_ppm=float(e*1e6),n_in=int(im.sum()),start=float(np.median(a['start'][im])),end=float(np.median(a['end'][im])),since_start_d=float(np.median(a['t'][im]-a['start'][im])),**forlab))
 vv=np.array(suffs);rng=np.random.default_rng(seed)
 if len(vv)==0:return rows,np.array([]),{}
 boot=np.array([from_suff(vv[rng.integers(0,len(vv),len(vv))].sum(0),s)[0]*1e6 for _ in range(2000)])
 return rows,boot,dict(n_transit_clusters=len(rows),n_boot=len(boot),depth_err_ppm=float(boot.std(ddof=1)),ci95_ppm=np.quantile(boot,[.025,.975]).tolist(),bootstrap_aggregate_ppm=float(from_suff(vv.sum(0),s)[0]*1e6))

def edge_diag(a,s,rows,seed):
 if not rows:return {}
 out={};sn=np.array([(r['depth_ppm']/r['err_ppm'])**2 for r in rows]);rng=np.random.default_rng(seed)
 for lab,(before,after) in {'A':(1.,.5),'B':(1.5,1.)}.items():
  frac=np.array([r['edge'+lab+'_fraction'] for r in rows]);share=float(np.dot(sn,frac)/sn.sum());controls=[];half=[]
  for _ in range(2000):
   shifted=np.array([r['tc'] for r in rows])+rng.uniform(-6,6,len(rows));fs=[]
   for tc,r in zip(shifted,rows):fs.append(float(tc<r['start']+before or tc>r['end']-after))
   controls.append(float(np.dot(sn,fs)/sn.sum()))
  for r in rows:
   tc=r['tc']+.5*(r['end']-r['start']);tc=r['start']+(tc-r['start'])%(r['end']-r['start']);half.append(float(tc<r['start']+before or tc>r['end']-after))
  x,oo=masks(a,s);im=x>.25 if s['key']=='toi4405' else x;ed=a['edge'+lab]
  out[lab]=dict(snr2_share=share,shifted_control_mean=float(np.mean(controls)),shifted_control_q025_975=np.quantile(controls,[.025,.975]).tolist(),half_orbit_shift_share=float(np.dot(sn,half)/sn.sum()),in_cadence_edge_fraction=float(ed[im].mean()),in_inverse_variance_edge_share=float(a['w'][im&ed].sum()/a['w'][im].sum()),definition='Fallback: sum per-transit (depth/error)^2 weighted by edge in-transit inverse variance; public pilot script unavailable. Shifted event midpoint ±6 d may leave orbit; outside is counted edge (documented conservative control).')
 return out

def run(key,idx):
 space();s=SS[key];raw=raw_for(key)
 if not raw:raise RuntimeError('No raw curves '+key)
 print(key,'RAW',[(d['sector'],len(d['t'])) for d in raw],flush=True)
 for d in raw:
  start,end,orb,source=bounds(d);d.update(start=start,end=end,orb=orb+1000*d['sector']);ap=(d['product']=='QLP' and 94<=d['sector']<=104)
  d['edgeA']=ap&((d['t']<start+1)|(d['t']>end-.5));d['edgeB']=ap&((d['t']<start+1.5)|(d['t']>end-1));d['orbit_source']=source
  # Exclude known transit cadences only from fitting; retain original raw counts for accounting.
 manifest=[{k:d[k] for k in ['sector','product','source','sha256','empty_error','orbit_source']}|{'n_cadences':len(d['t']),'cadence_s':float(np.median(d['cad'])*86400)} for d in raw]
 results=dict(key=key,signal=s,manifest=manifest,methods=json.loads(Path('evidence/method_sources.json').read_text()),cases={},complete=False)
 outpath=Path('results',key+'_holdout.json')
 nin0=sum((np.abs(phase(d['t'],s['P'],s['t0']))<s['duration']/48).sum() for d in raw)
 for j,lab in enumerate(['unmasked','A','B']):
  parts=[];nremove=0
  for d in raw:
   mm=np.ones(len(d['t']),bool) if lab=='unmasked' else ~d['edge'+lab]
   inn=np.abs(phase(d['t'],s['P'],s['t0']))<s['duration']/48;nremove+=int((inn&~mm).sum())
   dd=subset(d,mm);ans=flatten(dd,s,key);km=known(ans['t'],s,ans['cad']);ans=subset(ans,~km);parts.append(ans)
  a=cat(parts);r=stat(a,s);r['n_cadences']=len(a['t']);r['in_transit_raw_before']=int(nin0);r['in_transit_removed_mask']=nremove;r['in_transit_removed_fraction']=nremove/nin0;r['sectors']=sorted(set(map(int,a['sec'])));r['sector_beta']={str(d['sector']):d['beta'] for d in parts}
  rows,boot,bmeta=transit_boot(a,s,20261030+idx*100+j);r['transit_bootstrap']=bmeta;r['per_transit']=rows
  arrpath=Path('results',f'{key}_{lab}.npz');arrays=dict(a,bootstrap_ppm=boot)
  if lab in ['A','B']:
   ns,npms,nmeta=null(a,s,20261030+idx*1000+j);r['null']=nmeta;r['z_null']=r['depth_ppm']/nmeta['std_ppm'];r['null']['p_ge_observed']=(1+int((ns>=r['depth_ppm']).sum()))/(len(ns)+1);arrays.update(null_ppm=ns,null_ephemerides=npms)
   jack=[]
   if len(r['sectors'])>1:
    for sec in r['sectors']:
     aa={k:v[a['sec']!=sec] for k,v in a.items()};jr=stat(aa,s)
     try:
      jns,jnp,jmeta=null(aa,s,20261030+idx*1000+j*100+sec);jr.update(null=jmeta,z_null=jr['depth_ppm']/jmeta['std_ppm']);arrays['jack_null_'+str(sec)]=jns;arrays['jack_eph_'+str(sec)]=jnp
     except Exception as ex:jr.update(status='UNVERIFIED',error=repr(ex))
     jack.append(dict(omitted_sector=sec,**jr))
   r['jackknife']=jack;r['jackknife_status']='N/A single sector' if len(r['sectors'])==1 else 'computed'
  else:
   r['edge_attribution']=edge_diag(a,s,rows,20261030+idx)
   bins=np.array([0,.5,1,1.5,2,4,6,8,10,12,14,16]);profile=[]
   for lo,hi in zip(bins[:-1],bins[1:]):
    m=(a['t']-a['start']>=lo)&(a['t']-a['start']<hi)
    profile.append(dict(lo_d=lo,hi_d=hi,n_cadences=int(m.sum()),**stat(a,s,sel=m)))
   r['depth_by_orbit_age']=profile
  np.savez_compressed(arrpath,**arrays);results['cases'][lab]=r;outpath.write_text(json.dumps(results,indent=2));print(key,lab,{k:r.get(k) for k in ['depth_ppm','err_ppm','z_scatter','z_null','in_transit_removed_fraction']},flush=True)
 results['complete']=True;results['headline_bar']='FAIL' if results['cases']['A']['z_null']<3 else 'PASS';results['edge_bar']='FAIL' if results['cases']['unmasked']['edge_attribution']['B']['snr2_share']>.5 else 'PASS';outpath.write_text(json.dumps(results,indent=2))

