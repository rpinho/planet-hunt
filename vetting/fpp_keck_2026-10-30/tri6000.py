"""TOI-6000: single-cadence five-seed default-prior FPP with public Keck imaging."""
from pathlib import Path
import sys,os,json,hashlib,time,socket
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
import triceratops.triceratops as tr
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from net import space
W=Path('~/.worktrees/planet-hunt-r4-tic259')
OLD=W/'runs/round4_tic259';base=json.loads((OLD/'triceratops_base.json').read_text());pix=json.loads((OLD/'triceratops_pix.json').read_text())
s=next(x for x in json.loads(Path('signals.json').read_text()) if x['key']=='toi6000')
P,T,D=s['P'],s['t0'],s['duration']/24

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def phase(t):return (t-T+P/2)%P-P/2

def setup():
 field=Path('evidence/tic_field_259233660.json');j=json.loads(field.read_text());assert j['status']=='COMPLETE',j
 df=pd.DataFrame(j['data']);df['ID']=df.ID.astype(str);coords=SkyCoord(df.ra.to_numpy(float)*u.deg,df.dec.to_numpy(float)*u.deg);it=np.flatnonzero(df.ID=='259233660')[0];sep=coords[it].separation(coords).arcsec;order=np.argsort(sep);df=df.iloc[order].reset_index(drop=True);sep=sep[order]
 cols=['ID','Tmag','Vmag','GAIAmag','gaiabp','gaiarp','Jmag','Hmag','Kmag','ebv','ra','dec','mass','rad','Teff','logg','plx','disposition','duplicate_id'];df=df[cols].copy();df['sep (arcsec)']=np.round(sep,3);coords=SkyCoord(df.ra.to_numpy(float)*u.deg,df.dec.to_numpy(float)*u.deg);df['PA (E of N)']=np.round(coords[0].position_angle(coords).deg,3)
 target=tr.target.__new__(tr.target);target.ID=259233660;target.mission='TESS';target.search_radius=10;target.N_pix=22;target.stars=df;target.estimate_stellar_params(verbose=0)
 bg=Path('evidence/imaging/259233660_gaia_background.csv');bg.write_bytes((W/'259233660_gaia_background.csv').read_bytes());target.trilegal_fname=str(bg);target.trilegal_url=None;target.background_population_source='gaia'
 inv=json.loads(Path('results/inventory.json').read_text());paths=[Path(p) for p in inv['259233660'] if str(p).endswith('-s_lc.fits')];sectors=base['aperture_sectors'];aps=[];xy=[];col0=[];row0=[];images=[];apmeta=[]
 for sec in sectors:
  candidates=[p for p in paths if f'-s{sec:04d}-' in p.name];assert candidates,(sec,paths)
  p=candidates[0]
  with fits.open(p,memmap=False) as h:
   hdr=h['APERTURE'].header;im=np.array(h['APERTURE'].data);yy,xx=np.nonzero(im&2);xc=float(hdr['CRVAL1P']);yc=float(hdr['CRVAL2P']);ap=np.column_stack([xc+xx,yc+yy]);w=WCS(hdr);xx,yy=w.all_world2pix(df.ra.to_numpy(float),df.dec.to_numpy(float),0);pos=np.column_stack([xx+xc,yy+yc])
  aps.append(ap);xy.append(pos);col0.append(xc);row0.append(yc);images.append(im);apmeta.append(dict(sector=sec,source=str(p),sha256=digest(p),pixels=ap.tolist()))
 target.sectors=np.array(sectors);target.pix_coords=xy;target.col0s=col0;target.row0s=row0;target.TESS_images=images;target.calc_depths(tdepth=s['depth']*1e-6,all_ap_pixels=aps)
 comps=[]
 for row in base['stars']:
  hit=target.stars[target.stars.ID==str(row['ID'])];comps.append(dict(ID=row['ID'],old_fluxratio=row['fluxratio'],new_fluxratio=float(hit.iloc[0].fluxratio) if len(hit) else None,old_tdepth=row['tdepth'],new_tdepth=float(hit.iloc[0].tdepth) if len(hit) else None))
 meta=dict(field_sha256=digest(field),background_sha256=digest(bg),background_rows=len(pd.read_csv(bg)),apertures=apmeta,old_new=comps,old_star_count=base['n_stars_in_search'],new_star_count=len(target.stars),stars=target.stars.replace({np.nan:None}).to_dict('records'))
 return target,meta

def fold():
 # The original fold is explicitly SPOC 120-s only, unlike the mixed-cadence affected cases.
 z=np.load(OLD/'twomin_windows.npz');t,y,e,sec=[z[k] for k in ['t','y','dy','sector']];p=phase(t)
 cadence=120;sumiv=float(np.sum(1/e[abs(p)<D/2]**2));info=[dict(cadence_s=120,in_transit_inverse_variance=sumiv,n=len(t),sectors=sorted(set(map(int,sec))))]
 # Verify available raw SPOC sectors agree with the source convention and no other original cadence exists.
 oldfold=np.load(OLD/'fold_for_triceratops.npz');assert set(map(int,oldfold['sectors']))==set(map(int,sec));assert len(t)==len(oldfold['phase_d'])
 width=cadence/86400;edges=np.arange(p.min(),p.max()+width*1.01,width);bins=np.digitize(p,edges);w=1/e**2;xx=[];yy=[];ee=[]
 for k in np.unique(bins):
  m=bins==k;xx.append(np.average(p[m],weights=w[m]));yy.append(np.average(y[m],weights=w[m]));ee.append(1/np.sqrt(w[m].sum()))
 xx=np.array(xx);yy=np.array(yy);ee=np.array(ee);yy/=np.median(yy[abs(xx)>.1]);np.savez_compressed('results/tri6000_input.npz',phase_d=xx,flux=yy,err=ee)
 return xx,yy,ee,dict(selection=info,selected_cadence_s=120,source=str(OLD/'twomin_windows.npz'),source_sha256=digest(OLD/'twomin_windows.npz'),known_limitation='Only 120-s SPOC contributes to the note original FPP fold; no cadence mixing in this signal. Cadence competition checked against raw audit photometry separately.',bin_min=2,exposure_s=120,nsamples=1,period=P,epoch_BTJD=T,depth=s['depth']*1e-6,n_bins=len(xx))

def main():
 space();target,setupmeta=setup();xx,yy,ee,foldmeta=fold();allstars=target.stars.copy(deep=True)
 Path('results/tri6000_setup.json').write_text(json.dumps(dict(setup=setupmeta,fold=foldmeta,cleared_stars=pix['cleared_stars']),indent=2))
 # Runtime network forbidden: all priors and fields are local, only the shared limiter acquires public inputs.
 socket.create_connection=lambda *a,**kw:(_ for _ in ()).throw(RuntimeError('Offline FPP worker forbids network'))
 configs=({'No_imaging_baseline':([],[])} if '--baseline' in sys.argv else {'H_primary':(['evidence/imaging/keck_Hcont.csv'],['H']),'HJ_secondary':(['evidence/imaging/keck_Hcont.csv','evidence/imaging/keck_Jcont.csv'],['H','J'])})
 for config,(curves,filters) in configs.items():
  for clear in [False,True]:
   name=config+('_removed' if clear else '_retained');dest=Path('results/tri6000_'+name+'.json');target.stars=allstars.copy(deep=True)
   if clear:target.stars.loc[target.stars.ID.isin(list(map(str,pix['cleared_stars']))),'tdepth']=0
   cc=curves[0] if len(curves)==1 else curves or None;filt=filters[0] if len(filters)==1 else filters or 'TESS'
   result=dict(name=name,contrast_curves=curves,filters=filters,draws_per_run=200000,n_runs_planned=5,cleared_stars=pix['cleared_stars'] if clear else [],input=foldmeta,flux_err_scalar=float(np.median(ee)),runs=[],complete=False)
   if dest.exists():
    old=json.loads(dest.read_text())
    if old.get('complete'):print('EXISTS',name,flush=True);continue
    result['runs']=old['runs']
   for i in range(len(result['runs']),5):
    space();seed=202610300+100*clear+i;np.random.seed(seed);st=time.time();print('START',name,i,flush=True)
    target.calc_probs(time=xx,flux_0=yy,flux_err_0=float(np.median(ee)),P_orb=P,N=200000,parallel=False,exptime=120/86400,nsamples=1,verbose=0,contrast_curve_file=cc,filt=filt,molusc_file=None)
    assert not getattr(target,'FPP_degenerate',False)
    rr=dict(seed=seed,FPP=float(target.FPP),NFPP=float(target.NFPP),elapsed_s=time.time()-st,scenarios={f'{r.ID}:{r.scenario}':float(r.prob) for r in target.probs.itertuples()},lnZ=np.asarray(target.lnZ).tolist());result['runs'].append(rr);dest.write_text(json.dumps(result,indent=2));print(name,i,rr['FPP'],rr['NFPP'],rr['elapsed_s'],flush=True)
   result['complete']=True
   for k in ['FPP','NFPP']:
    ar=np.array([r[k] for r in result['runs']]);result[k+'_mean']=float(ar.mean());result[k+'_std']=float(ar.std(ddof=1))
   result['bar']='PASS' if not clear and result['FPP_mean']<.015 else 'FAIL' if not clear else 'N/A removed-neighbour result cannot validate';dest.write_text(json.dumps(result,indent=2))
if __name__=='__main__':main()
