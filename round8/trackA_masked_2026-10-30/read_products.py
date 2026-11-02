import bootstrap
from pathlib import Path
import json,hashlib
import numpy as np
import legacy_photometry as p

def hashfile(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def ready(prod):
 if prod['cached']:return Path(prod['cached']['path']).exists()
 dest=Path(prod['destination']);rec=Path(prod['destination']+'.receipt.json')
 if not dest.exists() or not rec.exists():return False
 try:j=json.loads(rec.read_text())
 except (OSError,ValueError):return False
 return j.get('returncode')==0 and j.get('curl',{}).get('http_code') in [200,206] and dest.stat().st_size>0
def read(prod,original_float32=False):
 c=prod['cached'];path=Path(c['path'] if c else prod['destination'])
 if c and c.get('multi'):
  with np.load(path) as z:
   i=c['index'];meta=json.loads(z[f's{i}_meta'].item());t=z[f's{i}_time'].copy();f=z[f's{i}_flux'].copy();err=z[f's{i}_ferr'].copy();cx=z[f's{i}_cx'].copy();cy=z[f's{i}_cy'].copy()
  cad=np.median(np.diff(t));sec=prod['sector'];d=dict(t=t,f=f,e=err,cad=np.full(len(t),cad),sec=np.full(len(t),sec,int),sector=sec,product=prod['product'],source=str(path),sha256=hashfile(path),empty_error=False,cx=cx,cy=cy)
 elif c and c.get('compact'):
  d=json.loads(path.with_suffix('.json').read_text())
  assert hashfile(path)==d['compact_sha256']
  with np.load(path) as z:d.update({k:z[k] for k in z.files})
 else:
  d=p.read(path)
  # The legacy hold-out estimator has no centroid dependence. Inverted original code uses rotations from flux only.
  from astropy.io import fits
  with fits.open(path,memmap=False) as h:
   tab=h[1].data;t=np.array(tab['TIME'],float);q=np.array(tab['QUALITY'],int);col=next(k for k in ['PDCSAP_FLUX','KSPSAP_FLUX','DET_FLUX'] if k in tab.names);f=np.array(tab[col],float);m=np.isfinite(t)&np.isfinite(f)&(q==0)
   for k,cols in [('cx',['MOM_CENTR1','SAP_X']),('cy',['MOM_CENTR2','SAP_Y'])]:
    col=next((k for k in cols if k in tab.names),None);d[k]=np.array(tab[col],float)[m].astype(np.float32) if col else np.zeros(len(d['t']))
 assert d['sector']==prod['sector'],('Sector mismatch',d['sector'],prod['sector'])
 assert ('SPOC-2min' if d['product']=='SPOC' else d['product'])==prod['product'],('Product mismatch',d['product'],prod['product'])
 if original_float32:
  d['f']=d['f'].astype(np.float32);d['e']=d['e'].astype(np.float32)
  if 'cx' not in d:d['cx']=np.zeros(len(d['t']),np.float32);d['cy']=np.zeros(len(d['t']),np.float32)
 return d
