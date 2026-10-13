"""Download every TESS light-curve product + 2-min TPFs for TIC 354944123."""
import sys, json
sys.path.insert(0, '.')
from planet_hunt import config
from pathlib import Path
import urllib.request, time
from astroquery.mast import Observations, conf
conf.timeout = 600
TIC = 354944123
OUT = config.CACHE / 'r6' / 'raw'
OUT.mkdir(parents=True, exist_ok=True)
obs = Observations.query_criteria(target_name=str(TIC), obs_collection=['TESS', 'HLSP'])
prods = Observations.get_product_list(obs)
df = prods.to_pandas()
keep = []
for _, r in df.iterrows():
    fn = str(r['productFilename'])
    if not fn.endswith('.fits'):
        continue
    if fn.endswith('_lc.fits') or fn.endswith('_fast-lc.fits') or fn.endswith('_tp.fits') and 'tess20' in fn \
       or fn.endswith('_fast-tp.fits') is False and fn.endswith('_tp.fits') and 'tess20' in fn \
       or ('hlsp_qlp' in fn and fn.endswith('_llc.fits')) or ('hlsp_tess-spoc' in fn and fn.endswith('_lc.fits')):
        keep.append((fn, r['dataURI']))
keep = sorted(set(keep))
print(len(keep), 'files')
for fn, uri in keep:
    p = OUT / fn
    if p.exists() and p.stat().st_size > 0:
        continue
    url = config.MAST_DOWNLOAD + uri
    for k in range(4):
        try:
            urllib.request.urlretrieve(url, p); break
        except Exception as e:
            print('retry', fn, e); time.sleep(5)
    print(fn, p.stat().st_size // 1024, 'kB', flush=True)
json.dump([k for k, _ in keep], open('runs/round6_toi4342/data/files.json', 'w'), indent=1)
