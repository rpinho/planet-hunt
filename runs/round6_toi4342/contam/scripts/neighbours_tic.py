import sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.utils.tap.core import TapPlus
from astroquery.mast import Catalogs
nbfile, out = sys.argv[1], sys.argv[2]
sid=6355915466181029376
nb = Table.read(nbfile); nb.rename_columns(nb.colnames,[c.lower() for c in nb.colnames])
t0 = nb[nb['source_id']==sid][0]
def tmag(G,c):
    if c is None or np.ma.is_masked(c) or not np.isfinite(c): return G-0.43
    return G-0.00522555*c**3+0.0891337*c**2-0.633923*c+0.0324473
T0 = tmag(float(t0['phot_g_mean_mag']), float(t0['bp_rp']))
c0 = SkyCoord(t0['ra']*u.deg, t0['dec']*u.deg)
rows=[]
for r in nb:
    if r['source_id']==sid: continue
    G=float(r['phot_g_mean_mag']) if not np.ma.is_masked(r['phot_g_mean_mag']) else np.nan
    c=None if np.ma.is_masked(r['bp_rp']) else float(r['bp_rp'])
    T=tmag(G,c); dT=T-T0
    f=10**(-0.4*dT)
    cc=SkyCoord(r['ra']*u.deg,r['dec']*u.deg)
    def fl(x): return None if np.ma.is_masked(x) else float(x)
    d={'source_id':int(r['source_id']),'sep_arcsec':round(float(r['sep']),2),'pa_deg':round(float(c0.position_angle(cc).deg),1),
       'G':fl(r['phot_g_mean_mag']),'bp_rp':c,'Tmag_est':round(T,3),'dTmag':round(dT,3),'parallax':fl(r['parallax']),'parallax_error':fl(r['parallax_error']),
       'pmra':fl(r['pmra']),'pmdec':fl(r['pmdec']),'ruwe':fl(r['ruwe']),
       'depth_needed_1150ppm':1150e-6*(1+f)/f,'depth_needed_284ppm':284e-6*(1+f)/f}
    d['mimic_1150']=d['depth_needed_1150ppm']<1; d['mimic_284']=d['depth_needed_284ppm']<1
    rows.append(d)
rows.sort(key=lambda d:d['sep_arcsec'])
# TIC
tic = Catalogs.query_object("TIC 354944123", radius=0.001, catalog="TIC")
tt = tic[tic['ID']=='354944123'][0]
ticd={k:(None if np.ma.is_masked(tt[k]) else (float(tt[k]) if isinstance(tt[k],(float,np.floating)) else str(tt[k]))) for k in ['ID','Tmag','e_Tmag','contratio','numcont','Teff','rad','mass','d','disposition','GAIA','ra','dec','pmRA','pmDEC','priority']}
ticn = Catalogs.query_region(SkyCoord(324.3895,-77.9792,unit='deg'), radius=150*u.arcsec, catalog="TIC")
ticnb=[{'ID':str(x['ID']),'Tmag':float(x['Tmag']),'dstArcSec':float(x['dstArcSec']),'GAIA':str(x['GAIA']),'objType':str(x['objType'])} for x in ticn]
json.dump({'T0':T0,'neighbours':rows,'tic':ticd,'tic_within150':ticnb},open(out,'w'),indent=1)
print('T0',T0, ticd)
print(len(rows),'Gaia nbrs; TIC within 150":',len(ticnb))
for d in rows: print(d['sep_arcsec'],d['pa_deg'],d['G'],d['bp_rp'],round(d['Tmag_est'],2),round(d['dTmag'],2),d['parallax'],d['pmra'],d['pmdec'],'%.3f'%d['depth_needed_1150ppm'],'%.3f'%d['depth_needed_284ppm'])
