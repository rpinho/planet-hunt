import sys, os, time, warnings
warnings.filterwarnings("ignore")
from astroquery.skyview import SkyView
from astropy.coordinates import SkyCoord
import astropy.units as u
outdir = sys.argv[1]
c = SkyCoord(324.3895, -77.9792, unit='deg')
surveys = sys.argv[2].split(',')
for s in surveys:
    fn = os.path.join(outdir, 'sv_'+s.replace(' ','_').replace('/','_')+'.fits')
    if os.path.exists(fn): print('have',s); continue
    try:
        hl = SkyView.get_images(position=c, survey=[s], pixels=300, radius=5*u.arcmin, cache=False)
        hl[0].writeto(fn, overwrite=True); print('ok', s)
    except Exception as e: print('FAIL', s, repr(e)[:200])
    time.sleep(1)
