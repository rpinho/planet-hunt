import sys, os, warnings
warnings.filterwarnings("ignore")
from astroquery.utils.tap.core import TapPlus
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
outdir=sys.argv[1]
ra0,de0=324.38949996727433,-77.97915951715359
tap = TapPlus(url="https://gea.esac.esa.int/tap-server/tap", verbose=False)
t=tap.launch_job(f"""SELECT source_id,ra,dec,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,bp_rp,parallax,pmra,pmdec,ruwe
 FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra0},{de0},{260/3600.}))""").get_results()
t.write(os.path.join(outdir,'gaia_260arcsec.ecsv'),overwrite=True); print('gaia',len(t))
v=Vizier(columns=['*','Date','JD'],row_limit=-1)
r=v.query_region(SkyCoord(ra0,de0,unit='deg'),radius=260*u.arcsec,catalog='II/246/out')[0]
r.write(os.path.join(outdir,'2mass_psc_260arcsec.ecsv'),overwrite=True); print('2mass',len(r)); print(r.colnames)
c=SkyCoord(r['RAJ2000'],r['DEJ2000'],unit='deg'); i=c.separation(SkyCoord(ra0,de0,unit='deg')).arcsec.argmin()
print(r[i])
