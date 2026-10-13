import sys, os, io, time, csv, urllib.request, warnings
warnings.filterwarnings("ignore")
from astropy.io.votable import parse_single_table
outdir = sys.argv[1]
def fetch(url):
    req = urllib.request.Request(url, headers={'User-Agent':'astro-amateur/1.0'})
    with urllib.request.urlopen(req, timeout=300) as r: return r.read()
# SkyMapper DR4 list
rows = list(csv.DictReader(io.StringIO(fetch("https://api.skymapper.nci.org.au/public/siap/dr4/query?POS=324.3895,-77.9792&SIZE=0.05&BAND=all&FORMAT=image/fits&RESPONSEFORMAT=CSV").decode())))
print(len(rows), 'SkyMapper images')
for r in rows: print(r['band'], r['exptime'], r['mjd_obs'], r['image_type'], r['unique_image_id'])
open(os.path.join(outdir,'skymapper_dr4_list.csv'),'w').write(fetch("https://api.skymapper.nci.org.au/public/siap/dr4/query?POS=324.3895,-77.9792&SIZE=0.05&BAND=all&FORMAT=image/fits&RESPONSEFORMAT=CSV").decode())
# 2MASS
vo = fetch("https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?POS=324.3895,-77.9792&SIZE=0.02&FORMAT=image/fits")
open(os.path.join(outdir,'2mass_sia.xml'),'wb').write(vo)
t = parse_single_table(io.BytesIO(vo)).to_table()
print(t.colnames)
print(t)
