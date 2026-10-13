import sys, os, io, csv, urllib.request, warnings
warnings.filterwarnings("ignore")
from astropy.io.votable import parse_single_table
outdir = sys.argv[1]
def get(url, fn):
    p=os.path.join(outdir,fn)
    if os.path.exists(p): print('have',fn); return
    req = urllib.request.Request(url, headers={'User-Agent':'astro-amateur/1.0'})
    with urllib.request.urlopen(req, timeout=300) as r: d=r.read()
    open(p,'wb').write(d); print('ok',fn,len(d))
t = parse_single_table(os.path.join(outdir,'2mass_sia.xml')).to_table()
for r in t:
    print(r['name'], r['crpix'], r['naxis'], r['date'], r['seesh'])
for r in t:
    nm=str(r['name'])
    if ('000716 s 072 0267' in nm):
        b=nm.split('All-Sky Release Survey ')[1][0]
        get(str(r['download']), f'2mass_{b}_000716s072_0267.fits.gz')
rows = list(csv.DictReader(open(os.path.join(outdir,'skymapper_dr4_list.csv'))))
rows.sort(key=lambda r: float(r['mjd_obs']))
print([(r['band'],r['exptime'],r['mjd_obs'],r['image_type']) for r in rows[:8]])
want=[]
first_r=[r for r in rows if r['band'] in ('r','i')][0]; want.append(first_r)
main_i=[r for r in rows if r['band']=='i' and r['image_type']=='main'][0]; want.append(main_i)
for r in want:
    url=r['get_fits'].replace('size=0.05,0.05','size=0.06,0.06')
    get(url, f"skymapper_{r['band']}_{r['unique_image_id']}.fits")
