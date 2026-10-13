import sys, os, time, urllib.request, urllib.parse
outdir = sys.argv[1]
ra, dec = 324.3895, -77.9792
size_arcmin = 6
def get(url, fn):
    p = os.path.join(outdir, fn)
    if os.path.exists(p) and os.path.getsize(p) > 5000: print('have', fn); return
    try:
        req = urllib.request.Request(url, headers={'User-Agent':'astro-amateur/1.0'})
        with urllib.request.urlopen(req, timeout=180) as r: data = r.read()
        open(p,'wb').write(data); print('ok', fn, len(data))
    except Exception as e: print('FAIL', fn, e)
    time.sleep(1)
# STScI DSS
for v in ['poss1_blue','poss1_red','poss2ukstu_blue','poss2ukstu_red','poss2ukstu_ir','quickv','phase2_gsc1','phase2_gsc2']:
    q = urllib.parse.urlencode({'v':v,'r':f'{ra}','d':f'{dec}','e':'J2000','h':size_arcmin,'w':size_arcmin,'f':'fits','c':'none','fov':'NONE','v3':''})
    get('https://archive.stsci.edu/cgi-bin/dss_search?'+q, f'stsci_{v}.fits')
