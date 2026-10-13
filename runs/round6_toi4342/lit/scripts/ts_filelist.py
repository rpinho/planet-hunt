import json,sys
d=json.load(open(sys.argv[1]))
for f in d['files']:
    if f['ftype']!='Light_Curve': continue
    n=f['fname']
    if n.endswith(('.dat','.tbl','.txt','_lightcurve.png','.plotcfg')):
        print(f['fid'], n)
