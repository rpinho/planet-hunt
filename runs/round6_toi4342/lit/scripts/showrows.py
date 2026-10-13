import csv,sys
rows=list(csv.DictReader(open(sys.argv[1],newline='')))
filt=sys.argv[2] if len(sys.argv)>2 else None
for r in rows:
    if filt and filt not in r.get('pl_refname',''): continue
    print('=====', r.get('pl_name'), r.get('pl_refname','')[:60])
    for k,v in r.items():
        if v not in ('',None): print(' ',k,'=',v)
