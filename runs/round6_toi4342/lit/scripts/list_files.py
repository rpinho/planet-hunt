import json,sys
d=json.load(open(sys.argv[1]))
print(json.dumps(d['files'][0],indent=1))
for f in d['files']:
    print(f.get('fid'), '|', f.get('ftype'), '|', f.get('fname'), '|', f.get('fsize'), '|', f.get('fdate'), '|', f.get('fuser'), '|', f.get('ftoi',''), '|', (f.get('fdesc') or '')[:120], '|', f.get('ftag'))
