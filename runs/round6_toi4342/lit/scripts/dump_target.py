import json,sys
d=json.load(open(sys.argv[1]))
for k in sys.argv[2:]:
    print('=====',k)
    print(json.dumps(d[k],indent=1))
