import json,sys
w=json.load(open(sys.argv[1]))
seen=set()
for r in w:
    key=(r['btjd_first'],r['btjd_last'])
    if key in seen or r['btjd_first']==0: continue
    seen.add(key)
    for name,(t0,P) in {'A':(1675.535,18.9083),'B':(1658.952,1.95769)}.items():
        mid=0.5*(r['btjd_first']+r['btjd_last'])
        n=round((mid-t0)/P); tc=t0+n*P
        print(r['file'][:45], name, 'epoch',n,'tc',round(tc,4),'tc-mid (h)',round((tc-mid)*24,2), 'half-span h', round(r['span_h']/2,2), ' dP needed per epoch to shift into window (min):', round(abs(tc-mid)*24*60/max(abs(n),1),2))
