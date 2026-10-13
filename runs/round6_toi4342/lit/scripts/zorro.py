import sys
fit=[];ann=[];mode=None
for line in open(sys.argv[1]):
    if line.startswith('# annulus'): mode='a'; continue
    if line.startswith('# fit'): mode='f'; continue
    if line.startswith('#'): continue
    p=[float(x) for x in line.split()]
    (ann if mode=='a' else fit).append(p)
for s in (0.1,0.15,0.2,0.3,0.5,0.8,1.0,1.17):
    f=min(fit,key=lambda r:abs(r[0]-s)); print('fit',s,f)
print('first/last annulus',ann[0],ann[-1])
