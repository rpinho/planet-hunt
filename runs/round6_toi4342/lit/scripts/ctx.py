import re,sys
t=open(sys.argv[1],encoding='utf-8').read()
pat=re.compile(sys.argv[2], re.I)
for m in pat.finditer(t):
    s=max(0,m.start()-160); e=min(len(t),m.end()+160)
    print('...', t[s:e].replace('\n',' '), '...\n')
