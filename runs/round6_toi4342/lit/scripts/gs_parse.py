import re,sys,html
t=open(sys.argv[1],encoding='utf-8',errors='replace').read()
m=re.search(r'About ([\d,]+) results|(\d+) results',t); print('COUNT:', m.group(0) if m else None)
for blk in re.split(r'<div class="gs_r gs_or gs_scl"', t)[1:]:
    ti=re.search(r'<h3 class="gs_rt".*?</h3>',blk,re.S)
    href=re.search(r'<h3 class="gs_rt"[^>]*>.*?<a [^>]*href="([^"]+)"',blk,re.S)
    au=re.search(r'<div class="gs_a">(.*?)</div>',blk,re.S)
    sn=re.search(r'<div class="gs_rs">(.*?)</div>',blk,re.S)
    f=lambda x: html.unescape(re.sub('<[^>]+>','',x.group(0) if x else '')).strip()
    print('TITLE:',f(ti)); print('  URL:', href.group(1) if href else ''); print('  BY:', html.unescape(re.sub('<[^>]+>','',au.group(1))) if au else ''); print('  SNIP:', html.unescape(re.sub('<[^>]+>','',sn.group(1)))[:400] if sn else ''); print()
