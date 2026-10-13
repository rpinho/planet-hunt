import sys,re,html
t=open(sys.argv[1],encoding='utf-8',errors='replace').read()
t=re.sub(r'(?is)<(script|style).*?</\1>','',t)
t=re.sub(r'(?i)<br\s*/?>|</tr>|</p>|</div>|</h\d>','\n',t)
t=re.sub(r'(?i)</t[dh]>',' | ',t)
t=re.sub(r'<[^>]+>','',t); t=html.unescape(t)
t=re.sub(r'[ \t]+',' ',t); t=re.sub(r'\n\s*\n+','\n',t)
print(t)
