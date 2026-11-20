import re
L=open('harvest.txt').read().split('\n')
cur=None
for l in L:
    m=re.match(r'^(NEW|SEEN)\s+(\S+)',l)
    if m: cur=(m.group(1),m.group(2)); continue
    if cur and cur[0]=='NEW' and l.startswith('   '):
        s=re.split(r'(?<=[.;])\s+',l.strip())
        hit=[x for x in s if re.search(r'GRH|ERH|Riemann [Hh]yp|generali[sz]ed Riemann',x)]
        print(cur[1]+': '+' | '.join(h[:240] for h in hit[:2]))
