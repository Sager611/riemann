import re, sys
# print GRH-bearing sentences of each NEW abstract in harvest.txt
L=open(sys.argv[1] if len(sys.argv)>1 else 'harvest.txt').read().split('\n')
cur=None
for l in L:
    m=re.match(r'^(NEW|SEEN)\s+(\S+)',l)
    if m: cur=(m.group(1),m.group(2)); continue
    if cur and l.startswith('   '):
        s=re.split(r'(?<=[.;])\s+',l.strip())
        hit=[x for x in s if re.search(r'GRH|ERH|Riemann [Hh]yp|generali[sz]ed Riemann|extended Riemann',x)]
        print(cur[0]+' '+cur[1]+': '+' | '.join(h[:260] for h in hit[:2]))
