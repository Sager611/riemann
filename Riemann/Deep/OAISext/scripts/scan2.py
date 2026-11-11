
import json,re,os
from collections import Counter
ROOT='/Users/adriansager/playground/riemann'
d=json.load(open(ROOT+'/.scratch/oaisext/closure.json'))
down=sorted(d['down'])
def p(m): return os.path.join(ROOT,*m.split('.'))+'.lean'
wb=re.compile(r'\bbeta\b')
ctx=Counter(); files=0; total=0
greek=0
samples={}
for m in down:
    src=open(p(m),encoding='utf-8').read()
    occ=list(wb.finditer(src))
    if occ: files+=1
    total+=len(occ)
    for o in occ:
        a=src[max(0,o.start()-25):o.start()]; b=src[o.end():o.end()+12]
        # classify by preceding token
        pre=re.search(r'([\w.]*)$',a).group(1)
        key=(pre[-25:] if pre else a[-3:])+'|'+re.match(r'[^\n]{0,6}',b).group(0)
        k2 = 'Hecke.beta' if pre.endswith('HeckeZeroSupremum.') else ('dot:'+pre if pre.endswith('.') else 'bare')
        ctx[k2]+=1
        samples.setdefault(k2,[]).append(m.split('.')[-1]+': ...'+a[-30:].replace('\n',' ')+'[beta]'+b.replace('\n',' '))
    greek+=len(re.findall(r'β',src))
print('files with beta:',files,'occurrences:',total,'greek β count:',greek)
for k,v in ctx.most_common(): print(k,v)
for k,v in samples.items():
    print('==',k)
    for s in v[:6]: print('   ',s)
# named args / binders
for m in down:
    src=open(p(m),encoding='utf-8').read()
    for pat in [r'\(beta\s*:=', r'\bbeta\s*:=', r'\(beta\s*:', r'fun\s+beta', r'⟨[^⟩]*\bbeta\b', r'\.beta\b']:
        for o in re.finditer(pat,src):
            print('PATTERN',pat,m.split('.')[-1],src[max(0,o.start()-40):o.end()+30].replace('\n',' '))

