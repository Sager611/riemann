# ugoaudit: run from .scratch/ugoaudit/ (TeX: curl -s -L -o 2512.11311.tgz https://arxiv.org/e-print/2512.11311; tar xzf into kl/); order: parse_tables2.py -> ugo class 6000 300 > class6000.log -> compare.py
import json, re, math, sys
from collections import Counter
tab=json.load(open('tables_pos.json'))
KL={}
for t,r in tab:
    KL.setdefault((t,r['D']),[]).append(r)
mine={}
for L in open('class6000.log'):
    if not L.startswith('TORS'): continue
    d=dict(kv.split('=') for kv in L.split()[1:])
    N=int(d['N']); D=int(d['D']); h=int(d['h']); hp=int(d['h+']); mu=int(d['mu']); w=int(d['wide2']); n2=int(d['narrow2'])
    if not w: print('narrow but not wide?!',L); continue
    cl=[2]*int(round(math.log2(h)))
    assert 2**len(cl)==h
    if n2: clp=[2]*(mu-1); assert hp==2**(mu-1)
    else: assert hp==2**mu, L; clp=sorted([2]*(mu-2)+[4])
    mine[('A' if N==1 else 'B',D)]=dict(cl=cl,clp=clp,f=int(d['f']),n=int(d['n']))
ok=True
for k in sorted(set(KL)|set(mine)):
    if k not in KL or k not in mine: print('MISSING',k,k in KL,k in mine); ok=False; continue
    r=KL[k][0]; m=mine[k]
    if sorted(r['cl'])!=m['cl'] or sorted(r['clp'])!=m['clp'] or r['f']!=m['f']:
        print('DIFF',k,r,m); ok=False
print('entries KL',len(KL),'mine',len(mine),'all match (Cl, Cl+, f):',ok)
print('max n in list:',max(m['n'] for m in mine.values()),'max D',max(k[1] for k in mine))
# compare with cert.py transcription
src=open('../../Riemann/Deep/KoppLag/cert.py').read()
ns={}
exec(src[src.index('KL_A = {'):src.index('ours_A =')],ns)
tA={D:(sorted(v[0]),sorted(v[1])) for D,v in ns['KL_A'].items()}
tB={D:sorted(v) for D,v in ns['KL_B'].items()}
pA={k[1]:(sorted(r[0]['cl']),sorted(r[0]['clp'])) for k,r in KL.items() if k[0]=='A'}
pB={k[1]:sorted(r[0]['cl']) for k,r in KL.items() if k[0]=='B'}
print('cert.py KL_A transcription == TeX parse:',tA==pA,' KL_B:',tB==pB)
# survivors
theirs={}
for L in open('../../Riemann/Deep/KoppLag/logs/survivors.txt'):
    if L.startswith('#'): continue
    d=dict(kv.split('=') for kv in L.split())
    theirs[(int(d['N']),int(d['n']))]=(int(d['h']),int(d['hplus']),d['wide2'])
ms={}
for L in open('class6000.log'):
    if not L.startswith('SURVCLASS'): continue
    d=dict(kv.split('=') for kv in L.split()[1:])
    ms[(int(d['N']),int(d['n']))]=(int(d['h']),int(d['h+']),'True' if d['wide2']=='1' else 'False')
print('survivors: theirs',len(theirs),'mine (n in (1000,6000], split test passed)',len(ms),'same set:',set(theirs)==set(ms),'same h,h+,wide2:',theirs==ms)
print('survivor h range',min(v[0] for v in ms.values()),max(v[0] for v in ms.values()),'max n',max(k[1] for k in ms))
