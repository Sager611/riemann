# ugoaudit: run from .scratch/ugoaudit/ (TeX: curl -s -L -o 2512.11311.tgz https://arxiv.org/e-print/2512.11311; tar xzf into kl/); order: parse_tables2.py -> ugo class 6000 300 > class6000.log -> compare.py
import re, json
T=open('kl/CLNO-Unit-generated-arxiv-apr-21-2026.tex').read().split('\n')
def grp(s):
    s=s.replace(' ','').replace('$','')
    if s in ('1',): return []
    out=[]
    for m in re.finditer(r'\\left\(\\Z/(\d+)\\Z\\right\)\^(\d+)|\\Z/(\d+)\\Z',s):
        if m.group(1): out+=[int(m.group(1))]*int(m.group(2))
        else: out.append(int(m.group(3)))
    return sorted(out)
def rows(a,b,ncl):
    out=[]; cur=None
    for i in range(a-1,b):
        L=T[i]
        if '&' not in L or 'Cl(' in L or 'f_' in L: continue
        cells=[c.strip() for c in L.replace('\\\\','').split('&')]
        cl=cells[0]; clp=cells[1] if ncl==2 else None
        Ds=[int(x) for x in re.findall(r'-?\d+', cells[ncl].replace('$',''))]
        f=int(re.findall(r'\d+',cells[ncl+1])[0])
        if cl.strip(): cur=(grp(cl), grp(clp) if ncl==2 else grp(cl))
        for D in Ds: out.append(dict(D=D, cl=cur[0], clp=cur[1], f=f, line=i+1))
    return out
A=rows(1629,1668,2); B=rows(1682,1704,1)
for nm,R in (('2A',A),('2B',B)):
    pos=[r for r in R if r['D']>0]
    print(nm,'entries',len(R),'positive',len(pos),'neg',[r['D'] for r in R if r['D']<0])
allpos=[('A',r) for r in A if r['D']>0]+[('B',r) for r in B if r['D']>0]
Ds=[r['D'] for _,r in allpos]
print('total positive entries',len(Ds),'distinct',len(set(Ds)),'max',max(Ds))
print('maximal f=1 distinct', len({r['D'] for _,r in allpos if r['f']==1}), 'nonmax', len({r['D'] for _,r in allpos if r['f']>1}))
narrow=[r['D'] for _,r in allpos if all(x==2 for x in r['clp'])]
print('one class per genus entries',len(narrow),'distinct',len(set(narrow)))
print('not narrow 2-torsion:',sorted({r['D'] for _,r in allpos if not all(x==2 for x in r['clp'])}), [ (r['D'],r['clp']) for _,r in allpos if not all(x==2 for x in r['clp'])])
# family check
import math
def fam(D):
    out=[]
    n=math.isqrt(D+4)
    if n*n==D+4 and n>=3: out.append(('-4',n))
    if D>4:
        n=math.isqrt(D-4)
        if n*n==D-4 and n>=1: out.append(('+4',n))
    return out
bad=[(t,r['D']) for t,r in allpos if not fam(r['D']) or (t=='A' and not any(x[0]=='-4' for x in fam(r['D']))) or (t=='B' and not any(x[0]=='+4' for x in fam(r['D'])))]
print('family mismatches',bad)
json.dump([(t,r) for t,r in allpos],open('tables_pos.json','w'))
