
import json, sys, time
from fractions import Fraction as F
sys.path.insert(0,'/Users/adriansager/playground/riemann/.scratch/dbn3/src/anc/certificate/code')
import check_geometry as cg
cat=json.load(open('catalogue.json'))['kernels']
import collections
print('my maxdepth histogram',sorted(collections.Counter(k['maxdepth'] for k in cat).items()))
def prove(b,maxdepth=80,maxnodes=200000):
    stack=[(b,0)]; leaves=0; nodes=0; md=0
    while stack:
        b,depth=stack.pop(); nodes+=1
        if min(b)>=0: leaves+=1; md=max(md,depth); continue
        if depth>=maxdepth or nodes>maxnodes: return False,leaves,nodes,md
        cand=[]
        for ax in range(3):
            a,c=cg.split(b,ax)
            score=(sum(v<0 for v in a)+sum(v<0 for v in c), -min(a+c), ax)
            cand.append((score,a,c))
        _,a,c=min(cand,key=lambda x:x[0])
        stack.append((c,depth+1)); stack.append((a,depth+1))
    return True,leaves,nodes,md
bad=[k for k in cat if k['maxdepth']>=29]
print('my kernels with depth>=29:',len(bad))
for k in cat:
    r,rho,lam=F(k['r']),F(k['rho']),F(k['lam'])
    p1=cg.prove(cg.bernstein(cg.poly(r,rho,lam)))
    if not p1[0]:
        p2=prove(cg.bernstein(cg.poly(r,rho,lam)))
        print(k['r'],k['rho'],'my depth',k['maxdepth'],'leaves',k['leaves'],'| Planat prove depth80:',p2)
