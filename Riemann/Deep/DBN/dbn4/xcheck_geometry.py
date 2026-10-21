
import json, sys, time
from fractions import Fraction as F
sys.path.insert(0,'/Users/adriansager/playground/riemann/.scratch/dbn3/src/anc/certificate/code')
import check_geometry as cg
cat=json.load(open('catalogue.json'))['kernels']
t=time.time(); ok=0; leaves=0; bad=[]
for k in cat:
    r,rho,lam,M=F(k['r']),F(k['rho']),F(k['lam']),F(k['self_cap'])
    a=cg.poly(r,rho,lam)
    passed,lv,nodes,depth=cg.prove(cg.bernstein(a))
    R=r*r; qa=M;qb=M*(2*R-2)+2;qc=M*(1+R)**2-2*(1+R)
    cands=[F(0),rho*rho]; vx=-qb/(2*qa)
    if 0<vx<rho*rho: cands.append(vx)
    capok=min(qa*w*w+qb*w+qc for w in cands)>=0
    if passed and capok and 0<rho<1 and 0<lam<=2: ok+=1; leaves+=lv
    else: bad.append((k['r'],k['rho'],k['lam'],passed,capok,nodes,depth))
print('Planat check_geometry on new kernels: pass',ok,'of',len(cat),'leaves',leaves,'%.1fs'%(time.time()-t))
print('failures',bad[:10])
