
import numpy as np, sys, json, time
sys.argv=[sys.argv[0]]+sys.argv[1:]
import importlib.util
spec=importlib.util.spec_from_file_location('dp','dual_pattern.py'); dp=importlib.util.module_from_spec(spec); spec.loader.exec_module(dp)
from scipy.optimize import linprog
import scipy.sparse as sp
theta=4/7; k0=7; X=float(sys.argv[1]); R=int(sys.argv[2]); lam=float(sys.argv[3])
h=lam/k0
t,mn,eta,it,p = dp.dual(lam,theta,h,X,R)
print('lam',lam,'X',X,'R',R,'t*',t)
K=len(eta); mu=h*(1+eta)
W=dp.words(R,k0)
# words with largest probability
idx=np.argsort(-p)[:12]
def wstr(w): return ''.join('1' if (w>>i)&1 else '.' for i in range(R))
print('top words (p):')
for j in idx: print('  %.4f %s'%(p[j], wstr(W[j])))
# Palm weights near multiples of lam
print('Palm atoms mu_k > 2h:')
for k in range(k0, min(K, int(6/h))):
    if mu[k-1] > 2*h: print('   u=%.3f (k=%d) mu=%.4f'%(k*h,k,mu[k-1]))
# gap distribution of the realizable process: P(next point at distance k | point at 0) from words with b0=1
gap=np.zeros(R)
for j,w in enumerate(W):
    if w&1:
        nxt=[i for i in range(1,R) if (w>>i)&1]
        if nxt: gap[nxt[0]]+=p[j]
        else: gap[0]+=p[j]  # no point within R
gap/=h
print('gap law (Palm): ', ', '.join('%.3f:%.3f'%(k*h,gap[k]) for k in range(1,R) if gap[k]>1e-3), ' none-within-R:', round(gap[0],4))
# consecutive gap pairs
pair={}
for j,w in enumerate(W):
    if w&1:
        nxt=[i for i in range(1,R) if (w>>i)&1]
        if len(nxt)>=2:
            key=(nxt[0], nxt[1]-nxt[0]); pair[key]=pair.get(key,0)+p[j]/h
print('top consecutive gap pairs (g1,g2):prob', sorted([(round(a*h,3),round(b*h,3),round(v,3)) for (a,b),v in pair.items() if v>0.01], key=lambda x:-x[2])[:12])
json.dump({'lam':lam,'X':X,'R':R,'t':t,'gap':gap.tolist(),'h':h}, open('logs/struct_R%d_X%d_lam%s.json'%(R,int(X),lam),'w'))
