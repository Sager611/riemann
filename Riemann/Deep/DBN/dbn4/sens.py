
import numpy as np, math, json, sys, time
from fractions import Fraction as F
sys.path.insert(0,'.')
from replay import load, kernels, replay
S=load(); K57=kernels(S)
SH=['9/20','19/40','49/100','1/2','51/100','21/40','107/200','113/200','23/40','3/5']
shifts=[float(F(x)) for x in SH]; modes=['r_'+x.replace('/','_') for x in SH]
H=[[float(F(S[p]['plan']['rows'][i]['height'])) for i in range(608)] for p in (0,1)]
ELL0=np.array([[[float(F(S[p]['rows'][i][m])) for m in modes] for i in range(608)] for p in (0,1)])
cat=json.load(open('catalogue.json'))['kernels']
# lookup: per shift, sorted rho list with lam,M (certified rational values) + originals
tab={s:[] for s in SH}
for k in cat: tab[k['r']].append((float(F(k['rho'])),float(F(k['lam'])),float(F(k['self_cap']))))
for k in S[0]['kernels']: tab[k['r']].append((float(F(k['rho'])),float(F(k['lambda'])),float(F(k['self_cap']))))
for s in SH: tab[s].sort()
def make_rate(ELL):
    def br(i,j,q):
        best=0
        for p in (0,1):
            h=H[p][i]
            for si,s in enumerate(SH):
                for (rho,lam,M) in tab[s]:
                    if q<=h*h*rho*rho:
                        # first eligible in sorted rho is the tightest; but lam decreasing in rho so check a few
                        c=2*lam/h*(ELL[p][i][si]-M/h)
                        if c>best: best=c
                        break
        return best
    return br
U0=replay(S,K57,rate_fn=make_rate(ELL0)); print('base (new catalogue, float) U=%.9f'%U0,flush=True)
def U_with(adj):
    E=ELL0.copy(); adj(E); return replay(S,K57,rate_fn=make_rate(E))
for lo,hi in [(0,152),(152,304),(304,456),(456,608),(0,608)]:
    def a(E,lo=lo,hi=hi): E[:,lo:hi,:]+=0.01
    U=U_with(a); print('floors +0.01 on cells %d-%d: dU=%.3e  (dU/dl=%.4f)'%(lo,hi,U-U0,(U-U0)/0.01),flush=True)
# REPORT3 torus margins (cell: margin) interpolated, applied to all shifts/providers
mc=[(0,0.0248),(5,0.0261),(25,0.0295),(60,0.0388),(100,0.0501),(120,0.0670),(200,0.1241),(240,0.1791),(300,0.2930),(400,0.2771),(500,0.2995),(594,0.2526),(607,0.2526)]
xs=[c for c,_ in mc]; ys=[m for _,m in mc]
marg=np.interp(np.arange(608),xs,ys)
for f in (0.25,0.5,1.0):
    def a(E,f=f): E+=f*marg[None,:,None]
    U=U_with(a); print('floors + %.2f x torus margin: U=%.6f'%(f,U),flush=True)
def a(E): E+=math.log(2)/4
U=U_with(a); print('X doubled (floors + log2/4): U=%.6f'%U)
