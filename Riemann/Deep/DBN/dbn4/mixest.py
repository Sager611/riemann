
import numpy as np, math, json, sys, time, gzip
from fractions import Fraction as F
sys.path.insert(0,'.')
from replay import load, kernels, replay
S=load(); K57=kernels(S)
SH=['9/20','19/40','49/100','1/2','51/100','21/40','107/200','113/200','23/40','3/5']
shifts=[float(F(x)) for x in SH]; modes=['r_'+x.replace('/','_') for x in SH]
H=[[float(F(S[p]['plan']['rows'][i]['height'])) for i in range(608)] for p in (0,1)]
ELL=np.array([[[float(F(S[p]['rows'][i][m])) for m in modes] for i in range(608)] for p in (0,1)])
cat=json.load(open('catalogue.json'))['kernels']
tab={s:[] for s in SH}
for k in cat: tab[k['r']].append((float(F(k['rho'])),float(F(k['lam'])),float(F(k['self_cap']))))
for k in S[0]['kernels']: tab[k['r']].append((float(F(k['rho'])),float(F(k['lambda'])),float(F(k['self_cap']))))
for s in SH: tab[s].sort()
def s_r(w,r): R=r*r; return 2*(1+R-w)/((1+R-w)**2+4*R*w)
zz=np.linspace(1e-7,1-1e-9,6001); DD=zz/(1-zz)
def parts(r,v,B):
    R=r*r; D=DD
    C=D+R+1+B; A=C*C+4*D*R-4*B
    D4=A*A-16*D*R*C*C; G=A*(C-2*B)-8*D*R*C
    return 2*G/D4   # J_{1,r}
def Tfun(v,B):
    Qs=(DD+v-B)**2+4*DD*B; return 4*(DD+v-B)/Qs   # 2K_y/y
thetas=np.linspace(0,1,41)
def best_mix(p,i,q):
    h=H[p][i]; rho=math.sqrt(q)/h*(1+1e-9)
    if rho>0.48: return (0,None)
    v=rho*rho
    best=(0,None)
    cache={}
    for b_frac in (1.0,):
        B=v*(1-1e-13)
        T=Tfun(v,B); Js=[parts(r,v,B) for r in shifts]
        Ms=[max(s_r(w,r) for w in np.linspace(0,v,40)) for r in shifts]
        for a in range(10):
            for b in range(a+1,10):
                for th in thetas[1:-1]:
                    lt=np.nanmin(T/((1-th)*Js[a]+th*Js[b]))
                    lt=min(lt,2.0)
                    c=2*lt/h*((1-th)*(ELL[p][i][a]-Ms[a]/h)+th*(ELL[p][i][b]-Ms[b]/h))
                    if c>best[0]: best=(c,(p,a,b,th,lt,rho))
    return best
Hh=json.load(gzip.open('histories/classical_dbn4.json.gz','rt'))
qcell={rec['cell']:float(F(rec['q_in'])) for rec in Hh['records'] if rec['sub']==0}
t=time.time()
mix={}
for i in range(0,345):
    if i not in qcell: break
    q=qcell[i]
    m=max(best_mix(0,i,q),best_mix(1,i,q),key=lambda x:x[0])
    mix[i]=m
    if i%25==0: print(i,'%.5f'%m[0],m[1],'%.0fs'%(time.time()-t),flush=True)
def single(i,q):
    best=0
    for p in (0,1):
        h=H[p][i]
        for si,s in enumerate(SH):
            for (rho,lam,M) in tab[s]:
                if q<=h*h*rho*rho:
                    c=2*lam/h*(ELL[p][i][si]-M/h); best=max(best,c); break
    return best
def br(i,j,q):
    s=single(i,q)
    m=mix.get(i,(0,None))[0]
    return max(s,m*(1-2e-6))
U0=replay(S,K57,rate_fn=lambda i,j,q:single(i,q))
U1=replay(S,K57,rate_fn=br)
print('single U=%.9f  with per-cell 2-shift mixtures U=%.9f  gain %.3e'%(U0,U1,U0-U1))
json.dump({str(k):[v[0],v[1]] for k,v in mix.items()},open('mix_cells.json','w'))
