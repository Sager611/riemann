
# zp_rs.py: |zeta'(rho)| for the first 2,001,052 zeros (Odlyzko zeros6), via flint (low) and vectorized Riemann-Siegel Z'(t) (high).
import gzip, time, sys, numpy as np, mpmath as mp
from flint import acb, ctx, acb_series
t0=time.time()
DL='/Users/adriansager/playground/riemann/.scratch/chaudit/dl/'
g=np.loadtxt(gzip.open(DL+'zeros6.gz','rt'))
print('zeros',len(g),g[0],g[-1],'load s',round(time.time()-t0,1),flush=True)
# Chebyshev fit of Psi(p)=cos(2pi(p^2-p-1/16))/cos(2pi p) on [0,1]
mp.mp.dps=50
def Psi(p):
    p=mp.mpf(p); return mp.cos(2*mp.pi*(p*p-p-mp.mpf(1)/16))/mp.cos(2*mp.pi*p)
K=90
nodes=[(mp.cos(mp.pi*(k+mp.mpf(1)/2)/K)+1)/2 for k in range(K)]
vals=[Psi(x) for x in nodes]
from numpy.polynomial import chebyshev as C
xs=np.array([float(2*x-1) for x in nodes]); ys=np.array([float(v) for v in vals])
cf=C.chebfit(xs,ys,60)
# check fit
pt=[0.1,0.25,0.2500001,0.5,0.75,0.9]
err=max(abs(C.chebval(2*p-1,cf)-float(Psi(mp.mpf(p)+mp.mpf('1e-30')))) for p in pt)
dcf=C.chebder(cf)*2.0   # d/dp = 2 d/dx
print('cheb fit err',err,flush=True)
def thp(t): return 0.5*np.log(t/(2*np.pi)) - 1/(48*t*t) - 7/(1920*t**4)
def th(t): return 0.5*t*np.log(t/(2*np.pi)) - 0.5*t - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
def Zp_rs(t):
    a=np.sqrt(t/(2*np.pi)); N=np.floor(a).astype(np.int64); p=a-N
    Nm=int(N.max()); n=np.arange(1,Nm+1,dtype=np.float64); ln=np.log(n); isq=1/np.sqrt(n)
    ph=th(t)[:,None]-t[:,None]*ln[None,:]
    w=(n[None,:]<=N[:,None])
    s=-2*np.sum(np.where(w, isq[None,:]*np.sin(ph)*(thp(t)[:,None]-ln[None,:]),0.0),axis=1)
    sg=np.where(N%2==1,1.0,-1.0)  # (-1)^(N-1)
    dadt=1/(4*np.pi*a)
    R=sg*(a**-0.5*C.chebval(2*p-1,dcf) - 0.5*a**-1.5*C.chebval(2*p-1,cf))*dadt
    return s+R
ctx.prec=80
def zp_flint(gam):
    r=acb(0.5,gam); return abs(complex(acb_series([r,1],prec=2).zeta()[1]))
NF=50000
out=np.empty(len(g))
for i in range(NF):
    out[i]=zp_flint(g[i])
print('flint done',NF,'s',round(time.time()-t0,1),flush=True)
B=2000
for j in range(NF,len(g),B):
    tt=g[j:j+B]; out[j:j+B]=np.abs(Zp_rs(tt))
print('rs done s',round(time.time()-t0,1),flush=True)
# validation: flint vs rs at sample indices
idx=list(range(NF-200,NF+200,7))+list(range(200000,2001052,90000))+[2001051]
rel=[]
for i in idx:
    f=zp_flint(g[i]); r=abs(Zp_rs(np.array([g[i]]))[0]); rel.append(abs(f-r)/f)
print('validation n=%d max rel %.3e median %.3e'%(len(rel),max(rel),np.median(rel)),flush=True)
# compare Odlyzko zero positions with flint zeros at a few indices
for n0 in [1,1000,100000,1000000,2001052]:
    z=acb.zeta_zero(n0); print('zero',n0,'odlyzko',g[n0-1],'flint',z.imag.mid(),flush=True)
np.save('/Users/adriansager/playground/riemann/.scratch/chaudit/zp.npy',np.vstack([g,out]))
print('saved; total s',round(time.time()-t0,1))
