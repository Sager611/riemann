
# independent float/mpmath cross-check of an augmented certificate (cert_aug.py input + EPS)
import sys, json, os, random
from fractions import Fraction as Fr
import mpmath as mp
import numpy as np
mp.mp.dps=40
d=json.load(open(sys.argv[1])); EPS=Fr(os.environ.get('EPS','0')); ccert=float(sys.argv[2])
lam=Fr(d['lam']); theta=Fr(d['theta']); D=d['D']; s=Fr(d['s']); U=Fr(d['U']); dg=Fr(d['dg'])
B=[Fr(b) for b in d['B']]; B[0]+=EPS; A=[Fr(a) for a in d['A']]; beta=float(Fr(d['beta']))
Bm=[mp.mpf(b.numerator)/b.denominator for b in B]; sm=mp.mpf(s.numerator)/s.denominator
k=mp.sqrt(2*mp.pi)
def f(a):
    y=k*a/sm; return sum(Bm[n]*mp.hermite(2*n,y) for n in range(D+1))*mp.exp(-y*y/2)
def fh(u):
    v=k*sm*u; return sm*sum(Bm[n]*(-1)**n*mp.hermite(2*n,v) for n in range(D+1))*mp.exp(-v*v/2)
th=mp.mpf(theta.numerator)/theta.denominator
J=f(0)+2*mp.quad(lambda a: f(a)*a,[0,th/2,th])
print('J (quadrature) =', mp.nstr(J,15), ' fhat(0)=', mp.nstr(fh(0),15))
for u in [0, 0.5, float(lam), 1.0, 1.5, float(U), 3.0]:
    num=2*mp.quad(lambda a: f(a)*mp.cos(2*mp.pi*a*u),[0,1,2,3,4,6,9])
    print('  u=%.4f  fhat Hermite %s  numeric FT %s'%(u, mp.nstr(fh(u),12), mp.nstr(num,12)))
# dense 2D random check in float (numpy Hermite via mpmath-converted coefficients -> use recurrence in float)
Bf=np.array([float(b) for b in B]); sf=float(s); kf=float(k)
def Hvals(n, x):
    H0=np.ones_like(x); H1=2*x; out=[H0,H1]
    for m in range(1,n): out.append(2*x*out[-1]-2*m*out[-2])
    return out
def fhf(u):
    u=np.asarray(u,float); v=kf*sf*u; Hs=Hvals(2*D,v)
    tot=np.zeros_like(v)
    for n in range(D+1): tot+=Bf[n]*(-1)**n*Hs[2*n]
    return sf*tot*np.exp(-v*v/2)
g=np.array([float(lam+i*dg) for i in range(len(A))]); Af=np.array([float(a) for a in A])
Ai=lambda x: np.interp(x,g,Af,right=Af[-1])
rng=np.random.default_rng(1)
worst=-9
for batch in range(20):
    x=rng.uniform(float(lam),3.5,50000); y=rng.uniform(float(lam),3.5,50000)
    phi=fhf(x)+fhf(x+y)+Ai(x)-Ai(y)+beta*(1-x)
    worst=max(worst,phi.max())
xx=np.linspace(float(lam),float(U),3001); 
for yv in [float(lam), float(lam)+1e-9, float(U), 5.0, 40.0]:
    worst=max(worst,(fhf(xx)+fhf(xx+yv)+Ai(xx)-Ai(yv)+beta*(1-xx)).max())
print('random/edge max Phi = %.10f vs c_cert %.10f : %s'%(worst, ccert, 'OK' if worst<=ccert else 'VIOLATION'))
uu=np.linspace(float(U),12,200001); print('max fhat on [U,12] = %.3e'%fhf(uu).max())
aa=np.linspace(0,12,400001); y=kf*aa/sf; Hs=Hvals(2*D,y); fa=sum(Bf[n]*Hs[2*n] for n in range(D+1))*np.exp(-y*y/2); print('min f on [0,12] = %.3e'%fa.min())
