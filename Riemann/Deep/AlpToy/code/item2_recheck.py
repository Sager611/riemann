#!/usr/bin/env python3
"""alptoy item 2 recheck by independent methods (exact rationals), K = Q(i), E: y^2 = x^3 - x (g2 = 4, g3 = 0), a <= A:
(a) sigma via DLMF 23.9.7 wp-recursion and log(sigma/z) = -sum c_n z^{2n}/((2n-1)2n)  vs  Weierstrass a_{m,n} values (main run);
(b) exp_E via formal logarithm log(t) = int omega dt, series reversion  vs  D-iteration values (main run)."""
import json, sys
from fractions import Fraction as Fr
from math import factorial
A=int(sys.argv[1]); main=json.load(open(sys.argv[2]))
def vp(x,p):
    x=Fr(x)
    if x==0: return None
    v=0; n,d=x.numerator,x.denominator
    while n%p==0: n//=p; v+=1
    while d%p==0: d//=p; v-=1
    return v
N=A//2+2
g2,g3=Fr(4),Fr(0)
c={2:g2/20,3:g3/28}
for n in range(4,N+1): c[n]=Fr(3,(2*n+1)*(n-3))*sum(c[m]*c[n-m] for m in range(2,n-1))
# L(z) = log(sigma/z) = -sum c_n z^{2n}/((2n-1)(2n)); series in z up to degree A
L=[Fr(0)]*(A+1)
for n in range(2,N+1):
    if 2*n<=A: L[2*n]=-c[n]/((2*n-1)*(2*n))
# exp(L) via E' = L' E
E=[Fr(0)]*(A+1); E[0]=Fr(1)
for k in range(1,A+1): E[k]=sum(i*L[i]*E[k-i] for i in range(1,k+1))/k
sig={a:factorial(a)*E[a-1] for a in range(1,A+1)}   # sigma = z*exp(L): [z^a] sigma = E[a-1]
# (b) formal group of y^2 = x^3 - x: w = t^3 + a4 t w^2 (a4=-1): w_n recursion over Z
a4=-1; M=A+6
w=[0]*(M+1); w[3]=1
w2=[0]*(M+1)
for n in range(4,M+1):
    w2=[sum(w[i]*w[k-i] for i in range(k+1)) for k in range(n)]+[0]*(M+1-n)
    w[n]=a4*w2[n-1]
w2=[sum(w[i]*w[k-i] for i in range(k+1)) for k in range(M+1)]
assert all(w[k]==(1 if k==3 else 0)+a4*(w2[k-1] if k>=1 else 0) for k in range(M+1))
# omega/dt = w'/(3t^2 + a4 w^2)  (a1=a2=a3=0): divide by t^2
num=[(k+1)*w[k+1] for k in range(M)]          # w'
den=[3*(1 if k==2 else 0)+a4*w2[k] for k in range(M+1)]
num=num[2:]; den=den[2:]
om=[Fr(0)]*(A+1)
for k in range(A+1):  # om = num/den, den[0]=3
    om[k]=(Fr(num[k]) - sum(om[i]*den[k-i] for i in range(k)))/den[0]
assert om[0]==1
lg=[Fr(0)]+[om[k-1]/k for k in range(1,A+1)]  # log(t) coefficients
# reversion: exp(z) = sum e_a z^a with log(exp(z)) = z; Lagrange-free iterative: solve degree by degree
e=[Fr(0)]*(A+1); e[1]=Fr(1)
powr={1:e[:]}  # powers of e truncated
def mulser(a,b):
    r=[Fr(0)]*(A+1)
    for i,x in enumerate(a):
        if x:
            for j in range(0,A+1-i):
                if b[j]: r[i+j]+=x*b[j]
    return r
for n in range(2,A+1):
    # coefficient of z^n in sum_k lg[k] e^k must vanish; e[n] enters only via k=1 term
    s=Fr(0); P=e[:]  # e^1 with e[n] unknown=0
    pk=e[:]
    for k in range(2,n+1):
        pk=mulser(pk,e)
        s+=lg[k]*pk[n]
    e[n]=-s
expv={a:factorial(a)*e[a] for a in range(1,A+1)}
bad=0; rows={(r['K'],r['p'],r['f']):r for r in main}
for p in [3,7,11]:
    Vs=rows[('Q(i)',p,'sigma')]['V']; Vt=rows[('Q(i)',p,'t')]['V']
    ds=[a for a in range(1,A+1) if vp(sig[a],p)!=Vs.get(str(a))]
    dt=[a for a in range(1,A+1) if vp(expv[a],p)!=Vt.get(str(a))]
    bad+=len(ds)+len(dt)
    print('p=%d: sigma mismatches %d (a<=%d), exp mismatches %d; e.g. a=%d: v(sigma)=%s v(exp)=%s, a p/(p^2-1)=%.2f a/(p^2-1)=%.2f'%(
        p,len(ds),A,len(dt),A-(A-1)%4,vp(sig[A-(A-1)%4],p),vp(expv[A-(A-1)%4],p),(A-(A-1)%4)*p/(p*p-1),(A-(A-1)%4)/(p*p-1)))
print('RECHECK', 'PASS' if bad==0 else 'FAIL', 'total mismatches', bad)
