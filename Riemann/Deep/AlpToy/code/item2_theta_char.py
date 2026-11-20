#!/usr/bin/env python3
"""alptoy item 2, third reading: theta functions NON-vanishing at 0 (half-period characteristics):
sigma_e(z) = exp(-eta_i z) sigma(z+w_i)/sigma(w_i) = sigma(z) sqrt(wp(z) - e_i)  (e_i = wp(w_i)), sigma_e(0) = 1, even.
D^a sigma_e (0) = a! [z^a] (sigma/z) * (z sqrt(wp - e)). Exact rationals; lattices with (g2,g3) = (4,0) [Z[i], e in {0,1,-1}]
and (0,4) [Z[omega], e = 1]."""
import sys
from fractions import Fraction as Fr
from math import factorial
A=int(sys.argv[1])
def vp(x,p):
    x=Fr(x)
    if x==0: return None
    v=0; n,d=x.numerator,x.denominator
    while n%p==0: n//=p; v+=1
    while d%p==0: d//=p; v-=1
    return v
def series(g2,g3,e):
    N=A//2+2; c={2:g2/20,3:g3/28}
    for n in range(4,N+1): c[n]=Fr(3,(2*n+1)*(n-3))*sum(c[m]*c[n-m] for m in range(2,n-1))
    L=[Fr(0)]*(A+1)
    for n in range(2,N+1):
        if 2*n<=A: L[2*n]=-c[n]/((2*n-1)*(2*n))
    E=[Fr(0)]*(A+1); E[0]=Fr(1)
    for k in range(1,A+1): E[k]=sum(i*L[i]*E[k-i] for i in range(1,k+1))/k
    U=[Fr(0)]*(A+1); U[2]=-e        # z^2 (wp - e) = 1 + U
    for n in range(2,N+1):
        if 2*n<=A: U[2*n]+=c[n]
    S=[Fr(0)]*(A+1); S[0]=Fr(1)
    for k in range(1,A+1): S[k]=(U[k]-sum(S[i]*S[k-i] for i in range(1,k)))/2
    P=[sum(E[i]*S[k-i] for i in range(k+1)) for k in range(A+1)]
    return P
for name,g2,g3,e,ps in [('Z[i] e=0',Fr(4),Fr(0),Fr(0),[3,7,11]),('Z[i] e=1',Fr(4),Fr(0),Fr(1),[3,7,11]),
                        ('Z[omega] e=1',Fr(0),Fr(4),Fr(1),[5,11,17])]:
    P=series(g2,g3,e)
    for p in ps:
        g=p/(p*p-1); r=1/(p*p-1)
        pts=[(a,vp(factorial(a)*P[a],p)) for a in range(2,A+1)]
        pts=[(a,v) for a,v in pts if v is not None]
        n=len(pts); sx=sum(a for a,_ in pts); sy=sum(v for _,v in pts); sxx=sum(a*a for a,_ in pts); sxy=sum(a*v for a,v in pts)
        sl=(n*sxy-sx*sy)/(n*sxx-sx*sx)
        print('%-13s p=%2d n=%3d slope=%.4f  gam=%.4f rho=%.4f  max(gam a - v)=%.2f  max(v - gam a)=%.2f  max(rho a - v)=%.2f'%(
            name,p,n,sl,g,r,max(g*a-v for a,v in pts),max(v-g*a for a,v in pts),max(r*a-v for a,v in pts)))
