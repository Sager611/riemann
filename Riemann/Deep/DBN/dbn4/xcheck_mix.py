
import json, gzip, random, sys
from fractions import Fraction as F
import numpy as np, sympy as sp
sys.path.insert(0,'.')
from certify_mix import build, homogenize, bernstein
H=json.load(gzip.open('histories/classical_dbn4_mix.json.gz','rt'))
MK=H['mixed_kernels']
random.seed(1)
# 1) sympy construction vs dict builder
D,v,B,R1,R2,l1,l2,s,p,u=sp.symbols('D v B R1 R2 l1 l2 s p u')
def parts(Rv):
    C=D+Rv+1+B; A=C**2+4*D*Rv-4*B
    return A**2-16*D*Rv*C**2, A*(C-2*B)-8*D*Rv*C
D41,G1=parts(R1); D42,G2=parts(R2); Qs=(D+v-B)**2+4*D*B
Pm=(D+v-B)*D41*D42-Qs/2*(l1*G1*D42+l2*G2*D41)
for m in random.sample(MK,4):
    vals={R1:sp.Rational(m['r1'])**2,R2:sp.Rational(m['r2'])**2,l1:sp.Rational(m['lam1']),l2:sp.Rational(m['lam2'])}
    rho2=sp.Rational(m['rho'])**2
    e=sp.Poly(sp.expand(Pm.subs(vals).subs({v:rho2*p,B:rho2*p*u})),D,p,u)
    dsym={k:F(int(sp.Rational(c).p),int(sp.Rational(c).q)) for k,c in e.terms()}
    ddict=build(m['r1'],m['r2'],m['rho'],m['lam1'],m['lam2'])
    print('kernel',m['id'],'sympy==dict:',dsym==ddict,len(dsym))
# 2) Bernstein representation evaluates to the same polynomial at random rational points
from math import comb
m=MK[7]; P=build(m['r1'],m['r2'],m['rho'],m['lam1'],m['lam2']); Pt=homogenize(P,9)
bb=bernstein(Pt,9)
from math import lcm
# recover scale: compare ratios at random points
def evalpow(a,z,pp,uu): return sum(c*z**i*pp**j*uu**k for (i,j,k),c in a.items())
def evalbern(b,z,pp,uu,n=9):
    N1=n+1; tot=F(0)
    for i in range(N1):
        bi=comb(n,i)*z**i*(1-z)**(n-i)
        for j in range(N1):
            bj=comb(n,j)*pp**j*(1-pp)**(n-j)
            for k in range(N1):
                tot+=b[(i*N1+j)*N1+k]*bi*bj*comb(n,k)*uu**k*(1-uu)**(n-k)
    return tot
pts=[(F(random.randint(1,99),100),F(random.randint(1,99),100),F(random.randint(1,99),100)) for _ in range(3)]
ratios=[evalbern(bb,*x)/evalpow(Pt,*x) for x in pts]
print('Bernstein/power ratio constant:',len(set(ratios))==1,float(ratios[0]))
# 3) direct float check of the physical kernel inequality for all mixed kernels (independent of polynomial algebra)
def Kp(d,b,y): return (y-b)/(d*d+(y-b)**2)+(y+b)/(d*d+(y+b)**2)
def J(d,b,r): return 0.5*(Kp(d+r,b,1.0)+Kp(d-r,b,1.0))
dg=np.concatenate([np.linspace(1e-4,5,4000),np.geomspace(5,1e5,400)])
worst=1e9
for m in MK:
    rho=float(F(m['rho'])); r1=float(F(m['r1'])); r2=float(F(m['r2'])); a1=float(F(m['lam1'])); a2=float(F(m['lam2']))
    for y in (rho,0.9*rho,0.5*rho,0.1*rho):
        for bf in (0,0.5,0.9,0.99,0.999999):
            b=bf*y
            lhs=2*Kp(dg,b,y); rhs=y*(a1*J(dg,b,r1)+a2*J(dg,b,r2))
            worst=min(worst,float(np.min(lhs/rhs)))
print('min over all mixed kernels of 2K/(y*sum lam J) on float grid = %.9f (must be >=1)'%worst)
# same for a sample of Phase-1 single kernels
cat=json.load(open('catalogue.json'))['kernels']; w2=1e9
for k in cat:
    rho=float(F(k['rho'])); r=float(F(k['r'])); lam=float(F(k['lam']))
    for y in (rho,0.5*rho):
        for bf in (0,0.9,0.999999):
            b=bf*y; w2=min(w2,float(np.min(2*Kp(dg,b,y)/(y*lam*J(dg,b,r)))))
print('min over all 786 single kernels = %.9f (must be >=1)'%w2)
