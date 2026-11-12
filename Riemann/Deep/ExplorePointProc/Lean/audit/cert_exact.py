
# explore40 certificate checker for the realizability-augmented small-gap criterion (python-flint; exact Sturm + Arb).
# Checks: (1) f>=0 on R  (p(w)>0 on [0,inf)), (2) fhat<=0 on [U,inf)  (q(w)<0 on [w0,inf), w0<2 pi s^2 U^2),
# (3) J_theta(f) and fhat(0) in Arb, (4) c_cert >= sup_{x,y>=lam} fhat(x)+fhat(x+y)+A(x)-A(y)+beta(1-x)
#     via lattice corners + bilinear-interpolation remainder 3 M2 delta^2/8, M2 >= sup_{[lam,2U]}|fhat''| (Arb balls),
# (5) margin J - fhat(0) - 2 c_cert > 0.
import sys, json, time
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq, arb, arb_poly, ctx
import numpy as np
ctx.prec = 1200
d=json.load(open(sys.argv[1])); m=int(sys.argv[2]) if len(sys.argv)>2 else 20
lam=Fr(d['lam']); theta=Fr(d['theta']); D=d['D']; s=Fr(d['s']); U=Fr(d['U']); dg=Fr(d['dg'])
B=[Fr(b) for b in d['B']]; A=[Fr(a) for a in d['A']]; beta=Fr(d['beta'])
import os
EPS=Fr(os.environ.get('EPS','0')); B[0]+=EPS   # f -> f + EPS e^{-pi a^2/s^2}: lifts p by EPS (exact), q by EPS
print('EPS added to B_0:', EPS)
assert U<=3*lam and beta>=0
def q2(x): return fmpq(x.numerator, x.denominator)
t0=time.time()
# Hermite H_n integer polys
H=[fmpq_poly([1]), fmpq_poly([0,2])]
for n in range(1,2*D): H.append(fmpq_poly([0,2])*H[n]-2*n*H[n-1])
P=fmpq_poly([0]); Q=fmpq_poly([0])
for n in range(D+1):
    P+= q2(B[n])*H[2*n]; Q+= q2(B[n]*(-1)**n)*H[2*n]
def even_to_w(R):
    c=[R[2*k] for k in range(R.degree()//2+1)]
    assert all(R[2*k+1]==0 for k in range((R.degree()+1)//2))
    return fmpq_poly(c)
p=even_to_w(P); q=even_to_w(Q)
def sturm_count(f, a):   # number of distinct real roots in (a, inf)
    seq=[f, f.derivative()]
    while seq[-1].degree()>0:
        r=-(seq[-2] % seq[-1])
        if r==0: break
        seq.append(r)
    def sgn(x): return (x>0)-(x<0)
    def changes(vals):
        v=[x for x in vals if x!=0]; return sum(1 for i in range(len(v)-1) if v[i]*v[i+1]<0)
    va=[sgn(g(a)) for g in seq]; vinf=[sgn(g[g.degree()]) for g in seq]
    return changes(va)-changes(vinf)
# (1)
p0=p(fmpq(0)); n1=sturm_count(p, fmpq(0))
print('(1) p(0)=%s >0: %s ; real roots of p in (0,inf): %d ; lead>0: %s'%(float(p0), p0>0, n1, p[p.degree()]>0))
ok1 = p0>0 and n1==0
# (2)
PI_LO=Fr(314159265,100000000)   # < pi
w0=2*PI_LO*s*s*U*U
qw0=q(q2(w0)); n2=sturm_count(q, q2(w0))
print('(2) w0=%.6f q(w0)=%.3e <0: %s ; roots of q in (w0,inf): %d ; lead<0: %s'%(float(w0), float(qw0), qw0<0, n2, q[q.degree()]<0))
ok2 = qw0<0 and n2==0
# (3) J and fhat(0)
pi=arb.pi(); sA=arb(q2(s)); th=arb(q2(theta))
X=pi*th*th/(sA*sA)
eX=(-X).exp()
acc=arb(0); fact=arb(1); pow2=arb(1)
for k in range(p.degree()+1):
    if k>0: fact*=k; pow2*=2
    # e_k(X)=sum_{j<=k} X^j/j!
    ek=arb(0); term=arb(1)
    for j in range(k+1):
        if j>0: term=term*X/j
        ek+=term
    acc+= arb(p[k])*pow2*fact*(1-eX*ek)
J = arb(p0) + (sA*sA/pi)*acc
fh0 = sA*arb(q(fmpq(0)))
print('(3) J=%s  fhat(0)=%s'%(J.str(12), fh0.str(12)))
# fhat evaluation helpers (Arb)
QA=arb_poly([arb(Q[i]) for i in range(Q.degree()+1)])
dQA=QA.derivative(); ddQA=dQA.derivative()
k_=(2*pi).sqrt()*sA
def fhat(u):   # u arb
    v=k_*u; return sA*QA(v)*(-(v*v)/2).exp()
def fhat2(u):  # second derivative
    v=k_*u; return sA*k_*k_*(ddQA(v)-2*v*dQA(v)+(v*v-1)*QA(v))*(-(v*v)/2).exp()
# (4) M2 >= sup_{[lam,2U]}|fhat''|: grid max + M3*eta/2, M3 >= sup|fhat'''| from Cramer |h_n|<=0.8163 and
#     h_n' = sqrt(n/2)h_{n-1} - sqrt((n+1)/2)h_{n+1}  =>  |h_m'''| <= 0.8163 (2(m+3))^{3/2}
cn=[arb(q2(B[n]))*(arb(2)**(2*n)*arb(2*n).fac()*pi.sqrt()).sqrt() for n in range(D+1)]
M3=arb(0)
for n in range(D+1): M3+= abs(cn[n])*arb(0.8163)*(arb(2*(2*n+3)))**arb(1.5)
M3=sA*k_**3*M3
eta=Fr(1,8000); lo=lam; ng2=int((2*U-lam)/eta)+1
mx=0.0
for i in range(ng2+1):
    v=abs(fhat2(arb(q2(lo+i*eta))))
    mx=max(mx,float(v.upper()))
M2=mx+float(M3.upper())*float(eta)/2
print('(4) max|fhat\'\'| on grid %.3f, M3 <= %.1f, M2 <= %.3f (eta=%s, %d pts)'%(mx,float(M3.upper()),M2,eta,ng2+1))
delta=dg/m; N=int((U-lam)/delta); assert lam+N*delta==U
def up(x): return float(x.upper())
F1=np.array([up(fhat(arb(q2(lam+i*delta)))) for i in range(N+1)])        # fhat(x_i), x_i = lam+i delta
F2=np.array([up(fhat(arb(q2(2*lam+k*delta)))) for k in range(2*N+1)])    # fhat(x_i+x_j) = fhat(2 lam+(i+j) delta)
# A at lattice points (linear interpolation between nodes at spacing dg = m delta), as floats + 1e-15
Al=np.empty(N+1)
for i in range(N+1):
    k,r=divmod(i,m)
    if r==0: Al[i]=float(A[k])
    else: Al[i]=float(A[k]+(A[k+1]-A[k])*Fr(r,m))
xs=np.array([float(lam+i*delta) for i in range(N+1)])
betaf=float(beta)
best=-1e9
for i in range(N+1):
    row=F1[i]+F2[i:i+N+1]+Al[i]-Al+betaf*(1-xs[i])
    best=max(best,row.max())
if os.environ.get('FULLARB')=='1':
    # rigorous corner pass entirely in Arb (no float64): Phi(x_i,y_j) upper endpoints
    F1a=[fhat(arb(q2(lam+i*delta))) for i in range(N+1)]
    F2a=[fhat(arb(q2(2*lam+k*delta))) for k in range(2*N+1)]
    Ala=[]
    for i in range(N+1):
        k,r=divmod(i,m)
        Ala.append(arb(q2(A[k])) if r==0 else arb(q2(A[k]+(A[k+1]-A[k])*Fr(r,m))))
    bA=arb(q2(beta)); lin=[Ala[i]+bA*(1-arb(q2(lam+i*delta))) for i in range(N+1)]
    bestA=None
    for i in range(N+1):
        base=F1a[i]+lin[i]
        for j in range(N+1):
            v=(base+F2a[i+j]-Ala[j]).upper()
            if bestA is None or v>bestA: bestA=v
    AU=Ala[-1]; R2A=None
    for i in range(N+1):
        v=(F1a[i]+lin[i]-AU).upper()
        if R2A is None or v>R2A: R2A=v
    print('(4A) full-Arb corner max R1 = %s, R2 lattice max = %s'%(bestA.str(12), R2A.str(12)))
rem2=3*M2*float(delta)**2/8; rem1=M2*float(delta)**2/8
R1=best+rem2
R2=(F1+Al-Al[-1]+betaf*(1-xs)).max()+rem1
R3=float(max(A[-1]-a for a in A)+beta*(1-U))
ccert=max(R1,R2,R3)+1e-12
if os.environ.get('FULLARB')=='1':
    cA=max(float(bestA)+rem2, float(R2A)+rem1, R3)
    print('(4A) full-Arb c_cert (upper, float of Arb upper + rem) = %.10f vs float-pass c_cert %.10f'%(cA, ccert))
    ccert=max(ccert, cA+1e-15)
print('(4) delta=%s N=%d: R1 = corners %.10f + rem %.3e = %.10f ; R2 (incl. rem %.3e) = %.10f ; R3 = %.10f ; c_cert = max(R1,R2,R3) + 1e-12 = %.10f'%(delta,N,best,rem2,R1,rem1,R2,R3,ccert))
margin = J - fh0 - 2*arb(ccert)
print('(5) margin J - fhat(0) - 2 c_cert = %s  -> %s'%(margin.str(10), 'PASS' if (ok1 and ok2 and margin>0) else 'FAIL'))
print('time %.1fs'%(time.time()-t0))
# ---- pplean: exact rationals for the Lean certificate interface ----
def arb_to_frac(x):
    m,e = x.man_exp(); m=int(m); e=int(e)
    return Fr(m)*(Fr(2)**e) if e>=0 else Fr(m, 2**(-e))
Jl = arb_to_frac(J.lower()); Ju = arb_to_frac(J.upper())
q0 = q(fmpq(0)); F0 = s*Fr(int(q0.p), int(q0.q))
p0f = Fr(int(p0.p), int(p0.q))
cc = Fr(ccert)
Amin = min(A); Amax = max(A)
R3e = A[-1]-Amin+beta*(1-U)
print('EXACT J_lo =', Jl, float(Jl)); print('EXACT J_hi =', Ju, float(Ju))
print('EXACT F0 =', F0, float(F0)); print('EXACT f(0)=p(0) =', p0f, float(p0f))
print('EXACT ccert =', cc, float(cc))
print('EXACT Amin =', Amin, float(Amin), 'index', A.index(Amin), 'Amax =', Amax, 'A[-1] =', A[-1])
print('EXACT R3 =', R3e, float(R3e))
print('EXACT beta =', beta, 'lam =', lam, 'U =', U, 'theta =', theta, 'dg =', dg, 'nodes', len(A))
print('EXACT margin_lo =', float(Jl - F0 - 2*cc))
import math
for k in (10, 12, 13, 14, 15):
    jlo = Fr(math.floor(Jl*10**k), 10**k); chi = Fr(math.ceil(cc*10**k), 10**k)
    print('DEC k=%d Jlo=%s chi=%s margin=%.15f R3<=chi:%s'%(k, jlo, chi, float(jlo-F0-2*chi), R3e<=chi))
