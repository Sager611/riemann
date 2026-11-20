#!/usr/bin/env python3
"""alptoy item 2b: d = 23 (≡ 7 mod 16, h = 3). For p inert in Q(sqrt-23), p <= 30, H_{-23} mod p = (x - j1)(x - j0)^2;
the double root lifts to a conjugate pair alpha in the unramified quadratic extension: R = Z_p[alpha]/(Q(alpha)),
H = (x - j1) Q(x) over Z_p. Curve E_alpha: y^2 + xy = x^3 - 36/(alpha-1728) x - 1/(alpha-1728) (j = alpha), good reduction
iff alpha !≡ 0, 1728. Valuation of u in R: v = v_p(Norm u)/2 (unramified). Same D-iteration as item2_cm_padic.py."""
import sys, json, time
from flint import fmpz_poly
AMAX=int(sys.argv[1]); OUT=sys.argv[2]
Hc=[12771880859375,-5151296875,3491750,1]
def red(f,M,n=None):
    c=f.coeffs()
    if n is not None: c=c[:n]
    return fmpz_poly([int(x)%M for x in c])
def vp(n,p):
    v=0
    while n%p==0: n//=p; v+=1
    return v
class Ctx:
    def __init__(s,p,M,q0,q1): s.p,s.M,s.q0,s.q1=p,M,q0,q1
    def mul(s,F,G,n):
        f0,f1=F; g0,g1=G; a=f0*g0; b=f1*g1; c=f0*g1+f1*g0
        return (red(a-s.q0*b,s.M,n),red(c-s.q1*b,s.M,n))
    def add(s,F,G,n=None): return (red(F[0]+G[0],s.M,n),red(F[1]+G[1],s.M,n))
    def sub(s,F,G,n=None): return (red(F[0]-G[0],s.M,n),red(F[1]-G[1],s.M,n))
    def sc(s,c,F): return (red(c[0]*F[0]-s.q0*c[1]*F[1],s.M),red(c[0]*F[1]+c[1]*F[0]-s.q1*c[1]*F[1],s.M))
    def norm(s,u0,u1): return (u0*u0-s.q1*u0*u1+s.q0*u1*u1)%s.M
    def einv(s,u0,u1):
        N=s.norm(u0,u1); Ni=pow(N,-1,s.M); return ((u0-s.q1*u1)*Ni%s.M,(-u1)*Ni%s.M)   # conj(u)=u0+u1*(−q1−α)
    def c0(s,F): 
        a=F[0].coeffs(); b=F[1].coeffs(); return (int(a[0]) if a else 0, int(b[0]) if b else 0)
    def inv(s,F,n):
        e0=s.einv(*s.c0(F)); G=(fmpz_poly([e0[0]]),fmpz_poly([e0[1]])); pr=1
        while pr<n:
            pr=min(2*pr,n); E=s.mul(F,G,pr); two=(fmpz_poly([2]),fmpz_poly([0]))
            G=s.mul(G,s.sub(two,E,pr),pr)
        return G
def der(f): c=f.coeffs(); return fmpz_poly([int(c[i])*i for i in range(1,len(c))])
def shift(f,k): c=f.coeffs(); return fmpz_poly([int(x) for x in c[k:]])
def run(p,prec):
    M=p**prec
    # simple root j1 (Hensel), quadratic cofactor Q
    rts=[x for x in range(p) if sum(c*x**i for i,c in enumerate(Hc))%p==0]
    simple=[x for x in rts if sum(i*c*x**(i-1) for i,c in enumerate(Hc) if i)%p]
    assert len(simple)==1, (p,rts)
    x=simple[0]
    for _ in range(14):
        fx=sum(c*pow(x,i,M) for i,c in enumerate(Hc))%M; dfx=sum(i*c*pow(x,i-1,M) for i,c in enumerate(Hc) if i)%M
        x=(x-fx*pow(dfx,-1,M))%M
    j1=x; q1=(Hc[2]+j1)%M; q0=(Hc[1]+j1*q1)%M   # H=(x-j1)(x^2+q1 x+q0)
    assert (Hc[0]+j1*q0)%M==0
    j0=[r for r in rts if r!=j1%p][0]
    C=Ctx(p,M,q0,q1)
    if j0%p in (0,1728%p): return dict(p=p,skip='double root j0=%d is 0 or 1728 mod p'%j0)
    # alpha - 1728 inverse in R; a4=-36 d, a6=-d
    d=C.einv((-1728)%M,1)
    a4=((-36*d[0])%M,(-36*d[1])%M); a6=((-d[0])%M,(-d[1])%M)
    L=AMAX+2; n=L+4
    one=(fmpz_poly([1]),fmpz_poly([0])); T=(fmpz_poly([0,1]),fmpz_poly([0]))
    t2=C.mul(T,T,n); t3=C.mul(t2,T,n); w=t3
    for it in range(14):  # a1=1,a2=0,a3=0
        w2=C.mul(w,w,n); w3=C.mul(w2,w,n)
        F=C.sub(C.sub(C.sub(w,t3,n),C.mul(T,w,n),n),C.add(C.sc(a4,C.mul(T,w2,n)),C.sc(a6,w3),n),n)
        if F[0].degree()<0 and F[1].degree()<0: break
        dF=C.sub(C.sub(one,T,n),C.add(C.sc(((2*a4[0])%M,(2*a4[1])%M),C.mul(T,w,n)),C.sc(((3*a6[0])%M,(3*a6[1])%M),w2),n),n)
        w=C.sub(w,C.mul(F,C.inv(dF,n),n),n)
    assert F[0].degree()<0 and F[1].degree()<0
    # omega/dt = (w - t w')/(w(-2 + t))  (a1=1, a3=0), p odd
    dw=(der(w[0]),der(w[1])); twp=C.mul(T,dw,n)
    u=C.sub(w,twp,n); u=(shift(u[0],3),shift(u[1],3)); v=(shift(w[0],3),shift(w[1],3))
    den=(fmpz_poly([-2,1]),fmpz_poly([0]))
    om=C.mul(C.mul(u,C.inv(v,L),L),C.inv(den,L),L)
    # cross-check with second formula omega/dt = w'/(3t^2 + a4 w^2 + w)
    den2=C.add(C.add((fmpz_poly([0,0,3]),fmpz_poly([0])),C.sc(a4,C.mul(w,w,n)),n),w,n)
    om2=C.mul((shift(dw[0],2),shift(dw[1],2)),C.inv((shift(den2[0],2),shift(den2[1],2)),L),L)
    dd=C.sub(om,om2,L); assert dd[0].degree()<0 and dd[1].degree()<0, 'omega mismatch'
    oi=C.inv(om,L); f=T; V={}; Dpz=None
    for j in range(AMAX+1):
        u0,u1=C.c0(f)
        if j>=1:
            N=C.norm(u0%M,u1%M)
            V[j]=None if N==0 else vp(N,p)/2
        if j==p:  # D^p t == 0 mod pO coefficientwise; Z_p[alpha] is non-maximal, so test v(u) >= 1 via v_p(Norm u) >= 2
            c0s=[int(c) for c in f[0].coeffs()]; c1s=[int(c) for c in f[1].coeffs()]; nn=max(len(c0s),len(c1s))
            c0s+= [0]*(nn-len(c0s)); c1s+=[0]*(nn-len(c1s))
            Dpz=all(C.norm(a%M,b%M)%(p*p)==0 for a,b in zip(c0s,c1s))
            Dp_naive=all(a%p==0 for a in c0s) and all(b%p==0 for b in c1s)
        f=C.mul(oi,(der(f[0]),der(f[1])),L-j-1)
    return dict(K='Q(sqrt-23)',p=p,kind='inert',f='t',prec=prec,Dp_zero=Dpz,Dp_componentwise=Dp_naive,note='j1=%d j0=%d mod p (double)'%(j1%p,j0),V=V)
rows=[]
for p in [17,19]:  # 5, 7, 11: no usable root (j = 0 / 1728 mod p, or triple/double root at 0, 1728), see logs/item2b_d23.log first run
    t0=time.time()
    try: r=run(p,AMAX//(p-1)+80)
    except AssertionError as e: r=dict(p=p,skip='assert %s'%(e,))
    print(p,{k:v for k,v in r.items() if k!='V'},'%.1fs'%(time.time()-t0),flush=True)
    if 'V' in r: rows.append(r)
json.dump(rows,open(OUT,'w'))
