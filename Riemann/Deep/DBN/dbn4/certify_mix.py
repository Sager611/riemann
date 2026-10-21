
"""Exact certifier for two-shift mixed translated probes (same observation height).
Claim: for 0<y<=rho*h, 0<=b<=y, all real d (own singular point excluded),
   2K_y(d,b) >= (y/h) [lam1 J_{h,r1}(d,b) + lam2 J_{h,r2}(d,b)].
Reduction (derive_mix.py, sympy-verified): with h=1,
   2K_y/y - lam1 J1 - lam2 J2 = 4 Pmix/(Q* D4_1 D4_2),
   Pmix = (D+v-B) D4_1 D4_2 - (Q*/2)(lam1 G1 D4_2 + lam2 G2 D4_1),  Q*, D4_k > 0.
Pt(z,p,u) = (1-z)^9 Pmix(z/(1-z), rho^2 p, rho^2 p u); Pt>=0 on [0,1]^3 via exact Bernstein subdivision.
"""
from fractions import Fraction as F
from math import comb, lcm
import time
def padd(a,b):
    c=dict(a)
    for k,v in b.items():
        c[k]=c.get(k,0)+v
        if c[k]==0: del c[k]
    return c
def pscale(a,x): return {k:v*x for k,v in a.items() if v*x!=0}
def pmul(a,b):
    c={}
    for ka,va in a.items():
        for kb,vb in b.items():
            k=(ka[0]+kb[0],ka[1]+kb[1],ka[2]+kb[2]); c[k]=c.get(k,0)+va*vb
    return {k:v for k,v in c.items() if v!=0}
def const(x): return {(0,0,0):F(x)} if x else {}
# variables (D,p,u)
Dv={(1,0,0):F(1)}; pv={(0,1,0):F(1)}; uv={(0,0,1):F(1)}
def build(r1,r2,rho,lam1,lam2):
    s=F(rho)**2; R1=F(r1)**2; R2=F(r2)**2; l1=F(lam1); l2=F(lam2)
    v=pscale(pv,s); B=pscale(pmul(pv,uv),s)
    def parts(R):
        C=padd(padd(Dv,const(R+1)),B)
        A=padd(padd(pmul(C,C),pscale(Dv,4*R)),pscale(B,-4))
        D4=padd(pmul(A,A),pscale(pmul(Dv,pmul(C,C)),-16*R))
        G=padd(pmul(A,padd(C,pscale(B,-2))),pscale(pmul(Dv,C),-8*R))
        return D4,G
    D41,G1=parts(R1); D42,G2=parts(R2)
    W=padd(padd(Dv,v),pscale(B,-1))
    Qs=padd(pmul(W,W),pscale(pmul(Dv,B),4))
    P=padd(pmul(W,pmul(D41,D42)),pscale(pmul(Qs,padd(pscale(pmul(G1,D42),l1),pscale(pmul(G2,D41),l2))),F(-1,2)))
    return P
def homogenize(P,n=9):
    """(1-z)^n P(z/(1-z),p,u): coefficient of D^i times z^i (1-z)^(n-i)"""
    assert max(k[0] for k in P)<=n
    out={}
    for (i,j,k),c in P.items():
        for m in range(n-i+1):   # (1-z)^(n-i) = sum_m C(n-i,m)(-z)^m
            key=(i+m,j,k); out[key]=out.get(key,0)+c*comb(n-i,m)*(-1)**m
    return {k:v for k,v in out.items() if v!=0}
def bernstein(a,n):
    N1=n+1
    assert all(max(k)<=n for k in a)
    b=[F(a.get((i,j,k),0)) for i in range(N1) for j in range(N1) for k in range(N1)]
    for ax in range(3):
        c=[F(0)]*(N1**3)
        for i in range(N1):
            for j in range(N1):
                for k in range(N1):
                    idx=[i,j,k]; t=idx[ax]; val=F(0)
                    for m in range(t+1):
                        tt=list(idx); tt[ax]=m
                        val+=b[(tt[0]*N1+tt[1])*N1+tt[2]]*F(comb(t,m),comb(n,m))
                    c[(i*N1+j)*N1+k]=val
        b=c
    den=1
    for x in b: den=lcm(den,x.denominator)
    return [int(x*den) for x in b]
def halves(b,ax,n):
    N1=n+1; L=[0]*(N1**3); R=[0]*(N1**3)
    for i in range(N1):
        for j in range(N1):
            for k in range(N1):
                idx=[i,j,k]; t=idx[ax]
                vl=0
                for m in range(t+1):
                    tt=list(idx); tt[ax]=m; vl+=comb(t,m)*b[(tt[0]*N1+tt[1])*N1+tt[2]]
                vr=0
                for m in range(t,N1):
                    tt=list(idx); tt[ax]=m; vr+=comb(n-t,m-t)*b[(tt[0]*N1+tt[1])*N1+tt[2]]
                L[(i*N1+j)*N1+k]=vl*2**(n-t); R[(i*N1+j)*N1+k]=vr*2**t
    return L,R
def prove(b,n,maxnodes=60000,maxdepth=70):
    stack=[(b,0)]; leaves=0; nodes=0; md=0
    while stack:
        bb,dep=stack.pop(); nodes+=1
        if min(bb)>=0: leaves+=1; md=max(md,dep); continue
        if dep>=maxdepth or nodes>maxnodes: return False,leaves,nodes,md
        best=None
        for ax in range(3):
            l,rr=halves(bb,ax,n)
            sc=(sum(x<0 for x in l)+sum(x<0 for x in rr),-min(l+rr),ax)
            if best is None or sc<best[0]: best=(sc,l,rr)
        stack.append((best[2],dep+1)); stack.append((best[1],dep+1))
    return True,leaves,nodes,md
def certify_mix(r1,r2,rho,lam1,lam2,maxnodes=60000):
    t=time.time()
    P=build(r1,r2,rho,lam1,lam2); Pt=homogenize(P,9)
    ok,leaves,nodes,md=prove(bernstein(Pt,9),9,maxnodes=maxnodes)
    return dict(r1=str(F(r1)),r2=str(F(r2)),rho=str(F(rho)),lam1=str(F(lam1)),lam2=str(F(lam2)),
                geometry=ok,leaves=leaves,nodes=nodes,maxdepth=md,seconds=round(time.time()-t,2))
if __name__=='__main__':
    import json,sys
    sys.path.insert(0,'.')
    # consistency tests: (r1, any r2, lam2=0) reduces to the single kernel times D4_2>0
    for (r,rho,lam,expect) in [('3/5','1/5','49903/25000',True),('21/40','9/20','199/100',True),('19/40','23/50','2',True),
                               ('3/5','1/5','1.9962',False),('21/40','9/20','1.9939',False)]:
        res=certify_mix(r,'1/2',rho,F(lam),0,maxnodes=3000)
        print('single-as-mix',r,rho,lam,'->',res['geometry'],'expected',expect,'leaves',res['leaves'],'nodes',res['nodes'],res['seconds'],'s',flush=True)
