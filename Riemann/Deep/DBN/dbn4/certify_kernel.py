
"""Independent exact certifier for fixed translated-probe kernels (r, rho, lam).

Claim certified: for 0 < y <= rho*h, 0 <= b <= y, all real d (singular point d=0,b=y excluded),
    2 K_y(d,b) >= lam * (y/h) * J_{h,r}(d,b).
Reduction (derive.py, sympy-verified): with h=1, D=d^2, v=y^2, B=b^2, R=r^2,
    2K_y/y - lam*J = 4 P /(Q* D4),  P=(D+v-B) D4 - (lam/2) Q* G,  Q*,D4 > 0.
So the claim <=> P >= 0 on D>=0, 0<=B<=v<=rho^2.  Substitute D=z/(1-z), v=rho^2 p, B=rho^2 p u,
Pt(z,p,u) = (1-z)^5 P  (degree <= 5 in each variable), and prove Pt >= 0 on [0,1]^3 by exact
integer Bernstein coefficients with de Casteljau subdivision.
"""
import sympy as sp, json, sys, time
from fractions import Fraction as F
from math import comb
z,p,u,Rs,s,lm=sp.symbols('z p u R s lam')
def build_generic():
    D=z/(1-z); v=s*p; B=s*p*u; R=Rs
    C=D+R+1+B; A=C**2+4*D*R-4*B
    D4=A**2-16*D*R*C**2; G=A*(C-2*B)-8*D*R*C
    Qs=(D+v-B)**2+4*D*B
    P=(D+v-B)*D4-lm/2*Qs*G
    Pt=sp.cancel(sp.together(P*(1-z)**5))
    Pt=sp.Poly(sp.expand(Pt),z,p,u)
    coeffs={}
    for mon,c in Pt.terms():
        cp=sp.Poly(c,Rs,s,lm)
        terms=[(e,F(int(sp.Rational(q).p),int(sp.Rational(q).q))) for e,q in cp.terms()]
        coeffs[mon]=(lambda T: (lambda R,S,L: sum((q*R**e[0]*S**e[1]*L**e[2] for e,q in T),F(0))))(terms)
    degs=Pt.degree(z),Pt.degree(p),Pt.degree(u)
    return coeffs,degs,Pt
COEFFS,DEGS,PT=build_generic()
assert max(DEGS)<=5, DEGS
N=5
IDX=[(i,j,k) for i in range(6) for j in range(6) for k in range(6)]
def power_coeffs(r,rho,lam):
    R=F(r)**2; S=F(rho)**2; L=F(lam)
    a=[F(0)]*216
    for (i,j,k),f in COEFFS.items():
        a[i*36+j*6+k]=F(f(R,S,L))
    return a
def to_bernstein(a):
    # b_k = sum_{j<=k} a_j C(k,j)/C(5,j) along each axis; scale by 10^3 to integers after lcm
    b=list(a)
    for ax in range(3):
        c=[F(0)]*216
        for (i,j,k) in IDX:
            idx=[i,j,k]; n=idx[ax]; val=F(0)
            for m in range(n+1):
                t=list(idx); t[ax]=m
                val+=b[t[0]*36+t[1]*6+t[2]]*F(comb(n,m),comb(5,m))
            c[i*36+j*6+k]=val
        b=c
    from math import lcm
    den=1
    for x in b: den=lcm(den,x.denominator)
    return [int(x*den) for x in b]
def halves(b,ax):
    L=[0]*216; Rr=[0]*216
    for (i,j,k) in IDX:
        idx=[i,j,k]; n=idx[ax]
        vl=0
        for m in range(n+1):
            t=list(idx); t[ax]=m; vl+=comb(n,m)*b[t[0]*36+t[1]*6+t[2]]
        vr=0
        for m in range(n,6):
            t=list(idx); t[ax]=m; vr+=comb(5-n,m-n)*b[t[0]*36+t[1]*6+t[2]]
        L[i*36+j*6+k]=vl*2**(5-n); Rr[i*36+j*6+k]=vr*2**n
    return L,Rr
def prove(b,maxnodes=200000,maxdepth=60):
    stack=[(b,0)]; leaves=0; nodes=0; md=0
    while stack:
        bb,dep=stack.pop(); nodes+=1
        if min(bb)>=0: leaves+=1; md=max(md,dep); continue
        if dep>=maxdepth or nodes>maxnodes: return False,leaves,nodes,md
        best=None
        for ax in range(3):
            l,rr=halves(bb,ax)
            sc=(sum(x<0 for x in l)+sum(x<0 for x in rr),ax)
            if best is None or sc<best[0]: best=(sc,l,rr)
        stack.append((best[2],dep+1)); stack.append((best[1],dep+1))
    return True,leaves,nodes,md
def s_cap_ok(M,r,rho):
    """s_r(w)<=M on [0,rho^2] <=> M w^2+(M(2R-2)+2)w+M(1+R)^2-2(1+R) >= 0 there (denominator>0)."""
    M=F(M);R=F(r)**2;V=F(rho)**2
    qa=M;qb=M*(2*R-2)+2;qc=M*(1+R)**2-2*(1+R)
    vals=[qc,qa*V*V+qb*V+qc]; vx=-qb/(2*qa)
    if 0<vx<V: vals.append(qa*vx*vx+qb*vx+qc)
    return min(vals)>=0
def certify(r,rho,lam,M=None,maxnodes=200000):
    t=time.time()
    ok,leaves,nodes,md=prove(to_bernstein(power_coeffs(r,rho,lam)),maxnodes=maxnodes)
    cap=None if M is None else s_cap_ok(M,r,rho)
    return dict(r=str(F(r)),rho=str(F(rho)),lam=str(F(lam)),self_cap=None if M is None else str(F(M)),
                geometry=ok,cap=cap,leaves=leaves,nodes=nodes,maxdepth=md,seconds=round(time.time()-t,3))
if __name__=='__main__':
    print('degrees',DEGS)
    base='/Users/adriansager/playground/riemann/.scratch/dbn3/src/anc/certificate/'
    S=json.load(open(base+'inputs/refined/SOURCE_RETUNED.json'))
    tot=0;t0=time.time()
    for idx,k in enumerate(S['kernels']):
        res=certify(k['r'],k['rho'],k['lambda'],k['self_cap'])
        tot+=res['leaves']
        if not (res['geometry'] and res['cap']): print('FAIL',idx,res)
    print('catalogue 57: total leaves',tot,'time %.1fs'%(time.time()-t0))
    # negative controls: lam slightly above numeric lam_max must FAIL
    for (r,rho,lam) in [('3/5','1/5','1.9962'),('21/40','9/20','1.9939'),('19/40','23/50','2.0000001'),('3/5','9/20','1.8902')]:
        res=certify(r,rho,F(lam),maxnodes=20000); print('NEG',r,rho,lam,'geometry',res['geometry'],'nodes',res['nodes'])
