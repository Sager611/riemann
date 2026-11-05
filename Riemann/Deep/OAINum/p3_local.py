
# p3_local.py -- finite checks of Section 2 (auxiliary polynomials, local admissibility) of
# "A Power Saving for Polynomial Differences at Prime Arguments"
import math, itertools, sys
from fractions import Fraction

def primes_upto(n):
    s=[True]*(n+1); s[0]=s[1]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]:
            for j in range(i*i,n+1,i): s[j]=False
    return [i for i in range(n+1) if s[i]]

def peval(h,x):
    r=0
    for c in reversed(h): r=r*x+c
    return r
def pmul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
def compose_affine(h,c0,c1):
    # h(c0 + c1 x) as coefficient list
    res=[0]; pw=[1]
    for c in h:
        res=[ (res[i] if i<len(res) else 0) + c*(pw[i] if i<len(pw) else 0) for i in range(max(len(res),len(pw)))]
        pw=pmul(pw,[c0,c1])
    while len(res)>1 and res[-1]==0: res.pop()
    return res
def deriv(h): return [i*h[i] for i in range(1,len(h))] or [0]
def vp(n,p):
    if n==0: return 10**9
    v=0
    while n%p==0: n//=p; v+=1
    return v

def unit_root_padic(h,p,J,prefer=None):
    # choose a unit root of h in Z_p, return its residue mod p^J; prefer exact integer roots
    if prefer is not None:
        for z in prefer:
            if z%p and peval(h,z)==0: return z%(p**J)
    # DFS on the root tree of unit roots mod p^j with lookahead to depth J
    def dfs(r,j):
        if j==J: return r
        for t in range(p):
            rr=r+t*p**j
            if peval(h,rr)%(p**(j+1))==0:
                out=dfs(rr,j+1)
                if out is not None: return out
        return None
    for r0 in range(1,p):
        if peval(h,r0)%p==0:
            out=dfs(r0,1)
            if out is not None: return out
    return None

def mult_at(h,z,p,J,exact):
    # m_p = order of vanishing of h(z+X) ; exact integer root -> exact; else simple-root polys only
    if exact:
        b=compose_affine(h,z,1)
        for j,c in enumerate(b):
            if c!=0: return j
    return 1

def check(hname,h,prefer,Lmax=48,Dmax=12):
    k=len(h)-1; a=h[-1]
    ps=primes_upto(Lmax*Dmax)
    J=12
    z={}; m={}
    for p in ps:
        zp=unit_root_padic(h,p,J+2,prefer)
        assert zp is not None,('no unit root',hname,p)
        exact = prefer is not None and any(zz%p and peval(h,zz)==0 and zz%(p**(J+2))==zp for zz in prefer)
        z[p]=zp; m[p]=mult_at(h,[zz for zz in prefer if zz%p and peval(h,zz)==0][0] if exact else zp,p,J,exact)
    def lam(n):
        r=1
        for p in ps:
            v=vp(n,p) if n%p==0 else 0
            r*=p**(v*m[p])
        return r
    def rl(l):
        if l==1: return 0
        # CRT r = z_p mod p^{v_p(l)}, r in (-l,0]
        r=0; mod=1
        for p in ps:
            if l%p: continue
            v=vp(l,p); pe=p**v
            # combine r mod mod with z_p mod pe
            t=((z[p]-r)*pow(mod,-1,pe))%pe
            r=r+mod*t; mod*=pe
        r%=l
        return r-l if r>0 else 0
    bad=0; cnt=0
    for l in range(1,Lmax+1):
        r=rl(l)
        assert -l<r<=0 and math.gcd(r,l)==1 or l==1
        num=compose_affine(h,r,l); L=lam(l)
        if any(c%L for c in num): bad+=1; print('NONINTEGRAL',hname,l); continue
        hl=[c//L for c in num]
        assert len(hl)-1==k and hl[-1]==a*l**k//L and hl[-1]>0
        for D in range(1,Dmax+1):
            rD=rl(l*D)
            i=(r-rD)
            if i%l: bad+=1; print('r mismatch',l,D); continue
            i//=l
            if not (0<=i<D): bad+=1; print('i range',l,D,i); continue
            # T_l(Dx-i) = T_{lD}(x)
            assert r+l*(-i)==rD
            lhs=compose_affine(hl,-i,D)
            numD=compose_affine(h,rD,l*D); LD=lam(l*D)
            if any(c%LD for c in numD): bad+=1; print('NONINTEGRAL lD',l,D); continue
            hlD=[c//LD for c in numD]
            lD=lam(D)
            rhs=[lD*c for c in hlD]
            if lhs!=rhs: bad+=1; print('IDENTITY FAIL',hname,l,D)
            if lD%D or lD>D**k: bad+=1; print('lambda(D) cond fail',D)
            cnt+=1
    return cnt,bad,{p:m[p] for p in ps[:6]}

def lifting_check(hname,h,prefer,Lmax=40):
    # Lemma local-lifting for small primes: uniform S_p and uniform distribution of values on the class
    res={}
    bad=0
    for p in [2,3,5,7]:
        Smax=0
        for l in range(1,Lmax+1):
            # recompute h_l via check-style construction (reuse by small local function)
            pass
        res[p]=None
    return res

if __name__=='__main__':
    tests=[('(x-1)^2',[1,-2,1],[1]),('x^2-1',[-1,0,1],[1,-1]),('x^3-x',[0,-1,0,1],[1,-1]),
           ('(x-1)(x-2)',[2,-3,1],[1]),('(x-1)^2(x+1)',[1,-1,-1,1],[1,-1]),('(x-1)^3',[-1,3,-3,1],[1]),
           ('(x^2-2)(x^2-17)(x^2-34)',None,None)]
    # build h8
    h8=pmul(pmul([-2,0,1],[-17,0,1]),[-34,0,1])
    for name,h,pref in tests:
        if h is None: h=h8; pref=None
        print(name,check(name,h,pref),flush=True)

