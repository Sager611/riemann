# Test of Section 3 (table construction): U_q field of degree q; omega_2; omega_q via V_q = U_q (x) K_q
import random, sys
from ffpoly import *
from unitroot import ddf
import numpy as np

def isprime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True

def aux_prime(p,q):
    l=12*q+1
    while True:
        if isprime(l) and l not in (2,3,p,q) and pow(p,(l-1)//q,l)!=1: return l
        l+=12*q

def prim_root_mod(l):
    fac=[]; n=l-1; d=2
    while d*d<=n:
        if n%d==0:
            fac.append(d)
            while n%d==0: n//=d
        d+=1
    if n>1: fac.append(n)
    for g in range(2,l):
        if all(pow(g,(l-1)//f,l)!=1 for f in fac): return g


def fold_phi(f, l, p):
    # reduce polynomial f (list) modulo Phi_l = 1+T+...+T^{l-1}
    a = np.zeros(l, dtype=np.int64)
    if f:
        arr = np.array(f, dtype=np.int64) % p
        n = len(arr)
        pad = (-n) % l
        if pad: arr = np.concatenate([arr, np.zeros(pad, dtype=np.int64)])
        a = arr.reshape(-1, l).sum(axis=0) % p
    top = a[l-1]
    r = (a[:l-1] - top) % p
    return norm([int(x) for x in r])

def fixed_algebra(p,q,l):
    # R = F_p[T]/Phi_l ; basis T^0..T^{l-2}; tau_a: T -> T^a
    Phi=[1]*l
    n=l-1
    g=prim_root_mod(l)
    a=pow(g,q,l)   # generator of H_q
    # matrix of tau_a
    cols=[]
    for j in range(n):
        img=[0]*(a*j+1); img[a*j]=1
        img=fold_phi(img,l,p)
        cols.append(img+[0]*(n-len(img)))
    A=(np.array(cols,dtype=np.int64).T - np.eye(n,dtype=np.int64))%p
    ker=nullspace(A,n,p)
    return Phi,[norm([int(x) for x in v]) for v in ker]

def minpoly_mod(x, modulus, p, maxdeg):
    # minimal polynomial of x in F_p[T]/modulus via linear dependence of powers
    pw=[[1]]
    for i in range(maxdeg+1):
        pw.append(pmod(pmul(pw[-1],x,p),modulus,p))
    n=deg(modulus)
    for d in range(1,maxdeg+1):
        M=np.array([ (v+[0]*(n-len(v)))[:n] for v in pw[:d+1]],dtype=np.int64).T%p
        ns=nullspace(M,d+1,p)
        if ns:
            c=[int(t) for t in ns[0]]
            return monic(norm(c),p)
    return None

def is_irreducible(f,p):
    D=ddf(monic(f,p),p)
    return list(D.keys())==[deg(f)]

def order_mod(x, modulus, p, n):
    # multiplicative order dividing n of x in F_p[T]/modulus (brute via factorization of n)
    fac=[]; m=n; d=2
    while d*d<=m:
        while m%d==0: fac.append(d); m//=d
        d+=1
    if m>1: fac.append(m)
    o=n
    for f in set(fac):
        while o%f==0 and ppowmod(x,o//f,modulus,p)==[1]:
            o//=f
    return o

random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 0)
results=[]
for p in [5,7,11,13,17,19,23,29,31,37,41,43,61,101,103]:
    for q in ([2,3,5,7] if p<=43 else [2,3,5]):
        if q>=p: continue
        l=aux_prime(p,q)
        Phi,basis=fixed_algebra(p,q,l)
        assert len(basis)==q, ('dim U_q',p,q,l,len(basis))
        # random element of U_q; its minimal polynomial over F_p
        ok=False
        for _ in range(5):
            x=[0]
            for v in basis: x=padd(x,pscal(v,random.randrange(p),p),p)
            mp=minpoly_mod(x,Phi,p,q)
            if mp and deg(mp)==q:
                assert is_irreducible(mp,p)
                ok=True; break
        assert ok, ('no generator found',p,q)
        if q==2:
            # solve v^p = -v in U_2: linear in the 2-dim space spanned by basis
            imgs=[]
            for v in basis:
                w=fold_phi(ppowmod(v,p,Phi,p),l,p)
                imgs.append(padd(w,v,p))   # v^p + v
            n=l-1
            M=np.array([ (w+[0]*(n-len(w))) for w in imgs],dtype=np.int64).T%p
            ns=nullspace(M,2,p)
            assert len(ns)==1
            v=[0]
            for c,b in zip(ns[0],basis): v=padd(v,pscal(b,int(c),p),p)
            c=fold_phi(pmul(v,v,p),l,p)
            assert deg(c)<=0
            cval=c[0] if c else 0
            s=0;t=p-1
            while t%2==0: t//=2; s+=1
            w2=pow(cval,t,p)
            assert pow(w2,2**(s-1),p)!=1 and pow(w2,2**s,p)==1, 'omega_2 order'
            results.append((p,q,l,'omega2 ok'))
        else:
            # K_q = F_p[Z]/phi, phi irreducible factor of Phi_q (found by brute-force factoring here)
            Phiq=[1]*q
            D=ddf(Phiq,p)
            (dq,Hd),=D.items()
            # split Hd into irreducible factors of degree dq (equal-degree; Cantor-Zassenhaus, testing aid only)
            def edf(f):
                if deg(f)==dq: return [f]
                while True:
                    a=norm([random.randrange(p) for _ in range(deg(f))])
                    if not a: continue
                    e=(p**dq-1)//2
                    g=pgcd(psub(ppowmod(a,e,f,p),[1],p),f,p)
                    if 0<deg(g)<deg(f): return edf(g)+edf(pdiv_exact(f,g,p))
            phi=sorted(edf(Hd))[0]
            Qq=p**dq
            sq=0; tq=Qq-1
            while tq%q==0: tq//=q; sq+=1
            # V_q elements: lists (length dq) of R_l elements (coeffs of Z^j)
            def vmul(x,y):
                r=[[] for _ in range(2*dq-1)]
                for i,a in enumerate(x):
                    if not a: continue
                    for j,b in enumerate(y):
                        if not b: continue
                        r[i+j]=padd(r[i+j],fold_phi(pmul(a,b,p),l,p),p)
                # reduce Z^k for k>=dq using phi (monic)
                for k in range(2*dq-2,dq-1,-1):
                    c=r[k]
                    if not c: continue
                    r[k]=[]
                    for j in range(dq):
                        coef=phi[j]
                        if coef: r[k-dq+j]=psub(r[k-dq+j],pscal(c,coef,p),p)
                return r[:dq]
            def vpow(x,e):
                res=[[1]]+[[] for _ in range(dq-1)]
                while e>0:
                    if e&1: res=vmul(res,x)
                    e>>=1
                    if e: x=vmul(x,x)
                return res
            # F_p-basis of V_q: basis_U[i] * Z^j
            Vb=[]
            for j in range(dq):
                for b in basis:
                    e=[[] for _ in range(dq)]; e[j]=b; Vb.append(e)
            nV=len(Vb); n=l-1
            def vec(x):
                # coordinates in the F_p basis: solve via U-basis coordinates of each Z^j coefficient
                out=[]
                for j in range(dq):
                    c=x[j]
                    # express c in basis of U_q (linear solve)
                    Bm=np.array([ (b+[0]*(n-len(b))) for b in basis],dtype=np.int64).T%p
                    sol=solve_lin(Bm,(c+[0]*(n-len(c))),p)
                    assert sol is not None, 'not in U_q'
                    out+= [int(t) for t in sol]
                return out
            zeta=[[] for _ in range(dq)]
            if dq>1: zeta[1]=[1]
            else: zeta[0]=[(-phi[0])%p]
            cols=[]
            for e in Vb:
                img=vpow(e,Qq)
                rhs=vmul(zeta,e)
                diff=[psub(a,b,p) for a,b in zip(img,rhs)]
                cols.append(vec(diff))
            M=np.array(cols,dtype=np.int64).T%p
            ns=nullspace(M,nV,p)
            assert ns, 'no eigenvector'
            u=[[] for _ in range(dq)]
            for c,e in zip(ns[0],Vb):
                if int(c): u=[padd(a,pscal(b,int(c),p),p) for a,b in zip(u,e)]
            w=vpow(u,(Qq**q-1)//q**(sq+1))
            om=vpow(w,q)
            # omega_q must lie in K_q: each Z^j coefficient is a constant of R_l
            assert all(deg(c)<=0 for c in om), 'omega_q not in K_q'
            # order check in K_q = F_p[Z]/phi
            omK=norm([ (c[0] if c else 0) for c in om])
            o=order_mod(omK,phi,p,Qq-1)
            assert o==q**sq, ('omega_q order',o,q**sq)
            results.append((p,q,l,dq,'omega_q ok')); print(results[-1],flush=True)
print(len(results),'table entries built and verified')
print(results[:6])
