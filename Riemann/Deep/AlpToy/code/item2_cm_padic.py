#!/usr/bin/env python3
"""alptoy item 2 (note item 2): 'D^p = 0 on a CM curve for p inert in K, D the invariant derivation, and
D^a(Omega theta(0)) is divisible by p^(a/p + a/p^3 + ... - O(log a))'.
D = d/dz (z = formal logarithm) acts on Z_p[[t]] as omega(t)^(-1) d/dt (omega = invariant differential, omega(0)=1).
Computed exactly mod M = p^prec:  F1: D^a t (0) = a! [z^a] exp_E(z)  (t = -x/y);  F2: D^a sigma(0), sigma = Weierstrass sigma
(Weierstrass' integer recursion a_{m,n}; sigma = sum a_{mn} (g2/2)^m (2g3)^n z^(4m+6n+1)/(4m+6n+1)!).
Also: D^p t == 0 mod p as a whole series (<=> D^p == 0 mod p as a derivation)."""
import sys, math, json, time
from flint import fmpz_poly
import mpmath as mp
AMAX=int(sys.argv[1]); OUT=sys.argv[2]
def red(f,M,n=None):
    c=f.coeffs()
    if n is not None: c=c[:n]
    return fmpz_poly([int(x)%M for x in c])
def mul(f,g,n,M): return red(f*g,M,n)
def inv(f,n,M):
    g=fmpz_poly([pow(int(f.coeffs()[0])%M,-1,M)]); pr=1
    while pr<n:
        pr=min(2*pr,n); e=mul(f,g,pr,M); g=mul(g,fmpz_poly([2])-e,pr,M)
    return g
def der(f): c=f.coeffs(); return fmpz_poly([int(c[i])*i for i in range(1,len(c))])
def coef(f,i): c=f.coeffs(); return int(c[i]) if i<len(c) else 0
def vp(n,p):
    v=0
    while n%p==0: n//=p; v+=1
    return v
T=fmpz_poly([0,1])
def formal_w(ai,n,M):
    a1,a2,a3,a4,a6=ai; t=T; t2=mul(t,t,n,M); t3=mul(t2,t,n,M)
    w=t3
    for it in range(14):
        w2=mul(w,w,n,M); w3=mul(w2,w,n,M)
        F=w-t3-a1*mul(t,w,n,M)-a2*mul(t2,w,n,M)-a3*w2-a4*mul(t,w2,n,M)-a6*w3
        F=red(F,M,n)
        if F.degree()<0: break
        dF=fmpz_poly([1])-a1*t-a2*t2-2*a3*w-2*a4*mul(t,w,n,M)-3*a6*w2
        w=red(w-mul(F,inv(red(dF,M,n),n,M),n,M),M,n)
    w2=mul(w,w,n,M); w3=mul(w2,w,n,M)
    F=red(w-t3-a1*mul(t,w,n,M)-a2*mul(t2,w,n,M)-a3*w2-a4*mul(t,w2,n,M)-a6*w3,M,n)
    assert F.degree()<0, 'w not converged'
    return w
def omega(ai,p,n,M):
    a1,a2,a3,a4,a6=ai; w=formal_w(ai,n+4,M); dw=der(w); t=T
    def shift(f,k): c=f.coeffs(); return fmpz_poly([int(x) for x in c[k:]])
    res=[]
    if p!=2:  # omega/dt = (w - t w')/(w (-2 + a1 t + a3 w))
        u=shift(red(w-mul(t,dw,n+4,M),M),3); v=shift(w,3)
        den=red(fmpz_poly([-2])+a1*t+a3*w,M,n)
        res.append(mul(mul(u,inv(red(v,M,n),n,M),n,M),inv(den,n,M),n,M))
    if p!=3:  # omega/dt = w'/(3t^2 + 2a2 t w + a4 w^2 + a1 w)
        den=red(3*mul(t,t,n+4,M)+2*a2*mul(t,w,n+4,M)+a4*mul(w,w,n+4,M)+a1*w,M)
        res.append(mul(shift(dw,2),inv(red(shift(den,2),M,n),n,M),n,M))
    if len(res)==2: assert red(res[0]-res[1],M).degree()<0, 'omega formulas disagree'
    om=res[0]; assert coef(om,0)%M==1
    return om
def Dvals(ai,p,amax,M):
    L=amax+2; om=omega(ai,p,L,M); oi=inv(om,L,M)
    f=T; vals=[]; Dp_zero=None
    for j in range(amax+1):
        vals.append(coef(f,0)%M)
        if j==p: Dp_zero=all(int(c)%p==0 for c in f.coeffs())
        f=mul(oi,der(f),L-j-1,M)
    return vals,Dp_zero,om
# Weierstrass sigma integers a_{m,n} (assert integrality)
def sigma_table(amax):
    A={}
    W=amax-1
    for w in range(0,W+1,2):
        for n in range(0,w//6+1):
            if (w-6*n)%4: continue
            m=(w-6*n)//4
            if m==0 and n==0: A[(0,0)]=1; continue
            x=3*(m+1)*A.get((m+1,n-1),0)
            y=16*(n+1)*A.get((m-2,n+1),0)-(2*m+3*n-1)*(4*m+6*n-1)*A.get((m-1,n),0)
            assert y%3==0, ('nonintegral',m,n)
            A[(m,n)]=x+y//3
    return A
def sigma_vals(A,g2h,g3d,amax,M):  # g2h=g2/2, g3d=2g3 as residues mod M
    out={}
    for (m,n),c in A.items():
        a=4*m+6*n+1
        if a<=amax: out[a]=(out.get(a,0)+c*pow(g2h,m,M)*pow(g3d,n,M))%M
    return out
def cinv(x,M): return pow(x%M,-1,M)
def frac_mod(num,den,M): return num*cinv(den,M)%M
def c4c6(ai):
    a1,a2,a3,a4,a6=ai
    b2=a1*a1+4*a2; b4=2*a4+a1*a3; b6=a3*a3+4*a6
    return b2*b2-24*b4, -b2**3+36*b2*b4-216*b6
mp.mp.dps=60
def H23():
    s=mp.sqrt(23)*1j; taus=[(-1+s)/2,(-1+s)/4,(1+s)/4]
    js=[1728*mp.kleinj(t) for t in taus]
    c=[1]
    for j in js:
        c=[ (c[i] if i<len(c) else 0) - j*(c[i-1] if i>=1 else 0) for i in range(len(c)+1)]
    return [int(mp.nint(mp.re(x))) for x in c], max(abs(mp.im(x)) for x in c)
Hc,imax=H23()
Hc=Hc[::-1]  # low->high
print('H_{-23} from kleinj (low->high):',Hc,'max imag %.1e'%imax)
assert Hc==[12771880859375,-5151296875,3491750,1]
def hens(poly,p,M):
    rts=[x for x in range(p) if sum(c*x**i for i,c in enumerate(poly))%p==0]
    out=[]
    for x in rts:
        dp=sum(i*c*x**(i-1) for i,c in enumerate(poly) if i)%p
        if dp==0: out.append((x,None)); continue
        for _ in range(12):
            fx=sum(c*pow(x,i,M) for i,c in enumerate(poly))%M; dfx=sum(i*c*pow(x,i-1,M) for i,c in enumerate(poly) if i)%M
            x=(x-fx*cinv(dfx,M))%M
        assert sum(c*pow(x,i,M) for i,c in enumerate(poly))%M==0
        out.append((x,True))
    return out
CASES=[('Q(i)','y^2=x^3-x',[0,0,0,-1,0],[3,7,11,19,23],[5,13]),
       ('Q(sqrt-3)','y^2+y=x^3',[0,0,1,0,0],[2,5,11,17,23,29],[7,13]),
       ('Q(sqrt-7)','49a1 y^2+xy=x^3-x^2-2x-1',[1,-1,0,-2,-1],[3,5,13,17,19],[11,2]),
       ('Q(sqrt-23)','E_j, j root of H_-23 in Z_p',None,[5,7,11,17,19],[])]
t0=time.time(); A=sigma_table(AMAX); print('sigma a_{m,n} integral: %d entries, %.1fs'%(len(A),time.time()-t0))
rows=[]
for K,model,ai0,inert,split in CASES:
    for p in inert+split:
        prec=AMAX//(p-1)+80; M=p**prec; kind='inert' if p in inert else 'split'
        ai=ai0; note=''
        if ai is None:
            rts=hens(Hc,p,M)
            good=[x for x,ok in rts if ok and x%p not in (0,1728%p)]
            if not good: print(K,p,'no usable root (roots mod p:',[x for x,_ in rts],')'); continue
            j=good[0]; d=cinv(j-1728,M); ai=[1,0,0,(-36*d)%M,(-d)%M]; note='j=%d mod %d (roots mod p: %s)'%(j%p,p,[x%p for x,_ in rts])
        t1=time.time()
        vals,Dpz,om=Dvals(ai,p,AMAX,M)
        c4,c6=c4c6(ai)
        sv=None
        if p>=5 or (p==3 and K=='Q(i)'):
            from fractions import Fraction as Fr
            if ai0 is not None:
                q2,q3=Fr(c4,24),Fr(c6,108)
                assert q2.denominator%p and q3.denominator%p
                g2h=q2.numerator*cinv(q2.denominator,M)%M; g3d=q3.numerator*cinv(q3.denominator,M)%M
            else:
                g2h=frac_mod(c4,24,M); g3d=frac_mod(c6,108,M)  # g2/2=c4/24, 2g3=c6/108
            sv=sigma_vals(A,g2h,g3d,AMAX,M)
        for fname,vv in [('t',{a:vals[a] for a in range(1,AMAX+1)}),('sigma',sv)]:
            if vv is None: continue
            V={a:(vp(x,p) if x else None) for a,x in vv.items()}
            rows.append(dict(K=K,p=p,kind=kind,f=fname,prec=prec,Dp_zero=Dpz,note=note,V=V))
        print('%-10s p=%2d %-5s D^p t==0 mod p: %-5s  prec=%d  %.1fs %s'%(K,p,kind,Dpz,prec,time.time()-t1,note),flush=True)
json.dump(rows,open(OUT,'w'))
