#!/usr/bin/env python3
# Lemma P check: for every discriminant Delta (all orders), the Hermitian form
#   Q'_Delta(b) = sum_{n,n'} b(n) conj b(n') sum_{c|(n,n')} (Delta/c) * W_Delta(n n'/c^2)
# (W = half the number of gamma in O_Delta with N gamma = M, no coprimality to the conductor)
# equals the sum over squarefree e | f of blocks for the order Delta/e^2 with b_e(m) = b(em)1[(m,f)=1],
# and each block is >= 0.  Variants: plain (D<0, blocks also computed independently from reduced forms),
# sectors (weight T_l(T/(2 sqrt M))), D>0 (weight w(arccosh(T/(2 sqrt M))), Gaussian w).
import math, cmath, random, signal, sys
signal.alarm(900)

def gcd(a,b):
    while b: a,b=b,a%b
    return abs(a)

def jacobi(a, n):
    assert n>0 and n%2==1
    a %= n; r = 1
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5): r=-r
        a,n = n,a
        if a%4==3 and n%4==3: r=-r
        a %= n
    return r if n==1 else 0

def isdisc(D): return D%4 in (0,1)

def conductor(D):
    f=1; Dc=D
    p=2
    while p*p<=abs(Dc):
        while Dc%(p*p)==0 and isdisc(Dc//(p*p)):
            Dc//= p*p; f*=p
        p+=1
    return f

def primes_of(n):
    ps=[];p=2
    while p*p<=n:
        if n%p==0:
            ps.append(p)
            while n%p==0: n//=p
        p+=1
    if n>1: ps.append(n)
    return ps

def divisors(n):
    ds=[1]
    for p in primes_of(n):
        ds=ds+[d*p for d in ds]
    return ds

def isqrt(n): return math.isqrt(n)

def wunits(D): return 6 if D==-3 else (4 if D==-4 else 2)

def forms(D):  # reduced primitive positive definite forms of disc D<0
    out=[]; a=1
    while 3*a*a<=-D:
        for b in range(-a+1,a+1):
            if (b*b-D)%(4*a)==0:
                c=(b*b-D)//(4*a)
                if c<a: continue
                if b<0 and a==c: continue
                if gcd(gcd(a,b),c)!=1: continue
                out.append((a,b,c))
        a+=1
    return out

def reps(form, mmax):
    a,b,c=form; D=b*b-4*a*c; cnt={}
    ymax=isqrt(4*a*mmax//(-D))+1
    for y in range(-ymax,ymax+1):
        # a x^2 + b x y + c y^2 <= mmax
        disc=b*b*y*y-4*a*(c*y*y-mmax)
        if disc<0: continue
        s=math.sqrt(disc)
        x0=math.floor((-b*y-s)/(2*a))-1; x1=math.ceil((-b*y+s)/(2*a))+1
        for x in range(x0,x1+1):
            m=a*x*x+b*x*y+c*y*y
            if 1<=m<=mmax: cnt[m]=cnt.get(m,0)+1
    return cnt

# weights on gamma = (T + y sqrt(Delta))/2, N gamma = M
def wsum_neg(D, M, kind, ell=0):
    # D<0: sum over (T,y) with T^2 - D y^2 = 4M, T = yD mod 2, of weight; returns half of it
    tot=0.0
    ymax=isqrt(4*M//(-D))
    for y in range(-ymax,ymax+1):
        t2=4*M+D*y*y
        if t2<0: continue
        t=isqrt(t2)
        if t*t!=t2: continue
        for T in ((t,-t) if t else (0,)):
            if (T-y*D)%2: continue
            if kind=='plain': tot+=1.0
            else:
                x=T/(2*math.sqrt(M)); x=max(-1.0,min(1.0,x))
                tot+=math.cos(ell*math.acos(x))
    return tot/2

SG=0.35  # Gaussian width for D>0
def wsum_pos(D, M):
    U=8*SG; tot=0.0
    ymax=int(2*math.sqrt(M)*math.sinh(U)/math.sqrt(D))+1
    for y in range(-ymax,ymax+1):
        t2=4*M+D*y*y
        t=isqrt(t2)
        if t*t!=t2 or t<=0: continue
        T=t
        if (T-y*D)%2: continue
        u=math.acosh(max(1.0,T/(2*math.sqrt(M))))
        tot+=math.exp(-u*u/(2*SG*SG))
    return tot

def form_value(D, b, ws, supp):
    # sum_{n,n'} b b' sum_c (D/c) ws(D, n n'/c^2)
    cache={}; tot=0j
    for n in supp:
        for n2 in supp:
            g=gcd(n,n2); s=0.0
            for c in divisors(g):
                j=jacobi(D,c)
                if j==0: continue
                M=n*n2//(c*c)
                if M not in cache: cache[M]=ws(D,M)
                s+=j*cache[M]
            tot+=b[n]*b[n2].conjugate()*s
    return tot

def blocks(D, b, ws, supp):
    f=conductor(D); out=[]
    fp=[p for p in primes_of(f) if p%2==1]
    es=[1]
    for p in fp: es=es+[e*p for e in es]
    for e in es:
        Dp=D//(e*e)
        be={}
        for n in supp:
            if n%e==0:
                m=n//e
                if gcd(m,f)==1: be[m]=b[n]
        if not be: out.append((e,0j)); continue
        out.append((e, form_value(Dp, be, ws, sorted(be))))
    return out

def blocks_forms(D, b, supp):
    # independent: (w'/2) sum_C |S_C|^2 for Delta' = D/e^2 via reduced forms
    f=conductor(D); out=[]
    fp=[p for p in primes_of(f) if p%2==1]
    es=[1]
    for p in fp: es=es+[e*p for e in es]
    for e in es:
        Dp=D//(e*e); wp=wunits(Dp)
        be={}
        for n in supp:
            if n%e==0:
                m=n//e
                if gcd(m,f)==1: be[m]=b[n]
        if not be: out.append((e,0.0)); continue
        mmax=max(be); tot=0.0
        for Q in forms(Dp):
            r=reps(Q,mmax)
            S=sum(be[m]*r.get(m,0) for m in be)/wp
            tot+=abs(S)**2
        out.append((e, wp/2*tot))
    return out

P=[3,5,7,11,13,17]
supp=[1]
for p in P: supp=supp+[n*p for n in supp]
supp=sorted(n for n in supp if n<=260)
print("support", len(supp), supp)

rng=random.Random(20261021)
maxdev=0.0; minblock=0.0; ntests=0; nnonfund=0; nwithcommon=0; extra_pos=0
mode=sys.argv[1] if len(sys.argv)>1 else 'plain'
if mode in ('plain','sector'):
    Ds=[D for D in range(-7,-1601,-1) if isdisc(D)]
else:
    Ds=[D for D in range(5,1201) if isdisc(D) and isqrt(D)**2!=D]
ells=[2,4,6] if mode=='sector' else [0]
for D in Ds:
    f=conductor(D)
    nonfund = f>1
    if not nonfund and rng.random()<0.85: continue   # sample fundamental ones sparsely
    for ell in ells:
        b={n: complex(rng.gauss(0,1),rng.gauss(0,1)) for n in supp}
        if mode=='plain': ws=lambda DD,M: wsum_neg(DD,M,'plain')
        elif mode=='sector': ws=lambda DD,M,ell=ell: wsum_neg(DD,M,'sector',ell)
        else: ws=wsum_pos
        L=form_value(D,b,ws,supp)
        B=blocks(D,b,ws,supp)
        R=sum(v for e,v in B)
        dev=abs(L-R)/(abs(L)+1e-9)
        maxdev=max(maxdev,dev)
        for e,v in B:
            minblock=min(minblock, v.real/(abs(L)+1e-9))
            if abs(v.imag)>1e-6*(abs(L)+1): print("imag block",D,e,v)
        if mode=='plain' and nonfund:
            BF=blocks_forms(D,b,supp)
            for (e,v),(e2,v2) in zip(B,BF):
                d=abs(v-v2)/(abs(v)+1e-9)
                maxdev=max(maxdev,d)
                if d>1e-9: print("FORMS MISMATCH",D,e,v,v2)
        ntests+=1; nnonfund+=nonfund
        if len(B)>1 and any(abs(v)>1e-9 for e,v in B[1:]): extra_pos+=1
print(mode, "tests", ntests, "non-fundamental", nnonfund, "with nonzero e>1 blocks", extra_pos)
print("max relative deviation (identity, and forms check for plain):", maxdev)
print("min block / |Q'| (should be >= -1e-12):", minblock)

