# T3: strong convergence needs a non-amenable limit.  Abelian family (Dirichlet characters mod q, Cayley graph of (Z/q)^*
# with generators p^{+-1}, p <= y) vs random 2k-regular Schreier graphs of S_n (random permutation reps of F_k).
import numpy as np, math, sys
def primes_upto(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
def primroot(q):
    f=[]; m=q-1; d=2
    while d*d<=m:
        if m%d==0:
            f.append(d)
            while m%d==0: m//=d
        d+=1
    if m>1: f.append(m)
    for g in range(2,q):
        if all(pow(g,(q-1)//p,q)!=1 for p in f): return g
def indtable(q):
    g=primroot(q); ind=np.zeros(q,dtype=np.int64); x=1
    for k in range(q-1):
        ind[x]=k; x=x*g%q
    return ind
P=primes_upto(4*10**6)
print("(a) fixed generator set {p <= y}, q -> infinity: max_{chi != chi0} 2 Re sum_p chi(p) / (2k)   [amenable limit: -> 1]")
for y in (30,100):
    gens=P[P<=y]; k=len(gens)
    line=[]
    for q in (1009,10007,100003,1000003):
        ind=indtable(q); v=np.zeros(q-1)
        for p in gens:
            if p!=q: v[ind[p%q]]+=1
        S=np.fft.fft(v)   # S[j] = sum_p exp(-2 pi i j ind(p)/(q-1)) = sum_p chi_j(p)-bar
        ev=2*np.real(S[1:])
        line.append(f"q={q}: {ev.max()/(2*k):.4f}")
    kes=math.sqrt(2*k-1)/k
    print(f"  y={y} k={k}: "+"  ".join(line)+f"   | Kesten (free group, non-amenable) 2sqrt(2k-1)/2k = {kes:.4f}")
print("(b) random permutation reps of F_k (n x n, k permutations, A = sum P_i + P_i^T): largest nontrivial eigenvalue / 2k")
rng=np.random.default_rng(1)
for k in (10,25):
    out=[]
    for n in (200,800,2000):
        A=np.zeros((n,n))
        for _ in range(k):
            s=rng.permutation(n); A[np.arange(n),s]+=1; A[s,np.arange(n)]+=1
        e=np.linalg.eigvalsh(A); out.append(f"n={n}: {max(abs(e[0]),e[-2])/(2*k):.4f}")
    print(f"  k={k}: "+"  ".join(out)+f"   | Kesten {math.sqrt(2*k-1)/k:.4f}")
print("(c) growing generator set y = q^theta: max_{chi!=chi0}|sum_{p<=y} chi(p)|/pi(y) vs Steinhaus union bound sqrt(log q / pi(y)) and Kesten sqrt(2/k)")
for q in (100003,1000003):
    ind=indtable(q)
    for th in (0.2,0.35,0.5,0.75,1.0,1.3):
        y=int(q**th)
        if y>4*10**6: continue
        gens=P[(P<=y)&(P!=q)]; k=len(gens)
        v=np.bincount(ind[gens%q],minlength=q-1).astype(float)
        S=np.fft.fft(v); m=np.abs(S[1:]).max()
        print(f"  q={q} theta={th:.2f} y={y} k={k}: max/k={m/k:.4f}  sqrt(log q/k)={math.sqrt(math.log(q)/k):.4f}  Kesten~sqrt(2/k)={math.sqrt(2/k):.4f}")
