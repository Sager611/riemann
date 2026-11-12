
import numpy as np, time
t0=time.time()
N=1000; n=np.arange(N+1,2*N+1); ln=np.log(n.astype(float))
def mobius(M):
    mu=np.ones(M+1,dtype=np.int64); mu[0]=0
    isp=np.ones(M+1,dtype=bool); isp[:2]=False
    for p in range(2,M+1):
        if isp[p]:
            isp[2*p::p]=False
            mu[p::p]*=-1
            mu[p*p::p*p]=0
    return mu
mu=mobius(2*N)
kap=(np.sqrt(10-2*np.sqrt(5))-2)/(np.sqrt(5)-1)
dh=np.array([0,1,kap,-kap,-1.0])[n%5]
chi4=np.array([0,1,0,-1.0])[n%4]
X=31
det=np.array([sum(mu[d] for d in range(1,X+1) if m%d==0) for m in n],dtype=float)
rng=np.random.default_rng(1)
fam={'zeta':(n/N)**-0.5,'mobius':mu[N+1:2*N+1].astype(float),'random':rng.choice([-1.0,1.0],N),
     'DH':dh,'chi4':chi4,'detector':det}
def Dabs2(a,t,chunk=2000):
    out=np.empty(len(t))
    for i in range(0,len(t),chunk):
        out[i:i+chunk]=np.abs(np.exp(-1j*np.outer(t[i:i+chunk],ln))@a)**2
    return out
def G2(a):
    c=np.zeros(4*N*N+1); 
    P=np.outer(n,n).ravel(); A=np.outer(a,a).ravel()
    np.add.at(c,P,A); return float(np.sum(c**2))
def greedy(t,mask):
    ts=t[mask]; cnt=0; last=-1e9
    for x in ts:
        if x-last>=1.0: cnt+=1; last=x
    return cnt
def cz(x,dt,lam):
    # dyadic CZ stopping intervals of avg(x) > lam, below the top level
    L=1
    while L<len(x): L*=2
    y=np.zeros(L); y[:len(x)]=x
    stops=[]; covered=np.zeros(L,dtype=bool)
    size=L//2
    while size>=1:
        av=y.reshape(-1,size).mean(axis=1)
        cov=covered.reshape(-1,size).any(axis=1)
        new=np.where((av>lam)&(~cov))[0]
        for j in new: stops.append(size); covered[j*size:(j+1)*size]=True
        size//=2
    return len(stops), (max(stops)*dt if stops else 0.0), sum(stops)*dt
sigs=[0.55,0.6,0.65,0.7,0.75]
for tau in (1.2,1.6):
    T=N**tau; dt=0.2; t=np.arange(0,T,dt)
    print(f"=== N={N} tau={tau} T={T:.0f}")
    print(" bounds (no constants, G=N): sigma: MVT  Huxley  GM  Montg N^(2-2s)")
    for s in sigs:
        mvt=(T+N)*N**(1-2*s); hux=N**(2-2*s)+T*N**(4-6*s); gm=N**(2-2*s)+N**(3.6-4*s)+T*N**(2.4-4*s)
        print(f"  {s:.2f}: {mvt:9.1f} {hux:9.1f} {gm:9.1f} {N**(2-2*s):7.1f}")
    for name,a in fam.items():
        G=float(np.sum(a*a)); g2=G2(a)
        x=Dabs2(a,t)
        mean=x.mean()/G; mx=x.max()/G; far=t>2*N
        line=f" {name:9s} G={G:6.0f} G2/G^2={g2/G**2:5.2f} mean|D|^2/G={mean:5.2f} max/G={mx:8.1f} max_far/G={x[far].max()/G:6.1f} |"
        for s in sigs:
            lam=G*N**(2*s-1)
            R=greedy(t,x>=lam); Rf=greedy(t[far],x[far]>=lam)
            ns,ml,tm=cz(x,dt,lam)
            Lstar=2*3*np.pi*N*G/(lam-G)
            line+=f" s{s}:R={R},far={Rf},CZ#={ns},maxI={ml:.1f}/L*={Lstar:.0f}"
        print(line,flush=True)
print("time",round(time.time()-t0,1))
