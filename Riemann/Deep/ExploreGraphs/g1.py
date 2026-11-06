# explore12 g1: quantum flower graphs with petals log p vs the explicit formula
import numpy as np, itertools, math, sys, time, json
from scipy.linalg import expm
from scipy.optimize import least_squares
from scipy.integrate import quad
np.random.seed(1)

def primes_upto(n):
    return [p for p in range(2,n+1) if all(p%q for q in range(2,int(p**0.5)+1))]

def factor_vec(N, P):
    a=[]
    for p in P:
        e=0
        while N%p==0: N//=p; e+=1
        a.append(e)
    return a if N==1 else None

def vonmangoldt(N):
    for p in primes_upto(N):
        if N%p==0:
            while N%p==0: N//=p
            return math.log(p) if N==1 else 0.0
    return 0.0

def flower_perm(E):
    # directed bond b=2j:(j,+) departs half-edge 2j, arrives half-edge 2j+1; b=2j+1 reverse
    D=2*E; Pi=np.zeros((D,D))
    for b in range(D):
        arr = b+1 if b%2==0 else b-1
        Pi[arr,b]=1.0
    return Pi

def coeff(S, E, alpha):
    """[z^alpha] Tr((S Z)^n), Z=diag(z_petal(b)), n=|alpha|"""
    n=sum(alpha); J=[j for j in range(E) if alpha[j]>0]
    D=2*E; tot=0j; cnt=0
    for th in itertools.product(range(n+1), repeat=len(J)):
        z=np.zeros(D,complex)
        ang=0.0
        for jj,j in enumerate(J):
            t=2*np.pi*th[jj]/(n+1); z[2*j]=z[2*j+1]=np.exp(1j*t); ang+=alpha[j]*t
        M=S*z[None,:]
        tot+=np.trace(np.linalg.matrix_power(M,n))*np.exp(-1j*ang); cnt+=1
    return tot/cnt

def smooth_list(P, Nmax):
    out=[]
    for N in range(2,Nmax+1):
        if factor_vec(N,P) is not None: out.append(N)
    return out

def weights(S, P, Nmax):
    E=len(P); W={}
    for N in smooth_list(P,Nmax):
        a=factor_vec(N,P); n=sum(a)
        c=coeff(S,E,a); W[N]=(math.log(N)/n*c.real, math.log(N)/n*c.imag, c)
    return W

# ---- spectrum via exact counting function
def spectrum(S, lens, K, dk=0.002, eps0=1e-9):
    D=len(lens); L=np.array(lens); TrL=L.sum()
    def args(k):
        mu=np.linalg.eigvals(S*np.exp(1j*k*L)[None,:])
        return np.mod(np.angle(mu),2*np.pi)
    Theta0=args(eps0).sum()-eps0*TrL
    def Ncount(k):
        return (k*TrL+Theta0-args(k).sum())/(2*np.pi)
    ks=np.arange(eps0,K+dk,dk)
    Nv=np.array([Ncount(k) for k in ks])
    Ni=np.round(Nv).astype(int)
    if np.max(np.abs(Nv-Ni))>1e-6: print('WARN count not integer', np.max(np.abs(Nv-Ni)))
    eig=[]
    for i in np.nonzero(np.diff(Ni))[0]:
        # bisect the jumps inside (ks[i], ks[i+1]]
        stack=[(ks[i],ks[i+1],Ni[i],Ni[i+1])]
        while stack:
            lo,hi,nlo,nhi=stack.pop()
            if nhi==nlo: continue
            if hi-lo<1e-12:
                eig+= [0.5*(lo+hi)]*(nhi-nlo); continue
            mid=0.5*(lo+hi); nm=int(round(Ncount(mid)))
            stack.append((lo,mid,nlo,nm)); stack.append((mid,hi,nm,nhi))
    m0=int(np.sum(np.abs(np.linalg.eigvals(S)-1)<1e-9))
    return np.sort(np.array(eig)), m0, TrL

# ---- test function: phi = Lambda_a * Lambda_a supported on [-X,X], h = FT
def make_test(X):
    a=X/2
    lam=lambda y: max(0.0,1-abs(y)/a)
    def phi(x):
        if abs(x)>=X: return 0.0
        lo=max(-a,x-a); hi=min(a,x+a)
        return quad(lambda y: lam(y)*lam(x-y), lo, hi, points=[0.0,x], epsabs=1e-14, epsrel=1e-13)[0]
    def h(k):  # entire; works for complex k
        u=np.asarray(k)*a/2
        s=np.where(np.abs(u)<1e-8, 1.0+0j, np.sin(u+0j)/np.where(np.abs(u)<1e-8,1,u))
        return (a*a/(2*np.pi))*s**4
    return phi,h,a

if __name__=='__main__':
    t0=time.time()
    out={}
    # ===== Part A: Kirchhoff flower, petals log2, log3, log5
    P=[2,3,5]; E=3; D=6
    sigma=np.ones((D,D))/E-np.eye(D)
    S=sigma@flower_perm(E)
    print('A: Kirchhoff flower petals log 2, log 3, log 5; S orthogonal:', np.allclose(S@S.T,np.eye(D)))
    W=weights(S,P,60)
    print('A: weights W_G(log N) [graph]  vs  W_zeta(log N) = -Lambda(N)/sqrt(N)')
    for N,(re,im,c) in W.items():
        print('   N=%3d  W_G=%+.6f  W_zeta=%+.6f  %s'%(N,re,-vonmangoldt(N)/math.sqrt(N),'composite' if vonmangoldt(N)==0 else ''))
    lens=[math.log(P[b//2]) for b in range(D)]
    eig,m0,TrL=spectrum(S,lens,K=400.0,dk=0.002)
    print('A: eigenvalues in (0,400]: %d, m0=%d, Weyl TrL*K/2pi=%.2f'%(len(eig),m0,TrL*400/(2*np.pi)))
    print('A: first eigenvalues', np.round(eig[:12],6))
    # trace formula check with X=log(12.5)
    X=math.log(12.5); phi,h,a=make_test(X)
    Wf=weights(S,P,12)
    for t in [0.0,3.0,14.134725,30.0]:
        ht=lambda k: h(k-t)+h(k+t)
        lhs=m0*ht(0).real+2*np.sum(ht(eig).real)
        rhs=TrL/(2*np.pi)*2*phi(0)+sum(Wf[N][0]*2*math.cos(t*math.log(N))*phi(math.log(N)) for N in Wf)/np.pi
        print('A: trace formula t=%7.3f  spectral=%.10f  orbits=%.10f  diff=%.2e'%(t,lhs,rhs,lhs-rhs))
    print('time',round(time.time()-t0,1))
