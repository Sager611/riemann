# explore12 g4: adversarial (a): window realisation of the Davenport-Heilbronn "prime side" by a directed flower
import numpy as np, math, time
from g1 import spectrum, make_test
from g2 import dh_coeffs

def perfect_power_base(N):
    for b in range(2,N+1):
        m=0; x=N
        while x%b==0: x//=b; m+=1
        if x==1 and m>=1: return b,m
    return N,1

def contraction_with_powersums(s, M, maxcopies=64):
    """real contraction B (normal) with Tr B^m = s[m-1], m=1..M; uses q copies of the roots of exp(-sum s_m y^m/(q m))"""
    for q in range(1,maxcopies+1):
        # truncated exp of -sum (s_m/q) y^m/m to degree M -> polynomial f(y) = prod(1 - lam y)
        a=np.zeros(M+1); a[0]=1.0
        g=np.zeros(M+1)
        for m in range(1,M+1): g[m]=-s[m-1]/q/m
        # f = exp(g) via recursion f' = g' f
        for k in range(1,M+1):
            a[k]=sum(j*g[j]*a[k-j] for j in range(1,k+1))/k
        # roots of y^M f(1/y) = sum a_k y^(M-k): these are the lam's
        lam=np.roots(a[::-1][::-1]) if M>0 else np.array([])
        lam=np.roots(a)  # numpy: a[0] y^M + ... + a[M] -> roots are lam
        if np.all(np.abs(lam)<=1-1e-9):
            blocks=[]
            used=np.zeros(len(lam),bool)
            for i,l in enumerate(lam):
                if used[i]: continue
                used[i]=True
                if abs(l.imag)<1e-12: blocks.append(np.array([[l.real]]))
                else:
                    j=[k for k in range(len(lam)) if not used[k] and abs(lam[k]-np.conj(l))<1e-9][0]; used[j]=True
                    blocks.append(np.array([[l.real,-l.imag],[l.imag,l.real]]))
            B1=np.zeros((0,0))
            for b in blocks:
                B1=np.block([[B1,np.zeros((B1.shape[0],b.shape[1]))],[np.zeros((b.shape[0],B1.shape[1])),b]])
            B=np.kron(np.eye(q),B1)
            return B,q,np.max(np.abs(lam))
    return None,None,None

def dilate(T):
    n=T.shape[0]; U,s,Vt=np.linalg.svd(T); J=np.nonzero(s<1-1e-12)[0]; d=len(J)
    S=np.zeros((n+d,n+d)); S[:n,:n]=T; S[:n,n:]=U[:,J]*np.sqrt(1-s[J]**2)[None,:]
    S[n:,:n]=(np.sqrt(1-s[J]**2)[:,None])*Vt[J,:]; S[n:,n:]=-np.diag(s[J]); return S,d

if __name__=='__main__':
    t0=time.time(); X=math.log(16.5); laux=2.9; Nmax=16
    a,Lam,kap=dh_coeffs(Nmax)
    W={N:-Lam[N]/math.sqrt(N) for N in range(2,Nmax+1)}
    classes={}
    for N in range(2,Nmax+1):
        b,m=perfect_power_base(N); classes.setdefault(b,{})[m]=W[N]
    blocks=[]; lens=[]
    print('DH window realisation, X=log 16.5, kappa=%.6f'%kap)
    for b,ws in sorted(classes.items()):
        M=int(math.floor(X/math.log(b)+1e-12))
        s=[ws.get(m,0.0)/math.log(b) for m in range(1,M+1)]
        if all(abs(x)<1e-14 for x in s): continue
        B,q,lmax=contraction_with_powersums(s,M)
        tr=[np.trace(np.linalg.matrix_power(B,m)) for m in range(1,M+1)]
        print('   class log %2d: M=%d targets %s  got %s  bonds=%d copies=%d max|lam|=%.3f'%(b,M,' '.join('%+.4f'%x for x in s),' '.join('%+.4f'%x for x in tr),B.shape[0],q,lmax))
        blocks.append(B); lens+=[math.log(b)]*B.shape[0]
    n=sum(B.shape[0] for B in blocks); T=np.zeros((n,n)); i=0
    for B in blocks:
        m=B.shape[0]; T[i:i+m,i:i+m]=B; i+=m
    S,d=dilate(T); lens=lens+[laux]*d
    print('   short bonds %d (length %.4f vs zeta psi(16)=%.4f), aux %d, max|SS^T-I|=%.1e, TrL=%.3f'%(n,sum(lens[:n]),math.log(720720),d,np.abs(S@S.T-np.eye(n+d)).max(),sum(lens)))
    K=400.0
    eig,m0,TrL=spectrum(S,lens,K=K,dk=0.002)
    print('   eigenvalues (0,%g]: %d (Weyl %.1f) m0=%d time %.1f'%(K,len(eig),TrL*K/(2*np.pi),m0,time.time()-t0))
    phi,h,aa=make_test(X)
    for t in [0.0,5.0,20.0,50.0]:
        ht=lambda k: h(k-t)+h(k+t)
        G=m0*ht(0).real+2*np.sum(ht(eig).real)
        R=TrL/(2*np.pi)*2*phi(0)+sum(W[N]*2*math.cos(t*math.log(N))*phi(math.log(N)) for N in W)/math.pi
        print('   t=%5.1f spectral %.10f  Weyl+DH-weights %.10f  diff %.1e'%(t,G,R,G-R))
    print('time',round(time.time()-t0,1))
