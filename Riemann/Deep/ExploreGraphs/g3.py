# explore12 g3: explicit window realisation F_X of zeta's prime side by a directed flower (momentum operator)
import numpy as np, math, time, sys, json
from scipy.special import digamma
from g1 import primes_upto, vonmangoldt, spectrum, make_test

def companion_block(r, c):
    # contraction c*R_r, R_r = cyclic shift on R^(r+1) restricted to the sum-zero subspace:
    # eigenvalues c*w, w^(r+1)=1, w!=1, so Tr B^m = c^m((r+1)[r+1|m]-1); all singular values = c
    Sh=np.roll(np.eye(r+1),1,axis=0)
    Q,_=np.linalg.qr(np.hstack([np.ones((r+1,1)),np.random.default_rng(0).normal(size=(r+1,r))]))
    Q=Q[:,1:]
    return c*(Q.T@Sh@Q)

def build_FX(X, laux):
    P=[p for p in primes_upto(int(math.exp(X))) ]
    blocks=[]; lens=[]
    for p in P:
        r=int(math.floor(X/math.log(p)+1e-12))
        if r<1: continue
        blocks.append(companion_block(r,p**-0.5)); lens+= [math.log(p)]*r
    n=sum(b.shape[0] for b in blocks)
    T=np.zeros((n,n)); i=0
    for b in blocks:
        m=b.shape[0]; T[i:i+m,i:i+m]=b; i+=m
    U,s,Vt=np.linalg.svd(T)
    J=np.nonzero(s<1-1e-12)[0]; d=len(J)
    S=np.zeros((n+d,n+d))
    S[:n,:n]=T
    S[:n,n:]=U[:,J]*np.sqrt(1-s[J]**2)[None,:]
    S[n:,:n]=(np.sqrt(1-s[J]**2)[:,None])*Vt[J,:]
    S[n:,n:]=-np.diag(s[J])
    lens=lens+[laux]*d
    return S, np.array(lens), P, T, n, d

if __name__=='__main__':
    t0=time.time()
    X=math.log(16.5); laux=2.9
    S,lens,P,T,n,d=build_FX(X,laux)
    print('F_X: X=log 16.5, short bonds n=%d (sum of short lengths %.6f, psi(16)=log 720720=%.6f), aux bonds d=%d of length %.2f'%(n,lens[:n].sum(),math.log(720720),d,laux))
    print('   S orthogonal: max|SS^T-I|=%.1e'%np.abs(S@S.T-np.eye(n+d)).max(),' TrL=%.4f'%lens.sum())
    # block power traces vs zeta targets
    i=0
    for p in P:
        r=int(math.floor(X/math.log(p)+1e-12)); B=T[i:i+r,i:i+r]; i+=r
        tr=[np.trace(np.linalg.matrix_power(B,m)) for m in range(1,r+3)]
        print('   p=%2d r=%d  Tr B^m (m=1..r+2):'%(p,r),' '.join('%+.6f'%x for x in tr),'  target -p^(-m/2):',' '.join('%+.6f'%(-p**(-m/2)) for m in range(1,r+3)))
    K=650.0
    eig,m0,TrL=spectrum(S,list(lens),K=K,dk=0.002)
    print('   eigenvalues in (0,%g]: %d (Weyl %.1f), m0=%d, time %.1f'%(K,len(eig),TrL*K/(2*np.pi),m0,time.time()-t0),flush=True)
    np.save('FX_eig.npy',eig)
    # zeta zeros
    import os
    if os.path.exists('zeta_zeros.npy'): zs=np.load('zeta_zeros.npy')
    else:
        import mpmath
        mpmath.mp.dps=20; zs=[]; j=1
        while True:
            g=float(mpmath.zetazero(j).imag)
            if g>K: break
            zs.append(g); j+=1
        zs=np.array(zs); np.save('zeta_zeros.npy',zs)
    print('   zeta zeros up to %g: %d  (time %.1f)'%(K,len(zs),time.time()-t0),flush=True)
    phi,h,a=make_test(X)
    rr=np.arange(-4000,4000,0.004)
    Om=(digamma(0.25+0.5j*rr).real-math.log(math.pi))/(2*math.pi)
    Ns=[N for N in range(2,17)]
    print('   t        graph-Weyl        zeros-pole-arch      diff        | zeta EF check: zeros  vs pole+arch-primes')
    for t in [0.0,7.0,14.134725,17.5,21.022040,23.0,25.010858,40.0,60.0,100.0]:
        ht=lambda k: h(k-t)+h(k+t)
        G=m0*ht(0).real+2*np.sum(ht(eig).real)-TrL/(2*np.pi)*2*phi(0)
        Z=2*np.sum(ht(zs).real)
        pole=2*ht(0.5j).real
        arch=np.sum(ht(rr).real*Om)*0.004
        primes=sum(vonmangoldt(N)/math.sqrt(N)*2*math.cos(t*math.log(N))*phi(math.log(N)) for N in Ns)/math.pi
        print('   %8.4f  %+.10f  %+.10f  %+.2e  | %+.10f %+.10f %+.1e'%(t,G,Z-pole-arch,G-(Z-pole-arch),Z,pole+arch-primes,Z-(pole+arch-primes)),flush=True)
    print('time',round(time.time()-t0,1))
