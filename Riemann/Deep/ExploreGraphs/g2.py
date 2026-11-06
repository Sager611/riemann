# explore12 g2: least-squares fit of a one-petal-per-prime flower (S in O(12)) to zeta / DH weights on N<=16
import numpy as np, math, time, sys, json
from scipy.linalg import expm
from scipy.optimize import least_squares
from g1 import *

def antisym(x, D):
    A=np.zeros((D,D)); iu=np.triu_indices(D,1); A[iu]=x; return A-A.T

def make_S(x, D, comp):
    R0=np.eye(D)
    if comp: R0[0,0]=-1.0
    return R0@expm(antisym(x,D))

def fit(P, targets, starts=60, seed=0, tol=1e-12, maxtries=None):
    E=len(P); D=2*E; npar=D*(D-1)//2
    Ns=sorted(targets); alph={N:factor_vec(N,P) for N in Ns}
    def res(x,comp):
        S=make_S(x,D,comp)
        return np.array([coeff(S,E,alph[N]).real-targets[N] for N in Ns])
    rng=np.random.default_rng(seed); best=(np.inf,None,None)
    for s in range(starts):
        comp=s%2; x0=rng.normal(0,1.5,npar)
        r=least_squares(res,x0,args=(comp,),method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
        c=np.sqrt(np.sum(r.fun**2))
        if c<best[0]: best=(c,r.x,comp)
        if c<tol: break
    return best

def zeta_targets(P, Nmax):
    T={}
    for N in smooth_list(P,Nmax):
        a=factor_vec(N,P); n=sum(a)
        T[N]= -vonmangoldt(N)/math.sqrt(N)*n/math.log(N)   # c-target (= -p^{-m/2} at prime powers, 0 at composites)
    return T

def dh_coeffs(Nmax):
    kap=(math.sqrt(10-2*math.sqrt(5))-2)/(math.sqrt(5)-1)
    chi={0:0,1:1,2:1j,3:-1j,4:-1}
    a=[0.0]+[((1-1j*kap)*chi[n%5]).real for n in range(1,Nmax+1)]
    Lam=[0.0]*(Nmax+1)
    for n in range(2,Nmax+1):
        s=a[n]*math.log(n)
        for d in range(2,n):
            if n%d==0: s-=Lam[d]*a[n//d]
        Lam[n]=s
    return a,Lam,kap

if __name__=='__main__':
    t0=time.time(); mode=sys.argv[1] if len(sys.argv)>1 else 'B'
    if mode=='C':
        P=[2,3,5,7,11,13]; E=6; D=12
        T=zeta_targets(P,16)
        c,x,comp=fit(P,T,starts=40,seed=5,tol=1e-13)
        print('C: zeta window N<=16, petals',P,': residual %.2e (comp %d)'%(c,comp))
        S=make_S(x,D,comp); np.save('C_S.npy',S)
        W=weights(S,P,40)
        for N in sorted(W):
            print('   N=%2d W_G=%+.10f W_zeta=%+.10f'%(N,W[N][0],-vonmangoldt(N)/math.sqrt(N)))
        a,Lam,kap=dh_coeffs(16)
        TD={}
        for N in smooth_list(P,16):
            al=factor_vec(N,P); TD[N]=-Lam[N]/math.sqrt(N)*sum(al)/math.log(N)
        cD,xD,compD=fit(P,TD,starts=40,seed=7,tol=1e-13)
        print('D: Davenport-Heilbronn (kappa=%.6f) window N<=16: residual %.2e'%(kap,cD))
        SD=make_S(xD,D,compD); np.save('D_S.npy',SD)
        WD=weights(SD,P,16)
        for N in sorted(WD):
            print('   N=%2d W_G=%+.10f W_DH=%+.10f'%(N,WD[N][0],-Lam[N]/math.sqrt(N)))
        print('time',round(time.time()-t0,1))
