# explore12 g5: frozen Riemann-Siegel cosine polynomial g_N(tau)=sum_{n<=N} n^-1/2 cos(tau log(N/n)):
# real-rootedness (graph RS-lookalikes are real-rooted by unitarity) vs the flat-density positivity barrier
import numpy as np, math
from g1 import vonmangoldt

def count_real_zeros(N, T=4000.0, dt=0.004):
    n=np.arange(1,N+1); w=n**-0.5; fr=np.log(N/n)
    tot=0; prev=None
    for t0 in np.arange(0,T,200.0):
        tau=np.arange(t0,min(T,t0+200.0),dt)
        g=(w[None,:]*np.cos(tau[:,None]*fr[None,:])).sum(1)
        s=np.sign(g)
        if prev is not None: tot+= int(prev!=s[0])
        tot+=int(np.sum(s[1:]!=s[:-1])); prev=s[-1]
    return tot, T*math.log(N)/math.pi

def barrier(N, m=600):
    """lambda_max of the prime form sum_{n<N} c_n (u(.-log n)+u(.+log n)) on L^2[-X/2,X/2], X=log N (Fejer/autocorrelation test functions); 
       real-rootedness of g_N requires lambda_max <= log N"""
    X=math.log(N); y=np.linspace(-X/2,X/2,m); dy=y[1]-y[0]
    Kmat=np.zeros((m,m))
    for k in range(2,N):
        c=vonmangoldt(k)/math.sqrt(k)
        if c==0: continue
        x=math.log(k)
        # (u * u~)(x) weight: sum_i u_i u(y_i - x) -> linear interpolation
        for sgn in (1,-1):
            pos=(y-sgn*x - y[0])/dy
            i0=np.floor(pos).astype(int); fr=pos-i0
            ok=(i0>=0)&(i0<m-1)
            for i in np.nonzero(ok)[0]:
                Kmat[i,i0[i]]+=c*(1-fr[i]); Kmat[i,i0[i]+1]+=c*fr[i]
    Kmat=0.5*(Kmat+Kmat.T)
    lam=np.linalg.eigvalsh(Kmat)[-1]
    fej=sum(vonmangoldt(k)/math.sqrt(k)*(1-math.log(k)/X) for k in range(2,N))
    return lam, fej, X

if __name__=='__main__':
    import sys
    print(' N   real zeros in [0,4000]   all-real count T logN/pi   deficit%   | real-rootedness needs log N >= B(N): B_opt=lam_max/2  B_Fejer   log N')
    for N in [2,3,4,5,8,10,20,30,40,50,55,60,65,70,100,200]:
        z,full=count_real_zeros(N)
        lam,fej,X=barrier(N)
        print('%4d   %8d              %10.1f        %5.1f%%   |  %8.4f   %8.4f   %7.4f  %s'%(N,z,full,100*(full-z)/full,lam/2,fej,X,'VIOLATED' if lam/2>X else ''),flush=True)
    for X in [math.log(16.5), math.log(100.5), math.log(1000.5)]:
        N=int(math.exp(X)); lam,fej,_=barrier(N+1)
        print('graph general-length bound at X=%.4f: TrL >= lam_max(K_X) = %.4f (Fejer %.4f); log-integer lengths: psi(e^X) = %.4f'%(X,lam,2*fej,sum(vonmangoldt(k) for k in range(2,N+1))))
