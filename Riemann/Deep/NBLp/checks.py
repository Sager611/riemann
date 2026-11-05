# /// script
# dependencies = ["mpmath", "numpy"]
# ///
# nblp checks: (1) N=10 exact infinite sum via periodicity (digamma), (2) direct-definition spot check of v(m),
# (3) constants C_Gamma, K_p, A_p, (4) illustration ||chi + H_eps||_p, H_eps(x) = -E_eps([1/x]).
import mpmath as mp, numpy as np, math, time, sys
mp.mp.dps = 30
t0=time.time()
def mobius_upto(n):
    mu=[1]*(n+1); isp=[True]*(n+1)
    for i in range(2,n+1):
        if isp[i]:
            for j in range(i,n+1,i):
                if j>i: isp[j]=False
                mu[j]*=-1
            for j in range(i*i,n+1,i*i): mu[j]=0
    mu[0]=0; return mu
mu=mobius_upto(1000)
from fractions import Fraction as Fr
# (1)+(2) N=10 exact
N=10; L=2520; p=mp.mpf('1.1')
a={k:Fr(mu[k]*(N-k),N) for k in range(1,N+1)}
def v_exact(m): return 1+sum(a[k]*Fr(m%k,k) for k in range(2,N+1))
S=mp.mpf(0)
for r in range(L):
    vr=v_exact(r)
    c=(mp.digamma(mp.mpf(1)+mp.mpf(1)/L)-mp.digamma(1))/L if r==0 else (mp.digamma(mp.mpf(r+1)/L)-mp.digamma(mp.mpf(r)/L))/L
    S+=abs(mp.mpf(vr.numerator)/vr.denominator)**p*c
print("N=10 exact p=1.1 norm:", mp.nstr(S**(1/p),15))
S2=mp.mpf(0)
for r in range(L):
    vr=v_exact(r)
    c=(mp.digamma(mp.mpf(1)+mp.mpf(1)/L)-mp.digamma(1))/L if r==0 else (mp.digamma(mp.mpf(r+1)/L)-mp.digamma(mp.mpf(r)/L))/L
    S2+=(mp.mpf(vr.numerator)/vr.denominator)**2*c
print("N=10 exact p=2 norm:", mp.nstr(mp.sqrt(S2),15))
# direct-definition spot check: v(m) via f_k(x) = {1/(kx)} - {1/x}/k at a point x in I_m, N=100
N2=100; a2={k:Fr(mu[k]*(N2-k),N2) for k in range(1,N2+1)}
import random; random.seed(1); bad=0
for m in [1,2,3,7,50,99,100,101,997,12345,10**6+3]:
    x=Fr(1,m)-Fr(1,10*m*(m+1))   # interior point of (1/(m+1),1/m]
    def frac(q): return q-(q.numerator//q.denominator)
    val=1+sum(a2[k]*(frac(1/(k*x))-frac(1/x)/k) for k in range(2,N2+1))
    val2=1+sum(a2[k]*Fr(m%k,k) for k in range(2,N2+1))
    bad+= (val!=val2)
print("direct-definition spot check N=100 mismatches:", bad)
# (3) constants
f=lambda t: (8+mp.log(2+abs(t)))*mp.sqrt(2+abs(t))/(mp.mpf(1)/4+t*t)
I=2*mp.quad(f,[0,1,10,100,1e4,1e6,mp.inf])
CG=mp.e**4*I
print("C_Gamma =", mp.nstr(CG,8))
for pp in [mp.mpf(8)/7, mp.mpf('1.1'), mp.mpf('1.5')]:
    Kp=CG/(2*mp.pi)+4*pp/(2-pp); r=2*pp/(2-pp); Ap=mp.sqrt(Kp)*(2/r)**(1/r)
    print("p=%s K_p=%s r=%s A_p=%s exponent 1-1/p=%s"%(mp.nstr(pp,6),mp.nstr(Kp,6),mp.nstr(r,6),mp.nstr(Ap,6),mp.nstr(1-1/pp,6)))
print("time", round(time.time()-t0,1)); sys.stdout.flush()
# (4) illustration: ||chi + H_eps||_{L^p(0,1)}^p = sum_m |1 - E_eps(m)|^p/(m(m+1)), E_eps(m)=sum_{j<=m} phi_eps(j) - m/zeta(1+eps)
M=10**7
spf_primes=None
sieve=np.ones(M+1,dtype=bool); sieve[:2]=False
for i in range(2,int(M**0.5)+1):
    if sieve[i]: sieve[i*i::i]=False
primes=np.nonzero(sieve)[0]
m=np.arange(1,M+1,dtype=np.float64); w=1.0/(m*(m+1.0))
for eps in [0.4,0.2,0.1,0.05]:
    phi=np.ones(M+1)
    fac=1.0-primes.astype(np.float64)**(-eps)
    for pr,fc in zip(primes,fac): phi[pr::pr]*=fc
    phi[0]=0.0
    E=np.cumsum(phi)[1:]-m/float(mp.zeta(1+eps))
    g=1.0-E   # value of chi + H_eps on I_m
    out=[]
    for pp in [1.1, 8/7, 2.0]:
        Sp=np.sum(np.abs(g)**pp*w); out.append("p=%.4f norm=%.6f"%(pp,Sp**(1/pp)))
    tailest=np.mean(np.abs(g[M//2:])**1.1)/M
    print("eps=%.3f: %s | tail(p=1.1) est %.2e | max|g| on [M/2,M]=%.3f"%(eps,"; ".join(out),tailest,np.max(np.abs(g[M//2:]))))
    sys.stdout.flush()
print("time", round(time.time()-t0,1))
