
# Finite-T test of DPC Prop 1 (TeX 262-388) with actual zeta zeros (chunked; self-aborts above 280 MB RSS):
# E_T[V^j] (numerical, Riemann-Siegel Z on a grid) vs  sum_r C(j,r) gamma_{j-r} <(K/pi)^r c,c>/<c,c>,
# gamma_m = Gaussian average of G^m (exact gamma-factor difference). Off-diagonal terms are e^{-12}-small for N=20.
import numpy as np, mpmath as mp, math, sys, time, resource
from scipy.special import loggamma
t0=time.time()
def rss_mb(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20   # macOS: bytes
def guard():
    if rss_mb()>280: sys.exit("ABORT: RSS %.0f MB > 280 MB"%rss_mb())
T=float(sys.argv[1]) if len(sys.argv)>1 else 2000.0
phi=float(sys.argv[2]) if len(sys.argv)>2 else 0.6
NPOL=int(sys.argv[3]) if len(sys.argv)>3 else 20
L=math.log(T); H=T/L; h=2*math.pi*phi/L
STEP=float(sys.argv[4]) if len(sys.argv)>4 else 0.002
m=int(round(h/2/STEP)); dl=h/2/m
mp.mp.dps=40
Psi=lambda p: mp.cos(2*mp.pi*(p*p-p-mp.mpf(1)/16))/mp.cos(2*mp.pi*p)
P=np.cos(np.pi*(np.arange(161)+0.5)/161)*0.5+0.5
C0=[];C1=[];C2=[]
for p in P:
    pp=mp.mpf(p)
    if abs(pp-0.25)<1e-9 or abs(pp-0.75)<1e-9: pp+=mp.mpf('1e-12')
    d=[mp.diff(Psi,pp,k) for k in (0,2,3,6)]
    C0.append(float(d[0])); C1.append(float(-d[2]/(96*mp.pi**2))); C2.append(float(d[1]/(64*mp.pi**2)+d[3]/(18432*mp.pi**4)))
cf=[np.polynomial.chebyshev.chebfit(2*P-1,np.array(C),60) for C in (C0,C1,C2)]
def theta(t): return t/2*np.log(t/(2*np.pi))-t/2-np.pi/8+1/(48*t)+7/(5760*t**3)
def Z(t):
    a=np.sqrt(t/(2*np.pi)); Nn=np.floor(a).astype(int); p=a-Nn; th=theta(t)
    s=np.zeros_like(t)
    for n in range(1,int(Nn.max())+1):
        s+=np.where(n<=Nn, np.cos(th-t*np.log(n))/np.sqrt(n), 0.0)
    x=2*p-1; u=np.sqrt(2*np.pi/t)
    R=((-1.0)**(Nn-1))*np.sqrt(u)*(np.polynomial.chebyshev.chebval(x,cf[0])+np.polynomial.chebyshev.chebval(x,cf[1])*u+np.polynomial.chebyshev.chebval(x,cf[2])*u*u)
    return 2*s+R
mp.mp.dps=20
tv=np.array([1000.3,1777.7,2500.1,3333.3,4100.9,4999.5,2222.2222,3000.0])
err=max(abs(Z(np.array([x]))[0]-float(mp.siegelz(x))) for x in tv)
print("RS vs mpmath.siegelz max abs err: %.2e"%err)
a,b=T/2,5*T/2
tg=np.arange(a-h/2, b+h/2+dl/2, dl); n_g=len(tg)
CH=200000
zg=np.empty(n_g); logF=np.empty(n_g); imF=np.empty(n_g)
for i in range(0,n_g,CH):
    tt=tg[i:i+CH]; zg[i:i+CH]=Z(tt); s=0.5+1j*tt; lg=-(s/2)*np.log(np.pi)+loggamma(s/2)
    logF[i:i+CH]=np.log((0.25+tt**2)/2)+lg.real; imF[i:i+CH]=lg.imag; guard()
sg=np.signbit(zg); ch=np.concatenate([[0],np.cumsum(sg[1:]!=sg[:-1])]).astype(np.int32); del sg
logxi=np.log(np.abs(zg))+logF; del zg
print("zeros found: %d ; RvM main term %.1f ; RSS %.0f MB"%(ch[-1], (theta(tg[-1])-theta(tg[0]))/np.pi, rss_mb()))
def mob(n):
    r=1;k=2
    while k*k<=n:
        if n%k==0:
            n//=k
            if n%k==0: return 0
            r=-r
        k+=1
    return -r if n>1 else r
def vmang(k):
    for p in range(2,k+1):
        if k%p==0:
            q=k
            while q%p==0: q//=p
            return math.log(p) if q==1 else 0.0
    return 0.0
idx=list(range(1,NPOL+1))
KM=np.zeros((NPOL,NPOL),complex)
for q in idx:
    for k in range(2,NPOL//q+1):
        lam=vmang(k)
        if lam>0: KM[k*q-1,q-1]=-2j*lam/(math.sqrt(k)*math.log(k))*math.sin(math.pi*phi*math.log(k)/L)
tests={'delta1':np.array([1.0]+[0.0]*(NPOL-1)),
       'mobius f=1':np.array([mob(n)/math.sqrt(n) for n in idx]),
       'mobius f=1-x':np.array([mob(n)/math.sqrt(n)*(1-math.log(n)/math.log(NPOL)) for n in idx])}
J=5
acc={k:np.zeros(J,complex) for k in tests}; accw={k:0.0 for k in tests}; accN={k:0.0 for k in tests}; acc2={k:0.0 for k in tests}
gsum=0.0; gG=np.zeros(J,complex); maxN=0
for i in range(m, n_g-m, CH):
    c_idx=np.arange(i, min(i+CH, n_g-m)); tc=tg[c_idx]
    X=(logxi[c_idx+m]-logxi[c_idx-m])/np.pi; N=(ch[c_idx+m]-ch[c_idx-m]).astype(float); maxN=max(maxN,N.max())
    V=X+1j*N
    G=((logF[c_idx+m]-logF[c_idx-m])+1j*(imF[c_idx+m]-imF[c_idx-m]))/np.pi
    g=np.exp(-(tc-1.5*T)**2/(2*H*H))
    gsum+=g.sum(); gG+=np.array([np.sum(g*G**mm) for mm in range(J)])
    for name,c in tests.items():
        Rt=np.zeros(len(tc),complex)
        for n in idx:
            if c[n-1]!=0: Rt+=c[n-1]*np.exp(-1j*tc*math.log(n))
        w=g*np.abs(Rt)**2
        accw[name]+=w.sum(); accN[name]+=np.sum(w*(N>=2)); acc2[name]+=np.sum(w*N)
        Vp=np.ones_like(V)
        for j in range(1,J):
            Vp=Vp*V; acc[name][j]+=np.sum(w*Vp)
    guard()
gam=gG/gsum
print("max N in a window: %d ; gamma_1 = %s (i phi = %sj) ; RSS %.0f MB"%(maxN, np.round(gam[1],5), phi, rss_mb()))
for name,c in tests.items():
    print("== %s   weighted E[N]=%.4f  P(N>=2)=%.4f"%(name, acc2[name]/accw[name], accN[name]/accw[name]))
    for j in range(1,J):
        num=acc[name][j]/accw[name]; pred=0
        for r in range(j+1):
            pred+=math.comb(j,r)*gam[j-r]*(np.conj(c)@np.linalg.matrix_power(KM/np.pi,r)@c)/(np.conj(c)@c)
        print("  j=%d numeric %s  predicted %s  |diff| %.1e"%(j,np.round(num,4),np.round(pred,4),abs(num-pred)))
print("t=%.1f s, max RSS %.0f MB"%(time.time()-t0, rss_mb()))
