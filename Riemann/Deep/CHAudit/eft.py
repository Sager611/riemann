
# eft.py: explicit-formula falsification tests of CH Cor. 1.4 (Moebius, sigma=0,1) and Thm 1.2 (psi, sigma=0) at small T,
# using exact M(n), m(n), psi(n) for n<=1e9 (sieve_out.txt) and zeros (Odlyzko zeros6) with complex zeta'(rho).
import numpy as np, gzip, time
from flint import acb, acb_series, ctx
t0=time.time()
S='/Users/adriansager/playground/riemann/.scratch/chaudit/'
g=np.loadtxt(gzip.open(S+'dl/zeros6.gz','rt'),max_rows=140000)
def th(t): return 0.5*t*np.log(t/(2*np.pi)) - 0.5*t - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
def thp(t): return 0.5*np.log(t/(2*np.pi)) - 1/(48*t*t) - 7/(1920*t**4)
import mpmath as mp
mp.mp.dps=40
from numpy.polynomial import chebyshev as C
def Psi(p):
    p=mp.mpf(p); return mp.cos(2*mp.pi*(p*p-p-mp.mpf(1)/16))/mp.cos(2*mp.pi*p)
K=90; nodes=[(mp.cos(mp.pi*(k+mp.mpf(1)/2)/K)+1)/2 for k in range(K)]
cf=C.chebfit(np.array([float(2*x-1) for x in nodes]),np.array([float(Psi(x)) for x in nodes]),60); dcf=C.chebder(cf)*2.0
def Zp_rs(t):
    a=np.sqrt(t/(2*np.pi)); N=np.floor(a).astype(np.int64); p=a-N
    n=np.arange(1,int(N.max())+1,dtype=np.float64); ln=np.log(n); isq=1/np.sqrt(n)
    ph=th(t)[:,None]-t[:,None]*ln[None,:]
    s=-2*np.sum(np.where(n[None,:]<=N[:,None], isq[None,:]*np.sin(ph)*(thp(t)[:,None]-ln[None,:]),0.0),axis=1)
    sg=np.where(N%2==1,1.0,-1.0)
    return s+sg*(a**-0.5*C.chebval(2*p-1,dcf)-0.5*a**-1.5*C.chebval(2*p-1,cf))/(4*np.pi*a)
ctx.prec=80
def zpf(gam): return complex(acb_series([acb(0.5,gam),1],prec=2).zeta()[1])
zc=np.empty(len(g),dtype=complex)
lo=g<=1000
for i in np.where(lo)[0]: zc[i]=zpf(g[i])
hi=np.where(~lo)[0]
for j in range(0,len(hi),2000):
    ii=hi[j:j+2000]; t=g[ii]; zc[ii]=-1j*np.exp(-1j*th(t))*Zp_rs(t)
chk=[700,1500,5000,20000,60000,100000,139000]
print('zeta-prime check (rel err RS-complex vs flint):',' '.join('%.1e'%(abs(zc[i]-zpf(g[i]))/abs(zpf(g[i]))) for i in chk),flush=True)
# exact data
dat=[l.split() for l in open(S+'sieve_out.txt') if l.startswith('S ')]
n=np.array([int(r[1]) for r in dat],dtype=np.float64); Mn=np.array([int(r[2]) for r in dat],dtype=np.float64)
mn=np.array([float(r[3]) for r in dat]); psin=np.array([float(r[4]) for r in dat])
print('samples',len(n),flush=True)
def zeta_min(T):
    vals=[abs(complex(acb(r,T).zeta())) for r in np.linspace(-1/64,1,400)]
    return min(vals)
for T in [1000.0, 10000.0, 100000.0]:
    gg=g[g<=T]; zz=zc[g<=T]
    if len(gg)<len(g[g<=T+1]) or (T-gg[-1])<0.05: print('T near zero!'); 
    zmin=zeta_min(T); d=np.pi/(2*T)
    rho=0.5+1j*gg
    c0=(1/np.tanh(d*rho)-np.tanh(d))/zz          # w_{d,0}(rho)/zeta'(rho)
    c1=(1/np.tanh(d*(rho-1)))/zz                 # w_{d,1}(rho)/zeta'(rho)
    u=rho/(1j*T); om=(1-u)/np.tan(np.pi*u)+(1/np.pi)*(np.pi/T)/np.tanh(np.pi/T)   # omega^+_{T,0}
    th1=1-(rho-1)/(1j*T)
    Cm=np.stack([c0,c1,om,th1],axis=1)
    sel=n>=np.e**2*T
    xs=np.concatenate([n[sel],n[sel]+1-1e-7]); M_=np.concatenate([Mn[sel],Mn[sel]]); m_=np.concatenate([mn[sel],mn[sel]]); p_=np.concatenate([psin[sel],psin[sel]])
    R=np.empty((len(xs),4),dtype=complex)
    for k in range(0,len(xs),16):
        L=np.log(xs[k:k+16]); E=np.exp(1j*np.outer(L,gg)); R[k:k+16]=E@Cm
    sq=np.sqrt(xs)
    P0=2*d*sq*R[:,0].real-2.0; P1=2*d*R[:,1].real/sq
    D0=M_-P0; B0=(np.pi/2)*xs/(T-1)+2
    D1=m_-P1; B1=(np.pi/2)/(T-1)+2/xs
    mainpsi=(np.pi/T)/np.tanh(np.pi/T)-np.log(2*np.pi)/xs-(2*np.pi/T)*R[:,2].imag/sq
    W=np.pi/T-(2*np.pi/T)*R[:,3].real/sq
    Dp=p_/xs-mainpsi
    small=d*xs<np.sqrt(xs)
    print('T=%g zeros=%d  min|zeta(r+iT)|,r in[-1/64,1]=%.4f (cond 1/that=%.2f <= log^2 x=%.0f)  x-range %.3g..%.3g (%d pts)'%(T,len(gg),zmin,1/zmin,np.log(xs.min())**2,xs.min(),xs.max(),len(xs)))
    print('  MOEBIUS s=0: max|M-P0|/bound=%.4f ; max|M-P0|/(delta x)=%.4f ; rms(M-P0)/rms(M) on delta x<sqrt x: %.4f'%(np.max(np.abs(D0)/B0),np.max(np.abs(D0)/(d*xs)),np.sqrt(np.mean(D0[small]**2)/np.mean(M_[small]**2)) if small.any() else float('nan')))
    print('  MOEBIUS s=1: max|m-P1|/bound=%.4f ; max|m-P1|/delta=%.4f'%(np.max(np.abs(D1)/B1),np.max(np.abs(D1)/d)))
    print('  PSI s=0: min width W*T/pi=%.4f ; max |psi/x-main|/W=%.4f ; max|psi/x-main|*T/pi=%.4f ; max|psi/x-1|*T/pi=%.2f'%(np.min(W)*T/np.pi,np.max(np.abs(Dp)/W),np.max(np.abs(Dp))*T/np.pi,np.max(np.abs(p_/xs-1))*T/np.pi),flush=True)
    # control: wrong sign of the zero sums must fail
    print('  CONTROL (zero sums with wrong sign): max|M-P0bad|/bound=%.3f  max|psi/x-mainbad|/W=%.3f'%(np.max(np.abs(M_-(-2*d*sq*R[:,0].real-2.0))/B0),np.max(np.abs(p_/xs-((np.pi/T)/np.tanh(np.pi/T)-np.log(2*np.pi)/xs+(2*np.pi/T)*R[:,2].imag/sq))/W)),flush=True)
print('total s',round(time.time()-t0,1))
