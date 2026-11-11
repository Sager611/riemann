
# platt_extrap.py: plausibility of Platt's values (lem:platt109) above 1e6 from exact data to 1.13e6 + HKO-type density.
import numpy as np, mpmath as mp
d=np.load('/Users/adriansager/playground/riemann/.scratch/chaudit/zp.npy'); g=d[0]; zp=d[1]
J=np.cumsum(1/zp)
def Lg(t): return np.log(t/(2*np.pi))
for t in [1e4,1e5,3e5,1e6,1.13e6]:
    i=np.searchsorted(g,t)-1; print('J(%g)=%.2f  J/((t/2pi)L^(1/4))=%.5f'%(t,J[i],J[i]/((t/(2*np.pi))*Lg(t)**0.25)))
# HKO constant G(3/2)^2 * a(1/2)
mp.mp.dps=20
G32=mp.barnesg(1.5); ap=mp.mpf(1)
import sympy
for p in sympy.primerange(2,200000):
    ap*=(1-mp.mpf(1)/p)**0.25*mp.hyp2f1(0.5,0.5,1,mp.mpf(1)/p)
print('HKO constant G(3/2)^2 a(1/2) ~ %s'%mp.nstr(G32**2*ap,8))
# fit J(t) ~ (t/2pi)(K1 L^{1/4} + K2 L^{-3/4}) on [t1,t2]
def fit(t1,t2):
    m=(g>=t1)&(g<=t2); t=g[m]; A=np.stack([(t/(2*np.pi))*Lg(t)**0.25,(t/(2*np.pi))*Lg(t)**-0.75],1)
    k,*_=np.linalg.lstsq(A,J[m],rcond=None); return k
def jdens(t,k): L=Lg(t); return (1/(2*np.pi))*(k[0]*(L**0.25+0.25*L**-0.75)+k[1]*(L**-0.75-0.75*L**-1.75))
def predict(T,t0,k):
    dl=np.pi/(2*T); m=g<=t0; z=dl*(0.5+1j*g[m])
    exact=dl*np.sum(np.abs(1/np.tanh(z))/zp[m])
    # integral t0..T of dl*|coth(dl(1/2+it))| j(t) dt, log-spaced trapezoid
    u=np.linspace(np.log(t0),np.log(T),200001); t=np.exp(u)
    f=dl*np.abs(1/np.tanh(dl*(0.5+1j*t)))*jdens(t,k)*t
    return exact+np.trapezoid(f,u)
paper={1e6+1:2.66161277991001,1893193.5:2.85417779533422,1e7+1:3.36904620179490,1e8+1:4.10963816503581,1e9:4.87936778767100,1e10+1:5.675256}
for (a,b) in [(1e4,2e5),(2e5,1.13e6)]:
    k=fit(a,b); t0=b
    print('fit on [%g,%g]: K1=%.5f K2=%.5f'%(a,b,k[0],k[1]))
    for T,v in paper.items():
        if T<=t0: continue
        p=predict(T,t0,k); print('   T=%-12g predicted %.5f  paper %.5f  diff %+.5f (%.3f%%)'%(T,p,v,p-v,100*(p-v)/v))
    k0=np.array([fit(a,b)[0]+fit(a,b)[1]/Lg(b),0.0])   # one-term model calibrated at t0
    for T,v in paper.items():
        if T<=t0: continue
        p=predict(T,t0,k0); print('   1-term: T=%-12g predicted %.5f  paper %.5f  diff %+.5f'%(T,p,v,p-v))
