import mpmath as mp, sys, time
D=int(sys.argv[2]); mp.mp.dps=D; pi=mp.pi
kind=sys.argv[1]; par=mp.mpf(sys.argv[3]); EMAX=int(sys.argv[4])
def theta(u):
    u=abs(u); z=mp.e**(2*u)
    return mp.fsum((4*pi*pi*n**4*mp.e**(mp.mpf(9)*u/2)-6*pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-pi*n*n*z) for n in range(1,14))
if kind=='riemann':   # Phi(v) e^{par v^2}
    K=lambda v: theta(v)*mp.e**(par*v*v); L=mp.mpf(3)
else:                 # gaussian-cos: (e^{-(v-b)^2}+e^{-(v+b)^2})/2, FT ~ cos(b t) e^{-t^2/4}; par=b; b=0: pure Gaussian
    K=lambda v: (mp.e**(-(v-par)**2)+mp.e**(-(v+par)**2))/2; L=par+mp.sqrt(2.3*D)+2
h=mp.mpf('0.04'); n=int(L/h)
nodes=[h*i for i in range(-n,n+1)]; kv=[K(v)*mp.e**(v/2) for v in nodes]; ev=[mp.e**v for v in nodes]
def P(x): return h*mp.fsum(k/(x+(1-x)*e) for k,e in zip(kv,ev))
def pmul(a,b):
    r=[mp.mpc(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
t0=time.time(); cache={}
for E in range(2,EMAX+1,2):
    # A_j(t) = (1/2+it)_j , B_k(t) = (1/2-it)_k as polynomials in t
    A=[[mp.mpc(1)]]; B=[[mp.mpc(1)]]
    for m in range(E):
        A.append(pmul(A[-1],[mp.mpf(1)/2+m, mp.mpc(0,1)])); B.append(pmul(B[-1],[mp.mpf(1)/2+m, mp.mpc(0,-1)]))
    Z=[mp.mpf(0)]*(E+1)
    for j in range(E+1):
        x=mp.mpf(j)/E; key=min(j,E-j)/mp.mpf(E)
        if key not in cache: cache[key]=P(key)
        w=mp.binomial(E,j)*cache[key]
        pr=pmul(A[j],B[E-j])
        for k,c in enumerate(pr): Z[k]+=w*mp.re(c)
    U=[Z[2*m]*(-1)**m for m in range(E//2+1)]   # t^2 = u  ->  poly in u (t^{2m} = u^m); sign from i^{2m} already in real parts
    U=[Z[2*m] for m in range(E//2+1)]
    rs=mp.polyroots(list(reversed(U)),maxsteps=500,extraprec=2*D,error=False)
    bad=[r for r in rs if abs(mp.im(r))>mp.mpf(10)**(-D//3) or mp.re(r)<0]
    tmin=min(mp.sqrt(mp.re(r)) for r in rs if mp.re(r)>0) if any(mp.re(r)>0 for r in rs) else None
    print(kind,mp.nstr(par,4),'E',E,'nonreal_t',len(bad),'tmin',mp.nstr(tmin,8) if tmin else None,'worst',mp.nstr(bad[0],6) if bad else '-','t',round(time.time()-t0),flush=True)

