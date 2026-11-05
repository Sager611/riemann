import mpmath as mp, sys, time
mp.mp.dps=50; pi=mp.pi
def raw(u):
    z=mp.e**(2*u)
    return mp.fsum((4*pi*pi*n**4*mp.e**(mp.mpf(9)*u/2)-6*pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-pi*n*n*z) for n in range(1,8))
def phi(u): return raw(abs(u))
I=mp.quad(phi,[-5,-2,-1,0,1,2,5]); print('intPhi',mp.nstr(I,15),'xi(1/2)=0.4971207781',flush=True)
cache={}
def P(t):
    t=mp.mpf(t); t=min(t,1-t)
    k=mp.nstr(t,40)
    if k not in cache: cache[k]=mp.quad(lambda v: phi(v)*mp.e**(v/2)/(t+(1-t)*mp.e**v),[-5,-2,-1,0,1,2,5])
    return cache[k]
t0=time.time()
for E in range(2,int(sys.argv[1])+1,2):
    c=[(-1)**j*mp.binomial(E,j)*P(mp.mpf(j)/E) for j in range(E+1)]
    rs=mp.polyroots(list(reversed(c)),maxsteps=400,extraprec=400,error=False)
    r=max(rs,key=lambda z: abs(abs(z)-1)); d=abs(abs(r)-1)
    print('E',E,'maxdev',mp.nstr(d,8),'worst root',mp.nstr(r,12),'t',round(time.time()-t0),flush=True)

