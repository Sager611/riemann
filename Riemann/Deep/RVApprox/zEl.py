import mpmath as mp, time, sys
LAM=mp.mpf(sys.argv[1]); EMAX=int(sys.argv[2])
mp.mp.dps=50; pi=mp.pi
def raw(u):
    z=mp.e**(2*u)
    return mp.fsum((4*pi*pi*n**4*mp.e**(mp.mpf(9)*u/2)-6*pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-pi*n*n*z) for n in range(1,8))
phi=lambda u: raw(abs(u))*mp.e**(LAM*u*u)
cache={}
def P(t):
    t=min(t,1-t); k=mp.nstr(t,40)
    if k not in cache: cache[k]=mp.quad(lambda v: phi(v)*mp.e**(v/2)/(t+(1-t)*mp.e**v),[-5,-2,-1,0,1,2,5])
    return cache[k]
def pmul(a,b):
    r=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
t0=time.time()
for E in range(2,EMAX+1,2):
    Z=[mp.mpf(0)]*(E+1)   # coefficients low->high in s
    for j in range(E+1):
        u=(-1)**j*mp.binomial(E,j)*P(mp.mpf(j)/E)
        poly=[mp.mpf(1)]
        for m in range(E):   # factor (E-j-m - s)
            poly=pmul(poly,[mp.mpf(E-j-m),mp.mpf(-1)])
        f=u/mp.factorial(E)
        for i,c in enumerate(poly): Z[i]+=f*c
    while abs(Z[-1])<mp.mpf(10)**(-40): Z.pop()
    rs=mp.polyroots(list(reversed(Z)),maxsteps=400,extraprec=400,error=False)
    w=max(rs,key=lambda z: abs(mp.re(z)-0.5))
    print('E',E,'deg',len(Z)-1,'max|Re s-1/2|',mp.nstr(abs(mp.re(w)-0.5),6),'root',mp.nstr(w,10),'t',round(time.time()-t0),flush=True)

