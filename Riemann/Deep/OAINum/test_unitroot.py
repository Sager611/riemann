import random, sys
from ffpoly import *
from unitroot import unitroot, ddf
def prim_root_unity(q,p):
    for g in range(2,p):
        z=pow(g,(p-1)//q,p)
        if z!=1: return z
def primary_gen(q,p):
    n=p-1; s=0; t=n
    while t%q==0: t//=q; s+=1
    for g in range(2,p):
        w=pow(g,t,p)
        if pow(w,q**(s-1),p)!=1: return w
def rand_sqfree(deg_,p):
    while True:
        f=[random.randrange(p) for _ in range(deg_)]+[1]
        if pgcd(f,pderiv(f,p),p)==[1]: return f
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 0)
tot=0; stats={}
for (p,q) in [(7,3),(13,3),(19,3),(31,3),(31,5),(11,5),(41,5),(29,7),(43,7),(10009,3),(10111,5),(97,3)]:
    zeta=prim_root_unity(q,p); omega=primary_gen(q,p)
    for trial in range(60):
        h=random.randrange(1,9)
        H=rand_sqfree(h,p)
        if h*h>=p: pass
        while True:
            b0=[random.randrange(p) for _ in range(h)]
            b0=norm(b0)
            if b0 and pgcd(b0,H,p)==[1]: break
        a=ppowmod(b0,q,H,p)
        b=unitroot(a,H,q,zeta,omega,p,stats)
        assert ppowmod(b,q,H,p)==pmod(a,H,p), (p,q,H,a,b)
        tot+=1
print('unitroot tests passed:',tot,'max c used in unit search:',stats.get('search_c_max'))
