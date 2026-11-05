# gammas.py: Jensen moments gamma_n of Chung (Zenodo 23977400) normalization
# gamma_n = (1/(2n-1)!!) int_0^inf u^{2n} Phi(u) du, Phi(u)=2e^{u/2} sum a_n(2a_n-3)e^{-a_n}, a_n=pi n^2 e^{2u}
# Run: uv run --no-project --with mpmath python3 gammas.py NMAX DPS OUT
import sys, json, mpmath as mp
NMAX=int(sys.argv[1]); DPS=int(sys.argv[2]); OUT=sys.argv[3]
mp.mp.dps=DPS
pi=mp.pi
def Phi(u):
    s=mp.mpf(0); e2u=mp.e**(2*u)
    for n in range(1,40):
        a=pi*n*n*e2u
        t=a*(2*a-3)*mp.e**(-a)
        s+=t
        if t < mp.mpf(10)**(-DPS-20): break
    return 2*mp.e**(u/2)*s
pts=[0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1,1.125,1.25,1.5,1.75,2,2.5,3,4]
g=[]
for n in range(NMAX+1):
    I=mp.quad(lambda u: u**(2*n)*Phi(u), pts)
    g.append(I/mp.fac2(2*n-1) if n>0 else I)
json.dump([mp.nstr(x, DPS-5) for x in g], open(OUT,'w'))
print('gamma0..3', [mp.nstr(x,20) for x in g[:4]])
print('mu2..7', [mp.nstr(g[n-1]/g[n],12) for n in range(2,8)])
paper=['23.25999270','24.71940992','26.06175622','27.31403609','28.49434404','29.61545695']
print('Lemma 1.1 moment ratios mu_2..mu_7 match paper to 1e-8:', 'PASS' if all(abs(g[n-1]/g[n]-mp.mpf(p))<mp.mpf('1e-8') for n,p in zip(range(2,8),paper)) else 'FAIL')
