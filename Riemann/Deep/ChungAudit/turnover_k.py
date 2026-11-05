# turnover_k.py: Lemma 4.0(2) turnover-lobe claim r*(z_m)=m (m=d-1,d-2) for shifts k>=1 on low anti-diagonals
# kappa_k(u)=Phi_{k-1}(u)/Phi_k(u), Phi_k(u)=(1/(k-1)!) int_u^inf ((t^2-u^2)/2)^{k-1} t Phi(t) dt  (k>=1)
import sys, json, mpmath as mp
from turnover_k0 import Phi_and_d, kappa0, hermite_pos_zeros
mp.mp.dps=30
def Phik(k,u):
    if k==0: return Phi_and_d(u)[0]
    f=lambda t: ((t*t-u*u)/2)**(k-1)*t*Phi_and_d(t)[0]
    return mp.quad(f,[u,u+mp.mpf('0.25'),u+mp.mpf('0.5'),u+1,u+2,u+4])/mp.factorial(k-1)
def kappak(k,u):
    if k==0: return kappa0(u)
    return Phik(k-1,u)/Phik(k,u)
def ystar_k(k,z):
    f=lambda y: 1/y**2+mp.mpf(1)/2-(2/z**2)*kappak(k,mp.sqrt(2)*y/z)
    lo=mp.mpf('0.05'); hi=mp.mpf(1)
    while f(hi)>0: hi*=2
    return mp.findroot(f,(lo,hi),solver='anderson')
g=[mp.mpf(s) for s in json.load(open('data/gammas_dps70.json'))]
for d,k in [(8,1),(7,2),(6,3),(5,4),(4,5),(9,1),(8,2),(10,2),(6,6)]:
    N=d-1; c=[mp.binomial(N,j)*g[k+1+j] for j in range(N+1)]
    xs=sorted([mp.re(x) for x in mp.polyroots(c[::-1],maxsteps=400,extraprec=200)],reverse=True)
    z=[mp.sqrt(-2*x) for x in xs]; yz=hermite_pos_zeros(2*d-1)
    rs=[]
    for m in [d-2,d-1]:
        ys=ystar_k(k,z[m-1]); rs.append((m,1+sum(1 for y in yz if y<ys),mp.nstr(ys,7)))
    ok=all(r==m for m,r,_ in rs)
    print('(d,k)=(%d,%d): (m, r*(z_m), y*) ='%(d,k), rs, ' claim r*=m:', 'PASS' if ok else 'FAIL')
