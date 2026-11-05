# audit_checks.py: one-shot PASS/FAIL audit of the quantitative claims of Zenodo 23977400 (Chung, Rev.3.0) on actual data.
# Run (from this directory): uv run --no-project --with mpmath --with numpy python3 audit_checks.py
# PASS = the paper's claim holds on the tested instance; FAIL = the claim is false there (witness printed).
import json, math, cmath, mpmath as mp, numpy as np
from turnover_k0 import Phi_and_d, kappa0, ystar, hermite_pos_zeros
from lobes_k0 import lobes
from saddle_scan import solve, consts
import omega_np as onp
mp.mp.dps=60
g=[mp.mpf(s) for s in json.load(open('data/gammas_dps70.json'))]
mu=lambda n: g[n-1]/g[n]; de=lambda n: 1-mu(n)/mu(n+1)
_cache={}
def zroots(d,k):
    """z_m = sqrt(-2 x_m), x_1>x_2>... roots of J^{d-1,k+1}(X)=sum C(d-1,j) gamma_{k+1+j} X^j"""
    if (d,k) in _cache: return _cache[(d,k)]
    N=d-1; c=[mp.binomial(N,j)*g[k+1+j] for j in range(N+1)]
    r=mp.polyroots(c[::-1],maxsteps=800,extraprec=300)
    assert max(abs(mp.im(x)) for x in r) < mp.mpf(10)**-30
    xs=sorted([mp.re(x) for x in r],reverse=True)
    _cache[(d,k)]=[mp.sqrt(-2*x) for x in xs]; return _cache[(d,k)]
def omega(d,k,m): return zroots(d,k)[m-1]*mp.sqrt(4*d+2*k-1)/2
def verdict(ok): return 'PASS' if ok else 'FAIL'
# C1 Thm 5.1/5.2 (pp.24-26): omega_{1,k}(d) >= Omega_LS(d,k) >= 12.2757
pairs=[(9,0),(8,1),(7,2),(6,3),(5,4),(4,5),(15,0),(30,0),(10,10)]
ok=True
for d,k in pairs:
    OL=mp.sqrt((2*d+k-mp.mpf(1)/2)*mu(k+2)/(1+(d-2)*mp.sqrt(de(k+2))))
    ok = ok and omega(d,k,1)>=OL>=mp.mpf('12.2757')
print('C1 Thm 5.1-5.2 frequency floor omega_1>=Omega_LS>=12.2757 on',pairs,':',verdict(ok))
# C2 Lemma 4.0 Part 2 (p.19): y_{d-2}<y*(z_{d-1})<y_{d-1} and y_{d-3}<y*(z_{d-2})<y_{d-2}
mp.mp.dps=40
d,k=9,0; z=zroots(d,k); yz=hermite_pos_zeros(2*d-1)
ys8=ystar(z[7]); ys7=ystar(z[6])
ok8=yz[6]<ys8<yz[7]; ok7=yz[5]<ys7<yz[6]
print('C2 Lemma 4.0(2) r*(z_{d-1})=d-1 at (9,0): z_8=%s y*=%s, need (y_7,y_8)=(%s,%s): %s'%(mp.nstr(z[7],12),mp.nstr(ys8,8),mp.nstr(yz[6],8),mp.nstr(yz[7],8),verdict(ok8)))
print('C2 Lemma 4.0(2) r*(z_{d-2})=d-2 at (9,0): z_7=%s y*=%s, need (y_6,y_7)=(%s,%s): %s'%(mp.nstr(z[6],12),mp.nstr(ys7,8),mp.nstr(yz[5],8),mp.nstr(yz[6],8),verdict(ok7)))
# C3 Lemma 4.0 Part 1 (p.18): L_j strictly increasing up to r* and strictly decreasing for j>r*
mp.mp.dps=30
L,_=lobes(d,z[7]); rstar=1+sum(1 for y in yz if y<ys8)
ok=all(L[j]<L[j+1] for j in range(rstar-1)) and all(L[j]>L[j+1] for j in range(rstar-1,d-1))
print('C3 Lemma 4.0(1) peak of L_j at r* at (9,0), z_8: r*=%d, L_7=%s, L_8=%s, L_9=%s: %s'%(rstar,mp.nstr(L[6],10),mp.nstr(L[7],10),mp.nstr(L[8],10),verdict(ok)))
# C4 Thm 4.1 Part 2 conclusion (p.20): Omega_1<...<Omega_{d-1} at m=d-1 (k=0, d=4..30 via omega_np at (9,0))
U=np.linspace(1e-6,4.0,400001); zf=float(z[7]); n=2*d-1
P,_,_=onp.Phi_derivs(U); gv=U*P*onp.He(n,zf*U/np.sqrt(2)); h=U[1]-U[0]
cum=np.concatenate([[0],np.cumsum((gv[1:]+gv[:-1])*h/2)]); W=cum[-1]-cum; kp=onp.kappa0p(U)
idx=[i for i in np.where(np.sign(W[1:])*np.sign(W[:-1])<0)[0] if U[i]>0.02]
b=[0]+[i+1 for i in idx]+[len(U)-1]
Om=[np.trapezoid(np.abs(W[b[r]:b[r+1]+1])*kp[b[r]:b[r+1]+1],U[b[r]:b[r+1]+1]) for r in range(len(b)-1)]
print('C4 Thm 4.1(2) conclusion Omega increasing at (9,0) m=8 (%d lobes):'%len(Om), verdict(all(Om[i]<Om[i+1] for i in range(len(Om)-1))))
# C5 Thm 5.3 Part 2 (p.28): s_m <= mu_{d+k} K3 with K3 < 1.8165 for all interior m
mp.mp.dps=60
for d,k in [(22,0),(23,0)]:
    zz=zroots(d,k)[d-4]; s=zz*zz/2; r=s/mu(d+k)
    print('C5 Thm 5.3(2) ceiling s_{d-3}/mu_{d+k} < 1.8165 at (%d,%d), m=%d: s=%s mu=%s ratio=%s: %s'%(d,k,d-3,mp.nstr(s,8),mp.nstr(mu(d+k),8),mp.nstr(r,6),verdict(r<mp.mpf('1.8165'))))
# C6 Thm 5.3 Part 3 (pp.28-30) disk |z|=eta|u_s|^2 <= 0.7945 (u_s = saddle of the printed equation, App. A.2/A.4)
for d,k,m in [(21,0,18),(22,0,19),(30,0,27)]:
    om=float(omega(d,k,m)); M=4*d+2*k-1
    u,res,Fu,S,Wv,V,gg,eta,nu2=solve(om,M,k); c=consts(u,Fu,S,Wv,V,gg,eta,nu2,om,M,k)
    print('C6 Thm 5.3(3) disk eta|u_s|^2<=0.7945 at (%d,%d,%d): omega=%.4f u_s=%.5f%+.5fi |F|=%.0e alpha=%.5f: %s'%(d,k,m,om,u.real,u.imag,res,c['alpha'],verdict(c['alpha']<=0.7945)))
# C7 Lemma 6.2b (p.34)/Lemma 6.4 (p.40): Delta Theta_0(u_s)>=0.4388, tau_s>=0.420; Zone I: tau_s>=0.486692
for d,k,m in [(9,0,1),(8,1,1),(4,5,1),(4,8,1)]:
    om=float(omega(d,k,m)); M=4*d+2*k-1
    u,res,Fu,S,Wv,V,gg,eta,nu2=solve(om,M,k); c=consts(u,Fu,S,Wv,V,gg,eta,nu2,om,M,k)
    ok=c['dTheta0']>=0.4388 and c['tau']>=0.420 and (c['sigma']>0.65 or (c['tau']>=0.486692 and c['dTheta0']>=0.4583))
    print('C7 Lemma 6.4 floors at (%d,%d,%d): omega=%.4f sigma=%.4f tau=%.4f dTheta0=%.4f: %s'%(d,k,m,om,c['sigma'],c['tau'],c['dTheta0'],verdict(ok)))
# C8 App. D (p.66): anchor u_s=0.470776360737+0.486692548117i solves the printed saddle equation F=0 (A.2)
from saddle_scan import F_and_Fu
F,Fu,*_=F_and_Fu(complex(0.470776360737099964,0.486692548116874293),12.272,35.0,0,None)
print('C8 App. D anchor is a root of printed F (A.2): F=%.4f%+.4fi: %s'%(F.real,F.imag,verdict(abs(F)<1e-8)))
# C9 consistency (not a refutation): exact B/A = 2x J''_{d,k}(x)/(d J_{d,k}(x)) >= (1-0.853) Lambda_pos(u_s) at (9,0)
d,k=9,0; ok=True
for m in range(1,d-2):
    zz=zroots(d,k)[m-1]; x=-zz*zz/2; c9=[mp.binomial(d,j)*g[k+j] for j in range(d+1)]; Pp=lambda t: mp.polyval(c9[::-1],t)
    BA=2*x*mp.diff(Pp,x,2)/(d*Pp(x))
    om=float(omega(d,k,m)); u,res,Fu,S,Wv,V,gg,eta,nu2=solve(om,4*d-1,0); c=consts(u,Fu,S,Wv,V,gg,eta,nu2,om,4*d-1,0)
    th=math.atan2(c['tau'],c['sigma']); Lp=(c['sigma']**2+c['tau']**2)*math.sin(c['Phia']+0.09756+th)/math.sin(c['Phia']+0.09756-th)
    ok = ok and float(BA) >= (1-0.853)*Lp
print('C9 Thm 6.3 consequence B/A >= (1-0.853)Lambda_pos at (9,0), m=1..6:', verdict(ok))
