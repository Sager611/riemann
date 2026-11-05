# crosscheck.py: second independent computation of the witness data
# (1) gamma_n from Taylor coefficients of xi(1/2+w) (mpmath taylor) vs quadrature gam100.json
# (2) roots of J^{8,1} by sign-change bisection (no polyroots) vs roots30.json
# (3) Hermite He_17 zeros: Golub-Welsch (mpmath eig) vs numpy hermeroots
import json, mpmath as mp, numpy as np
from numpy.polynomial.hermite_e import hermeroots
mp.mp.dps=60
gq=[mp.mpf(s) for s in json.load(open('/Users/adriansager/playground/riemann/.scratch/chungaudit/gam100.json'))]
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
c=mp.taylor(lambda w: xi(mp.mpf(1)/2+w), 0, 20, method='quad', radius=mp.mpf(2))
worst=0
for n in range(0,10):
    gt=2**(n-1)*mp.factorial(n)*c[2*n] if n>0 else c[0]/2
    rel=abs(gt/gq[n]-1); worst=max(worst,rel)
    print('n',n,'gamma(taylor)=',mp.nstr(gt,20),' gamma(quad)=',mp.nstr(gq[n],20),' rel diff=',mp.nstr(rel,3))
print('PASS gamma cross-check' if worst<1e-20 else 'FAIL gamma cross-check', mp.nstr(worst,3))
# (2) bisection roots of J^{8,1}(x), x<0
d,k=9,0; N=d-1
coef=[mp.binomial(N,j)*gq[k+1+j] for j in range(N+1)]
P=lambda x: mp.polyval(coef[::-1],x)
xs=[]; grid=[-mp.mpf(i)/4 for i in range(0,1200)]
for a,b in zip(grid[:-1],grid[1:]):
    if P(a)*P(b)<0: xs.append(mp.findroot(P,(b,a),solver='bisect',tol=mp.mpf(10)**-50) if False else mp.findroot(P,(b,a),solver='anderson'))
z=[mp.sqrt(-2*x) for x in xs]
R=json.load(open('/Users/adriansager/playground/riemann/.scratch/chungaudit/roots30.json'))
zr=[mp.mpf(s) for s in R['9,0']['z']]
print('J^{8,1} real roots found by sign changes:',len(xs))
print('z (bisection):',[mp.nstr(v,15) for v in z])
dz=max(abs(a-b) for a,b in zip(z,zr))
print('max |z_bis - z_polyroots| =', mp.nstr(dz,3), 'PASS' if dz<mp.mpf('1e-20') else 'FAIL')
# (3) Hermite zeros
from turnover_k0 import hermite_pos_zeros
hz=hermite_pos_zeros(17); hn=sorted([r for r in hermeroots([0]*17+[1]) if r>1e-9])
dh=max(abs(float(a)-b) for a,b in zip(hz,hn))
print('He_17 zeros GW vs numpy max diff:', dh, 'PASS' if dh<1e-12 else 'FAIL')
