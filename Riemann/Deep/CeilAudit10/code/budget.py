import mpmath as mp
mp.mp.dps = 40
d = {}
for line in open('Riemann/Deep/ExplorePositivity/logs/d2n_zeta.txt'):
    k, v = line.split()[:2]; d[int(k)] = mp.mpf(v)
rho = [d[2*n]**2/(d[2*n-2]*d[2*n+2]) for n in range(1, 100)]
P = mp.mpf(1); out = []
for m in range(2, 100):
    P *= rho[m-2]
    kap = lambda k: mp.sqrt(mp.rf(1.5, k))
    b = (kap(0)/kap(2))*(kap(2*m)/kap(2*m-2))*P/m
    if m in (2,3,4,5,6,8,10,20,50,99): out.append('m=%d: bound mu_m/mu_{m-1} <= %s * mu1/mu0 (prod rho = %s)' % (m, mp.nstr(b,4), mp.nstr(P,4)))
print('\n'.join(out))
# MP lambda=3/4 Lemma D cumulative cap with mu0=mu1 (N -> inf), using only c log-convexity
lam = mp.mpf(3)/4; cap = mp.mpf(1); o2 = []
for m in range(2, 11):
    cap *= 2*lam*(2*lam+1)*(2*m-1)/((2*lam+2*m-2)*(2*lam+2*m-1))
    o2.append('%d:%s' % (m, mp.nstr(cap, 3)))
print('MP lam=3/4, mu0=mu1=1: mu_m <=', ', '.join(o2))
