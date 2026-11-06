# t1_verify.py (explore21): independent checks of the D = -15 Epstein zeros found by t1_lattice.py (python-flint):
# (1) coefficient identity r_Q(n) = (1*chi_-15)(n) -/+ (chi_-3*chi_5)(n) for n <= 3000 (Q2 = 2x^2+xy+2y^2: minus; Q1 = x^2+xy+4y^2: plus);
# (2) the zeros re-checked with mpmath (its own zeta and Dirichlet L) at 30 digits; (3) the sigma > 1 zero of Z_Q1 by the
#     absolutely convergent Euler-product ratio R = L(chi_-3)L(chi_5)/(zeta L(chi_-15)) = -1 (primes < 2e6).
import mpmath as mp
import numpy as np
def chi_m15(n): return [0,1,1,0,1,0,0,-1,1,0,0,-1,0,-1,-1][n % 15]
def chi_m3(n): return [0,1,-1][n % 3]
def chi_5(n): return [0,1,-1,-1,1][n % 5]
Nmax = 3000
r1 = np.zeros(Nmax+1, int); r2 = np.zeros(Nmax+1, int)
R = 80
for m in range(-R, R+1):
    for n in range(-R, R+1):
        q1 = m*m + m*n + 4*n*n; q2 = 2*m*m + m*n + 2*n*n
        if 0 < q1 <= Nmax: r1[q1] += 1
        if 0 < q2 <= Nmax: r2[q2] += 1
bad1 = bad2 = 0
for n in range(1, Nmax+1):
    a = sum(chi_m15(d) for d in range(1, n+1) if n % d == 0)
    b = sum(chi_m3(d)*chi_5(n//d) for d in range(1, n+1) if n % d == 0)
    if r1[n] != a + b: bad1 += 1
    if r2[n] != a - b: bad2 += 1
print('coefficient identity n <= %d: Q1 mismatches %d, Q2 mismatches %d' % (Nmax, bad1, bad2))
mp.mp.dps = 30
c15 = [chi_m15(n) for n in range(15)]; c3 = [chi_m3(n) for n in range(3)]; c5 = [chi_5(n) for n in range(5)]
def ZQ(s, sg):
    return mp.zeta(s)*mp.dirichlet(s, c15) + sg*mp.dirichlet(s, c3)*mp.dirichlet(s, c5)
for (sg, name, s0) in ((-1, 'Q2', mp.mpc('0.7580708934', '24.4828213943')), (1, 'Q1', mp.mpc('0.6955863428', '20.3459678815')),
                       (1, 'Q1', mp.mpc('0.8000109882', '12.0385986287')), (1, 'Q1', mp.mpc('1.0259726162', '61.4227983088'))):
    r = mp.findroot(lambda s: ZQ(s, sg), s0)
    print('%s zero (mpmath, 30 digits): %s   |Z| = %s' % (name, mp.nstr(r, 20), mp.nstr(abs(ZQ(r, sg)), 3)))
    if name == 'Q1' and mp.re(r) > 1:
        # Euler-product ratio at the root, primes < 2e6
        s = r
        sv = np.ones(2*10**6+1, bool); sv[:2] = False
        for i in range(2, int((2*10**6)**0.5)+1):
            if sv[i]: sv[i*i::i] = False
        P = np.nonzero(sv)[0]
        sc = complex(s)
        lr = 0j
        for p in P:
            if p == 3 or p == 5: lr += np.log((1 - p**(-sc))/(1 + p**(-sc))); continue
            if p % 15 in (2, 8): lr += 2*np.log((1 - p**(-sc))/(1 + p**(-sc)))
        print('   Euler-product ratio R(s) over primes < 2e6: %s  (zero of Z_Q1 = zeta_K (1 + R) needs R = -1)' % (np.exp(lr)))
