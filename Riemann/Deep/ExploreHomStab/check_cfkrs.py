# explore9: independent check of cfkrs.py -- evaluate the CFKRS recipe sum directly at small shifts
#   sum_eps prod_{eps_j=-1} X(1/2+a_j) R(eps a),  R = A * prod_{i<=j} zeta(1+b_i+b_j),  X(1/2+a) = e^{-aL} g(a)
# with the same truncated Euler product as cfkrs.py (primes <= P0 exactly, p^{-2} tail by prime sums).
# Usage: python3 check_cfkrs.py k L delta [P0]
import sys
from itertools import product
from math import factorial
import numpy as np, mpmath as mp
from scipy.special import gammaincc, exp1, gammaln
k = int(sys.argv[1]); Lval = mp.mpf(sys.argv[2]); delta = mp.mpf(sys.argv[3])
P0 = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
P1 = 10**7
mp.mp.dps = 60
def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
pr = sieve(P1)
small = [int(p) for p in pr if 2 < p <= P0]
big = pr[pr > P0].astype(float)
lb = np.log(big); M = 90
tau = []
for m in range(M):
    t = np.sum(np.exp(m * np.log(lb) - gammaln(m + 1) - 2 * lb))
    x0 = np.log(P1)
    t += exp1(x0) if m == 0 else gammaincc(m, x0) / m
    tau.append(mp.mpf(t))
def T2(w):
    return mp.fsum(tau[m] * (-w) ** m for m in range(M))
def hmon(m):
    return [e for e in product(range(m + 1), repeat=k) if sum(e) == m]
def c2(b):  # sum_{p>P0} p^{-2} c2(p^{-b}) via monomials
    def mon(e): return T2(mp.fsum(e[j] * b[j] for j in range(k)))
    h4 = mp.fsum(mon(e) for e in hmon(4)); h2 = mp.fsum(mon(e) for e in hmon(2))
    h2sq = mp.fsum(mon(tuple(x + y for x, y in zip(e1, e2))) for e1 in hmon(2) for e2 in hmon(2))
    S2 = mp.fsum(mon(tuple((2 if t == i else 0) + (2 if t == j else 0) for t in range(k))) for i in range(k) for j in range(i, k))
    return h4 - h2 - h2sq / 2 - S2 / 2
lsmall = [mp.log(p) for p in small]
def A(b):
    v = mp.mpf(1)
    for i in range(k):
        for j in range(i, k):
            v *= 1 - mp.power(2, -1 - b[i] - b[j])
    for p, lp in zip(small, lsmall):
        x = [mp.exp(-(mp.mpf(1) / 2 + bi) * lp) for bi in b]
        Pm = mp.mpf(1); Pp = mp.mpf(1)
        for xi in x:
            Pm /= (1 - xi); Pp /= (1 + xi)
        R = 1 + mp.mpf(p) / (p + 1) * ((Pm + Pp) / 2 - 1)
        C = mp.mpf(1)
        for i in range(k):
            for j in range(i, k):
                C *= 1 - x[i] * x[j]
        v *= R * C
    return v * mp.exp(c2(b))
def g(a):
    return mp.gamma((mp.mpf(1) / 2 - a) / 2) / mp.gamma((mp.mpf(1) / 2 + a) / 2)
def Xf(a):
    return mp.exp(-a * Lval) * g(a)
alpha = [delta * mp.sqrt(j + 1) for j in range(k)]
tot = mp.mpf(0)
for eps in product([1, -1], repeat=k):
    b = [e * a for e, a in zip(eps, alpha)]
    Rv = A(b)
    for i in range(k):
        for j in range(i, k):
            Rv *= mp.zeta(1 + b[i] + b[j])
    pref = mp.mpf(1)
    for e, a in zip(eps, alpha):
        if e == -1: pref *= Xf(a)
    tot += pref * Rv
print(k, float(Lval), float(delta), mp.nstr(tot, 20))
