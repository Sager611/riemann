#!/usr/bin/env python3
"""indep_small.py (unitaudit) -- independent exact recheck of the small-Q claims of UnitCover sec 6 (numpy closure,
own code: S <- S | S[p^-1 x mod q]).  Claims: first 17 primes (2..59) cover every q <= 2000 and 16 do not;
16 least primes > 2000 cover, 15 do not; Q = 100/300/1000: first 11/13/15, least primes > Q: 9/12/14."""
import numpy as np, math, time
t0 = time.time()
def primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return [int(x) for x in np.nonzero(s)[0]]
PR = primes(30000)
def least_m(q, seq):
    """least m with <seq[:m] minus divisors of q> subset products = (Z/qZ)^x; None if seq too short"""
    units = np.array([math.gcd(x, q) == 1 for x in range(q)])
    ph = int(units.sum())
    if ph == 1: return 0
    S = np.zeros(q, bool); S[1] = True
    ar = np.arange(q, dtype=np.int64)
    for i, p in enumerate(seq):
        if q % p == 0: continue
        pinv = pow(p, -1, q)
        S = S | S[(ar * pinv) % q]          # x in pS  <=>  p^-1 x in S
        assert not (S & ~units).any()
        if int(S.sum()) == ph: return i + 1
    return None
def sweep(Qs, seq):
    worst = 0; arg = []
    for q in range(2, Qs + 1):
        m = least_m(q, seq)
        if m is None: return None, None
        if m > worst: worst, arg = m, [q]
        elif m == worst: arg.append(q)
    return worst, arg
first = PR[:80]
for Qs in (100, 300, 1000, 2000):
    above = [p for p in PR if p > Qs][:80]
    w1, a1 = sweep(Qs, first); w2, a2 = sweep(Qs, above)
    print("Q=%d: first primes need m=%d (last prime %d; hardest %s); primes > Q need m=%d (last %d; hardest %s)"
          % (Qs, w1, first[w1 - 1], a1[:6], w2, above[w2 - 1], a2[:6]), flush=True)
print("log2 phi(1999) = %.4f; log2 phi(q) max over q<=2000 = %.4f" % (math.log2(1998), math.log2(1998)))
print("done %.1f s" % (time.time() - t0))
