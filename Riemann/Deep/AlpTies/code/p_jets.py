#!/usr/bin/env python3
# Probe B (toy analogue of MomAt, NOT a faithful instance): jets of the CM kernel wp for L = Omega*Z[i] (y^2 = 4x^3 - 4x,
# g2 = 4, g3 = 0) in the integral normalisation: D^a(wp - z^-2)(0) = a! * c_{a/2}. MomAt asks p-adic integrality of kernel jets
# (up to a bounded table normalisation) at a SPLIT level prime; the inert gain a/p + a/p^3 + ... is a different phenomenon.
# Exact rationals; reports min over a of v_p(a! c) - (expected main term) for split p = 5, 13, 17 and inert p = 3, 7, 11.
from fractions import Fraction as Fr
from math import log
def vp(x, p):
    if x == 0: return None
    n, d = x.numerator, x.denominator; v = 0
    while n % p == 0: n //= p; v += 1
    while d % p == 0: d //= p; v -= 1
    return v
M = 300
c = [Fr(0)] * (M + 1); c[1] = Fr(4, 20); c[2] = Fr(0)
for n in range(3, M + 1):
    c[n] = Fr(3, (2 * n + 3) * (n - 2)) * sum(c[k] * c[n - 1 - k] for k in range(1, n - 1))
fact = [1] * (2 * M + 1)
for i in range(1, 2 * M + 1): fact[i] = fact[i - 1] * i
def odd_main(a, p):
    s, q = 0.0, p
    while q <= a: s += a / q; q *= p * p
    return s
for p, kind in ((5, 'split'), (13, 'split'), (17, 'split'), (3, 'inert'), (7, 'inert'), (11, 'inert')):
    worst, worst_a, cnt = None, None, 0
    for n in range(1, M + 1):
        if c[n] == 0: continue
        a = 2 * n; v = vp(fact[a] * c[n], p)
        main = odd_main(a, p) if kind == 'inert' else 0.0
        dev = (v - main) / max(1.0, log(a, p)); cnt += 1
        if worst is None or dev < worst: worst, worst_a = dev, a
    print(f"p={p} {kind}: {cnt} jets a<= {2*M}; min over a of (v_p(a! c_(a/2)) - main)/log_p a = {worst:.3f} at a={worst_a}; main = {'sum a/p^odd' if kind=='inert' else '0'}")
