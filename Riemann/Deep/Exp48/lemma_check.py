#!/usr/bin/env python3
"""exp48 lemma_check.py: test 2605.26020 Lemma property-of-f, L(f) | 12 u E(D), on every D = D0 f^2 with E(D) | 8 and
|D| <= 10^10 from our unconditional scan (logs/scan_E_divides_8_upto_1e10.txt). Run: python3 lemma_check.py"""
from math import gcd
L = [tuple(map(int, l.split())) for l in open('logs/scan_E_divides_8_upto_1e10.txt')]
def kron(D, p):
    if p == 2:
        return 0 if D % 2 == 0 else (1 if D % 8 in (1, 7) else -1)
    r = D % p
    return 0 if r == 0 else (1 if pow(r, (p - 1) // 2, p) == 1 else -1)
def factor(n):
    f = {}; p = 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def fund_part(D):  # D = D0 f^2 with D0 fundamental
    n = -D; fa = factor(n); m = 1; f = 1
    for p, e in fa.items():
        m *= p ** (e % 2); f *= p ** (e // 2)
    D0 = -m
    if D0 % 4 != 1:  # D0 = -m with m = 1, 2 mod 4 needs factor 4
        D0 *= 4; assert f % 2 == 0, D; f //= 2
    assert D0 * f * f == D
    return D0, f
bad = 0; nonfund = 0
for (D, fund, h, E, n2, n4) in L:
    D0, f = fund_part(D)
    assert (f == 1) == bool(fund), (D, D0, f, fund)
    if f == 1: continue
    nonfund += 1
    u = 6 if D0 == -3 else (4 if D0 == -4 else 2)
    Lf = 1
    for p, k in factor(f).items():
        Lp = p ** (k - 1) * (p - kron(D0, p)); Lf = Lf * Lp // gcd(Lf, Lp)
    if (12 * u * E) % Lf: bad += 1; print('VIOLATION', D, D0, f, E, Lf)
print(f"non-fundamental D with E(D) | 8, |D| <= 1e10: {nonfund}; violations of L(f) | 12 u E(D): {bad}")
