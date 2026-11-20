#!/usr/bin/env python3
"""Test C: normalisation of the core inequality ASK <= MAIN + Err in the one-character toy case.
FirstPrimeNeSecondPrime.lean: ASK = sum over INERT p, odd i, of 2 log p (a/p^i - w)^+ (+ ramified d term);
MAIN = a (log a - 1 - gamma + ell(chC) - mean log w), ell = L'/L(1, chi) (TieArith/Ident).
With w = 1, a = y: sum_n Lambda(n)(1 - chi(n))(y/n - 1) - y(log y - 1 - gamma + L'/L(1,chi)) = Gz(y) - Gchi(y),
which under GRH is O(y^(1/2)). A sign error in MAIN would leave a linear term 2 y L'/L(1,chi)."""
import numpy as np, mpmath as mp, math, sys
Y = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
def kron(D, n):  # Kronecker symbol (D/n) for fundamental D<0, n>0
    from sympy import jacobi_symbol
    if n % 2 == 0:
        if D % 2 == 0: return 0
        v = 1 if D % 8 in (1, 7) else -1
        e = 0
        while n % 2 == 0: n //= 2; e += 1
        return v ** e * kron(D, n) if n > 1 else v ** e
    return jacobi_symbol(D % n, n) if n > 1 else 1
def chi_table(D, N):
    q = abs(D); t = np.array([kron(D, r) if math.gcd(r, q) == 1 else 0 for r in range(q)], dtype=np.int8)
    return t
sieve = np.ones(Y + 1, dtype=bool); sieve[:2] = False
for p in range(2, int(Y ** 0.5) + 1):
    if sieve[p]: sieve[p * p::p] = False
primes = np.nonzero(sieve)[0]
mp.mp.dps = 60  # the Hurwitz poles at s = 1 cancel in the character sum; evaluate at s = 1 + 1e-30
gamma = float(mp.euler)
for d in (7, 23, 71, 151, 199):
    D = -d; t = chi_table(D, Y)
    # L'/L(1,chi) via Hurwitz zeta: L(s) = q^-s sum_r chi(r) zeta(s, r/q)
    L = lambda s: sum(int(t[r]) * mp.zeta(s, mp.mpf(r) / d) for r in range(1, d)) / mp.mpf(d) ** s
    s1 = 1 + mp.mpf(10) ** -30
    lpl = float(mp.diff(L, s1) / L(s1))
    out = []
    for y in (10 ** 4, 10 ** 5, 10 ** 6, Y):
        lhs = 0.0; sz = 0.0; schi = 0.0
        for p in primes[primes <= y]:
            lp = math.log(p); c = int(t[p % d]); pk = int(p); i = 1
            while pk <= y:
                ci = c ** i if c != 0 else 0
                w = y / pk - 1.0
                lhs += lp * (1 - ci) * w; sz += lp * w; schi += lp * ci * w
                pk *= int(p); i += 1
        main = y * (math.log(y) - 1 - gamma + lpl)
        wrong = y * (math.log(y) - 1 - gamma - lpl)
        gz = sz - y * math.log(y) + (1 + gamma) * y
        gchi = schi + lpl * y
        out.append((y, round((lhs - main) / math.sqrt(y), 3), round((gz - gchi) / math.sqrt(y), 3), round((lhs - wrong) / math.sqrt(y), 1)))
    print(f"d={d:4d} L'/L(1,chi)={lpl:+.4f}  (y, (ASK-MAIN)/sqrt y, (Gz-Gchi)/sqrt y, (ASK-MAIN_wrongsign)/sqrt y): {out}")
