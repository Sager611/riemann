#!/usr/bin/env python3
"""lsanity.py -- sanity check of Theorem 4 (|L'/L(1,chi)| <= (1/(1-B)) log log q + 10.53, chi primitive mod q):
all nonprincipal characters of prime moduli 3 <= q <= 200 (all primitive).
L(1,chi) = -(1/q) sum_a chi(a) psi(a/q),  L'(1,chi) = -log q L(1,chi) - (1/q) sum_a chi(a) gamma_1(a/q)
(Laurent expansion zeta(s,a) = 1/(s-1) - psi(a) - gamma_1(a)(s-1) + ...; mpmath digamma and stieltjes).
Run: uv run --no-project --with mpmath python3 lsanity.py"""
import mpmath as mp, math, time
mp.mp.dps = 20
B = 0.8733819666132076; INVB = 1 / (1 - B)
t0 = time.time()
def primroot(p):
    fac = []; m = p - 1; d = 2
    while d * d <= m:
        if m % d == 0:
            fac.append(d)
            while m % d == 0: m //= d
        d += 1
    if m > 1: fac.append(m)
    for g in range(2, p):
        if all(pow(g, (p - 1) // f, p) != 1 for f in fac): return g
# self-test against a known value: L(1, chi_{-4}) = pi/4
chi4 = {1: 1, 3: -1}
L4 = -sum(chi4[a] * mp.digamma(mp.mpf(a) / 4) for a in (1, 3)) / 4
assert abs(L4 - mp.pi / 4) < 1e-15
worst = (0, None); worstratio = (0, None); n = 0
for q in [p for p in range(3, 201) if all(p % d for d in range(2, int(p ** 0.5) + 1))]:
    g = primroot(q); ind = {}; v = 1
    for k in range(q - 1):
        ind[v] = k; v = v * g % q
    ps = {a: mp.digamma(mp.mpf(a) / q) for a in range(1, q)}
    g1 = {a: mp.stieltjes(1, mp.mpf(a) / q) for a in range(1, q)}
    for j in range(1, q - 1):
        chi = {a: mp.exp(2j * mp.pi * j * ind[a] / (q - 1)) for a in range(1, q)}
        L = -sum(chi[a] * ps[a] for a in chi) / q
        Ld = -mp.log(q) * L - sum(chi[a] * g1[a] for a in chi) / q
        r = abs(Ld / L); n += 1
        bound = INVB * math.log(math.log(q)) + 10.53
        if r > worst[0]: worst = (float(r), (q, j))
        if r / bound > worstratio[0]: worstratio = (float(r / bound), (q, j, round(bound, 3)))
print("prime q in [3,200], %d nonprincipal chi: max |L'/L(1,chi)| = %.4f at (q, j) = %s; max ratio to the bound = %.4f at %s (%.0f s)"
      % (n, worst[0], worst[1], worstratio[0], worstratio[1], time.time() - t0))
