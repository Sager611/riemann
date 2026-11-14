#!/usr/bin/env python3
# values.py -- oddraudit: independent recomputation (mpmath interval arithmetic iv, 80 digits; gamma hard-coded to 60 digits)
# of every constant of PROOF.md §1 and of the eps-cost enumeration (exponent sets and combinations with cost <= delta0).
#   uv run --no-project --with mpmath python3 values.py
import math, itertools
from mpmath import mp, iv, mpf
mp.dps = 80; iv.dps = 80
G60 = '0.577215664901532860606512090082402431042159335939923598805767'
assert abs(mp.euler - mpf(G60)) < mpf(10)**-59
gam = iv.mpf([mpf(G60) - mpf(10)**-58, mpf(G60) + mpf(10)**-58])
L2 = iv.log(2)
def sg(p, a): return (p**(a + 1) - 1)//(p - 1)
fac = {3: 4, 5: 3, 7: 2, 11: 1, 13: 1, 17: 1, 19: 1, 23: 1, 29: 1, 31: 1, 37: 1, 41: 1, 43: 1, 47: 1, 53: 1, 59: 1, 61: 1}
def val(f):
    n = 1; s = 1
    for p, a in f.items(): n *= p**a; s *= sg(p, a)
    return n, s
N, sN = val(fac); c0 = 67*N
def marg(n, s): return iv.log(iv.mpf(s)/iv.mpf(n)) - (gam - L2) - iv.log(iv.log(iv.log(iv.mpf(n))))
def st(x, d=12): return f'[{mp.nstr(mpf(x.a), d)}, {mp.nstr(mpf(x.b), d)}]'
print('sigma(N) =', sN, ' sigma(N)/N =', mp.nstr(mpf(sN)/N, 15))
d0 = marg(N, sN); print('delta0 = log-margin(N) =', st(d0), ' N violates:', d0.a > 0)
print('f21(N) = sigma(N)/(N loglog N) =', st(iv.mpf(sN)/N/iv.log(iv.log(iv.mpf(N))), 10), ' e^g/2 =', st(iv.exp(gam)/2, 10), ' f(N)/2 =', st(iv.mpf(sN)/N/iv.log(iv.log(iv.mpf(N)))/2, 6))
mc = marg(c0, val({**fac, 67: 1})[1]); print('log-margin(c0) =', st(mc), ' c0 satisfies:', mc.b < 0)
eps = iv.log(iv.mpf(68)/67)/iv.log(67)
tc = iv.log(iv.mpf(c0)); bound = 1/(tc*iv.log(tc))
print('eps* =', st(eps), ' 1/(tc0 log tc0) =', st(bound), ' eps* < bound:', eps.b < bound.a)
c73 = (1 + eps)*iv.log(73) - iv.log(74); print('c_73(1) =', st(c73), ' > delta0:', c73.a > d0.b)
def logg(p, a): return iv.log(sg(p, a)) - a*(1 + eps)*iv.log(p)
P = [p for p in range(3, 1000, 2) if all(p % q for q in range(3, int(p**0.5) + 1, 2))]
# maximality of N's exponents at eps* (strict except the 67 tie)
okmax = True
for p in P:
    a0 = fac.get(p, 0)
    up = logg(p, a0 + 1) - logg(p, a0)
    if p == 67: assert abs(mpf(up.a)) < mpf(10)**-70 and abs(mpf(up.b)) < mpf(10)**-70
    elif not up.b < 0: okmax = False
    if a0 and not (logg(p, a0) - logg(p, a0 - 1)).a > 0: okmax = False
print('N maximises g_p at eps* for odd p < 1000 (tie only at 67):', okmax)
allowed = {}
for p in [q for q in P if q <= 71]:
    a0 = fac.get(p, 0)
    allowed[p] = [(a, logg(p, a0) - logg(p, a)) for a in range(0, 60) if (logg(p, a0) - logg(p, a)).a <= d0.b]
print('exponents within delta0:', {p: [(a, mp.nstr(mpf(c.a), 6)) for a, c in v] for p, v in allowed.items() if len(v) > 1})
# all combinations (product over primes with >1 option; others fixed)
multi = [p for p in allowed if len(allowed[p]) > 1]
vecs = []
for choice in itertools.product(*[allowed[p] for p in multi]):
    cost = sum((c for _, c in choice), iv.mpf(0))
    if cost.a <= d0.b:
        f = {p: fac.get(p, 0) for p in allowed}
        for p, (a, _) in zip(multi, choice): f[p] = a
        f = {p: a for p, a in f.items() if a}
        vecs.append((f, cost))
print(f'{len(vecs)} vectors with total cost <= delta0')
for f, cost in sorted(vecs, key=lambda t: val(t[0])[0]):
    n, s = val(f); m = marg(n, s)
    where = 'gap' if N < n < c0 else ('=N' if n == N else ('<N' if n < N else ('=c0' if n == c0 else '>c0')))
    print(f'   n/N = {mp.nstr(mpf(n)/N, 8):>10} {where:4s} cost {mp.nstr(mpf(cost.a), 6):>12}  log-margin {st(m, 8)}  (OR) holds: {m.b < 0}')
