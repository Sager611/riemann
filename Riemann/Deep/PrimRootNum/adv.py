#!/usr/bin/env python3
"""
adv.py -- adversarial sanity checks (not part of the proof).
 1. Exact g*(p) for primes with p - 1 highly composite (p = k P_m + 1, smallest k, m = 1..45), Sophie Germain primes p = 2l + 1,
    Fermat primes, and p = 2^a 3^b + 1; compared with Theorem A (kappa_w table) and Theorem C.
 2. Numerical check of the Rosser-Schoenfeld inequalities R1-R5 at every prime x in [286, 10^7] (R5: x >= 41).
Run: uv run --no-project --with numpy --with sympy python3 adv.py
"""
import math, os, sys
import numpy as np
from sympy import isprime, primerange, factorint
B = 0.8733819666132076; b = 1 - B
HERE = os.path.dirname(os.path.abspath(__file__))
kap = {}
for line in open(os.path.join(HERE, "logs", "kappa_table.txt")):
    if line.startswith("#"):
        continue
    f = line.split(); kap[int(f[0])] = float(f[6])

def gstar(p):
    ls = list(factorint(p - 1).keys())
    for r in primerange(2, 10 ** 6):
        if r == p:
            continue
        if all(pow(r, (p - 1) // l, p) != 1 for l in ls):
            return r, len(ls)
    return None, len(ls)

rows = []
P = 1
for m, q in enumerate(primerange(2, 400), start=1):
    P *= q
    if m > 45:
        break
    k = 1
    while not isprime(k * P + 1):
        k += 1
    p = k * P + 1
    g, w = gstar(p)
    L = math.log(p)
    boundA = (kap[w] * L) ** (1 / b)
    rows.append(("primorial m=%d" % m, p, w, g, boundA))
for l in [5, 11, 23, 29, 41, 53, 83, 89, 113, 131, 173, 179, 191, 233, 239, 251, 281, 293, 10 ** 12 + 39, 2 ** 89 - 1]:
    if isprime(l) and isprime(2 * l + 1):
        p = 2 * l + 1; g, w = gstar(p); rows.append(("Sophie Germain", p, w, g, (kap[w] * math.log(p)) ** (1 / b)))
for p in [3, 5, 17, 257, 65537]:
    g, w = gstar(p); rows.append(("Fermat", p, w, g, (kap[w] * math.log(p)) ** (1 / b)))
for a in range(1, 60, 7):
    for c in range(0, 40, 9):
        p = 2 ** a * 3 ** c + 1
        if isprime(p) and p > 3:
            g, w = gstar(p); rows.append(("2^a3^b+1", p, w, g, (kap[w] * math.log(p)) ** (1 / b)))
worst = 0
for (kind, p, w, g, bd) in rows:
    r = math.log(g) / math.log(bd)
    worst = max(worst, r)
    assert g <= bd
print("%d adversarial primes (log10 p up to %.0f): all g*(p) <= (kappa_w log p)^(1/(1-B)); max log g*/log bound = %.4f" % (len(rows), max(math.log10(r[1]) for r in rows), worst))
for (kind, p, w, g, bd) in rows[::5]:
    print("  %-16s log10 p = %7.2f  omega = %2d  g* = %4d  bound = 10^%.1f" % (kind, math.log10(p), w, g, math.log10(bd)))
# 2. Rosser-Schoenfeld checks
N = 10 ** 7
s = np.ones(N + 1, dtype=bool); s[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if s[i]:
        s[i * i::i] = False
pr = np.nonzero(s)[0].astype(np.float64)
S = np.cumsum(1.0 / pr); lp = np.log(pr); V = np.exp(np.cumsum(np.log1p(-1.0 / pr))); th = np.cumsum(lp)
M = 0.2614972128476427837554; g = 0.5772156649015328606
m = pr >= 286
x = pr[m]; lx = np.log(x)
e1 = S[m] - np.log(lx) - M
print("R1/R2 at primes 286 <= x <= 1e7: max |E| log^2 x = %.4f (< 1/2 required)" % np.max(np.abs(e1) * lx ** 2))
v = V[m] * np.exp(g) * lx
print("R3/R4: max |V e^g log x - 1| log^2 x = %.4f (< 1/2 required)" % np.max(np.abs(v - 1) * lx ** 2))
m5 = pr >= 41
print("R5: min theta(x)/(x(1 - 1/log x)) at primes 41 <= x <= 1e7 = %.5f (> 1 required)" % np.min(th[m5] / (pr[m5] * (1 - 1 / lp[m5]))))
