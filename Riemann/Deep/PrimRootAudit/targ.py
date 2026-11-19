#!/usr/bin/env python3
"""targ.py (ppraudit): targeted primes where kappa_w is tight: p = 1 + k P_m (least k, m = 1..80, p up to ~10^170),
Fermat primes, p = 2q + 1 and p = 4q + 1 (q prime) near 10^20, 10^40, 10^60. Exact g*(p) (p - 1 fully factored by construction),
compared with Theorem A (kappa_w from the certificate's table), B' and D.  Pure Python."""
import math, random
INV_B = 1 / (1 - 0.8733819666132076)
kap = {}
for line in open("../PrimRootNum/logs/kappa_table.txt"):
    if not line.startswith("#"):
        f = line.split(); kap[int(f[0])] = float(f[6])
def mr(n, bases=(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53)):
    if n < 2: return False
    for p in bases:
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in bases:
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True
PR = [q for q in range(2, 3000) if mr(q)]
def factor_small(k):
    fs, d = set(), 2
    while d * d <= k:
        while k % d == 0: fs.add(d); k //= d
        d += 1
    if k > 1: fs.add(k)
    return fs
def gstar(p, fs):
    for r in PR:
        if all(pow(r, (p - 1) // l, p) != 1 for l in fs): return r
    raise SystemExit("no root found")
worst = [0, None]
def test(p, fs, tag):
    w = len(fs); g = gstar(p, fs); L = math.log(p)
    rA = math.log(g) / (INV_B * math.log(kap[w] * L)) if w <= 10000 else None
    rD = math.log(g) / (3.4 * INV_B * math.log(1.956 * L))
    if rA > worst[0]: worst[0], worst[1] = rA, (tag, w, g)
    return w, g, rA, rD
out = []
P = 1
for m in range(1, 81):
    P *= PR[m - 1]
    k = 1
    while not mr(1 + k * P): k += 1
    p = 1 + k * P
    fs = set(PR[:m]) | factor_small(k)
    w, g, rA, rD = test(p, fs, "1+%d*P_%d" % (k, m))
    if m in (1, 2, 3, 5, 8, 10, 20, 27, 40, 52, 60, 70, 80): out.append("p = 1 + %d P_%d (10^%.1f): omega = %d, g* = %d, logratio A %.4f, D %.4f" % (k, m, math.log10(p), w, g, rA, rD))
for F in (3, 5, 17, 257, 65537):
    w, g, rA, rD = test(F, {2}, "Fermat %d" % F); out.append("Fermat p = %d: g* = %d, logratio A %.4f" % (F, g, rA))
random.seed(1)
for e in (20, 40, 60):
    for c in (2, 4):
        q = random.randrange(10 ** e, 2 * 10 ** e) | 1
        while not (mr(q) and mr(c * q + 1)): q += 2
        p = c * q + 1; w, g, rA, rD = test(p, {2, q}, "%dq+1" % c)
        out.append("p = %dq + 1 ~ 10^%d (q prime): omega = %d, g* = %d, logratio A %.4f" % (c, e, w, g, rA))
print("\n".join(out)); print("worst log g*/log (kappa_w log p)^(1/b) over all targeted primes: %.4f at %s" % (worst[0], worst[1]))
