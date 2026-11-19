#!/usr/bin/env python3
"""bf.py (ppraudit): brute-force g*(p) for small p, two independent ways, to validate gstar.py. Usage: python3 bf.py N1 N2"""
import sys
N1, N2 = int(sys.argv[1]), int(sys.argv[2])
def isprime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True
PR = [q for q in range(2, 5000) if isprime(q)]
def order(r, p):
    x, k = r % p, 1
    while x != 1:
        x = x * r % p; k += 1
    return k
def gstar_order(p):
    for r in PR:
        if r % p and order(r, p) == p - 1: return r
def gstar_pow(p):
    n, fs, d = p - 1, [], 2
    while d * d <= n:
        if n % d == 0:
            fs.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: fs.append(n)
    for r in PR:
        if r % p and all(pow(r, (p - 1) // l, p) != 1 for l in fs): return r
a = [(p, gstar_order(p)) for p in range(3, N1 + 1) if isprime(p)]
b = [(p, gstar_pow(p)) for p in range(3, N1 + 1) if isprime(p)]
assert a == b, "order vs pow disagree"
rec = 0; out = []
for p in range(3, N2 + 1):
    if isprime(p):
        g = gstar_pow(p)
        if g > rec: rec = g; out.append((p, g))
if len(sys.argv) > 3:
    got = dict(tuple(map(int, l.split())) for l in open(sys.argv[3]))
    bad = [p for p in range(3, N2 + 1) if isprime(p) and got.get(p) != gstar_pow(p)]
    print("bf: per-prime comparison with %s for p <= %d: %d mismatches (%d entries)" % (sys.argv[3], N2, len(bad), len(got)))
print("bf: order == pow for p <= %d (%d primes); records p <= %d:" % (N1, len(a), N2), ", ".join("%d:%d" % t for t in out))
