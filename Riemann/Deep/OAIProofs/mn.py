"""Characters of symmetric groups by the Murnaghan-Nakayama rule (beta-set form),
written independently for the audit. Exact integers only."""
from functools import lru_cache
from math import factorial
from collections import Counter

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k):
            yield (k,) + rest

def beta(lam, L):
    lam = tuple(lam) + (0,)*(L-len(lam))
    return tuple(sorted((lam[i] + (L-1-i) for i in range(L)), reverse=True))

@lru_cache(maxsize=None)
def _mn(bset, mu):
    # bset: tuple of distinct nonnegative ints (beta-set), mu: tuple of parts still to remove
    if not mu:
        return 1
    r = mu[0]; rest = mu[1:]
    s = set(bset); tot = 0
    for b in bset:
        c = b - r
        if c >= 0 and c not in s:
            sign = (-1) ** sum(1 for x in bset if c < x < b)
            nb = tuple(sorted((s - {b}) | {c}, reverse=True))
            tot += sign * _mn(nb, rest)
    return tot

def chi(lam, mu):
    """character of S^lam at cycle type mu"""
    L = max(len(lam), 1) + 0
    return _mn(beta(lam, L), tuple(sorted(mu, reverse=True)))

def zee(mu):
    c = Counter(mu); z = 1
    for k, m in c.items():
        z *= k**m * factorial(m)
    return z

def class_size(mu):
    return factorial(sum(mu)) // zee(mu)

def dim(lam):
    n = sum(lam)
    return chi(lam, (1,)*n)

def conj(lam):
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()

def dominates(a, b):
    """a dominance-greater-or-equal b (same size)"""
    sa = sb = 0
    for i in range(max(len(a), len(b))):
        sa += a[i] if i < len(a) else 0
        sb += b[i] if i < len(b) else 0
        if sa < sb: return False
    return True
