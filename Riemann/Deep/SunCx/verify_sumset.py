"""Counterexample to Wang-Sun, arXiv 2210.12044 (EJC 33(2) 2026, #P2.64), Conjecture 1.1 (L) = OEIS A357130 Conj. 2(ii).
L(A_1..A_n) = {a_1+...+a_n : a_i in A_i, a_i != a_{i+1} (0<i<n)}; conjectured |L| >= min{p(G), sum|A_i| - 2n + 1 + (n mod 2)}.
Run: python3 verify.py   (stdlib only, ~1 min)"""
import itertools
def sums(As, mod=None, cyc=False):
    S = set(); n = len(As)
    def go(i, prev, first, s):
        if i == n:
            S.add(s % mod if mod else s); return
        for a in As[i]:
            if i > 0 and a == prev: continue
            if cyc and i == n - 1 and a == first: continue
            go(i + 1, a, a if i == 0 else first, s + a)
    go(0, None, None, 0)
    return S
# 1. In Z (p(G) = infinity): bound 3+2+3-6+1+1 = 4, but L = {2,3,4}.
A = [{0, 1, 2}, {0, 2}, {0, 1, 2}]
L = sums(A); bound = sum(map(len, A)) - 2 * 3 + 1 + 1
assert L == {2, 3, 4} and bound == 4
print("Z: L =", sorted(L), "|L| =", len(L), "< bound", bound)
# 2. Exhaustive in Z/p: count violating n-tuples of subsets (|A_i| >= 2).
def count(p, n, cyc):
    subs = [frozenset(c) for r in range(2, p + 1) for c in itertools.combinations(range(p), r)]
    bad = 0
    for As in itertools.product(subs, repeat=n):
        s = sum(map(len, As)); par = n % 2
        b = min(p, s - 2 * n + 1 + par if not cyc else s - 2 * n + (-1) ** n * (1 + par))
        if b > 0 and len(sums(As, p, cyc)) < b: bad += 1
    return bad
for p, n in [(5, 3), (7, 3), (5, 4)]:
    print("Z/%d n=%d: violations of (L) %d, of (C) %d" % (p, n, count(p, n, False), count(p, n, True)))
# expected: Z/5 n=3: 20, 0; Z/7 n=3: 203, 0; Z/5 n=4: 0, 0

