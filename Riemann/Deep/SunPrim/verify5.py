#!/usr/bin/env python3
"""Z.-W. Sun, "On permutations of {1,...,n} and related topics", arXiv 1811.10503, Conjecture 4.9: let n >= 1, F a field with
p(F) > n+1 (p(F) = characteristic, or +infinity in characteristic 0), A a finite subset of F with |A| >= n + delta_(n,3). Then
S(A) = { sum_(k=1..n) k a_k : a_1..a_n distinct elements of A } has |S(A)| >= min{ p(F), (|A|-n) n(n+1)/2 + n(n^2-1)/6 + 1 }.
Counterexample: n = 2, A = {0, 1, 2}: S(A) = {a + 2b : a != b} has 4 elements, the bound is min{p(F), 5} = 5 for every
field with p(F) >= 5 (including Q). More generally any 3-term progression in such a field."""
from fractions import Fraction as Fr

def S(A, n2=True, mod=None):
    vals = set()
    for a in A:
        for b in A:
            if a != b:
                v = a + 2 * b
                vals.add(v % mod if mod else v)
    return vals

bound = (3 - 2) * 2 * 3 // 2 + 2 * (4 - 1) // 6 + 1          # = 5
assert bound == 5
assert len(S([Fr(0), Fr(1), Fr(2)])) == 4                    # Q
for p in (5, 7, 11, 13, 101):
    assert len(S([0, 1, 2], mod=p)) == 4 < min(p, bound)
for p in (7, 11, 13):
    for d in range(1, p):
        assert len(S([0, d, 2 * d % p], mod=p)) == 4          # every 3-term progression
print("1811.10503 Conj. 4.9 fails for n = 2, A = {0,1,2}: |S(A)| = 4 < 5 over Q and F_p (p >= 5); also every 3-term progression")
