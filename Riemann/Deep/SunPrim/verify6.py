#!/usr/bin/env python3
"""Z.-W. Sun, MathOverflow question 316143 (2018-11-25, unanswered), Conjecture: let lambda_1..lambda_n (n >= 3) be positive
integers with gcd 1 and lambda_(k+1) - lambda_k in {0, 1}; F a field with p(F) > n+1; A finite with |A| >= n + delta_(n,3). Then
|{sum lambda_k a_k : a_k distinct in A}| >= min{ p(F), (lambda_1+...+lambda_n)(|A|-n) + sum_(k<=n/2) (n+1-2k)(lambda_(n+1-k)-lambda_k) + 1 }.
The bound is the length of the integer interval [min, max] for A = {1..m}; the conjecture fails already for such A over Q."""
from itertools import permutations
from fractions import Fraction as Fr

def S(lam, A, mod=None):
    out = set()
    for t in permutations(A, len(lam)):
        v = sum(l * a for l, a in zip(lam, t))
        out.add(v % mod if mod else v)
    return out

def bound(lam, m):
    n = len(lam)
    return sum(lam) * (m - n) + sum((n + 1 - 2 * k) * (lam[n - k] - lam[k - 1]) for k in range(1, n // 2 + 1)) + 1

lam, A = (1, 1, 2), [0, 1, 2, 3]
assert bound(lam, 4) == 7
assert sorted(S(lam, A)) == [3, 4, 5, 7, 8, 9]            # over Z (= Q): 6 is missing, 6 < 7
for p in (5, 7, 11, 13):
    assert len(S(lam, A, p)) < min(p, 7)
for lam, m in (((2, 3, 3), 4), ((2, 3, 3), 7), ((3, 3, 4), 5), ((1, 1, 2, 3), 6)):
    s = len(S(lam, list(range(m))))
    print(lam, "A = {0..%d}:" % (m - 1), s, "<", bound(lam, m)) if s < bound(lam, m) else None
print("MO 316143 conjecture fails: lambda = (1,1,2), A = {0,1,2,3}: sums {3,4,5,7,8,9} (6 missing), 6 < 7, over Q and F_p, p = 5, 7, 11, 13")
