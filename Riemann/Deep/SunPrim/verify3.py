#!/usr/bin/env python3
"""Two more false generalizations stated by Z.-W. Sun in OEIS comments (checked 2026-10-22; no correction note there).
A228766 (comment): if a_1..a_n are n > 2 distinct elements of a finite additive abelian group G with n odd or |G| not
  divisible by n, then some circular permutation b_1..b_n has b_1+b_2, ..., b_n+b_1 pairwise distinct.
A229827 (ii): F any field with |F| > 7, A a finite subset with |A| = n > 2, P a polynomial (deg < p-1 in characteristic p)
  not of the form c - x: some circular permutation has P(a_1)+a_2, ..., P(a_n)+a_1 pairwise distinct."""
from itertools import permutations, combinations
from fractions import Fraction as Fr

def circ_ok(A, f, mod=None):
    A = list(A)
    for p in permutations(A[1:]):
        b = (A[0],) + p
        s = [f(b[i], b[(i + 1) % len(b)]) for i in range(len(b))]
        if mod:
            s = [x % mod for x in s]
        if len(set(s)) == len(s):
            return b
    return None

assert circ_ok([0, 1, 3, 4], lambda x, y: x + y, 6) is None      # n = 4 does not divide 6
assert circ_ok([0, 1, 5, 6], lambda x, y: x + y, 10) is None     # n = 4 does not divide 10
fails = {}
for m in range(3, 13):
    for n in range(3, m + 1):
        if n % 2 == 1 or m % n != 0:
            c = [A for A in combinations(range(m), n) if circ_ok(A, lambda x, y: x + y, m) is None]
            if c:
                fails[(m, n)] = len(c)
print("A228766 general form fails; failing n-subsets of Z/m (m <= 12) by (m, n):", fails)

two_x = lambda x, y: 2 * x + y   # P(x) = 2x
assert circ_ok([0, 1, 2], two_x, 11) is None
assert circ_ok([Fr(0), Fr(1), Fr(2)], two_x) is None
print("A229827(ii) fails: P(x) = 2x, A = {0,1,2} in F_11 and in Q (every 3-term progression: both circular orders repeat a sum)")
