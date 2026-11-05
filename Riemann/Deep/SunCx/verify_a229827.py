"""Counterexamples to Z.-W. Sun's OEIS A229827 Conjecture (ii) (2013): for any field F with |F| > 7, finite A in F with
|A| = n > 2, and P(x) not of the form c - x (deg P < p - 1 in characteristic p), some circular ordering a_1..a_n of A has
the n sums P(a_1)+a_2, ..., P(a_n)+a_1 pairwise distinct.  Run: python3 verify_a229827.py (stdlib only)."""
from fractions import Fraction
import itertools
def ok_orders(A, P, mod=None):
    A = list(A); first = A[0]; good = []
    for rest in itertools.permutations(A[1:]):
        c = (first,) + rest
        s = [P(c[i]) + c[(i + 1) % len(c)] for i in range(len(c))]
        if mod: s = [x % mod for x in s]
        if len(set(s)) == len(s): good.append(c)
    return good
# 1. F = Q (or F_p, p >= 5), A = {0,1,2}, P(x) = 2x: both circular orders repeat a sum.
assert ok_orders([Fraction(0), Fraction(1), Fraction(2)], lambda x: 2 * x) == []
for p in [11, 13, 101]:
    assert ok_orders([0, 1, 2], lambda x: 2 * x, p) == []
print("Q and F_11, F_13, F_101: A = {0,1,2}, P(x) = 2x has no valid circular order")
# 2. F = F_11, A = {0,3,6}, P(x) = x^2.
assert ok_orders([0, 3, 6], lambda x: x * x, 11) == []
print("F_11: A = {0,3,6}, P(x) = x^2 has no valid circular order")

