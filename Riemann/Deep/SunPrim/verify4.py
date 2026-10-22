#!/usr/bin/env python3
"""Two more counterexamples (checked 2026-10-22).
Sun, arXiv 1309.1679, Conjecture 3.10(ii) (TeX label Conj3.12): for any set A of n > 3 distinct nonzero reals there is a circular permutation with
  the n adjacent differences a_i - a_(i+1) pairwise distinct AND the n adjacent products a_i a_(i+1) pairwise distinct,
  unless |A| = 4 and A = {+-s, +-t}.
OEIS A228728 (Sun, comment): for distinct reals a_1 < ... < a_n, if some cyclic arrangement has the n cyclic distances
  |b_i - b_(i+1)| pairwise distinct, then some arrangement with c_1 = a_1 and c_n = a_n has the n cyclic distances
  |c_1 - c_2|, ..., |c_(n-1) - c_n|, |c_n - c_1| pairwise distinct ("we are unable to find a counterexample")."""
from itertools import permutations
from fractions import Fraction as Fr

def orders(A):
    A = list(A)
    for p in permutations(A[1:]):
        yield (A[0],) + p

def ok312ii(A):
    for b in orders(A):
        n = len(b)
        D = [b[i] - b[(i + 1) % n] for i in range(n)]
        P = [b[i] * b[(i + 1) % n] for i in range(n)]
        if len(set(D)) == n and len(set(P)) == n:
            return b
    return None

for A in ([2, 3, 4, 6], [-8, -2, 1, 4], [Fr(-3), Fr(-1), Fr(-1, 3), Fr(1, 3), Fr(1), Fr(3)]):
    assert ok312ii(A) is None
print("3.10(ii) fails for {2,3,4,6}, {-8,-2,1,4} (n = 4, not of the form {+-s,+-t}) and {+-1/3, +-1, +-3} (n = 6)")

A = [0, 1, 2, 3, 4, 5, 8]
b = (0, 3, 4, 2, 8, 1, 5)
d = [abs(b[i] - b[(i + 1) % 7]) for i in range(7)]
assert len(set(d)) == 7          # antecedent holds: distances 3,1,2,6,7,4,5
for p in permutations([1, 2, 3, 4, 5]):
    c = (0,) + p + (8,)
    e = [abs(c[i] - c[(i + 1) % 7]) for i in range(7)]
    assert len(set(e)) < 7       # conclusion fails for all 120 arrangements with c_1 = 0, c_7 = 8
print("A228728 conjecture fails for {0,1,2,3,4,5,8}: antecedent", b, "distances", d, "; no arrangement with ends 0 and 8")
