#!/usr/bin/env python3
"""Counterexamples to (printed) Conjectures 3.2, 3.3(ii), 3.4(i) (TeX labels Conj3.3, Conj3.4, Conj3.5) of Z.-W. Sun, "Some new problems in additive combinatorics",
arXiv 1309.1679 (v9, 2020; TeX l. 512-516, 525-534, 541-555). Exhaustive searches in cyclic groups Z/m. Runs in a few seconds.

3.2   G finite abelian, n distinct elements a_1..a_n, n does not divide |G| (or n even and the Sylow 2-subgroup cyclic):
      some permutation (b_1..b_n) with b_1 = a_1 has the n-1 differences b_i - b_(i+1) pairwise distinct.
3.3(ii) n > 2, A an n-subset of G with |G| odd and n not dividing |G|: some numbering has the n circular differences
      a_1-a_2, ..., a_n-a_1 pairwise distinct.
3.4(i) |A| = n > 3, G finite with |G| not divisible by 3: some numbering has the n circular sums a_i + 2a_(i+1)
      (indices mod n) pairwise distinct."""
from itertools import permutations, combinations

def path_diff_ok(m, A, a1):
    rest = [x for x in A if x != a1]
    for p in permutations(rest):
        b = (a1,) + p
        d = [(b[i] - b[i + 1]) % m for i in range(len(b) - 1)]
        if len(set(d)) == len(d):
            return b
    return None

def circ_ok(m, A, f):
    A = list(A)
    for p in permutations(A[1:]):
        b = (A[0],) + p
        s = [f(b[i], b[(i + 1) % len(b)]) % m for i in range(len(b))]
        if len(set(s)) == len(s):
            return b
    return None

# 3.2: Z/6 (4 does not divide 6; Sylow 2-subgroup Z/2 cyclic), a_1 = 0, others 1, 3, 5; and Z/9 with 0, 1, 4, 7
assert path_diff_ok(6, [0, 1, 3, 5], 0) is None
assert path_diff_ok(9, [0, 1, 4, 7], 0) is None
print("3.2 fails: G = Z/6, (a1..a4) = (0,1,3,5); G = Z/9, (0,1,4,7)")

# 3.3(ii): Z/5 (odd order, 3 does not divide 5), A = {0,1,2}; every 3-term progression fails when 3 does not divide |G|
assert circ_ok(5, [0, 1, 2], lambda x, y: x - y) is None
for m in (5, 7, 11, 13):
    for d in range(1, m):
        assert circ_ok(m, [0, d, (2 * d) % m], lambda x, y: x - y) is None
print("3.3(ii) fails: G = Z/5, A = {0,1,2}; every 3-term progression in Z/5, Z/7, Z/11, Z/13")

# 3.4(i): Z/7 (7 not divisible by 3), n = 5, A = {0,1,2,3,4}
assert circ_ok(7, [0, 1, 2, 3, 4], lambda x, y: x + 2 * y) is None
fails = {}
for m in (4, 5, 7, 8, 10, 11):
    for n in range(4, m + 1):
        c = sum(1 for A in combinations(range(m), n) if circ_ok(m, A, lambda x, y: x + 2 * y) is None)
        if c:
            fails[(m, n)] = c
print("3.4(i) fails: G = Z/7, A = {0,1,2,3,4}; failing n-subsets by (m, n):", fails)
