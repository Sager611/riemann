#!/usr/bin/env python3
"""f2_class_bound.py -- exact arithmetic behind 'an F2-linear class has at most 24 lines' (REPORT.md section 5).

L/2L = F_2^24 with the 98280 classes T2 of minimal vectors (one line each). For y != 0,
fhat(y) = sum_{t in T2} (-1)^{b(t,y)} is 4072 if y in T2 and -24 otherwise (two-valued; checked here by Parseval
from the counts).  Let V be T2-free, dim d, W = V^perp (dim 24-d), a = |T2 cap W|.  T2 cap V = 0 forces
a = 3*2^(15-d) - 24, and every coset x+V (x not in V) then holds 2^(d-24) * (98304 + 4096 s) lines with |s| <= a,
which is at most 24 for every d.  So a coset of a minimum-raising sublattice 2L <= M <= L carries at most one frame.
"""
from fractions import Fraction as Fr
N2 = 98280
lhs = N2**2 + N2 * 4072**2 + (2**24 - 1 - N2) * 24**2
assert lhs == 2**24 * N2, 'Parseval'
for d in range(0, 13):
    a = Fr(3 * 2**(15 - d) - 24)
    if a < 0: print('d=%d impossible' % d); continue
    count_max = Fr(2)**(d - 24) * (98304 + 4096 * a)
    print('dim V = %2d : |T2 cap V^perp| = %s, max lines in a coset = %s' % (d, a, count_max))
    assert count_max == 24
print('dim V = 13 :', 3 * 2**(15 - 13) - 24, '< 0, impossible')
