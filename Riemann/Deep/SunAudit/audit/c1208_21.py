import sys
from decimal import Decimal as D, getcontext
getcontext().prec = 60
sys.path.insert(0, '.')
from lib import *
P = primes_upto(200000)
def diff(alpha, n):
    a = D(alpha)
    lhs = sum(D(p) ** a for p in P[:n]) / n
    rhs = D(P[n - 1]) ** a / (a + 1)
    return lhs - rhs   # conjecture requires < 0
import math
def thr(alpha):
    c = math.ceil(alpha)
    return 2 * c * c + c + 6
for alpha in ['0.0001', '0.001', '0.002', '0.0023', '0.0024', '0.003', '0.01', '0.1', '0.5', '1']:
    n = thr(float(alpha))
    print('alpha', alpha, 'n', n, 'p_n', P[n-1], 'LHS-RHS =', '%.6e' % diff(alpha, n))
# bisection for the critical alpha at n = 9
lo, hi = D('0.0001'), D('0.01')
for _ in range(60):
    mid = (lo + hi) / 2
    if diff(mid, 9) > 0: lo = mid
    else: hi = mid
print('n=9 fails for 0 < alpha <', lo)
# scan alpha grid for other failures at n in [thr, thr+400]
getcontext().prec = 30
fails = []
import itertools
grid = [k / 1000 for k in range(1, 50)] + [k / 100 for k in range(5, 601)]
for al in grid:
    n0 = thr(al)
    a = D(al)
    S = sum(D(p) ** a for p in P[:n0 - 1])
    for n in range(n0, n0 + 400):
        S += D(P[n - 1]) ** a
        if S / n >= D(P[n - 1]) ** a / (a + 1):
            fails.append((al, n))
fl = {}
for al, n in fails: fl.setdefault(n, []).append(al)
print('grid failures (n: alphas):', {n: (min(v), max(v), len(v)) for n, v in fl.items()})
