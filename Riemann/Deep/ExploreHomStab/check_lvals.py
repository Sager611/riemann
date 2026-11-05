# explore9: check lvals.py against mpmath.dirichlet for a few small d.
import sys
import numpy as np, mpmath as mp
mp.mp.dps = 20

def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0

def chi8d(d, n):
    if n % 2 == 0: return 0
    two = 1 if n % 8 in (1, 7) else -1
    return two * jacobi(d, n)

z = np.load(sys.argv[1]); dd, LL = z['d'], z['L']
for d in [1, 3, 5, 7, 11, 101, 251, 1001, 1999]:
    k = np.searchsorted(dd, d)
    if k >= len(dd) or dd[k] != d: continue
    q = 8 * d
    chi = [chi8d(d, n) for n in range(q)]
    ref = mp.dirichlet(0.5, chi)
    print(d, LL[k], mp.nstr(ref, 15), float(LL[k] - ref))
