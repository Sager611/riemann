import sys
sys.path.insert(0, '.')
from lib import *
from c1809b_det import det_bareiss
from math import gcd, isqrt
N = int(sys.argv[1])
def phi(n):
    r = n; m = n; d = 2
    while d * d <= m:
        if m % d == 0:
            while m % d == 0: m //= d
            r -= r // d
        d += 1
    if m > 1: r -= r // m
    return r
def sum2sq(n):
    a = 0
    while a * a <= n:
        b = isqrt(n - a * a)
        if a * a + b * b == n: return True
        a += 1
    return False
def cd(c, d, n):
    rng = range(2, n - 1)
    return det_bareiss([[jacobi(j*j + c*j*k + d*k*k, n) for k in rng] for j in rng])
def issq(x): return x >= 0 and isqrt(x) ** 2 == x
fails = []
for n in range(3, N, 2):
    if n % 4 == 1 and not sum2sq(n):
        if cd(3, 2, n) != 0: fails.append(('5.1i', n))
    if n % 4 == 3:
        v = cd(3, 2, n)
        if v % (phi(n) // 2): fails.append(('5.1ii', n, v))
        if n % 8 == 3:
            q = v * 2 // phi(n) if v % (phi(n)//2) == 0 else None
            if q is None or not issq(q): fails.append(('5.1iii', n, v))
    if is_prime(n) and n % 24 in (13, 19) and cd(2, 2, n) != 0: fails.append(('5.2i', n))
    if is_prime(n) and n % 24 in (17, 23) and cd(2, 2, n) % n: fails.append(('5.2ii', n))
    if n % 8 == 5 and (cd(4, 2, n) != 0 or cd(8, 8, n) != 0): fails.append(('5.3i', n))
    if n % 12 == 5 and cd(3, 3, n) != 0: fails.append(('5.3ii', n))
    if n % 4 == 1 and jacobi(n, 7) == -1 and (cd(42, -7, n) != 0 or cd(21, 112, n) != 0): fails.append(('5.4', n))
    if n > 3:
        v = cd(2, 3, n)
        if v % n: fails.append(('5.5i-a', n))
        if n % 12 not in (1, 11) and v % (n * n): fails.append(('5.5i-b', n, v))
    if n > 7 and cd(6, 15, n) % n: fails.append(('5.5ii', n))
    if n % 20 in (13, 17) and sum2sq(n) and cd(5, 5, n) != 0: fails.append(('5.6i', n))
    if n > 1 and (n % 20 in (11, 19) or (n % 60 == 9 and n > 69)):
        if jacobi(cd(5, 5, n), n) != 0: fails.append(('5.6ii', n))
    if n % 12 == 5 and sum2sq(n) and cd(10, 9, n) != 0: fails.append(('5.7i', n))
    if is_prime(n) and n % 12 == 11 and cd(10, 9, n) % n: fails.append(('5.7ii', n))
    if n % 24 in (13, 17) and sum2sq(n) and cd(8, 18, n) != 0: fails.append(('5.8i', n))
    if is_prime(n) and n % 24 == 19 and cd(8, 18, n) % (n * n): fails.append(('5.8ii', n))
    if is_prime(n) and n % 24 == 23 and cd(8, 18, n) % n: fails.append(('5.8iii', n))
print('Section 5 failures, odd n <', N, ':', fails)
