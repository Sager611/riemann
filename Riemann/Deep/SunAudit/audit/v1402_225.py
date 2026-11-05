import sys
sys.path.insert(0, '.')
from lib import *
from math import isqrt
P = primes_upto(400000)
def true_count(x):
    c = 0
    for p in P:
        if p > x: break
        if p == 2: c += 1
        elif p % 4 == 1: c += 2
        elif p * p <= x: c += 1
    return c
def sun_count(x):
    r = isqrt(x)
    return sum(1 for p in P if p <= r) + sum(1 for p in P if r < p <= x and p % 4 != 3)
# brute-force true count by enumerating Gaussian primes up to associates for small x
def brute(x):
    # count ideals (a+bi) with a>0,b>=0 (one generator per associate class), norm<=x, prime
    cnt = 0
    for a in range(1, isqrt(x) + 1):
        for b in range(0, isqrt(x) + 1):
            N = a*a + b*b
            if N > x: continue
            if b == 0:
                if is_prime(a) and a % 4 == 3: cnt += 1
            elif is_prime(N): cnt += 1
    return cnt
for x in [2, 5, 9, 10, 25, 50, 100, 1000]:
    print(x, 'true', true_count(x), 'brute', brute(x), 'sun-formula', sun_count(x))
def ok(n, cnt):
    return any(is_prime(cnt(k*n)) and cnt(k*n) % 4 == 1 for k in range(1, n + 1))
print('literal (true count) failures n in [5,300]:', [n for n in range(5, 301) if not ok(n, true_count)])
print('Sun-formula failures n in [5,300]:', [n for n in range(5, 301) if not ok(n, sun_count)])
print('n=5 true counts at kn:', [true_count(5*k) for k in range(1, 6)], 'sun:', [sun_count(5*k) for k in range(1, 6)])
