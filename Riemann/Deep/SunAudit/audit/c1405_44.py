import sys
sys.path.insert(0, '.')
from lib import *
def phi(n):
    r = n; m = n; d = 2
    while d * d <= m:
        if m % d == 0:
            while m % d == 0: m //= d
            r -= r // d
        d += 1
    if m > 1: r -= r // m
    return r
def two_pr(p):
    if not is_prime(p) or p == 2: return False
    m = p - 1; d = 2; fs = []
    while d * d <= m:
        if m % d == 0:
            fs.append(d)
            while m % d == 0: m //= d
        d += 1
    if m > 1: fs.append(m)
    return all(pow(2, (p - 1) // q, p) != 1 for q in fs)
N = 3000
PH = [0] + [phi(i) for i in range(1, N + 1)]
bad_print = []; bad_oeis = []
for n in range(2, N + 1):
    if not any(m % 2 == 0 and two_pr(PH[n - m] + PH[m // 2] - 1) for m in range(1, n)): bad_print.append(n)
    if not any(two_pr(PH[k] + PH[n - k] // 2 - 1) for k in range(1, n)): bad_oeis.append(n)
print('printed phi(m/2): n without rep:', bad_print)
print('OEIS phi(m)/2: n without rep:', bad_oeis)
# n=8 detail, printed reading
print([(m, 8 - m, PH[8 - m] + PH[m // 2] - 1) for m in range(2, 8, 2)])
