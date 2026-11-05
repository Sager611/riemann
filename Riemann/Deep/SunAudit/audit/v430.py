import sys, math
sys.path.insert(0, '.')
from lib import *
PP = primes_upto(6 * 10**6)
def pm(m, n):
    x = n
    for _ in range(m): x = PP[x - 1]
    return x
bad = []
for m in (1, 2, 3, 4):
    n = 1
    while True:
        try:
            a1, a2, b1, b2 = pm(m, n), pm(m, n + 1), pm(m + 1, n), pm(m + 1, n + 1)
        except IndexError:
            break
        r_m = math.log(a2) / (n + 1) - math.log(a1) / n
        r_m1 = math.log(b2) / (n + 1) - math.log(b1) / n
        if not (r_m1 < r_m < 0): bad.append((m, n))
        n += 1
    print('m=%d checked n<%d' % (m, n), flush=True)
print('4.30(iii) failures:', bad)
print('exact: 11^2*3^3 =', 11**2 * 3**3, ' 5^5 =', 5**5)
