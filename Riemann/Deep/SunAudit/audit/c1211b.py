import sys
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from lib import *
import math
P = primes_upto(200000)
# #72 (Conj 4.2(i)) m=1: sums_{i=j}^k 1/(p_i - 1) pairwise distinct; fractional-part coincidences with j,s>=2 (or k=1) only the 3 listed
K = 120
def segs(f, K, jmin_rule):
    pref = [Fr(0)]
    for i in range(1, K + 1): pref.append(pref[-1] + f(i))
    out = {}
    for k in range(1, K + 1):
        for j in range(1, k + 1):
            if jmin_rule and j < min(2, k): continue
            v = pref[k] - pref[j - 1]
            out.setdefault(v - (v.numerator // v.denominator), []).append((j, k, v))
    return out
for m in (1, 2):
    d = segs(lambda i: Fr(1, (P[i-1] - 1)**m), K, True)
    coinc = [v for v in d.values() if len(v) > 1]
    print('4.2(i) m=%d K=%d coincidences of fractional parts:' % (m, K), [[(a, b) for a, b, _ in c] for c in coinc])
    # pairwise distinct values (all j>=1)
    d2 = {}
    pref = [Fr(0)]
    for i in range(1, K + 1): pref.append(pref[-1] + Fr(1, (P[i-1]-1)**m))
    vals = [pref[k] - pref[j-1] for k in range(1, K + 1) for j in range(1, k + 1)]
    print('   all values distinct:', len(vals) == len(set(vals)))
for m in (1, 2):
    d = segs(lambda i: Fr(1, (P[i-1] + 1)**m), K, False)
    coinc = [v for v in d.values() if len(v) > 1]
    print('4.2(ii) m=%d coincidences:' % m, [[(a, b) for a, b, _ in c] for c in coinc])
for dd in (2, 3, 4, 5, 6):
    d = segs(lambda i: Fr(1, P[i-1] + dd), 80, False)
    coinc = [v for v in d.values() if len(v) > 1]
    print('4.2(iii) d=%d coincidences:' % dd, [[(a, b) for a, b, _ in c] for c in coinc][:5])
# #73 (Conj 4.3(i)) harmonic segments
d = segs(lambda i: Fr(1, i), 150, True)
coinc = [v for v in d.values() if len(v) > 1]
print('4.3(i) harmonic coincidences (K=150):', [[(a, b) for a, b, _ in c] for c in coinc])
# #100 (Conj 4.30(iii)) higher-order primes
PP = primes_upto(3 * 10**6)
def pm(m, n):
    x = n
    for _ in range(m): x = PP[x - 1]
    return x
bad = []
for m in (1, 2, 3):
    for n in range(1, 400):
        try:
            a1, a2 = pm(m, n), pm(m, n + 1)
            b1, b2 = pm(m + 1, n), pm(m + 1, n + 1)
        except IndexError:
            break
        r_m = math.log(a2) / (n + 1) - math.log(a1) / n
        r_m1 = math.log(b2) / (n + 1) - math.log(b1) / n
        if not (r_m1 < r_m < 0): bad.append((m, n, r_m1, r_m))
print('4.30(iii) failures (m<=3, n<400):', bad[:10])
