from fractions import Fraction as F
from math import isqrt
a = F('0.37338196661320752')
def t0(q):
    if q % 2 == 0: return max(F(10**8, q), F(75*10**6, q) + 200)
    return max(F(10**8, q), F(375*10**5, q) + 200)
def D(q):
    r = a / t0(q); v = (1 + r*r)**2 / (4*r*r); return v.numerator // v.denominator
def sqfree(m):
    if m < 1: return False
    k = 2
    while k*k <= m:
        if m % (k*k) == 0: return False
        k += 1
    return True
N = 400000
# squarefree sieve
sf = bytearray([1])*(N+1); sf[0] = 0
for k in range(2, isqrt(N)+1):
    for j in range(k*k, N+1, k*k): sf[j] = 0
mods = [q for q in range(3, N+1) if (q % 2 == 1 and sf[q]) or (q % 4 == 0 and q % 8 != 0 and sf[q//4] and (q//4) % 4 == 3) or (q % 8 == 0 and q % 16 != 0 and sf[q//8] and (q//8) % 2 == 1)]
print('real primitive moduli <= 400000:', len(mods), 'largest', mods[-1])
worst = min(mods, key=D); print('min D over all', worst, D(worst), 't0', float(t0(worst)))
w2 = min([q for q in mods if q <= 125000], key=D); print('min D over q<=125000', w2, D(w2))
for q in (3, 4, 5, 8, 125003, 200003, 312497, 399997, mods[-1]):
    if q in set(mods[-50:]) or q in (3,4,5,8) or q in mods: print(q, D(q), 'old 1e8/q:', (lambda r: ((1+r*r)**2/(4*r*r)).numerator//((1+r*r)**2/(4*r*r)).denominator)(a/F(10**8,q)))
