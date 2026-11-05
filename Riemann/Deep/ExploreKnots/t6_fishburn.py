# t6_fishburn.py — sign patterns (node N1): coefficients of the Habiro elements at q = 1 - x for the trefoil
# (Kontsevich–Zagier F: Fishburn numbers) and for Hikami's T(2,5) series, whose L-series is of Davenport–Heilbronn type.
D = 30
def mul(a, b):
    c = [0]*(D+1)
    for i, x in enumerate(a):
        if x:
            for j in range(D+1-i):
                if b[j]: c[i+j] += x*b[j]
    return c
def binom_poly(j):  # (1-x)^j truncated
    from math import comb
    return [comb(j, k)*(-1)**k for k in range(min(j, D)+1)] + [0]*(D-min(j, D))
def one_minus_qj(j):  # 1 - (1-x)^j, valuation 1
    p = binom_poly(j); return [1-p[0]] + [-c for c in p[1:]]
ONE = [1] + [0]*D
poch = [ONE]
for k in range(1, D+1): poch.append(mul(poch[-1], one_minus_qj(k)))
def qbinom(n, k):  # Gaussian binomial as polynomial in q, then substitute q = 1-x (exact, integer)
    # build in q first: coefficients of [n choose k]_q
    from functools import lru_cache
    @lru_cache(None)
    def gb(n, k):
        if k < 0 or k > n: return (0,)
        if k == 0 or k == n: return (1,)
        a = gb(n-1, k-1); b = gb(n-1, k)  # [n,k] = [n-1,k-1] + q^k [n-1,k]
        r = [0]*(k*(n-k)+1)
        for i, c in enumerate(a): r[i] += c
        for i, c in enumerate(b): r[i+k] += c
        return tuple(r)
    poly = gb(n, k); out = [0]*(D+1)
    for e, c in enumerate(poly):
        if c:
            p = binom_poly(e)
            for i in range(D+1): out[i] += c*p[i]
    return out
def qpow(e): return binom_poly(e)
F1 = [0]*(D+1)
for k in range(D+1): F1 = [x+y for x, y in zip(F1, poch[k])]
X = [0]*(D+1)
for k2 in range(D+1):
    for k1 in range(k2+1):
        t = mul(mul(poch[k2], qpow(k1*k1+k1)), qbinom(k2, k1))
        X = [x+y for x, y in zip(X, t)]
print('trefoil F(1-x):', F1[:14])
print('  all positive up to x^%d:' % D, all(c > 0 for c in F1))
print('T(2,5) X(1-x):', X[:14])
print('  all positive up to x^%d:' % D, all(c > 0 for c in X), ' signs:', ''.join('+' if c > 0 else ('-' if c < 0 else '0') for c in X))
