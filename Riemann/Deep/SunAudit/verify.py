#!/usr/bin/env python3
"""sunaudit verify.py -- stdlib-only re-check of every counterexample in REPORT.md.
Run: python3 verify.py   (about 2 s, single thread, < 50 MB). Each check prints PASS/FAIL."""
import itertools
from decimal import Decimal as D, getcontext

def isprime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0: return False
        d += 2
    return True

def primes(k):
    out, x = [], 2
    while len(out) < k:
        if isprime(x): out.append(x)
        x += 1
    return out

def det(M):  # Leibniz, exact integers (small matrices only)
    n = len(M); s = 0
    for perm in itertools.permutations(range(n)):
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j])
        t = 1
        for i in range(n): t *= M[i][perm[i]]
        s += (-1) ** inv * t
    return s

def leg(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

results = []
def check(name, ok, info):
    results.append(ok)
    print(('PASS ' if ok else 'FAIL ') + name + ': ' + info)

# 1. arXiv 1811.10503, Conj 4.9: |S(A)| >= min{p(F), (|A|-n)n(n+1)/2 + n(n^2-1)/6 + 1}
def S(A, n, p=None):
    return {(sum((k + 1) * a for k, a in enumerate(t)) % p) if p else sum((k + 1) * a for k, a in enumerate(t))
            for t in itertools.permutations(A, n)}
def bound(sz, n, p=None):
    b = (sz - n) * n * (n + 1) // 2 + n * (n * n - 1) // 6 + 1
    return min(p, b) if p else b
for A, n, p in [((0, 1, 2), 2, None), ((0, 1, 2), 2, 5), ((0, 1, 2, 5), 3, 11), ((0, 1, 2, 5), 4, 11), ((0, 1, 3, 4), 4, 13)]:
    s = S(A, n, p)
    check('1811.10503 C4.9 %s n=%d A=%s' % ('Q' if p is None else 'F_%d' % p, n, A), len(s) < bound(len(A), n, p),
          '|S(A)|=%d < bound %d, S(A)=%s' % (len(s), bound(len(A), n, p), sorted(s)))

# 2. arXiv 1208.2683, Conj 2.1 at alpha = 1/1000, n = 9 (threshold 2*1^2+1+6 = 9)
getcontext().prec = 50
a = D(1) / D(1000)
P9 = primes(9)
lhs = sum(D(q) ** a for q in P9) / 9
rhs = D(P9[-1]) ** a / (1 + a)
check('1208.2683 C2.1 alpha=0.001 n=9', lhs > rhs, 'LHS-RHS = %.6E (must be < 0 by the conjecture)' % (lhs - rhs))
prodP = 1
for q in P9: prodP *= q
check('1208.2683 C2.1 first-order sign', D(prodP) * D(9).exp() > D(23) ** 9,
      'theta(23)/9 > log 23 - 1  <=>  223092870*e^9 = %.6E > 23^9 = %d' % (D(prodP) * D(9).exp(), 23 ** 9))

# 3. arXiv 1208.2683, Conj 3.1 (u_0=0,u_1=1,u_{n+1}=A u_n - B u_{n-1}): (A,B)=(3,-26), A>2, A^2>4B
A_, B_ = 3, -26
u = [0, 1]
for k in range(1, 8): u.append(A_ * u[k] - B_ * u[k - 1])
# r_5 > r_6 required;  r_5 <= r_6  <=>  u6^(2*5*7) <= u5^(6*7) * u7^(5*6)
check('1208.2683 C3.1 (A,B)=(3,-26) n=5', u[6] ** 70 < u[5] ** 42 * u[7] ** 30,
      'u5,u6,u7=%d,%d,%d; u6^70 < u5^42*u7^30 so r_5 < r_6 (ratio not decreasing)' % (u[5], u[6], u[7]))

# 4. arXiv 1211.1588, Conj 4.30(iii) at m=1, n=2: p_5^(1/3)/p_3^(1/2) < p_3^(1/3)/p_2^(1/2) required
check('1211.1588 C4.30(iii) m=1 n=2', 11 ** 2 * 3 ** 3 > 5 ** 5,
      '11^(1/3)/5^(1/2) > 5^(1/3)/3^(1/2)  <=>  11^2*3^3 = %d > 5^5 = %d' % (11 ** 2 * 3 ** 3, 5 ** 5))

# 5. arXiv 1809.07766, Conj 6.8 (6.15): det[floor(i^2 j^2/n)] != 0 iff n = 9 or prime > 7, n = 3 mod 4
M7 = [[(i * i * j * j) // 7 for j in range(1, 4)] for i in range(1, 4)]
check('1809.07766 C6.8(6.15) n=7', det(M7) != 0, 'matrix %s has det %d != 0, but 7 is not 9 nor a prime > 7' % (M7, det(M7)))

# 6. arXiv 1402.6641, Conj 2.25: count nonzero prime ideals of Z[i] (one generator a+bi, a>0, b>=0, per ideal)
def gauss_ideals(x):
    c = 0
    for a_ in range(1, int(x ** 0.5) + 2):
        for b_ in range(0, int(x ** 0.5) + 2):
            N = a_ * a_ + b_ * b_
            if N > x: continue
            if b_ == 0:
                if isprime(a_) and a_ % 4 == 3: c += 1
            elif isprime(N): c += 1
    return c
cnts = [gauss_ideals(5 * k) for k in range(1, 6)]
check('1402.6641 C2.25 n=5', not any(isprime(c) and c % 4 == 1 for c in cnts),
      'ideal counts with norm <= 5k, k=1..5: %s (none is a prime = 1 mod 4); the paper formula gives %s' % (cnts, [2, 3, 4, 5, 5]))

# 7. arXiv 2405.03626, Conj 3.9 (p = 3 mod 4): det[x+((j+k)/p)+((j-k)/p)+(j/p)y+(k/p)z+(jk/p)w]_{0<=j,k<=(p-3)/2}
p = 7
def M39(x, y, z, w):
    return [[x + leg(j + k, p) + leg(j - k, p) + leg(j, p) * y + leg(k, p) * z + leg(j * k, p) * w for k in range(3)] for j in range(3)]
d = det(M39(0, 0, 1, 0))
check('2405.03626 C3.9 p=7 at (x,y,z,w)=(0,0,1,0)', d != 0, 'det = %d, formula (-1)^((h+1)/2) p^((p-7)/4) (p-2 m_p w) x gives 0' % d)

# 8. arXiv 1405.0290, Conj 4.4(ii) as printed: p = phi(k) + phi(m/2) - 1 prime with 2 a primitive root, n = k+m = 8
def phi(n):
    return sum(1 for t in range(1, n + 1) if __import__('math').gcd(t, n) == 1)
def two_primroot(q):
    if not isprime(q) or q == 2: return False
    return all(pow(2, (q - 1) // r, q) != 1 for r in range(2, q) if (q - 1) % r == 0 and isprime(r))
vals = [phi(8 - m) + phi(m // 2) - 1 for m in range(2, 8, 2)]
check('1405.0290 C4.4(ii) printed, n=8', not any(two_primroot(v) for v in vals), 'm even -> p in %s; p=2 has no primitive root 2' % vals)

# 9. arXiv 1211.1588, Conj 2.6: 2C contains {4006, 4007, ...}, C = {x>0: 2x-3, 2x+3 prime}
C = [x for x in range(1, 4006) if isprime(2 * x - 3) and isprime(2 * x + 3)]
Cs = set(C)
check('1211.1588 C2.6 2C at 4006', not any((4006 - x) in Cs for x in C), '4006 is not a sum of two elements of C')

# Leads (expected MISREAD): the scouts' claimed failures do not occur
def s_t(p):
    n = (p - 1) // 2
    sq = [k * k % p for k in range(1, n + 1)]
    s = sum(1 for i in range(n) for j in range(i + 1, n) if sq[i] > sq[j])
    t = sum(1 for j in range(1, n + 1) for k in range(j + 1, n + 1) if 2 * ((k * k - j * j) % p) > p)
    return s, t
ok = all((sum(s_t(p)) - sum(1 for k in range(1, p) if 4 * k < p and leg(k, p) == 1)) % 2 == 0 for p in (5, 13, 17, 29))
check('lead 1809.07766 C6.1 holds at p=5,13,17,29 (MISREAD)', ok, 'p=11 is not 1 mod 4')
def c8(p):
    n = (p - 1) // 2
    v = [pow(j, 8, p) for j in range(1, n + 1)]
    return sum(1 for i in range(n) for j in range(i + 1, n) if v[i] > v[j])
ok = all((c8(p) - sum(1 for k in range(1, p) if 4 * k < p and leg(k, p) == 1)) % 2 == 0 for p in (17, 41, 73))
check('lead 1809.07766 C6.6 (6.12) holds at p=17,41,73 (MISREAD)', ok, 'p = 1 mod 8 case')
P = primes(10)
check('lead 1402.6641 C3.5(i): 6 is in the set (MISREAD)', 2 * P[3] - P[1] == P[4], 'p_2, p_4, p_5 = 3, 7, 11 form an AP and 2+4 = 6')

print('ALL PASS' if all(results) else 'SOME FAIL')
