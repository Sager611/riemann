#!/usr/bin/env python3
# genus_orders.py -- numerical check of Lemma G (real ring class characters of ALL quadratic orders).
# For every discriminant Delta (negative: -NNEG <= Delta <= -3; positive non-square: 5 <= Delta <= NPOS), and every
# fundamental discriminant d (d = 1 included) with |d| <= 2|Delta|, test whether psi_d : a -> chi_d(N a) is trivial on all
# values n = N(gamma) > 0, gcd(n, Delta) = 1, of the principal form of discriminant Delta (n <= B), i.e. whether psi_d is a
# (narrow, if Delta > 0) ring class character of O_Delta.  Lemma G predicts:
#   (i)  the number of passing d equals 2 |G_Delta[2]| = 2^mu  (each real character <-> exactly two d: d and fund(d*D_K)),
#   (ii) every passing d has odd part dividing Delta, and 2 | d implies 4 | Delta,
#   (iii) for each pair {d, d'} = {d, fund(d D_K)}: |d d'| <= |Delta| and d d' divides Delta (Lemma G(ii): |d d'| = |d_K| N(conductor)).
# |G[2]| for Delta < 0 is counted directly (ambiguous reduced forms); 2^mu is Gauss's genus count (both signs).
import math, signal, sys
signal.alarm(900)
NNEG = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
NPOS = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
MUT = sys.argv[3] if len(sys.argv) > 3 else ''  # mutation 'wide': for Delta > 0 also use |n| for norms n < 0 (wide classes); should FAIL somewhere

def isqf(m):
    m = abs(m)
    p = 2
    while p * p <= m:
        if m % (p * p) == 0: return False
        p += 1
    return True

def is_fund(d):
    if d == 1: return True
    if d % 4 == 1 or d % 4 == -3:
        return d % 4 in (1, -3) and isqf(d)
    if d % 4 == 0:
        m = d // 4
        return (m % 4 in (2, 3)) and isqf(m)
    return False

def fund_part(n):  # fundamental discriminant d with n in d * squares (n != 0)
    s = -1 if n < 0 else 1
    m = abs(n); core = 1; p = 2
    while p * p <= m:
        while m % (p * p) == 0: m //= p * p
        if m % p == 0: core *= p; m //= p
        p += 1
    core *= m
    c = s * core
    return c if c % 4 == 1 else 4 * c

def kron(a, n):  # Kronecker symbol (a/n), n > 0
    if n == 0: return 1 if abs(a) == 1 else 0
    res = 1
    while n % 2 == 0:
        n //= 2
        if a % 2 == 0: return 0
        if a % 8 in (3, 5): res = -res
    a %= n
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: res = -res
        a %= n
    return res if n == 1 else 0

def conductor(D):
    f = 1; Dc = D; p = 2
    while p * p <= abs(Dc):
        while Dc % (p * p) == 0 and (Dc // (p * p)) % 4 in (0, 1): Dc //= p * p; f *= p
        p += 1
    return Dc, f

def mu_genus(D):
    r = len([p for p in range(3, abs(D) + 1, 2) if abs(D) % p == 0 and all(p % q for q in range(3, int(p ** .5) + 1, 2))])
    if D % 4 != 0: return r
    m = D // 4
    if m % 4 == 1: return r
    if m % 4 in (2, 3): return r + 1
    if m % 8 == 4: return r + 1
    return r + 2

def ambiguous_count(D):  # D < 0: number of classes of order <= 2 (reduced forms with b = 0, b = a or a = c)
    A = -D; cnt = 0; a = 1
    while 3 * a * a <= A:
        for b in range(-a + 1, a + 1):
            num = b * b + A
            if num % (4 * a): continue
            c = num // (4 * a)
            if c < a or (b < 0 and a == c): continue
            if math.gcd(math.gcd(a, b), c) != 1: continue
            if b == 0 or b == a or a == c: cnt += 1
        a += 1
    return cnt

def principal_values(D, B):
    b0 = D % 2; c0 = (b0 - D) // 4  # x^2 + b0 x y + c0 y^2, discriminant b0^2 - 4 c0 = D
    vals = set()
    if D < 0:
        Y = int(2 * math.sqrt(B / -D)) + 1
        for y in range(0, Y + 1):
            r = B + D * y * y / 4.0
            if r < 0: continue
            R = math.sqrt(r) + 1
            for x in range(int(-b0 * y / 2 - R) - 1, int(-b0 * y / 2 + R) + 2):
                n = x * x + b0 * x * y + c0 * y * y
                if 0 < n <= B and math.gcd(n, D) == 1: vals.add(n)
    else:
        for y in range(0, 400):
            lo = math.sqrt(D) * y / 2.0; hi = math.sqrt(B + D * y * y / 4.0)
            for sgn in (1, -1):
                for u in range(int(lo) - 1, int(hi) + 2):  # u ~ |x + b0 y/2|
                    for x in (sgn * u - (b0 * y) // 2, sgn * u - (b0 * y) // 2 - 1, sgn * u - (b0 * y) // 2 + 1):
                        n = x * x + b0 * x * y + c0 * y * y
                        if MUT == 'wide' and -B <= n < 0 and math.gcd(n, D) == 1: vals.add(-n)
                        if 0 < n <= B and math.gcd(n, D) == 1: vals.add(n)
    return sorted(vals)

fund_cache = [d for d in range(-4 * max(NNEG, NPOS) - 8, 4 * max(NNEG, NPOS) + 9) if d != 0 and is_fund(d)]

def run(D):
    DK, f = conductor(D)
    B = max(4000, 30 * abs(D))
    vals = principal_values(D, B)
    vals = vals[:600]
    cands = [d for d in fund_cache if abs(d) <= 2 * abs(D)]
    passing = []
    for d in cands:
        ok = True
        for n in vals:
            if kron(d, n) != 1: ok = False; break
        if ok: passing.append(d)
    mu = mu_genus(D)
    bad = []
    if len(passing) != 2 ** mu: bad.append('count %d vs 2^mu=%d' % (len(passing), 2 ** mu))
    if D < 0:
        g2 = ambiguous_count(D)
        if 2 * g2 != 2 ** mu: bad.append('ambig %d vs 2^(mu-1)=%d' % (g2, 2 ** (mu - 1)))
    for d in passing:
        dodd = abs(d)
        while dodd % 2 == 0: dodd //= 2
        if D % dodd != 0: bad.append('odd part of %d does not divide' % d)
        if d % 2 == 0 and D % 4 != 0: bad.append('even %d, 4 does not divide' % d)
        dp = fund_part(d * DK) if d != DK else 1
        if d == 1: dp = DK
        if dp not in passing: bad.append('partner of %d (%d) missing' % (d, dp))
        if abs(d * dp) > abs(D): bad.append('|d d\'| = %d > |D| for d=%d' % (abs(d * dp), d))
        if D % (d * dp) != 0: bad.append('d d\' = %d does not divide D for d=%d' % (d * dp, d))
    return DK, f, len(vals), passing, mu, bad

tot = 0; nonfund = 0; fails = 0; maxpairs = 0
for D in list(range(-3, -NNEG - 1, -1)) + list(range(5, NPOS + 1)):
    if D % 4 not in (0, 1): continue
    if D > 0 and int(math.isqrt(D)) ** 2 == D: continue
    DK, f, nv, passing, mu, bad = run(D)
    tot += 1; nonfund += (f > 1); maxpairs = max(maxpairs, len(passing))
    if bad:
        fails += 1
        if fails <= 30: print('FAIL', D, 'DK', DK, 'f', f, 'nvals', nv, 'passing', passing, 'mu', mu, bad)
print('tested %d discriminants (%d non-fundamental), range [-%d,-3] and [5,%d] non-square; failures: %d; max #passing d: %d'
      % (tot, nonfund, NNEG, NPOS, fails, maxpairs))
