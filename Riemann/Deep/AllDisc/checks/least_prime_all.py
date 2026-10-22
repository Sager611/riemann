#!/usr/bin/env python3
# least_prime_all.py -- sanity check for the target statement on ALL discriminants (not a proof).
# For every non-square discriminant Delta with 5 <= |Delta| <= NMAX (both signs), and every class of primitive forms of
# discriminant Delta (Delta < 0: positive definite, reduced forms; Delta > 0: proper equivalence = cycles of reduced
# indefinite forms under rho, i.e. narrow classes), find the least prime p not dividing Delta represented by the class
# (p is represented by Q iff Q ~ (p, +-b, c) with b^2 = Delta mod 4p). Report P(Delta) := max over classes, and
# e(Delta) := log P(Delta) / log|Delta|, separately for fundamental and non-fundamental Delta.
import math, signal, sys
signal.alarm(1200)
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
PLIM = 4 * 10**6

sieve = bytearray([1]) * (PLIM + 1); sieve[0] = sieve[1] = 0
for i in range(2, int(PLIM ** 0.5) + 1):
    if sieve[i]: sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
PRIMES = [i for i in range(2, PLIM + 1) if sieve[i]]

def kron(a, n):
    if n == 2: return 0 if a % 2 == 0 else (1 if a % 8 in (1, 7) else -1)
    return pow(a % n, (n - 1) // 2, n) if a % n else 0

def sqrt_mod(a, p):  # p odd prime, a QR mod p
    a %= p
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p); m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r

def is_fund(d):
    def sqf(m):
        m = abs(m); p = 2
        while p * p <= m:
            if m % (p * p) == 0: return False
            p += 1
        return True
    if d % 4 == 1: return sqf(d)
    if d % 4 == 0: return (d // 4) % 4 in (2, 3) and sqf(d // 4)
    return False

def red_neg(a, b, c):  # positive definite
    while True:
        if a > c: a, b, c = c, -b, a
        if b > a or b <= -a:
            k = (a - b) // (2 * a); b2 = b + 2 * k * a; c = (b2 * b2 - D_) // (4 * a); b = b2; continue
        if a == c and b < 0: b = -b
        return (a, b, c)

def rr(b, a, D, sD):  # Cohen's r(b, a)
    A = abs(a); m = 2 * A
    if A > sD:
        r = b % m
        if r > A: r -= m
    else:
        lo = sD - m  # need sD - 2A < r < sD
        r = b % m
        k = math.floor((lo - r) / m) + 1
        r += k * m
        while r <= lo: r += m
        while r >= sD: r -= m
    return r

def is_red_pos(a, b, c, sD):
    return 0 < b < sD and sD - b < 2 * abs(a) < sD + b

def rho(a, b, c, D, sD):
    b2 = rr(-b, c, D, sD)
    return (c, b2, (b2 * b2 - D) // (4 * c))

def red_pos(a, b, c, D, sD):
    for _ in range(10000):
        if is_red_pos(a, b, c, sD): return (a, b, c)
        a, b, c = rho(a, b, c, D, sD)
    raise RuntimeError('no reduction')

def classes(D):
    """map reduced primitive form -> class id; returns (dict, h)"""
    cid = {}
    if D < 0:
        A = -D; h = 0; a = 1
        while 3 * a * a <= A:
            for b in range(-a + 1, a + 1):
                if (b * b + A) % (4 * a): continue
                c = (b * b + A) // (4 * a)
                if c < a or (b < 0 and a == c) or math.gcd(math.gcd(a, b), c) != 1: continue
                cid[(a, b, c)] = h; h += 1
            a += 1
        return cid, h
    sD = math.sqrt(D); reds = []
    for b in range(1, int(sD) + 1):
        if b >= sD or (b * b - D) % 4: continue
        N = (D - b * b) // 4
        for a in range(1, N + 1):
            if N % a: continue
            for sa in (a, -a):
                c = -N // sa
                if is_red_pos(sa, b, c, sD) and math.gcd(math.gcd(a, b), abs(c)) == 1: reds.append((sa, b, c))
    h = 0
    for f in reds:
        if f in cid: continue
        g = f
        while True:
            cid[g] = h
            g = rho(*g, D, sD)
            if g == f: break
            if g in cid and cid[g] != h: raise RuntimeError('cycle clash')
        h += 1
    return cid, h

# self-test of the class counts against known values (h for Delta < 0, narrow h+ for Delta > 0)
KNOWN = {-23: 3, -47: 5, -56: 4, -71: 7, -108: 3, -36: 2, -12: 1, -16: 1, -27: 1, -28: 1, -400: 4, -1984: 12,
         5: 1, 8: 1, 12: 2, 13: 1, 20: 1, 21: 2, 24: 2, 28: 2, 32: 2, 60: 4, 136: 4, 145: 4, 401: 5, 316: 6}
for D, hk in KNOWN.items():
    D_ = D
    h = classes(D)[1]
    print('self-test h(%d) = %d (expected %d) %s' % (D, h, hk, 'ok' if h == hk else 'MISMATCH'))

res = {}
for D in list(range(-5, -NMAX - 1, -1)) + list(range(5, NMAX + 1)):
    if D % 4 not in (0, 1): continue
    if D > 0 and math.isqrt(D) ** 2 == D: continue
    D_ = D
    cid, h = classes(D)
    sD = math.sqrt(abs(D))
    least = [0] * h; left = h; P = 0
    for p in PRIMES:
        if D % p == 0: continue
        if p == 2:
            if D % 8 != 1: continue
            b = 1
        else:
            if kron(D, p) != 1: continue
            r = sqrt_mod(D, p)
            b = r if (r - D) % 2 == 0 else r - p
        c = (b * b - D) // (4 * p)
        for bb in (b, -b):
            f = red_neg(p, bb, c) if D < 0 else red_pos(p, bb, c, D, sD)
            k = cid[f]
            if not least[k]: least[k] = p; left -= 1; P = p
        if left == 0: break
    if left: print('NOT ALL CLASSES HIT', D); continue
    res[D] = (h, P, is_fund(D))

for sgn in (-1, 1):
    for fund in (True, False):
        L = [(math.log(P) / math.log(abs(D)), D, h, P) for D, (h, P, fu) in res.items() if fu == fund and D * sgn > 0 and abs(D) >= 200]
        L.sort()
        if not L: continue
        n = len(L); med = L[n // 2][0]
        top = L[-1]
        rat = max(P / abs(D) for (e, D, h, P) in L)
        print('Delta %s 0, %-15s n=%4d  median e=%.3f  max e=%.3f at Delta=%d (h=%d, P=%d)  max P/|Delta|=%.2f'
              % ('<' if sgn < 0 else '>', 'fundamental' if fund else 'non-fundamental', n, med, top[0], top[1], top[2], top[3], rat))
big = [(D, h, P) for D, (h, P, fu) in res.items() if not fu and abs(D) >= 200]
big.sort(key=lambda t: -t[2] / abs(t[0]))
print('largest P/|Delta| among non-fundamental:', big[:8])
print('discriminants tested: %d (non-fundamental %d)' % (len(res), sum(1 for v in res.values() if not v[2])))
