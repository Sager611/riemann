"""Theorem D check: deterministic construction of a generator of the r-Sylow subgroup of k^x for an explicitly given
finite field k = F_p[t]/(g) containing a given primitive r-th root of unity zeta_r (DERIVATION section 4).
Construction (deterministic given zeta_r; only F_p-linear algebra, exponentiation, and one small prime):
  K_r = F_p(zeta_r) inside k, d0 = [K_r:F_p], D = [k:F_p], e = v_r(D/d0), t = v_r(|k|-1).
  (i)  e >= 1, and r odd or p^(d0 r^(e-1)) = 1 mod 4: solve x^(p^(d0 r^(e-1))) = zeta_r x (kernel of Frob^f - zeta_r),
       g = x^((|k|-1)/r^t).
  (ii) e = 0: auxiliary prime l (l = 1 mod 12r, p not an r-th power mod l), Gauss-period field U_r, V = U_r (x) k
       (a field: gcd(r, D) = 1), u in V with u^|k| = zeta_r u, w = u^((|k|^r-1)/r^(t+1)), g = w^r  (OpenAI's table step
       with K_q replaced by k).  r = 2 and |k| = 3 mod 4: g = -1.
  (iii) r = 2, p = 3 mod 4, e = 1 (k_1 = F_{p^2}): u with u^p = -u, i = u/c (c^2 = -u^2 in F_p), least prime
       l = 1 mod 4 with (l/p) = -1, l = a^2 + b^2, x = a + b i (norm l: a non-square of F_{p^2}), g = x^((|k|-1)/2^t).
The test harness picks g irreducible and zeta_r at random (allowed: they are the INPUT), runs the construction, and checks
that g has order exactly r^t.  Pure Python.  Run: python3 primary.py
"""
import random, sys
random.seed(12345)

def is_prime(n):
    if n < 2: return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0: return n == sp
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def vr(n, r):
    s = 0
    while n % r == 0: n //= r; s += 1
    return s

# ---------- polynomials over F_p (lists, low degree first) ----------
def trim(a):
    while a and a[-1] == 0: a.pop()
    return a
def pmul(a, b, p):
    if not a or not b: return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] = (r[i + j] + x * y) % p
    return trim(r)
def pmod(a, g, p):
    a = a[:]; dg = len(g) - 1; inv = pow(g[-1], p - 2, p)
    while len(a) - 1 >= dg and a:
        c = a[-1] * inv % p; sh = len(a) - 1 - dg
        for i in range(len(g)): a[sh + i] = (a[sh + i] - c * g[i]) % p
        trim(a)
    return a
def psub(a, b, p):
    n = max(len(a), len(b)); r = [((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p for i in range(n)]
    return trim(r)
def pgcd(a, b, p):
    while b: a, b = b, pmod(a, b, p)
    return a

class Field:
    def __init__(s, p, g):
        s.p, s.g, s.D = p, g, len(g) - 1; s.Q = p ** s.D
    def mul(s, a, b): return pmod(pmul(a, b, s.p), s.g, s.p)
    def pw(s, a, e):
        r = [1]; a = a[:]
        while e:
            if e & 1: r = s.mul(r, a)
            a = s.mul(a, a); e >>= 1
        return r
    def inv(s, a): return s.pw(a, s.Q - 2)
    def vec(s, a): return [(a[i] if i < len(a) else 0) for i in range(s.D)]
    def el(s, v): return trim([x % s.p for x in v])

def irreducible(p, D):
    while True:
        g = [random.randrange(p) for _ in range(D)] + [1]
        x = [0, 1]; ok = True; F = Field(p, g); y = x
        for i in range(1, D // 2 + 1):
            y = F.pw(y, p)
            if len(pgcd(g[:], psub(y, x, p), p)) > 1: ok = False; break
        if ok and g[0] != 0: return g

# ---------- linear algebra over F_p ----------
def kernel(M, p):
    """basis of {v : M v = 0}, M given as list of rows."""
    m = [row[:] for row in M]; R = len(m); C = len(m[0]); piv = []; r = 0
    for c in range(C):
        pr = next((i for i in range(r, R) if m[i][c] % p), None)
        if pr is None: continue
        m[r], m[pr] = m[pr], m[r]; iv = pow(m[r][c], p - 2, p); m[r] = [x * iv % p for x in m[r]]
        for i in range(R):
            if i != r and m[i][c] % p:
                f = m[i][c]; m[i] = [(m[i][k] - f * m[r][k]) % p for k in range(C)]
        piv.append(c); r += 1
        if r == R: break
    free = [c for c in range(C) if c not in piv]; basis = []
    for fc in free:
        v = [0] * C; v[fc] = 1
        for i, c in enumerate(piv): v[c] = (-m[i][fc]) % p
        basis.append(v)
    return basis

def frob_matrix(F, f):
    """matrix (columns) of x -> x^(p^f) on k in the basis 1, t, ..., t^(D-1)."""
    cols = [F.vec(F.pw([0] * j + [1], F.p ** f)) for j in range(F.D)]
    return cols

def mult_matrix(F, z):
    return [F.vec(F.mul(z, [0] * j + [1])) for j in range(F.D)]

def cols_to_rows(cols):
    return [list(r) for r in zip(*cols)]

def subfield_degree(F, z):
    """d0 = degree over F_p of F_p(z) = least f with z^(p^f) = z."""
    f, y = 1, F.pw(z, F.p)
    while y != z: y = F.pw(y, F.p); f += 1
    return f

# ---------- Gauss periods (as in table.py) ----------
def prime_factors(n):
    f, d = [], 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: f.append(n)
    return f
def structure_constants(l, q, p):
    g = next(g for g in range(2, l) if all(pow(g, (l - 1) // r, l) != 1 for r in prime_factors(l - 1)))
    cos = [[] for _ in range(q)]; x = 1
    for k in range(l - 1): cos[k % q].append(x); x = x * g % l
    idx = {e: j for j, c in enumerate(cos) for e in c}
    a = [[[0] * q for _ in range(q)] for _ in range(q)]
    for i in range(q):
        for j in range(q):
            cnt = [0] * q; ones = 0
            for x in cos[i]:
                for y in cos[j]:
                    e = (x + y) % l
                    if e == 0: ones += 1
                    else: cnt[idx[e]] += 1
            a[i][j] = [(cnt[k] // len(cos[k]) - ones) % p for k in range(q)]
    return a
def aux_prime(p, q):
    l = 12 * q + 1
    while not (l != p and is_prime(l) and pow(p, (l - 1) // q, l) != 1): l += 12 * q
    return l

def construct(F, r, zeta):
    """returns (g, case, extra) with g of order r^t in k^x."""
    p, D, Q = F.p, F.D, F.Q
    t = vr(Q - 1, r)
    d0 = subfield_degree(F, zeta)
    e = vr(D // d0, r)
    if r == 2 and Q % 4 == 3:
        return [p - 1], 'r=2,|k|=3 mod 4', None
    if e >= 1 and (r != 2 or pow(p, d0 * r ** (e - 1), 4) == 1):
        f = d0 * r ** (e - 1)
        A = frob_matrix(F, f); Z = mult_matrix(F, zeta)
        Mcols = [[(A[j][i] - Z[j][i]) % p for i in range(D)] for j in range(D)]
        ker = kernel(cols_to_rows(Mcols), p)
        x = F.el(ker[0])
        return F.pw(x, (Q - 1) // r ** t), 'eigenvector e=%d' % e, None
    if r == 2 and e >= 1:   # then p = 3 mod 4 and e = 1 (e >= 2 gives |k_{e-1}| = 1 mod 4)
        assert e == 1 and p % 4 == 3
        A = frob_matrix(F, 1)
        Mcols = [[(A[j][i] + (1 if i == j else 0)) % p for i in range(D)] for j in range(D)]   # x^p = -x
        u = F.el(kernel(cols_to_rows(Mcols), p)[0])
        u2 = F.mul(u, u); assert len(u2) <= 1
        m = (-u2[0]) % p if u2 else 0
        c = pow(m, (p + 1) // 4, p); assert c * c % p == m
        i_el = F.mul(u, [pow(c, p - 2, p)])
        assert F.mul(i_el, i_el) == [p - 1]
        l = 5
        while not (is_prime(l) and l % 4 == 1 and pow(l, (p - 1) // 2, p) == p - 1): l += 4
        a = next(a for a in range(1, l) if int(round((l - a * a) ** 0.5)) ** 2 == l - a * a and l - a * a > 0)
        b = int(round((l - a * a) ** 0.5))
        xb = F.vec(F.mul(i_el, [b % p]))
        x = F.el([((a % p) if j == 0 else 0) + xb[j] for j in range(D)])     # x = a + b*i
        return F.pw(x, (Q - 1) // 2 ** t), 'r=2,p=3 mod 4,e=1 (l=%d=%d^2+%d^2)' % (l, a, b), l
    # e = 0: OpenAI table step with K_q replaced by k
    l = aux_prime(p, r)
    a = structure_constants(l, r, p)
    def vmul(x, y):
        res = [[] for _ in range(r)]
        for i in range(r):
            if not x[i]: continue
            for j in range(r):
                if not y[j]: continue
                c = F.mul(x[i], y[j])
                for k in range(r):
                    if a[i][j][k]: res[k] = F.el([(u + v) % p for u, v in zip(F.vec(res[k]), F.vec(F.mul(c, [a[i][j][k]])))])
        return res
    one = [[p - 1] for _ in range(r)]
    def vpow(x, n):
        res = one; x = [c[:] for c in x]
        while n:
            if n & 1: res = vmul(res, x)
            x = vmul(x, x); n >>= 1
        return res
    # F_p-linear map u -> u^Q - zeta u on V (dimension rD), basis eta_j * t^m
    cols = []
    for j in range(r):
        for m_ in range(D):
            b = [[] for _ in range(r)]; b[j] = [0] * m_ + [1]
            img = vpow(b, Q)
            zb = [F.mul(zeta, c) if c else [] for c in b]
            col = []
            for k in range(r): col += [(x - y) % p for x, y in zip(F.vec(img[k]), F.vec(zb[k]))]
            cols.append(col)
    ker = kernel(cols_to_rows(cols), p)
    vec = ker[0]
    u = [F.el(vec[k * D:(k + 1) * D]) for k in range(r)]
    w = vpow(u, (Q ** r - 1) // r ** (t + 1))
    om = vpow(w, r)
    assert all(F.vec(c) == F.vec(om[0]) for c in om), 'omega not in k'
    g = F.el([(-x) % p for x in F.vec(om[0])])
    return g, 'auxiliary prime l=%d (e=0)' % l, l

def order_is(F, g, r, t):
    if F.pw(g, r ** t) != [1]: return False
    return t == 0 or F.pw(g, r ** (t - 1)) != [1]

def random_zeta(F, r):
    while True:
        x = [random.randrange(F.p) for _ in range(F.D)]; x = F.el(x)
        if not x: continue
        z = F.pw(x, (F.Q - 1) // r)
        if z != [1]: return z

cases = []
for p in [1009, 1019, 1031, 1033, 10007, 10009, 65537, 1000003, 1000033]:
    for r in [2, 3, 5, 7]:
        if r == p: continue
        dr = 1
        while pow(p, dr, r) != 1 % r if r > 1 else False: dr += 1
        d0 = 1
        while (p ** d0 - 1) % r: d0 += 1
        for mult in [1, r, r * r, 2, 3, 2 * r]:
            D = d0 * mult
            if D > 14 or D < 1: continue
            cases.append((p, r, D))
cases = sorted(set(cases))
random.shuffle(cases)
cases = cases[:70]
stats = {}
for (p, r, D) in sorted(cases):
    g = irreducible(p, D); F = Field(p, g)
    zeta = random_zeta(F, r)
    t = vr(F.Q - 1, r)
    out, case, extra = construct(F, r, zeta)
    ok = order_is(F, out, r, t)
    key = case.split(' (')[0].split('=%')[0]
    key = case if 'eigenvector' not in case else 'eigenvector'
    key = 'auxiliary prime' if 'auxiliary' in case else key
    key = 'r=2,p=3 mod 4,e=1' if case.startswith('r=2,p=3') else key
    stats[key] = stats.get(key, 0) + 1
    print('p=%-8d r=%d D=%-2d t=v_r(|k|-1)=%-2d  %-40s  order r^t: %s' % (p, r, D, t, case, ok))
    assert ok
print('cases by branch:', stats)
print('ALL THEOREM D CHECKS PASSED (%d fields)' % len(cases))
