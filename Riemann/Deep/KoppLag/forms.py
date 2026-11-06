"""Exact arithmetic for indefinite binary quadratic forms (a, b, c), discriminant D = b^2 - 4ac > 0 nonsquare.
Reduced (Gauss/Lagrange):  0 < b < sqrt(D),  sqrt(D) - b < 2|a| < sqrt(D) + b.
Proper classes of primitive forms  <->  narrow class group Cl+(O_D); wide group Cl(O_D) = Cl+ / <J>, J = class of (-1, b, c).
All integer arithmetic is exact (Python ints).  No floating point anywhere in this module.
"""
from math import isqrt, gcd

def egcd(a, b):
    """return (g, x, y) with x*a + y*b = g = gcd(a, b) >= 0"""
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return a, x0, y0

def disc(f):
    a, b, c = f
    return b * b - 4 * a * c

def is_reduced(f, s):
    """s = isqrt(D), D nonsquare.  b < sqrt D <=> b <= s;  2|a| + b > sqrt D <=> 2|a| + b > s;  2|a| - b < sqrt D <=> 2|a| - b <= s."""
    a, b, c = f
    A = abs(a)
    return 0 < b <= s and 2 * A + b > s and 2 * A - b <= s

def rho(f, D, s):
    """Cohen, CCANT Def. 5.6.4: rho(a,b,c) = (c, r, (r^2 - D)/(4c)), r = -b mod 2|c|, r in (-|c|,|c|] if |c| > sqrt D, else r in (sqrt D - 2|c|, sqrt D)."""
    a, b, c = f
    C = abs(c)
    m = 2 * C
    if C > s:                       # |c| > sqrt(D)  (|c| integer, D nonsquare)
        r = (-b) % m
        if r > C:
            r -= m
    else:
        r = s - ((s + b) % m)       # largest r <= s with r = -b (mod 2|c|)
    num = r * r - D
    assert num % (4 * c) == 0
    return (c, r, num // (4 * c))

def reduce_form(f, D, s, maxit=10**6):
    it = 0
    while not is_reduced(f, s):
        f = rho(f, D, s)
        it += 1
        if it > maxit:
            raise RuntimeError('reduction did not terminate')
    return f

def make_pos(f):
    """properly equivalent form with a > 0 (input: a*c < 0 or a > 0)."""
    a, b, c = f
    if a > 0:
        return f
    if c > 0:
        return (c, -b, a)           # (x,y) -> (-y, x), det 1
    raise ValueError('need a > 0 or c > 0')

def compose(f1, f2, D):
    """Cohen, CCANT Algorithm 5.4.7 (composition of primitive forms, a1, a2 > 0). Output not reduced."""
    a1, b1, c1 = make_pos(f1)
    a2, b2, c2 = make_pos(f2)
    if a1 > a2:
        a1, b1, c1, a2, b2, c2 = a2, b2, c2, a1, b1, c1
    s2 = (b1 + b2) // 2
    n = b2 - s2
    if a2 % a1 == 0:
        y1, d = 0, a1
    else:
        d, u, v = egcd(a2, a1)      # u*a2 + v*a1 = d
        y1 = u
    if s2 % d == 0:
        y2, x2, d1 = -1, 0, d
    else:
        d1, x2, y2 = egcd(s2, d)    # x2*s2 + y2*d = d1
        y2 = -y2
    v1, v2 = a1 // d1, a2 // d1
    r = (y1 * y2 * n - x2 * c2) % v1
    b3 = b2 + 2 * v2 * r
    a3 = v1 * v2
    num = b3 * b3 - D
    assert num % (4 * a3) == 0, (f1, f2)
    return (a3, b3, num // (4 * a3))

def reduced_forms(D, spf=None):
    """all reduced primitive forms of discriminant D (exact enumeration over b, divisors a of (D-b^2)/4)."""
    s = isqrt(D)
    assert s * s != D and D % 4 in (0, 1)
    out = []
    b = s if (s - D) % 2 == 0 else s - 1
    while b > 0:
        M = (D - b * b) // 4
        lo = (s - b) // 2 + 1          # 2a + b > s  <=> a > (s-b)/2
        hi = (s + b) // 2              # 2a - b <= s <=> a <= (s+b)/2
        for a in divisors(M, spf):
            if lo <= a <= hi:
                c = M // a
                if gcd(gcd(a, b), c) == 1:
                    out.append((a, b, -c))
                    out.append((-a, b, c))
        b -= 2
    return out

def factor(m, spf=None):
    fs = {}
    if spf is not None and m < len(spf):
        while m > 1:
            p = int(spf[m])
            while m % p == 0:
                m //= p
                fs[p] = fs.get(p, 0) + 1
        return fs
    p = 2
    while p * p <= m:
        while m % p == 0:
            m //= p
            fs[p] = fs.get(p, 0) + 1
        p += 1 if p == 2 else 2
    if m > 1:
        fs[m] = fs.get(m, 0) + 1
    return fs

def divisors(m, spf=None):
    ds = [1]
    for p, e in factor(m, spf).items():
        ds = [d * p ** k for d in ds for k in range(e + 1)]
    return ds

def class_data(D, spf=None):
    """Exact narrow/wide class group data of the order of discriminant D > 0 (nonsquare).
    returns dict: hplus, h, cyc (form -> cycle id), reps, pid (principal cycle), jid (cycle of (-1,b,c)),
    narrow2 (Cl+ = Cl+[2]), wide2 (Cl = Cl[2]), checks agreeing composition test and inverse test."""
    s = isqrt(D)
    forms = reduced_forms(D, spf)
    fset = set(forms)
    cyc = {}
    reps = []
    for f in forms:
        if f in cyc:
            continue
        cid = len(reps)
        reps.append(f)
        g = f
        while True:
            cyc[g] = cid
            g = rho(g, D, s)
            assert g in fset, ('rho left the reduced set', D, g)
            if g == f:
                break
            assert g not in cyc, ('cycles overlap', D, g)
    # principal form (1, b0, c0) reduced, and J = (-1, b0, -c0)
    b0 = s if (s - D) % 2 == 0 else s - 1
    P = (1, b0, (b0 * b0 - D) // 4)
    J = (-1, b0, -((b0 * b0 - D) // 4))
    assert is_reduced(P, s) and is_reduced(J, s)
    pid, jid = cyc[P], cyc[J]
    hplus = len(reps)
    h = hplus if pid == jid else hplus // 2
    assert h * (1 if pid == jid else 2) == hplus
    narrow2 = wide2 = True
    for f in reps:
        sq = reduce_form(compose(f, f, D), D, s)
        q = cyc[sq]
        a, b, c = f
        inv_in = cyc[(c, b, a)] == cyc[f]          # [f]^{-1} = [(c,b,a)]
        jinv_in = cyc[(-c, b, -a)] == cyc[f]       # J [f]^{-1} = [(-c,b,-a)]
        # agreement between composition test and inverse test
        assert (q == pid) == inv_in, ('narrow test disagreement', D, f)
        assert (q in (pid, jid)) == (inv_in or jinv_in), ('wide test disagreement', D, f)
        if q != pid:
            narrow2 = False
        if q not in (pid, jid):
            wide2 = False
    return dict(D=D, hplus=hplus, h=h, nforms=len(forms), cyc=cyc, reps=reps, pid=pid, jid=jid,
                narrow2=narrow2, wide2=wide2)

def class_order(f, D, s, cyc, pid, maxord=10**6):
    """order of the narrow class of reduced f (exact)."""
    g = f
    k = 1
    while cyc[g] != pid:
        g = reduce_form(compose(g, f, D), D, s)
        k += 1
        if k > maxord:
            raise RuntimeError
    return k

def group_structure_2(D, data):
    """invariants of Cl+ and Cl when they are 2-groups (from counts of elements of order 2^k). returns (narrow, wide) as tuples of cyclic orders, or None if not a 2-group."""
    s = isqrt(D)
    cyc, reps, pid, jid = data['cyc'], data['reps'], data['pid'], data['jid']
    orders = [class_order(f, D, s, cyc, pid) for f in reps]
    def inv_from_orders(ords):
        if any(o & (o - 1) for o in ords):
            return None
        n = len(ords)
        # |G[2^k]| for k = 0.. ; invariants of abelian 2-group
        cnt = {}
        k = 0
        res = []
        sizes = []
        while True:
            c = sum(1 for o in ords if (2 ** k) % o == 0)
            sizes.append(c)
            if c == n:
                break
            k += 1
        # number of cyclic factors of order >= 2^k  = log2(|G[2^k]| / |G[2^(k-1)]|)
        import math
        r = [round(math.log2(sizes[i] / sizes[i - 1])) for i in range(1, len(sizes))]
        # r[k-1] = number of invariants >= 2^k
        inv = []
        for k in range(1, len(r) + 1):
            nk = r[k - 1] - (r[k] if k < len(r) else 0)
            inv += [2 ** k] * nk
        return tuple(sorted(inv))
    narrow = inv_from_orders(orders)
    # wide: orders in Cl+/<J>: order of class C modulo J = least k with C^k in {1, J}
    if pid == jid:
        wide = narrow
    else:
        word = []
        for f in reps:
            g = f
            k = 1
            while cyc[g] not in (pid, jid):
                g = reduce_form(compose(g, f, D), D, s)
                k += 1
            word.append(k)
        # each wide class has 2 narrow preimages: take every narrow class once, halve multiset
        word.sort()
        wide = inv_from_orders(word[::2]) if len(word) % 2 == 0 else None
    return narrow, wide
