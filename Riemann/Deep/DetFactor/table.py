"""Executable check of the 'table' construction in OpenAI, Deterministic Polynomial Factorization over Prime
Fields (4 Oct 2026), section 'Auxiliary primes and primary generators' (01-table.tex), for q = 2 and q = 3,
plus one F_{p^2} example (p = 2 mod 3: K_3 = F_{p^2}, where no element of F_p is a cubic non-residue).
Pure Python (no packages).  Run: python3 table.py
  - U_q = H_q-invariants of R_l = F_p[T]/(1+T+...+T^(l-1)): basis = the q Gauss periods eta_j (orbit sums).
  - Lemma table:degreeq: U_q is a field of degree q iff p is not a q-th power mod l (also a negative control).
  - q = 2: v with v^p = -v, c = v^2 is a non-square in F_p (it equals l, by quadratic reciprocity).
  - q = 3: K_3 = F_p(zeta_3), V_3 = U_3 (x) K_3, u with u^Q = zeta_3 u, w = u^((Q^3-1)/3^(s+1)), omega = w^3:
    omega lies in K_3 and has order exactly 3^s, s = v_3(Q-1)  (Q = |K_3|).
"""
import sys

def is_prime(n):
    if n < 2: return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0: return n == sp
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):   # deterministic for n < 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def aux_prime(p, q, want_aux=True):
    l = 12 * q + 1
    while True:
        if l != p and is_prime(l) and ((pow(p, (l - 1) // q, l) != 1) == want_aux): return l
        l += 12 * q

def gauss_periods(l, q):
    g = next(g for g in range(2, l) if all(pow(g, (l - 1) // r, l) != 1 for r in prime_factors(l - 1)))
    cos = [[] for _ in range(q)]
    x = 1
    for k in range(l - 1):
        cos[k % q].append(x); x = x * g % l      # coset j = {g^k : k = j mod q}
    return g, cos

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
    """eta_i * eta_j = sum_k a[i][j][k] eta_k in R_l (uses 1 = -sum_k eta_k)."""
    g, cos = gauss_periods(l, q)
    idx = {}
    for j, c in enumerate(cos):
        for e in c: idx[e] = j
    a = [[[0] * q for _ in range(q)] for _ in range(q)]
    for i in range(q):
        for j in range(q):
            cnt = [0] * q; ones = 0
            for x in cos[i]:
                for y in cos[j]:
                    e = (x + y) % l
                    if e == 0: ones += 1
                    else: cnt[idx[e]] += 1
            # the product is invariant: coefficient of T^e constant on cosets; T^0 = 1 = -(eta_0+...+eta_{q-1})
            per = [cnt[k] // len(cos[k]) for k in range(q)]
            assert all(cnt[k] % len(cos[k]) == 0 for k in range(q))
            a[i][j] = [(per[k] - ones) % p for k in range(q)]
    return g, cos, a

# ---- K = F_p or F_p[z]/(z^2+z+1); elements as tuples (a0, a1) meaning a0 + a1 z ----
class K:
    def __init__(s, p, quad, zeta=None):
        s.p, s.quad, s.zeta_val = p, quad, zeta
    def add(s, x, y): return ((x[0] + y[0]) % s.p, (x[1] + y[1]) % s.p)
    def sub(s, x, y): return ((x[0] - y[0]) % s.p, (x[1] - y[1]) % s.p)
    def mul(s, x, y):
        p = s.p
        if not s.quad: return (x[0] * y[0] % p, 0)
        a0, a1 = x; b0, b1 = y          # z^2 = -z - 1
        c0 = a0 * b0 - a1 * b1; c1 = a0 * b1 + a1 * b0 - a1 * b1
        return (c0 % p, c1 % p)
    def inv(s, x):
        return s.pow(x, s.size() - 2)
    def pow(s, x, e):
        r = (1, 0)
        while e:
            if e & 1: r = s.mul(r, x)
            x = s.mul(x, x); e >>= 1
        return r
    def size(s): return s.p ** 2 if s.quad else s.p
    def zeta(s): return (0, 1) if s.quad else (s.zeta_val, 0)
    def zero(s): return (0, 0)
    def one(s): return (1, 0)

class V:
    """U_q (x) K with basis eta_0..eta_{q-1}; multiplication from structure constants a."""
    def __init__(s, Kf, a, q): s.K, s.a, s.q = Kf, a, q
    def mul(s, x, y):
        K, q, a = s.K, s.q, s.a
        r = [K.zero()] * q
        for i in range(q):
            if x[i] == (0, 0): continue
            for j in range(q):
                if y[j] == (0, 0): continue
                c = K.mul(x[i], y[j])
                for k in range(q):
                    if a[i][j][k]: r[k] = K.add(r[k], K.mul(c, (a[i][j][k], 0)))
        return r
    def one(s):
        m1 = (s.K.p - 1) % s.K.p
        return [(m1, 0)] * s.q      # 1 = -(eta_0 + ... + eta_{q-1})
    def pow(s, x, e):
        r = s.one()
        while e:
            if e & 1: r = s.mul(r, x)
            x = s.mul(x, x); e >>= 1
        return r
    def is_scalar(s, x):   # lies in K (x) 1 iff all coordinates equal
        return all(c == x[0] for c in x)

def nullvec(Kf, M):
    """nonzero kernel vector of square matrix M over K (Gaussian elimination)."""
    n = len(M); A = [row[:] for row in M]; piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, n) if A[i][c] != (0, 0)), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]
        iv = Kf.inv(A[r][c]); A[r] = [Kf.mul(iv, v) for v in A[r]]
        for i in range(n):
            if i != r and A[i][c] != (0, 0):
                f = A[i][c]; A[i] = [Kf.sub(A[i][k], Kf.mul(f, A[r][k])) for k in range(n)]
        piv.append(c); r += 1
    free = [c for c in range(n) if c not in piv]
    assert free, "no kernel"
    fc = free[0]; v = [Kf.zero()] * n; v[fc] = Kf.one()
    for i, c in enumerate(piv): v[c] = Kf.sub(Kf.zero(), A[i][fc])
    return v

def v3(n):
    s = 0
    while n % 3 == 0: n //= 3; s += 1
    return s

def sqrt_mod(a, p, nonres):
    """Tonelli-Shanks with a supplied non-residue."""
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    s, t = 0, p - 1
    while t % 2 == 0: t //= 2; s += 1
    c = pow(nonres, t, p); x = pow(a, (t + 1) // 2, p); b = pow(a, t, p); m = s
    while b != 1:
        i, bb = 0, b
        while bb != 1: bb = bb * bb % p; i += 1
        cc = pow(c, 1 << (m - i - 1), p); x = x * cc % p; c = cc * cc % p; b = b * c % p; m = i
    return x

def period_poly_has_root(l, q, p):
    """U_q split (not a field) iff the minimal polynomial of eta_0 has a root in F_p; test via p-th power
    Frobenius: eta_0^p = eta_0 iff eta_0 in F_p.  Frobenius permutes eta_j by p's coset index."""
    g, cos = gauss_periods(l, q)
    return p % l in set(cos[0])     # p in H_q  <=>  Frobenius fixes every eta_j  <=>  U_q = F_p^q

def run(p, verbose=True):
    out = []
    # ---- q = 2 ----
    l2 = aux_prime(p, 2)
    g, cos, a = structure_constants(l2, 2, p)
    Kp = K(p, False, 0)
    Vq = V(Kp, a, 2)
    # Frobenius matrix of x -> x^p on U_2 (F_p-linear); solve v^p = -v
    cols = [Vq.pow([(1, 0) if k == j else (0, 0) for k in range(2)], p) for j in range(2)]
    M = [[Kp.add(cols[j][i], ((1 if i == j else 0), 0)) for j in range(2)] for i in range(2)]
    vv = nullvec(Kp, M)
    c = Vq.mul(vv, vv)
    assert Vq.is_scalar(c)
    cval = (-c[0][0]) % p            # scalar s*1 = -s*(eta_0+eta_1): coordinate = -s
    assert pow(cval, (p - 1) // 2, p) == p - 1, "c is not a non-square"
    out.append("q=2: l_2=%d, c = v^2 = %d is a non-square mod p (c == l_2 mod p: %s)" % (l2, cval, cval == l2 % p))
    nonres = cval
    # ---- q = 3 ----
    l3 = aux_prime(p, 3)
    assert not period_poly_has_root(l3, 3, p)
    lbad = aux_prime(p, 3, want_aux=False)          # negative control: p IS a cube mod lbad
    assert period_poly_has_root(lbad, 3, p)
    g, cos, a = structure_constants(l3, 3, p)
    # independent check of Lemma table:degreeq by arithmetic: eta_0^p != eta_0 but eta_0^(p^3) = eta_0 in U_3,
    # and in the control algebra (p a cube mod lbad) eta_0^p = eta_0 (U_3 split).
    U3 = V(K(p, False, 0), a, 3); e0 = [(1, 0), (0, 0), (0, 0)]
    assert U3.pow(e0, p) != e0 and U3.pow(e0, p ** 3) == e0
    _, _, abad = structure_constants(lbad, 3, p); U3b = V(K(p, False, 0), abad, 3)
    assert U3b.pow(e0, p) == e0
    if p % 3 == 1:
        r = sqrt_mod(p - 3, p, nonres)               # zeta_3 = (-1 + sqrt(-3))/2
        z = (-1 + r) * pow(2, p - 2, p) % p
        assert (z * z + z + 1) % p == 0
        Kf = K(p, False, z)
    else:
        Kf = K(p, True)
    Vq = V(Kf, a, 3)
    Q = Kf.size()
    zeta = Kf.zeta()
    cols = [Vq.pow([Kf.one() if k == j else Kf.zero() for k in range(3)], Q) for j in range(3)]
    M = [[Kf.sub(cols[j][i], zeta if i == j else Kf.zero()) for j in range(3)] for i in range(3)]
    u = nullvec(Kf, M)
    assert Vq.pow(u, Q) == [Kf.mul(zeta, x) for x in u]
    s = v3(Q - 1)
    w = Vq.pow(u, (Q ** 3 - 1) // 3 ** (s + 1))
    om = Vq.pow(w, 3)
    assert Vq.is_scalar(om), "omega not in K_3"
    omK = Kf.sub(Kf.zero(), om[0])                  # coordinate -> scalar
    assert Kf.pow(omK, 3 ** s) == Kf.one() and (s == 0 or Kf.pow(omK, 3 ** (s - 1)) != Kf.one())
    if s > 0:
        assert Kf.pow(omK, (Q - 1) // 3) != Kf.one()   # cubic non-residue in K_3
    out.append("q=3: l_3=%d (control l=%d splits U_3), K_3=F_p^%d, s=v_3(Q-1)=%d, omega_3=%s has order 3^%d, in K_3: yes"
               % (l3, lbad, 2 if Kf.quad else 1, s, omK, s))
    # least small cubic non-residue a + b*z in the basis {1, z} of K_3 (empirical, for the F_{p^2} case)
    if Kf.quad and s > 0:
        best = None
        for h in range(0, 50):
            for a0 in range(0, h + 1):
                b0 = h - a0
                x = (a0 % p, b0 % p)
                if x == (0, 0): continue
                if Kf.pow(x, (Q - 1) // 3) != Kf.one():
                    best = (a0, b0); break
            if best: break
        out.append("     least cubic non-residue a + b*zeta_3 in F_{p^2} by height a+b: %s" % (best,))
        # quadratic non-residue in F_{p^2} via the norm: x non-square iff N(x) = x^{p+1} non-square in F_p
        for a0 in range(0, 50):
            x = (a0, 1)
            nx = Kf.pow(x, p + 1)
            assert nx[1] == 0
            if pow(nx[0], (p - 1) // 2, p) == p - 1:
                assert Kf.pow(x, (Q - 1) // 2) != Kf.one()
                out.append("     least non-square of the form a + zeta_3 in F_{p^2}: a = %d (norm test agrees)" % a0); break
    if verbose:
        print("p = %d (p mod 3 = %d)" % (p, p % 3))
        for line in out: print("   " + line)

if __name__ == "__main__":
    ps = [7, 11, 13, 1000003, 1000033, 10**9 + 7, 10**9 + 9, 998244353, (1 << 61) - 1, 10**18 + 9]
    ps = [p for p in ps if is_prime(p)]
    for p in ps: run(p)
    print("ALL TABLE CHECKS PASSED (%d primes)" % len(ps))
