"""leafcert.py - exact certificates for translated-probe kernel inequalities (standard library only).

Claim for a kernel with shifts r_1..r_m, weights lam_k >= 0 and radius 0 < rho < 1 (after scaling h = 1):
    for 0 < y <= rho, 0 <= b <= y and every real d, except the point d = 0, b = y (the zero's own pair):
        2 K_y(d,b) >= y * sum_k lam_k J_{r_k}(d,b),
    K_y(d,b) = (y-b)/(d^2+(y-b)^2) + (y+b)/(d^2+(y+b)^2),   J_r(d,b) = (K_1(d+r,b) + K_1(d-r,b))/2.
With D = d^2, v = y^2, B = b^2, R_k = r_k^2 (derive_identities.py checks every step in sympy):
    2K_y/y - sum_k lam_k J_k = 4 P / (Q* prod_k D4_k),
    P = (D+v-B) prod_k D4_k - (Q*/2) sum_k lam_k G_k prod_{j != k} D4_j,
    Q* = (d^2+(y-b)^2)(d^2+(y+b)^2)            > 0 except at d = 0, b = y,
    D4_k = prod over both signs of ((d +- r_k)^2+(1-b)^2)((d +- r_k)^2+(1+b)^2)   > 0 because b <= rho < 1.
Substituting D = z/(1-z), v = rho^2 p, B = rho^2 p u maps [0,1) x (0,1] x [0,1] onto the whole domain, and
Pt(z,p,u) = (1-z)^n P with n = 1 + 4m is a polynomial of degree <= n in each variable. So the claim follows from
Pt >= 0 on [0,1]^3.
A certificate is a finite list of dyadic boxes. check_certificate() verifies that the boxes tile [0,1]^3
(exact total volume 1 and pairwise disjoint interiors) and, on every box, recomputes all Bernstein coefficients
of Pt directly from the power basis (exact Taylor shift and basis change, integers only) and requires them >= 0.
Every scaling used is a positive integer, so signs are preserved.
"""
from fractions import Fraction as F
from math import comb, lcm


# ---------------- sparse polynomials in three variables: dict {(i,j,k): Fraction} ----------------
def const(x):
    return {(0, 0, 0): F(x)} if x else {}


def add(a, b):
    c = dict(a)
    for k, v in b.items():
        c[k] = c.get(k, 0) + v
        if not c[k]:
            del c[k]
    return c


def scale(a, x):
    x = F(x)
    return {k: v * x for k, v in a.items() if v * x}


def mul(a, b):
    c = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = (ka[0] + kb[0], ka[1] + kb[1], ka[2] + kb[2])
            c[k] = c.get(k, 0) + va * vb
    return {k: v for k, v in c.items() if v}


def prod(ps):
    out = const(1)
    for p in ps:
        out = mul(out, p)
    return out


X = {(1, 0, 0): F(1)}
Y = {(0, 1, 0): F(1)}
Zv = {(0, 0, 1): F(1)}


def degree_n(m):
    return 1 + 4 * m


# ---------------- construction A: homogenised building blocks, as in Planat's check_geometry.poly --------
def poly_A(shifts, rho):
    """shifts = [(r, lam), ...]. Variables (z, p, u); u1 = 1 - z. Every block is (1-z)^k times its D-form."""
    rho = F(rho)
    z, p, q = X, Y, Zv
    u1 = add(const(1), scale(z, -1))
    v = scale(p, rho * rho)
    B = mul(v, q)
    W = add(z, mul(u1, add(v, scale(B, -1))))                      # (1-z)(D+v-B)
    Q = add(mul(W, W), scale(mul(mul(z, u1), B), 4))               # (1-z)^2 Q*
    D4s, Gs = [], []
    for r, lam in shifts:
        R = F(r) ** 2
        C = add(z, mul(u1, add(const(R + 1), B)))                    # (1-z) C
        A = add(add(mul(C, C), scale(mul(z, u1), 4 * R)), scale(mul(B, mul(u1, u1)), -4))  # (1-z)^2 A
        D4s.append(add(mul(A, A), scale(mul(mul(z, u1), mul(C, C)), -16 * R)))              # (1-z)^4 D4
        Gs.append(add(mul(A, add(C, scale(mul(B, u1), -2))), scale(mul(mul(z, u1), C), -8 * R)))  # (1-z)^3 G
    first = mul(W, prod(D4s))
    second = {}
    for k, (r, lam) in enumerate(shifts):
        second = add(second, scale(mul(Gs[k], prod(D4s[j] for j in range(len(shifts)) if j != k)), lam))
    return add(first, scale(mul(Q, second), F(-1, 2)))


# ---------------- construction B: D-form polynomial, then homogenised by binomial expansion -------------
def poly_B(shifts, rho):
    rho = F(rho)
    s = rho * rho
    D, p, u = X, Y, Zv
    v = scale(p, s)
    B = scale(mul(p, u), s)
    W = add(add(D, v), scale(B, -1))
    Qs = add(mul(W, W), scale(mul(D, B), 4))
    D4s, Gs = [], []
    for r, lam in shifts:
        R = F(r) ** 2
        C = add(add(D, const(R + 1)), B)
        A = add(add(mul(C, C), scale(D, 4 * R)), scale(B, -4))
        D4s.append(add(mul(A, A), scale(mul(D, mul(C, C)), -16 * R)))
        Gs.append(add(mul(A, add(C, scale(B, -2))), scale(mul(D, C), -8 * R)))
    P = mul(W, prod(D4s))
    second = {}
    for k, (r, lam) in enumerate(shifts):
        second = add(second, scale(mul(Gs[k], prod(D4s[j] for j in range(len(shifts)) if j != k)), lam))
    P = add(P, scale(mul(Qs, second), F(-1, 2)))
    n = degree_n(len(shifts))
    assert max(k[0] for k in P) <= n
    out = {}
    for (i, j, k), c in P.items():          # D^i -> z^i (1-z)^(n-i)
        for t in range(n - i + 1):
            key = (i + t, j, k)
            out[key] = out.get(key, 0) + c * comb(n - i, t) * (-1) ** t
    return {k: v for k, v in out.items() if v}


def degrees(a):
    return [max(k[i] for k in a) for i in range(3)]


# ---------------- exact point check of the identity against the definitions of K and J -------------------
def Kf(y, d, b):
    return (y - b) / (d * d + (y - b) ** 2) + (y + b) / (d * d + (y + b) ** 2)


def Jf(r, d, b):
    return (Kf(F(1), d + r, b) + Kf(F(1), d - r, b)) / 2


def evalp(a, z, p, u):
    return sum((c * z ** i * p ** j * u ** k for (i, j, k), c in a.items()), F(0))


def point_identity(a, shifts, rho, d, y, b):
    """Exact: 2K_y/y - sum lam J == 4 Pt/((1-z)^n Q* prod D4) at the point (d, y, b)."""
    rho = F(rho)
    n = degree_n(len(shifts))
    lhs = 2 * Kf(y, d, b) / y - sum((F(lam) * Jf(F(r), d, b) for r, lam in shifts), F(0))
    Qs = (d * d + (y - b) ** 2) * (d * d + (y + b) ** 2)
    D4 = F(1)
    for r, lam in shifts:
        r = F(r)
        for e in (r, -r):
            D4 *= ((d + e) ** 2 + (1 - b) ** 2) * ((d + e) ** 2 + (1 + b) ** 2)
    D = d * d
    z = D / (1 + D)
    val = evalp(a, z, y * y / (rho * rho), (b * b) / (y * y)) * (1 + D) ** n   # = P(D, v, B)
    return lhs == 4 * val / (Qs * D4)


# ---------------- Bernstein machinery (integers) --------------------------------------------------------
def to_int(a):
    den = 1
    for v in a.values():
        den = lcm(den, v.denominator)
    return {k: int(v * den) for k, v in a.items()}, den


def dense(aint, n):
    N = n + 1
    c = [0] * (N ** 3)
    for (i, j, k), v in aint.items():
        assert max(i, j, k) <= n
        c[(i * N + j) * N + k] = v
    return c


def _lines(n, ax):
    N = n + 1
    for s in range(N):
        for t in range(N):
            if ax == 0:
                yield [(m * N + s) * N + t for m in range(N)]
            elif ax == 1:
                yield [(s * N + m) * N + t for m in range(N)]
            else:
                yield [(s * N + t) * N + m for m in range(N)]


_LINES = {}


def lines(n, ax):
    key = (n, ax)
    if key not in _LINES:
        _LINES[key] = list(_lines(n, ax))
    return _LINES[key]


def _weights(n):
    L = lcm(*[comb(n, m) for m in range(n + 1)])
    return [[comb(i, m) * (L // comb(n, m)) for m in range(i + 1)] for i in range(n + 1)]


_W = {}


def bern_on_box(c, n, box):
    """Bernstein coefficients (times a positive integer) on box = ((kz,dz),(kp,dp),(ku,du)), i.e. the product
    of the intervals [k/2^d, (k+1)/2^d], computed directly from dense power coefficients c."""
    if n not in _W:
        _W[n] = _weights(n)
    W = _W[n]
    c = list(c)
    for ax in range(3):
        k, d = box[ax]
        for idx in lines(n, ax):
            e = [c[idx[m]] << (d * (n - m)) for m in range(n + 1)]     # x = (k + t)/2^d, times 2^(d n)
            if k:
                for i in range(n):                                      # Taylor shift t -> t + k
                    for j in range(n - 1, i - 1, -1):
                        e[j] += k * e[j + 1]
            for m in range(n + 1):
                c[idx[m]] = e[m]
    for ax in range(3):
        for idx in lines(n, ax):
            e = [c[i] for i in idx]
            for i in range(n + 1):
                c[idx[i]] = sum(W[i][m] * e[m] for m in range(i + 1))
    return c


# ---------------- de Casteljau prover (Planat's split rule) that records the leaf boxes ----------------
def split(b, ax, n):
    left = [0] * len(b)
    right = [0] * len(b)
    for idx in lines(n, ax):
        e = [b[i] for i in idx]
        for t in range(n + 1):
            vl = sum(comb(t, j) * e[j] for j in range(t + 1))
            vr = sum(comb(n - t, j - t) * e[j] for j in range(t, n + 1))
            left[idx[t]] = vl << (n - t)
            right[idx[t]] = vr << t
    return left, right


def child(box, ax, bit):
    """Half of a dyadic box along axis ax (bit 0: lower half, bit 1: upper half)."""
    k, d = box[ax]
    nb = (2 * k + bit, d + 1)
    if ax == 0:
        return (nb, box[1], box[2])
    if ax == 1:
        return (box[0], nb, box[2])
    return (box[0], box[1], nb)


def prove_boxes(b, n, maxdepth=80, maxnodes=400000):
    stack = [(b, ((0, 0), (0, 0), (0, 0)), 0)]
    leaves, nodes, md = [], 0, 0
    while stack:
        bb, box, dep = stack.pop()
        nodes += 1
        if min(bb) >= 0:
            leaves.append(box)
            md = max(md, dep)
            continue
        if dep >= maxdepth or nodes > maxnodes:
            return None, nodes, md
        cand = []
        for ax in range(3):
            l, r = split(bb, ax, n)
            cand.append(((sum(v < 0 for v in l) + sum(v < 0 for v in r), -min(l + r), ax), l, r, ax))
        _, l, r, ax = min(cand, key=lambda x: x[0])
        stack.append((r, child(box, ax, 1), dep + 1))
        stack.append((l, child(box, ax, 0), dep + 1))
    return leaves, nodes, md


# ---------------- certificate checker ------------------------------------------------------------------
def tiles_unit_cube(boxes):
    """Exact: boxes lie in [0,1]^3, total volume 1, pairwise disjoint interiors => they cover [0,1]^3."""
    vol = F(0)
    iv = []
    for box in boxes:
        e = []
        for k, d in box:
            if not (type(k) is int and type(d) is int and d >= 0 and 0 <= k < 2 ** d):
                return False, 'box outside the cube'
            e.append((F(k, 2 ** d), F(k + 1, 2 ** d)))
        vol += F(1, 2 ** sum(d for k, d in box))
        iv.append(e)
    if vol != 1:
        return False, 'volume %s' % vol
    for a in range(len(iv)):
        for b in range(a + 1, len(iv)):
            if all(max(iv[a][x][0], iv[b][x][0]) < min(iv[a][x][1], iv[b][x][1]) for x in range(3)):
                return False, 'overlap %d %d' % (a, b)
    return True, 'tiles'


def check_certificate(a, n, boxes):
    """a: exact polynomial Pt (dict). Returns (ok, reason, minimum scaled coefficient over all leaves)."""
    if not boxes:
        return False, 'empty', None
    if max(max(k) for k in a) > n:
        return False, 'degree', None
    ok, why = tiles_unit_cube(boxes)
    if not ok:
        return False, why, None
    aint, den = to_int(a)
    c = dense(aint, n)
    mn = None
    for box in boxes:
        bc = bern_on_box(c, n, box)
        m = min(bc)
        mn = m if mn is None else min(mn, m)
        if m < 0:
            return False, 'negative Bernstein coefficient on %s' % (box,), mn
    return True, 'ok', mn


def encode_boxes(boxes):
    return ' '.join('%d.%d,%d.%d,%d.%d' % (b[0][0], b[0][1], b[1][0], b[1][1], b[2][0], b[2][1]) for b in boxes)


def decode_boxes(s):
    out = []
    for tok in s.split():
        parts = [tuple(int(x) for x in t.split('.')) for t in tok.split(',')]
        out.append(tuple(parts))
    return out


# ---------------- own-pair cap: s_r(w) = J_r(0, sqrt w) <= M on [0, V] ---------------------------------
def cap_ok(M, r, V):
    """s_r(w) = 2(1+R-w)/((1+R-w)^2+4Rw) <= M on [0,V]  <=>  M w^2 + (M(2R-2)+2) w + M(1+R)^2 - 2(1+R) >= 0,
    because the denominator (1+R-w)^2 + 4Rw is > 0. Convex quadratic: endpoints and interior vertex."""
    M, R, V = F(M), F(r) ** 2, F(V)
    if not (M > 0 and V >= 0):
        return False
    qa, qb, qc = M, M * (2 * R - 2) + 2, M * (1 + R) ** 2 - 2 * (1 + R)
    vals = [qc, qa * V * V + qb * V + qc]
    vx = -qb / (2 * qa)
    if 0 < vx < V:
        vals.append(qa * vx * vx + qb * vx + qc)
    return min(vals) >= 0


# ---------------- certificate files --------------------------------------------------------------------
def load_certs(path):
    """Lines: kind<TAB>key<TAB>boxes ; key = 'r|rho|lam' (single) or 'r1|r2|rho|lam1|lam2' (mixture)."""
    certs = {}
    try:
        text = open(path).read()
    except OSError as e:
        raise SystemExit('cannot read certificate file %s: %s' % (path, e))
    for line in text.splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        kind, key, boxes = line.rstrip('\n').split('\t')
        certs[(kind, key)] = boxes
    return certs


class Certificates:
    """Check kernels and mixtures used by a history against certificate boxes; the polynomial is rebuilt from
    the parameters in the (hashed) source file or history, never taken from the certificate file."""

    def __init__(self, paths):
        self.raw = {}
        for p in paths:
            self.raw.update(load_certs(p))
        self.done = {}

    def _check(self, kind, key, shifts, rho):
        if (kind, key) in self.done:
            return self.done[(kind, key)]
        boxes = self.raw.get((kind, key))
        ok = False
        if boxes is not None and 0 < F(rho) < 1 and all(F(lam) >= 0 for r, lam in shifts):
            a = poly_A([(F(r), F(lam)) for r, lam in shifts], F(rho))
            try:
                parsed = decode_boxes(boxes)
            except ValueError:
                parsed = []                     # malformed certificate: fail closed
            ok = check_certificate(a, degree_n(len(shifts)), parsed)[0]
        self.done[(kind, key)] = ok
        return ok

    def kernel_ok(self, k):
        key = '%s|%s|%s' % (F(k['r']), F(k['rho']), F(k['lambda']))
        return self._check('single', key, [(k['r'], k['lambda'])], k['rho'])

    def mixture_ok(self, m):
        key = '%s|%s|%s|%s|%s' % (F(m['r1']), F(m['r2']), F(m['rho']), F(m['lam1']), F(m['lam2']))
        return self._check('mix', key, [(m['r1'], m['lam1']), (m['r2'], m['lam2'])], m['rho'])
