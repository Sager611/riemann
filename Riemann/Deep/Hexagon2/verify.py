#!/usr/bin/env python3
"""Re-check the finite claims of

  B. Liu, "Selected character moments, zero-free half-planes for Hecke
  L-functions over number fields, and Quasi-Riemann Hypothesis",
  Hexagon preprint hexagon:2611.00082v1 (hexagonmath.org, 4 Nov 2026;
  source number_field_new_v2-1.tex), boundary B_2 = (777 - sqrt 1081)/852.

"l.N" = line N of number_field_new_v2-1.tex.  Stdlib + mpmath only.

  * Exact arithmetic: fractions.Fraction and the class QZ for Q(sqrt 1081)
    (signs decided exactly by squaring).
  * Rectangle inequalities (Prop. sc:rational-certificate, Lemma
    sc:double-contact, App. A): the multiplied function is a polynomial; we
    recover it exactly by tensor Lagrange interpolation over Q or Q(sqrt 1081),
    confirm it at extra random rational points, and run the paper's own
    rational Bernstein bisection (delta at even depth, x at odd depth).  A PASS
    is an exact proof.  We also report our (depth, rectangles, k) next to the
    paper's table entries.
  * Exponent bookkeeping of Sections 3-4 (affine identities and local
    valuation tables): exact, exhaustive over valuations <= 60 where tables
    are claimed for all valuations.
  * Counter-instance search: (i) the Gauss-mean Lemma mom:gauss-mean and the
    correlation facts it uses, tested over Z with Jacobi symbols (all odd
    moduli, natural zeros kept); (ii) random search on the restoration
    optimisation of Thm op:restoration.

Run:  python3 verify.py      (single thread; about 20 s and 30 MB on an M-series Mac,
      Python 3.14 + mpmath 1.3.0; 99 checks; runtime and memory are printed at the end)
"""
import cmath
import math
import random
import resource
import sys
import time
from fractions import Fraction as Fr

import mpmath
from mpmath import mp, mpf

mp.dps = 60
random.seed(20261105)
T0 = time.time()
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("  [" + detail + "]") if detail else ""))
    sys.stdout.flush()


# ----------------------------------------------------------------------------
# Q(sqrt 1081)
# ----------------------------------------------------------------------------
NSQ = 1081
ZLO = Fr(3287856444554719302046, 10**20)   # paper's enclosure (App. A, l.4600)
ZHI = Fr(3287856444554719302047, 10**20)


class QZ:
    __slots__ = ("a", "b")

    def __init__(s, a, b=0):
        if isinstance(a, QZ):
            s.a, s.b = a.a, a.b + Fr(b)
        else:
            s.a = Fr(a); s.b = Fr(b)

    @staticmethod
    def c(o):
        return o if isinstance(o, QZ) else QZ(o)

    def __add__(s, o):
        o = QZ.c(o); return QZ(s.a + o.a, s.b + o.b)
    __radd__ = __add__

    def __sub__(s, o):
        o = QZ.c(o); return QZ(s.a - o.a, s.b - o.b)

    def __rsub__(s, o):
        return QZ.c(o) - s

    def __mul__(s, o):
        o = QZ.c(o); return QZ(s.a * o.a + NSQ * s.b * o.b, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__

    def __neg__(s):
        return QZ(-s.a, -s.b)

    def __truediv__(s, o):
        o = QZ.c(o)
        n = o.a * o.a - NSQ * o.b * o.b
        return s * QZ(o.a / n, -o.b / n)

    def __rtruediv__(s, o):
        return QZ.c(o) / s

    def __pow__(s, k):
        r = QZ(1)
        for _ in range(k):
            r = r * s
        return r

    def sign(s):
        a, b = s.a, s.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a >= 0 and b >= 0:
            return 1
        if a <= 0 and b <= 0:
            return -1
        if a > 0:  # b < 0
            return 1 if a * a > NSQ * b * b else -1
        return 1 if NSQ * b * b > a * a else -1  # a < 0 < b

    def __eq__(s, o):
        o = QZ.c(o); return s.a == o.a and s.b == o.b

    def __hash__(s):
        return hash((s.a, s.b))

    def lower(s):
        return s.a + s.b * (ZLO if s.b >= 0 else ZHI)

    def upper(s):
        return s.a + s.b * (ZHI if s.b >= 0 else ZLO)

    def mpf(s):
        return mpf(s.a.numerator) / s.a.denominator + mpf(s.b.numerator) / s.b.denominator * mpmath.sqrt(NSQ)


def lowr(v):
    return v.lower() if isinstance(v, QZ) else Fr(v)


def uppr(v):
    return v.upper() if isinstance(v, QZ) else Fr(v)


def sgn(v):
    return v.sign() if isinstance(v, QZ) else ((v > 0) - (v < 0))


def tompf(v):
    return v.mpf() if isinstance(v, QZ) else mpf(v.numerator) / v.denominator


Z0 = QZ(0, 1)
B2 = (777 - Z0) / 852
D2 = (61 - Z0) / 66                       # delta_2 (l.3876)
b2 = Fr(-11, 71) + Z0 * Fr(595, 76751)    # b_2
BRAT = Fr(8733821, 10**7)
bRAT = Fr(499781491, 5 * 10**9)
P1 = Fr(1124661, 250000)
KAP1 = Fr(374887, 500000)
LST = Fr(21, 100)


# ----------------------------------------------------------------------------
# Section 6 model (l.3195-3460)
# ----------------------------------------------------------------------------
def geom(B, b):
    l = Fr(11, 3) - 4 * B
    lx = (1 - b - l) / 2
    ly = (1 + b - l) / 2
    h = (1 + b + 3 * l) / 2
    T = 2 * ly - h
    C0 = Fr(-1, 4) + b / 6 + 5 * l / 4
    return dict(l=l, lx=lx, ly=ly, h=h, T=T, C0=C0)


def model(B, b, dl, x, p1=P1):
    g = geom(B, b)
    l, ly, h, T, C0 = g["l"], g["ly"], g["h"], g["T"], g["C0"]
    th0 = dl * x / p1
    A = Fr(5, 6) - dl
    md = (Fr(5, 6) - dl * x) / (Fr(5, 2) - dl - 2 * dl * x)
    rd, td = 2 * md, 3 * md
    Fd = x + 2 * (1 - x) * md
    Rd = 1 - dl * Fd
    c = 2 / p1
    Dc = 3 - (1 + 2 * c) * x
    Pc = (2 - 2 * c * x) * (1 - x)
    J = A * Dc + dl * Pc
    Rs = 1 - dl + A * dl * Pc / (2 * J)
    E = lambda R: C0 + dl / 2 + dl * x * (l + (h / 2 - ly) / p1) - 3 * h * (1 - R) / 4
    S = lambda R: (3 * R + 1) / 4 + th0 / 2 - Fr(17, 50)
    N = lambda R: (3 * R + 1) / 4 - th0 / 2
    Dsh = lambda R: R + (Fr(1, 2) - th0) * Fr(3, 40) - N(R)
    JN = Fr(5, 2) - 5 / p1
    gam2 = lambda m: dl * x * JN / 10 * (1 - 2 * m)
    GN = dl * x * JN / 10 * (2 * ly - h)
    GL = th0 * (Fr(70, 100) * h - ly)
    GR = (dl / 2 - th0) * (ly - Fr(53, 100) * h)
    return dict(g, th0=th0, A=A, md=md, rd=rd, td=td, Fd=Fd, Rd=Rd, J=J, Rs=Rs, E=E, S=S, N=N, Dsh=Dsh,
                gam2=gam2, GN=GN, GL=GL, GR=GR, Dc=Dc, Pc=Pc, c=c)


def gact(B, b, dl, x):
    """Active margin g = G_N - E(R_d), closed form (sc:active-function, l.3881)."""
    G = geom(B, b)
    md = (Fr(5, 6) - dl * x) / (Fr(5, 2) - dl - 2 * dl * x)
    Rd = 1 - dl * (x + 2 * (1 - x) * md)
    return (Fr(1, 4) - b / 6 - 5 * G["l"] / 4 - dl / 2 - dl * x * G["l"] + 3 * G["h"] * (1 - Rd) / 4
            + dl * x * G["T"] / 4)


# ----------------------------------------------------------------------------
# exact bivariate polynomials, interpolation, Bernstein
# ----------------------------------------------------------------------------
def lagrange_basis(nodes):
    """Coefficient lists (low->high) of the Lagrange basis polynomials."""
    res = []
    for i, xi in enumerate(nodes):
        coef = [Fr(1)]
        den = Fr(1)
        for j, xj in enumerate(nodes):
            if j == i:
                continue
            coef = [Fr(0)] + coef            # multiply by t
            for k in range(len(coef) - 1):
                coef[k] -= xj * coef[k + 1]
            den *= (xi - xj)
        res.append([cf / den for cf in coef])
    return res


def zero_like(v):
    return QZ(0) if isinstance(v, QZ) else Fr(0)


def interp2(F, nd, nx, dnodes=None, xnodes=None):
    dn = dnodes or [Fr(2 * i + 1, 17) for i in range(nd + 1)]
    xn = xnodes or [Fr(3 * j + 1, 19) for j in range(nx + 1)]
    Ld, Lx = lagrange_basis(dn), lagrange_basis(xn)
    vals = [[F(di, xj) for xj in xn] for di in dn]
    z = zero_like(vals[0][0])
    C = {}
    for a in range(nd + 1):
        for bb in range(nx + 1):
            s = z
            for i in range(nd + 1):
                if Ld[i][a] == 0:
                    continue
                for j in range(nx + 1):
                    if Lx[j][bb] == 0:
                        continue
                    s = s + vals[i][j] * (Ld[i][a] * Lx[j][bb])
            if not (s == 0 if not isinstance(s, QZ) else (s.a == 0 and s.b == 0)):
                C[(a, bb)] = s
    return C


def peval(C, d, x):
    s = Fr(0)
    for (i, j), c in C.items():
        s = s + c * (d ** i) * (x ** j)
    return s


def isnull(v):
    return (v.a == 0 and v.b == 0) if isinstance(v, QZ) else v == 0


def poly_of(F, nd=6, nx=6, ntest=4):
    """Exact polynomial equal to F, with a confirmation at random points."""
    C = interp2(F, nd, nx)
    for _ in range(ntest):
        d = Fr(random.randint(1, 997), random.randint(1000, 1999))
        x = Fr(random.randint(1, 499), random.randint(1000, 1999))
        if not isnull(peval(C, d, x) - F(d, x)):
            return None
    return C


def comb(n, k):
    return math.comb(n, k)


def shift_poly(C, a, b, c, d):
    """Coefficients of P(a+(b-a)v, c+(d-c)w) in (v,w)."""
    sd, sx = b - a, d - c
    out = {}
    for (i, j), cf in C.items():
        for k in range(i + 1):
            ck = comb(i, k) * a ** (i - k) * sd ** k
            for l in range(j + 1):
                t = cf * (ck * comb(j, l) * c ** (j - l) * sx ** l)
                out[(k, l)] = out.get((k, l), zero_like(cf)) + t
    return out


def bernstein(C, rect, M=None, N=None):
    (a, b), (c, d) = rect
    if M is None:
        M = max([i for i, _ in C] + [0])
        N = max([j for _, j in C] + [0])
    S = shift_poly(C, a, b, c, d)
    out = []
    for k in range(M + 1):
        for l in range(N + 1):
            s = Fr(0)
            for (i, j), cf in S.items():
                if i <= k and j <= l:
                    s = s + cf * (Fr(comb(k, i), comb(M, i)) * Fr(comb(l, j), comb(N, j)))
            out.append(s)
    return out


def bnb_bern(C, rect, maxdepth=26):
    """Paper's rule: accept if all Bernstein coefficients > 0 (rational lower
    bounds for Q(sqrt1081)); else bisect delta at even depth, x at odd depth.
    Returns (ok, maxdepth_used, n_accepted, min accepted coefficient)."""
    M = max([i for i, _ in C] + [0]); N = max([j for _, j in C] + [0])
    stack = [(rect, 0)]
    nacc, dmax, mmin = 0, 0, None
    while stack:
        (r, dep) = stack.pop()
        bc = bernstein(C, r, M, N)
        lows = [lowr(v) for v in bc]
        if min(lows) > 0:
            nacc += 1; dmax = max(dmax, dep)
            mm = min(lows); mmin = mm if mmin is None else min(mmin, mm)
            continue
        if dep >= maxdepth:
            return False, dep, nacc, mmin
        (a, b), (c, d) = r
        if dep % 2 == 0:
            m = (a + b) / 2; stack += [(((a, m), (c, d)), dep + 1), (((m, b), (c, d)), dep + 1)]
        else:
            m = (c + d) / 2; stack += [(((a, b), (c, m)), dep + 1), (((a, b), (m, d)), dep + 1)]
    return True, dmax, nacc, mmin


def upper_bern(C, rect):
    return max(uppr(v) for v in bernstein(C, rect))


# ----------------------------------------------------------------------------
print("== A. Endpoint algebra (l.53-90, 3129-3140, 3869-3995, 4741-4790)")
# A1 quadratic
ok = (1278 * B2 * B2 - 2331 * B2 + 1061) == 0 and (B2 - (777 + Z0) / 852).sign() < 0
check("A1 B_2=(777-sqrt1081)/852 is the smaller root of 1278B^2-2331B+1061", ok,
      mpmath.nstr(B2.mpf(), 40))
dig = "0.8733819666132075199290267786423718694"
check("A2 printed digits of B_2 (l.88)", mpmath.nstr(B2.mpf(), 37).startswith(dig[:37]) or
      str(B2.mpf())[:38] == dig[:38], str(B2.mpf())[:40])

# A3 derive B_2 from scratch: g affine in b; coefficient of b at x=1/2
def gb_coef(dl):
    return gact(Fr(0), Fr(1), dl, Fr(1, 2)) - gact(Fr(0), Fr(0), dl, Fr(1, 2))
ok = all(gb_coef(d) * 24 * (5 - 4 * d) == -(33 * d * d - 61 * d + 20)
         for d in [Fr(1, 7), Fr(2, 5), Fr(3, 11), Fr(5, 9), Fr(13, 20)])
check("A3 d/db g(B,b,delta,1/2) = -(33delta^2-61delta+20)/(24(5-4delta)) (sc:b-derivative)", ok)
ok = (33 * D2 * D2 - 61 * D2 + 20) == 0
check("A4 delta_2=(61-sqrt1081)/66 is the smaller zero (l.3876)", ok, mpmath.nstr(D2.mpf(), 25))
l2 = (5 - 6 * D2) / (3 * (3 + 4 * D2))
Bder = (Fr(11, 3) - l2) / 4
# g(B,b,delta2,1/2) is b-free and equals 5/12-delta/2-(3/4+delta)l : check at two (B,b)
def gQZ(B, b, dl, x):
    l = Fr(11, 3) - 4 * B
    h = (1 + b + 3 * l) / 2; ly = (1 + b - l) / 2; T = 2 * ly - h
    md = (Fr(5, 6) - dl * x) / (Fr(5, 2) - dl - 2 * dl * x)
    Rd = 1 - dl * (x + 2 * (1 - x) * md)
    return Fr(1, 4) - b / 6 - 5 * l / 4 - dl / 2 - dl * x * l + 3 * h * (1 - Rd) / 4 + dl * x * T / 4
ok = all(gQZ(QZ(B), QZ(b), D2, Fr(1, 2)) == (Fr(5, 12) - D2 / 2 - (Fr(3, 4) + D2) * (Fr(11, 3) - 4 * B))
         for B, b in [(Fr(87, 100), Fr(1, 10)), (Fr(9, 10), Fr(0)), (Fr(1, 3), Fr(7, 5))])
check("A5 at delta_2, g(B,b,delta_2,1/2)=5/12-delta_2/2-(3/4+delta_2)l is independent of b", ok)
check("A6 g=0 there gives l=(5-6delta_2)/(3(3+4delta_2)) and B=(11/3-l)/4 = B_2 exactly", Bder == B2,
      "from first principles")
# A7 factorisation and b_2 (double contact in delta)
coef = Fr(178, 213) - Z0 * Fr(5179, 460506)
ok = all((5 - 4 * d) * gQZ(B2, b2, QZ(d), Fr(1, 2)) == coef * (QZ(d) - D2) ** 2
         for d in [Fr(0), Fr(1, 3), Fr(1, 2), Fr(5, 7), Fr(1), Fr(13, 20)])
check("A7 (5-4delta) g(B_2,b_2,delta,1/2) = (178/213-5179 sqrt1081/460506)(delta-delta_2)^2 (sc:contact-factorization)",
      ok, "degree-2 identity at 6 points")
# b_2 is forced: with b free, the double root needs d/ddelta g = 0 at delta_2
def gx(b, d):
    return (5 - 4 * d) * gQZ(B2, QZ(b) if not isinstance(b, QZ) else b, QZ(d) if not isinstance(d, QZ) else d, Fr(1, 2))
# g is affine in b: G(b,delta) = G0(delta)+b G1(delta); (5-4delta)G0,(5-4delta)G1 are quadratics -> derivative at delta_2
def quad_coeffs(f):
    v = [f(Fr(t)) for t in (0, 1, 2)]
    c0 = v[0]; c2 = (v[2] - 2 * v[1] + v[0]) / 2; c1 = v[1] - v[0] - c2
    return c0, c1, c2
G0 = quad_coeffs(lambda d: gx(Fr(0), d)); G1 = quad_coeffs(lambda d: gx(Fr(1), d) - gx(Fr(0), d))
der = lambda c, t: c[1] + 2 * c[2] * t
val = lambda c, t: c[0] + c[1] * t + c[2] * t * t
bstar = -der(G0, D2) / der(G1, D2)
check("A8 b_2 = -11/71+595 sqrt1081/76751 is the unique b with a double root at delta_2", bstar == b2 and
      val(G0, D2) + b2 * val(G1, D2) == 0, mpmath.nstr(b2.mpf(), 20))
# A9 dB g formula
ok = True
for _ in range(6):
    d = Fr(random.randint(1, 99), 100); x = Fr(random.randint(0, 50), 100)
    lhs = gact(Fr(1), bRAT, d, x) - gact(Fr(0), bRAT, d, x)
    Dl = 5 - 2 * d - 4 * d * x
    ok &= lhs == (-26 * d * d * x * x + 14 * d * d * x + 5 * d * x - 25 * d + 25) / Dl
check("A9 g is affine in B with dB g=(-26d^2x^2+14d^2x+5dx-25d+25)/(5-2d-4dx), free of B,b (sc:boundary-derivative)", ok)
# A10 balanced dual
ok = all((3 + 4 * D2) == (221 - 2 * Z0) / 33 for _ in [0]) and \
    all(-gQZ(QZ(B), QZ(Fr(1, 10)), D2, Fr(1, 2)) == (221 - 2 * Z0) / 33 * (B2 - B) for B in [Fr(87, 100), Fr(9, 10)])
check("A10 E_*-G_* = (221-2sqrt1081)/33 (B_2-B) at delta_2 for every b (sc:balanced-dual)", ok)
# A11 Bebel's PR103 fraction T_103 vs B_2
T103 = Fr(2804966284017490358030969257914123810640482798488248385236764587735806455694206308975383836842958115341512355256943573632654341,
          3211614609922119834694179379888743028634953588197852962423516288191911449231164741438582626063161096341458934143644085157872271)
dT = QZ(T103) - B2
ok = dT.sign() > 0 and (QZ(Fr(3, 10**254)) - dT).sign() > 0
mp.dps = 320
check("A11 0 < T_103 - B_2 < 3e-254 (l.191; T_103 = Bebel PR103 Lean bound)", ok, mpmath.nstr(dT.mpf(), 4))
mp.dps = 60
check("A12 B_2 < B_rat = 0.8733821 < B_tw(prev.) and B_2 > 0.87 (principal box)", (QZ(BRAT) - B2).sign() > 0 and
      BRAT < Fr(874058503, 10**9) and (B2 - Fr(87, 100)).sign() > 0)

# ----------------------------------------------------------------------------
print("== B. Optimality algebra (Thm sc:optimality, l.3998-4160) and unequal-factor model (l.4167-4265)")
d = D2
F = lambda dl: (25 - 18 * dl) / (30 - 24 * dl)
R2 = 1 - d * F(d)
tstar = (5 - 3 * d) / (5 - 4 * d)
ok = (1 - d * (Fr(1, 2) + tstar / 3)) == (Fr(1, 6) + (Fr(5, 6) - d) * tstar) == R2
mdx = (Fr(5, 6) - d * Fr(1, 2)) / (Fr(5, 2) - d - d)
ok &= (3 * mdx == tstar) and (Fr(1, 2) + mdx == F(d))
check("B1 t_*=(5-3d)/(5-4d), F=(25-18d)/(30-24d), R=1-dF: the two splits cross and Lemma sc:main-count attains it", ok)
check("B2 delta_2-5/12>0 (second split); 1-5delta_2/6 > R(delta_2); (R-1)/4+0.37/2>0",
      (d - Fr(5, 12)).sign() > 0 and (1 - 5 * d / 6 - R2).sign() > 0 and ((R2 - 1) / 4 + Fr(37, 200)).sign() > 0)
mu = Fr(37, 100)
okid = True
for _ in range(5):
    M, b, l = (Fr(random.randint(1, 999), 1000) for _ in range(3))
    lx = (M - b) / 2; ly = (M + b) / 2; h = 1 - lx + l
    L0 = Fr(5, 6) + M / 12 - l / 6; LK = Fr(1, 3) + 7 * M / 12 + l / 3; G0_ = 2 * b / 3 - M / 4
    rhs = (11 * (4 - 3 * mu) / 6 - (Fr(22, 3) - 25 * mu / 4) * L0 - (Fr(2, 3) + mu / 4) * LK
           - 3 * mu / 4 * (L0 + G0_))
    okid &= (h - ly - mu * h) == rhs
lowb = 11 * (4 - 3 * mu) / 6 - (8 - 21 * mu / 4) * Fr(8741, 10000)
check("B3 identity (sc:free-short-exclusion) and its bound 41671/12000000 at B<=0.8741", okid and lowb == Fr(41671, 12000000))
C = Fr(1, 8)
okH = True
for _ in range(5):
    M, b, l = (QZ(Fr(random.randint(1, 999), 1000)) for _ in range(3))
    lx = (M - b) / 2; ly = (M + b) / 2; h = 1 - lx + l; T = 2 * ly - h
    H1 = (1 + d) / 2 - ly / 2 + h * ((3 * R2 + 1) / 4 - Fr(1, 6)) + (d - 1) / 2 * l - C * d * T
    H2 = 1 + (2 * C + Fr(1, 2)) * d - (Fr(1, 2) + 2 * C * d) * M + (2 * C + Fr(1, 2)) * d * l
    okH &= (H1 == H2)
    Dp = 16 * C * d + 2 * d + 3
    u = (44 * C * d + 7 * d + 4) / (3 * Dp); v = (4 * C * d - d + 2) / (3 * Dp); w = 1 / Dp
    L0 = Fr(5, 6) + M / 12 - l / 6; LK = Fr(1, 3) + 7 * M / 12 + l / 3
    okH &= (u * L0 + v * LK + w * H2) == B2
okH &= (u.sign() > 0 and v.sign() > 0 and w.sign() > 0 and u + v + w == 1 and Dp == 4 * d + 3)
check("B4 (sc:free-high-plane) two forms agree; u,v,w>0, sum 1, uL_0+vL_K+wH = B_2 identically (sc:free-positive-dual)", okH)
okN = True
Dm = 19 - (12 * C + 9) * d
um = 3 * (2 - (4 * C + 3) * d) / Dm; vm = QZ(13) / Dm; km = (18 * C * d + 7 * d + 4) / (2 * Dm)
for _ in range(5):
    M, b, l = (QZ(Fr(random.randint(1, 999), 1000)) for _ in range(3))
    lx = (M - b) / 2; ly = (M + b) / 2; h = 1 - lx + l; T = 2 * ly - h
    H = 1 + (2 * C + Fr(1, 2)) * d - (Fr(1, 2) + 2 * C * d) * M + (2 * C + Fr(1, 2)) * d * l
    H0 = H + C * d * T
    LK = Fr(1, 3) + 7 * M / 12 + l / 3
    okN &= (um * LK + vm * H0 + km * T) == (52 + b * (44 * C * d + 7 * d + 4)) / (4 * Dm)
okN &= um.sign() > 0 and vm.sign() > 0 and km.sign() > 0 and um + vm == 1
check("B5 (sc:free-negative-dual) identity, positivity; 13/D_- = 0.8949386587587.. > B_2", okN and
      str((13 / Dm).mpf())[:15] == "0.8949386587587" and (13 / Dm - B2).sign() > 0, str((13 / Dm).mpf())[:16])
mD = (1 + Z0) / 90; rD = 2 * mD; tD = 3 * mD
ok = (mD == mdx) and (R2 == (281 - Z0) / 396)
Cf = lambda v: (1 - v) / 2
okm = True
for n in [QZ(Fr(1, 10)), QZ(Fr(2, 5)), QZ(Fr(7, 9))]:
    okm &= (R2 + 1 - d * mD - d * Cf(mD + n)) / 2 == (3 * R2 + 1 - d * Cf(2 * n)) / 4
for y1, y2 in [(mD, mD), (mD, QZ(Fr(1, 2))), (QZ(Fr(45, 100)), QZ(Fr(49, 100)))]:
    s = lambda y: d * (y + mD) / 4 if (y - mD).sign() >= 0 else d * y / 2
    okm &= (R2 + 2 * s(y1) + 2 * s(y2) + d * Cf(y1 + y2)) == 1
okm &= (R2 + d * rD + d * (1 - rD) / 2 == 1) and (R2 + d * tD == (1 + 5 * tD) / 6) and (R2 + 3 * d / 8 - 1).sign() < 0
okm &= (d - Fr(1, 2)).sign() < 0 and (d - KAP1).sign() < 0
check("B6 contact data m_D=(1+sqrt1081)/90, R_2=(281-sqrt1081)/396; neutrality A_mix=A_sep; saturation identities", ok and okm)

# ----------------------------------------------------------------------------
print("== C. Rational certificate and contact inequalities (Prop sc:rational-certificate, Lemma sc:double-contact, App. A)")
DOM = ((Fr(11, 50), Fr(13, 20)), (Fr(6, 25), Fr(1, 2)))
OMG = ((Fr(1, 8), Fr(5, 6)), (Fr(0), Fr(1, 2)))
D1 = ((Fr(1, 8), Fr(5, 6)), (Fr(0), Fr(6, 25)))
D2r = ((Fr(1, 8), Fr(11, 50)), (Fr(6, 25), Fr(1, 2)))
D3 = ((Fr(13, 20), Fr(5, 6)), (Fr(6, 25), Fr(1, 2)))
DD0 = ((Fr(2, 5), Fr(23, 50)), (Fr(6, 25), Fr(1, 2)))
DDm = ((Fr(11, 50), Fr(2, 5)), (Fr(6, 25), Fr(1, 2)))
DDp = ((Fr(23, 50), Fr(13, 20)), (Fr(6, 25), Fr(1, 2)))
E8, E5 = Fr(1, 10**8), Fr(1, 10**5)
Dl = lambda d, x: 5 - 2 * d - 4 * d * x


def entries():
    """(name, f(m), rect, multiplier, paper (d,L,k) rational, algebraic)."""
    return [
        ("G_L-E(R_d)-1e-8", lambda m: m["GL"] - m["E"](m["Rd"]) - E8, DOM, "D", (1, 2, 587138129), (1, 2, 586639331)),
        ("G_N-E(R_d)-1e-8", lambda m: m["GN"] - m["E"](m["Rd"]) - E8, DOM, "D", (13, 14, 194222), None),
        ("G_R-E(R_d)-1e-8", lambda m: m["GR"] - m["E"](m["Rd"]) - E8, DOM, "D", (1, 2, 1302783884), (1, 2, 1302190830)),
        ("S(R_d)-1/5", lambda m: m["S"](m["Rd"]) - Fr(1, 5), DOM, "D", (0, 1, 24051095581), (0, 1, 24051095581)),
        ("1-S(R_d)", lambda m: 1 - m["S"](m["Rd"]), DOM, "D", (0, 1, 416191181453), (0, 1, 416191181453)),
        ("1/250-E(R_d)", lambda m: Fr(1, 250) - m["E"](m["Rd"]), DOM, "D", (3, 4, 107812329), (3, 4, 107411143)),
        ("(p1-2dx)L-R_d", lambda m: (P1 - 2 * m["th0"] * P1) * LST - m["Rd"], DOM, "D", (0, 1, 81252174270), (0, 1, 81252174270)),
        ("-E(R_d)-hDsh(R_d)-1e-5", lambda m: -m["E"](m["Rd"]) - m["h"] * m["Dsh"](m["Rd"]) - E5, DOM, "D", (0, 1, 2018042319), (0, 1, 2017445082)),
        ("S+Dsh", lambda m: m["S"](m["Rd"]) + m["Dsh"](m["Rd"]), DOM, "D", (0, 1, 90495844449), (0, 1, 90495844449)),
        ("2-S-Dsh", lambda m: 2 - m["S"](m["Rd"]) - m["Dsh"](m["Rd"]), DOM, "D", (0, 1, 1013256914932), (0, 1, 1013256914932)),
        ("m_d-7/20", lambda m: m["md"] - Fr(7, 20), DOM, "D", (0, 1, 1073092960), (0, 1, 1073092960)),
        ("49/100-m_d", lambda m: Fr(49, 100) - m["md"], DOM, "D", (0, 1, 35084007848), (0, 1, 35084007848)),
        ("r_d-3/5", lambda m: m["rd"] - Fr(3, 5), DOM, "D", (0, 1, 96884964434), (0, 1, 96884964434)),
        ("1-r_d", lambda m: 1 - m["rd"], DOM, "D", (0, 1, 84314446897), (0, 1, 84314446897)),
        ("t_d-1", lambda m: m["td"] - 1, DOM, "D", (0, 1, 50588668138), (0, 1, 50588668138)),
        ("3/2-t_d", lambda m: Fr(3, 2) - m["td"], DOM, "D", (0, 1, 126471670345), (0, 1, 126471670345)),
        ("49/50+x/50-F_d", None, DOM, "D", (0, 1, 36638459651), (0, 1, 36638459651)),
        ("1/50-gamma_2(3/10)", lambda m: Fr(1, 50) - m["gam2"](Fr(3, 10)), DOM, "D", (0, 1, 1075495777), (0, 1, 1075495777)),
        ("J", lambda m: m["J"], OMG, "1", (0, 1, 648092328858), (0, 1, 648092328858)),
        ("S(R_s)", lambda m: m["S"](m["Rs"]), OMG, "J", (0, 1, 22188021055), (0, 1, 22188021055)),
        ("1-S(R_s)", lambda m: 1 - m["S"](m["Rs"]), OMG, "J", (0, 1, 250692959516), (0, 1, 250692959516)),
        ("(p1-2dx)L-R_s", lambda m: (P1 - 2 * m["th0"] * P1) * LST - m["Rs"], OMG, "J", (0, 1, 13927405530), (0, 1, 13927405530)),
        ("-E(R_s)-1e-5 on D1", lambda m: -m["E"](m["Rs"]) - E5, D1, "J", (9, 11, 15034272), (9, 11, 14627855)),
        ("-E(R_s)-1e-5 on D2", lambda m: -m["E"](m["Rs"]) - E5, D2r, "J", (0, 1, 463076646), (0, 1, 462618241)),
        ("-E(R_s)-1e-5 on D3", lambda m: -m["E"](m["Rs"]) - E5, D3, "J", (0, 1, 1764325203), (0, 1, 1763978740)),
    ]


def run_entry(name, f, rect, mult, pap, B, b, tag):
    def ff(d, x):
        dd, xx = (QZ(d), QZ(x)) if isinstance(B, QZ) else (d, x)
        m = model(B, b, dd, xx)
        if f is None:   # 49/50 + x/50 - F_d (x is a coordinate; F_d is the detector saving)
            val = Fr(49, 50) + xx / 50 - m["Fd"]
        else:
            val = f(m)
        mul = {"D": Dl(dd, xx), "J": m["J"], "1": 1}[mult]
        return val * mul
    C = poly_of(ff)
    if C is None:
        check("C %s [%s]: multiplied function is a polynomial" % (name, tag), False)
        return
    ok, dep, nacc, mmin = bnb_bern(C, rect)
    if mult == "D":
        Cm = poly_of(lambda d, x: Dl(d, x)); Mx = upper_bern(Cm, rect)
    elif mult == "J":
        Cm = poly_of(lambda d, x: model(B, b, QZ(d) if isinstance(B, QZ) else d,
                                        QZ(x) if isinstance(B, QZ) else x)["J"]); Mx = upper_bern(Cm, rect)
    else:
        Mx = Fr(1)
    k = math.floor(mmin / Mx * 10**12) if ok else None
    same = pap is not None and ok and (dep, nacc) == (pap[0], pap[1]) and k == pap[2]
    check("C %s [%s] > 0 on %s" % (name, tag, [[float(t) for t in iv_] for iv_ in rect]), ok,
          "ours (d,L,k)=(%s,%s,%s); paper %s%s" % (dep, nacc, k, pap, "; identical" if same else ""))


for nm, f, rect, mult, pr, pa in entries():
    run_entry(nm, f, rect, mult, pr, BRAT, bRAT, "B_rat")
for nm, f, rect, mult, pr, pa in entries():
    if pa is not None:
        run_entry(nm, f, rect, mult, pa, B2, QZ(b2), "B_2")

# contact entries (App. A second table)
def g2(d, x):
    return gQZ(B2, b2, QZ(d), QZ(x))
Cc = poly_of(lambda d, x: (5 - 4 * d) * gQZ(B2, b2, QZ(d), Fr(1, 2)) / (QZ(d) - D2) ** 2, 4, 0)
ok, dep, nacc, mmin = bnb_bern(Cc, ((DOM[0][0], DOM[0][1]), (Fr(0), Fr(1))))
k = math.floor(mmin / (5 - 4 * DOM[0][0]) * 10**12)
check("C contact g(B_2,d,1/2)/(d-d_2)^2 > 0.11308683 on [11/50,13/20]", ok and k == 113086830076,
      "k=%s (paper 113086830076)" % k)
Ct = poly_of(lambda d, x: (5 - 4 * d) * Dl(d, x) * (gQZ(B2, b2, QZ(d), QZ(x)) - gQZ(B2, b2, QZ(d), Fr(1, 2))) / (Fr(1, 2) - x))
ok, dep, nacc, mmin = bnb_bern(Ct, DD0)
Mx = upper_bern(poly_of(lambda d, x: (5 - 4 * d) * Dl(d, x)), DD0)
k = math.floor(mmin / Mx * 10**12)
check("C transverse (g(x)-g(1/2))/(1/2-x) > 0.00908267 on D_0=[2/5,23/50]x[6/25,1/2]", ok,
      "(d,L,k)=(%s,%s,%s), paper (0,1,9082672972)" % (dep, nacc, k))
for nm, rect, kp in [("D_-", DDm, 72884214), ("D_+", DDp, 147316073)]:
    Cg = poly_of(lambda d, x: Dl(d, x) * gQZ(B2, b2, QZ(d), QZ(x)))
    ok, dep, nacc, mmin = bnb_bern(Cg, rect)
    Mx = upper_bern(poly_of(lambda d, x: Dl(d, x)), rect)
    k = math.floor(mmin / Mx * 10**12)
    check("C g(B_2,b_2,.) > 0 on %s" % nm, ok, "(d,L,k)=(%s,%s,%s), paper (0,1,%s)" % (dep, nacc, k, kp))
Cb = poly_of(lambda d, x: -26 * d * d * x * x + 14 * d * d * x + 5 * d * x - 25 * d + 25)
ok, dep, nacc, mmin = bnb_bern(Cb, DOM)
k = math.floor(mmin / upper_bern(poly_of(lambda d, x: Dl(d, x)), DOM) * 10**12)
check("C dB g > 2.3723473 on D", ok and k >= 2372347300000, "k=%s (paper 2372347314201)" % k)
# uniqueness of the contact: on D_0, g(d,x) >= g(d,1/2) >= 0 with equality only at (d_2,1/2); also d_2 in D_0
check("C contact point (delta_2,1/2) lies in D_0 interior in delta; D_0 u D_- u D_+ = D",
      (D2 - Fr(2, 5)).sign() > 0 and (Fr(23, 50) - D2).sign() > 0)
# the cover of [1/8,5/6]x[0,1/2]
check("C cover D u D1 u D2 u D3 = [1/8,5/6]x[0,1/2] (sc:complement-cover)",
      DOM[0] == (D2r[0][1], D3[0][0]) and D1[1][1] == DOM[1][0] == D2r[1][0] and D1[0] == OMG[0])
# denominators
check("C Delta >= 12/5 on D; J >= 4373305/6747966 on Omega (App. A)",
      min(bernstein(poly_of(lambda d, x: Dl(d, x)), DOM)) >= Fr(12, 5) and
      min(bernstein(poly_of(lambda d, x: model(BRAT, bRAT, d, x)["J"]), OMG)) >= Fr(4373305, 6747966))

# ----------------------------------------------------------------------------
print("== D. Constant inequalities (App. A third table, l.4660-4700; Lemma sc:outside l.3640-3760)")
def consts(B, b):
    G = geom(B, b); l, lx, ly, h, T, C0 = G["l"], G["lx"], G["ly"], G["h"], G["T"], G["C0"]
    m0 = 1 - ly / h
    Qsm = l * (Fr(1, 2) + 2 * l / h) - 11 * h / 54
    return [
        ("b", b, 99956298200, 99956298225), ("l-1/6", l - Fr(1, 6), 6471600000, 6472133547),
        ("1/5-l", Fr(1, 5) - l, 26861733333, 26861199786), ("lx-l", lx - l, 190314450900, 190313650566),
        ("ly-l-11b/6", ly - l - 11 * b / 6, 107017535733, 107016735378), ("p1 l-3/4", P1 * l - Fr(3, 4), 28887424510, 28889824749),
        ("T", T, 117132482433, 117131148578), ("1/2-ly", Fr(1, 2) - ly, 36590984233, 36591250994),
        ("B-87/100", B - Fr(87, 100), 3382100000, 3381966613),
        ("l/(h+1e-9)-L", l / (h + Fr(1, 10**9)) - LST, 3833958930, 3834406521),
        ("2/9-l/(h-1/40)", Fr(2, 9) - l / (h - Fr(1, 40)), 1575535322, 1575080419),
        ("m0-3/10", m0 - Fr(3, 10), 127667918389, 127668813572), ("47/100-m0", Fr(47, 100) - m0, 42332081610, 42331186427),
        ("ly/h-57/100", ly / h - Fr(57, 100), 2332081610, 2331186427), ("59/100-ly/h", Fr(59, 100) - ly / h, 17667918389, 17668813572),
        ("-Q_sm", -Qsm, 4321129754, 4320642840), ("-1/250-C0-(1-h)/16", -Fr(1, 250) - C0 - (1 - h) / 16, 1023130452, 1022513534),
        ("-1/1000-C0-(1/2+3l/2)/50", -Fr(1, 1000) - C0 - (Fr(1, 2) + 3 * l / 2) / 50, 723635633, 722952688),
        ("-1/20-13h/75+ly/2-1/50", -Fr(1, 20) - 13 * h / 75 + ly / 2 - Fr(1, 50), 21359012706, 21358740601),
        ("h/600", h / 600, 1349475915, 1349477249), ("ly/20", ly / 20, 23170450788, 23170437450),
        ("73/1000-13/(40p1)", Fr(73, 1000) - 13 / (40 * P1), 756008254, 756008254),
        ("1/5-(73/1000)(37/40)-1/50-1/10", Fr(1, 5) - Fr(73, 1000) * Fr(37, 40) - Fr(1, 50) - Fr(1, 10), 12475000000, 12475000000),
    ]
bad = []
for (nm, v, kr, _), (_, w, _, ka) in zip(consts(BRAT, bRAT), consts(B2, QZ(b2))):
    k1 = math.floor(v * 10**12)
    k2 = int(mpmath.floor(tompf(w) * 10**12))
    if not (v > 0 and sgn(w) > 0 and k1 == kr and k2 == ka):
        bad.append((nm, k1, kr, k2, ka))
check("D1 all 23 constant entries positive at both parameter sets, k = floor(1e12 value) reproduced", not bad, str(bad)[:300])
G = geom(BRAT, bRAT)
check("D2 l,h,l_y,T,m_0 at B_rat match (sc:rational-parameters/-lengths)",
      G["l"] == Fr(1298537, 7500000) and G["h"] == Fr(8096855491, 10**10) and G["ly"] == Fr(13902270473, 3 * 10**10)
      and G["T"] == Fr(3513974473, 3 * 10**10) and 1 - G["ly"] / G["h"] == Fr(10388296000, 24290566473))
# direct majorant at k=0,G=b/12 equals B (both parameter sets)
okd = True
for B, b in [(QZ(BRAT), QZ(bRAT)), (B2, QZ(b2))]:
    l = Fr(11, 3) - 4 * B; lx = (1 - b - l) / 2; ly = (1 + b - l) / 2; M = lx + ly
    kk = [(1 + 3 * l - 2 * M) / 4, M + l - 1, (3 * l - 1) / 2]
    okd &= all(sgn(t) <= 0 for t in kk) and (ly - 11 * b / 6).sign() > 0
    okd &= (Fr(5, 6) + lx / 6 - l / 6 + b / 12) == B
check("D3 direct majorant (sc:direct-majorant): k=0, G=b/12, B_low = 11/12-l/4 = B exactly", okd)
# small-amplitude identity
okp = True
for dd in [Fr(i, 7) for i in range(5)]:
    for x in [Fr(j, 5) for j in range(6)]:
        A = Fr(5, 6) - dd; Dq = 3 - Fr(17, 9) * x; Pq = (2 - Fr(8, 9) * x) * (1 - x); J = A * Dq + dd * Pq
        v = Fr(1, 2) - x
        okp &= (2 * J - 3 * A * Pq - Fr(44, 27) * x * J ==
                Fr(250, 729) * x * x + dd / 243 * (352 * v ** 3 + 948 * v ** 2 + 984 * v + 199))
check("D4 small-amplitude polynomial identity (sc:small-amplitude-identity), 5x6 grid > degrees (1,3)", okp)
okr = True
for _ in range(200):
    dd = Fr(random.randint(1, 125), 1000); x = Fr(random.randint(0, 500), 1000)
    A = Fr(5, 6) - dd; Dq = 3 - Fr(17, 9) * x; Pq = (2 - Fr(8, 9) * x) * (1 - x); J = A * Dq + dd * Pq
    R = 1 - dd + A * dd * Pq / (2 * J)
    okr &= J > 0 and R <= 1 - 2 * dd / 3 - 22 * dd * x / 81
check("D5 hence R <= 1-2d/3-22dx/81 (200 random points, J>0); 'C0+(1-h)/16 < -1/250' margin 0.001023130452 is in D1", okr)
# row slope and principal box
check("D6 row slope 1/5-0.073*37/40-0.02 = 0.112475 > 1/10 (sc:row-slope)", Fr(1, 5) - Fr(73, 1000) * Fr(37, 40) - Fr(2, 100) == Fr(112475, 10**6))
s0, w0, z0 = Fr(87, 100), Fr(19, 20), Fr(33, 200)
D_, R_, V_, Y_ = s0, 6 * s0 + 6 * z0 - 4, 6 * z0, w0
tab = [V_ + Y_, D_ + Y_, D_ + V_, D_ + V_ + Y_, D_ + Y_ + V_ - 1, D_ + 2 * Y_ + V_ - 1, R_, D_ + R_, R_ + 2 * Y_, R_ + 2 * Y_ + V_, R_ + 1, R_ + Y_ + 1]
paper = [Fr(194, 100), Fr(182, 100), Fr(186, 100), Fr(281, 100), Fr(181, 100), Fr(276, 100), Fr(221, 100), Fr(308, 100), Fr(411, 100), Fr(510, 100), Fr(321, 100), Fr(416, 100)]
err = [s0, w0, 6 * z0, w0 + 6 * z0 - 1, 5 * s0 + 6 * z0 - 4]
check("D7 principal box (Lemma principal-box): 12 decay exponents and 5 error exponents as printed, min 1.81 > 1, min error 0.87",
      tab == paper and err == [Fr(87, 100), Fr(95, 100), Fr(99, 100), Fr(94, 100), Fr(134, 100)])
# principal factorisation at w=1, z=1/6 (exact, symbolic in D,R with q numeric)
okf = True
for q in [5, 7, 49, 101]:
    for (Dv, Rv) in [(Fr(1, 3), Fr(1, 7)), (Fr(-2, 9), Fr(3, 11))]:
        Yv = Vv = Fr(1, q)
        lhs = (-Vv * Yv + Dv * Yv + Dv * Vv * (1 - Yv) - Dv * (q - 1) * Yv * Vv * (1 - Yv) + Rv - Dv * Rv
               - Rv * Yv ** 2 * (1 - Vv) - Rv * (1 - Yv) / q)
        Hp = (1 - Fr(1, q * q)) * (1 - Dv + Dv / q - Rv / q) / ((1 - Rv) * (1 - Dv))
        okf &= lhs == (1 - Rv) * (1 - Dv) * (Hp - 1)
check("D8 (principal-polynomial) at w=1,z=1/6 gives (principal-factorization) exactly", okf)
# Appendix B: Sturm sequence of P(B) for B_flat
def pdiv(a, b):
    a = a[:]; q = [Fr(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        c = a[-1] / b[-1]; k = len(a) - len(b); q[k] = c
        for i in range(len(b)):
            a[i + k] -= c * b[i]
        a.pop()
    while a and a[-1] == 0:
        a.pop()
    return q, a
P = [Fr(26231), Fr(-167541), Fr(384534), Fr(-376596), Fr(133488)]
dP = [i * P[i] for i in range(1, 5)]
seq = [P, dP]
while len(seq[-1]) > 1:
    _, r = pdiv(seq[-2], seq[-1]); seq.append([-t for t in r])
ev = lambda p, t: sum(c * t ** i for i, c in enumerate(p))
V = lambda t: sum(1 for i in range(len(seq) - 1) if ev(seq[i], t) * ev(seq[i + 1], t) < 0)
signs = lambda t: "".join("+" if ev(p, t) > 0 else "-" for p in seq)
rows = [(Fr(8748, 10**4), "+--++", 2), (Fr(87488664703996, 10**14), "+--++", 2), (Fr(87488664704, 10**11), "---++", 1), (Fr(8749, 10**4), "---++", 1)]
check("D9 App. B Sturm table for B_flat (signs and V(t)); B_flat < 0.87488664704 < (1+kappa_1)/2 = 0.874887",
      all(signs(t) == s and V(t) == v for t, s, v in rows) and Fr(87488664704, 10**11) < (1 + KAP1) / 2 == Fr(874887, 10**6)
      and P1 == 6 * KAP1 and (1 + KAP1) / 2 - Fr(87488664704, 10**11) == Fr(35296, 10**11))

# ----------------------------------------------------------------------------
print("== E. Section 6 derivations (l.3195-3470)")
oke = True
for _ in range(8):
    d = Fr(random.randint(1, 830), 1000); x = Fr(random.randint(0, 500), 1000)
    m = model(BRAT, bRAT, d, x)
    FI = lambda r: x + (1 - x) * r
    F2 = lambda mm: x + 2 * (1 - x) * mm
    oke &= FI(m["rd"]) == F2(m["md"]) == m["Fd"] and m["Rd"] == Fr(1, 6) + m["A"] * m["td"]
    oke &= m["GN"] - m["E"](m["Rd"]) == gact(BRAT, bRAT, d, x)        # p_1 cancels
    oke &= m["GN"] - m["E"](m["Rd"]) == (lambda mm: mm["GN"] - mm["E"](mm["Rd"]))(model(BRAT, bRAT, d, x, p1=Fr(9, 2)))
    # old count crossing
    c = m["c"]; t = 1 + d * m["Pc"] / (2 * m["J"])
    mm_ = (x * (1 - c) + (1 - x) * t) / m["Dc"]
    oke &= (1 - d * (c * x + (2 - 2 * c * x) * mm_)) == m["Rs"] == Fr(1, 6) + (Fr(5, 6) - d) * t
    oke &= FI(t - mm_) == c * x + (2 - 2 * c * x) * mm_
    # numerator gain and contour derivatives
    for mmm in [Fr(3, 10), Fr(2, 5), Fr(47, 100)]:
        oke &= m["gam2"](mmm) == d * x / 2 * ((1 - 2 * mmm) / 2 - (1 - 2 * mmm) / P1)
    h, ly = m["h"], m["ly"]; m0 = 1 - ly / h; th0 = m["th0"]
    Am = lambda mm: h * m["gam2"](mm) + th0 * h * (m0 - mm)
    Ap = lambda mm: h * m["gam2"](mm) + (d / 2 - th0) * h * (mm - m0)
    oke &= (Am(1) - Am(0)) == -h * d * x / 2 and (Ap(1) - Ap(0)) == h * d * (1 - x) / 2
    oke &= Am(m0) == Ap(m0) == m["GN"]
    oke &= th0 * h * (m0 - Fr(3, 10)) == m["GL"] and (d / 2 - th0) * h * (Fr(47, 100) - m0) == m["GR"]
check("E1 detector balance r_d=2m_d, R_d=1/6+A t_d; g=G_N-E(R_d) free of p_1; old count R_s crossing; gains A_+-, G_N,G_L,G_R", oke)
# reference lemma algebra N(R)
okn = True
for _ in range(6):
    R = Fr(random.randint(1, 999), 1000); g_ = Fr(random.randint(1, 400), 1000); mm = Fr(random.randint(0, 500), 1000)
    th0 = g_ / P1; z = (1 - 2 * mm) / P1
    okn &= (3 * R + 1) / 4 - g_ * z / 2 - th0 * mm == (3 * R + 1) / 4 - th0 / 2
    mst = (1 - P1 * LST) / 2
    okn &= R + (Fr(1, 2) - th0) * mst - ((3 * R + 1) / 4 - th0 / 2) == (R - (P1 - 2 * g_) * LST) / 4
check("E2 reference Lemma sc:reference: N(R) identity and trivial-range identity (R-(p1-2g)L)/4", okn)

# ----------------------------------------------------------------------------
print("== F. Sections 3-4 exponent bookkeeping (exact)")
# F1 restoration optimisation: random search for a counter-instance
worst = -1e9
for _ in range(4000):
    K = random.uniform(0.5001, 0.75); m = random.uniform(0, 3); d = random.uniform(0, 3); q = random.uniform(0, 2)
    z = random.uniform(0, 1); pz = random.uniform(0, 1) * z
    A = d + z
    for u in [m * i / 400 for i in range(401)]:
        s = (m - u) / 6
        e1 = s + K * z + min(max(d, (1 - K) * u + K * d), u + q)
        b1 = max((1 - K) * m + K * A, (m + q) / 6 + 5 * A / 6)
        e2 = s + min(max(A, (1 - K) * u + K * A), u + q + pz)
        b2_ = max((1 - K) * m + K * A, (m + q + pz) / 6 + 5 * A / 6)
        worst = max(worst, e1 - b1, e2 - b2_)
check("F1 Thm op:restoration suprema (whole-supremum/block-supremum): no counter-instance in 4000x401 samples", worst <= 1e-12, "max excess %.2e" % worst)
# F2 recursion of Prop op:squarefree
okq = True
for K in [Fr(501, 1000), Fr(6, 10), Fr(3, 4), Fr(9, 10)]:
    q = Fr(2)
    for _ in range(6):
        a = q + 2 * K - 1
        okq &= a / K >= 1                       # L = H^{a/K} >= H
        qn = (1 - K) / K * (q + 2 * K - 1)
        okq &= qn - (1 - K) == (1 - K) / K * (q - (1 - K))
        q = qn
    okq &= (1 - K) / K < 1
    # reflection budget exponents
    for qq in [1 - K, Fr(2)]:
        okq &= (1 - qq - 2 * K) <= -K and K - 1 < 0
check("F2 Prop op:squarefree recursion q->(1-K)(q+2K-1)/K contracts to 1-K; L>=H; reflection budget signs", okq)
# F3 25-cell table and natural-near exponent bounds
okt = all((((i - j) % 6) == [[0, 5, 4, 3, 2], [1, 0, 5, 4, 3], [2, 1, 0, 5, 4], [3, 2, 1, 0, 5], [4, 3, 2, 1, 0]][i - 1][j - 1])
          and ((i == j) == ((i - j) % 6 == 0)) for i in range(1, 6) for j in range(1, 6))
for _ in range(2000):
    Kf = random.uniform(0.5001, 0.999)
    Gg = math.exp(random.uniform(0, 8)); Ff = Gg * random.uniform(0, 1); GL = Gg * random.uniform(0, 1); GR = Gg * random.uniform(0, 1)
    okt &= Ff ** (Kf - 0.5) * math.sqrt(GL * GR) / Gg ** (2 * Kf) <= Gg ** (0.5 - Kf) * (1 + 1e-12) <= 1 + 1e-12
    okt &= math.sqrt(Ff * GL * GR) / Gg ** 2 <= 1 + 1e-12
check("F3 25-cell digit table = (i-j) mod 6, off-diagonal nonzero; natural-near nonlinear/linear factor bounds (2000 random)", okt)
# F4 local valuation tables for all valuations <= 60
ok4 = True
for i in range(1, 61):
    for j in range(1, i + 1):
        r = 1 if (i - j) % 6 else 0
        avail = i + j - 2 - r
        pos = Fr(max(0, 3 * i - 5 * j - r), 6)
        ok4 &= pos <= avail
        if 3 * i - 5 * j - r > 0:
            ok4 &= 6 * avail - (3 * i - 5 * j - r) == 3 * i + 11 * j - 12 - 5 * r
# second-sign table (s + lambda = 1, 0<=t<=s<=1) via the increment formula
for s in [Fr(k, 10) for k in range(11)]:
    lam = 1 - s
    for t in [Fr(k, 10) for k in range(int(s * 10) + 1)]:
        for i in range(1, 61):
            if i % 6:
                inc_u = 1 + i - 1 - i - lam * i - s * (i + 0) + t * (1 + 0)
                inc_n = 1 + i - 0 - i - lam * i - s * (i + 1) + t * (0 + 1)
                ok4 &= inc_u == -i + t and inc_n == 1 - i - s + t and inc_u <= 0 and inc_n <= 0
            else:
                inc = 1 + i - 0 - i - lam * i - s * i
                ok4 &= inc == 1 - i and inc <= 0
        for j in range(6, 61, 6):
            for i in range(j + 1, 61):
                for istar in (i, j):
                    inc = 1 + j - Fr(i + j, 2) - lam * istar - s * j
                    tgt = 1 - j - (Fr(1, 2) + lam) * (i - j) if istar == i else 1 - j - Fr(i - j, 2)
                    ok4 &= inc == tgt and inc <= 0
# F_2 table >= 2 cbar/3
for i in range(1, 61):
    if i % 6:
        ok4 &= Fr(7 * i, 6) >= Fr(2 * i, 3)
        ok4 &= Fr(7 * i - 5, 6) + (Fr(1, 3) if i == 1 else 0) >= Fr(2 * i, 3)
        ok4 &= (2 * i - Fr(5 * i, 6) - 1 + 1) == Fr(7 * i, 6)
        ok4 &= (2 * i - Fr(5 * i, 6) - 1 + Fr(1, 6) + Fr(2 * (1 if i == 1 else 0), 6)) == Fr(7 * i - 5, 6) + (Fr(1, 3) if i == 1 else 0)
    else:
        ok4 &= Fr(7 * i, 6) - 1 >= Fr(2 * i, 3) and (2 * i - Fr(5 * i, 6) - 1) == Fr(7 * i, 6) - 1
for j in range(6, 61, 6):
    for i in range(j + 1, 61):
        ok4 &= (i + j) - Fr(5 * j, 6) - 1 == i + Fr(j, 6) - 1 and i + Fr(j, 6) - 1 >= Fr(i + j, 3)
check("F4 local tables for all valuations <= 60: first-brace (3i+11j-12-5r), second-sign table, F_2 >= 2cbar/3", ok4)
# F5 pool errors and parent pool
ok5 = True
for s in [Fr(k, 20) for k in range(21)]:
    for t in [Fr(k, 20) for k in range(int(s * 20) + 1)]:
        lam = 1 - s
        sig = Fr(3)  # sigma (units: sigma/3 = 1)
        P_ = sig / 3
        main = (1 - s) * 2 * sig + (1 - lam - 2 * s) * 0 + (t - 1) * 0 - sig / 3
        ok5 &= main == (Fr(5, 6) - s) * 2 * sig
        for (i, e, wa, wao) in [(1, 1, P_, P_), (6, 0, 6 * P_, 0), (7, 1, 7 * P_, P_)]:
            err = (1 - s) * sig + (1 - lam - 2 * s) * wa + (t - 1) * wao - 0
            ok5 &= err - main == sig / 3 * (-2 + 3 * s - s * i + (t - 1) * e) and err - main <= 0
for chi in [Fr(k, 50) for k in range(1, 50)]:
    for B0 in [Fr(k, 20) for k in range(1, 21)]:
        ok5 &= -chi - 2 * chi * B0 < 0 and -5 - 6 * chi + 9 * chi * B0 < 0 and -6 - 7 * chi + 10 * chi * B0 < 0
        # parent first/second increments
        for c1 in [Fr(1), Fr(3)]:
            for d1 in [Fr(0), Fr(1)]:
                if c1 >= d1:
                    ok5 &= -chi * c1 + chi * B0 * (c1 - d1) == -chi * ((1 - B0) * c1 + B0 * d1) <= 0
        for i in range(1, 30):
            for tu in (0, 1):
                ok5 &= 1 - chi * B0 * tu - (1 + chi * (1 - B0)) * i <= 0
check("F5 auxiliary pool: main (5/6-s)Gamma, errors minus main = (sigma/3)(-2+3s-si+(t-1)e) <= 0; parent pool/first/second increments", ok5)
# F6 marked-product identities (sec mom:marked), exact affine identities
names = "d z w wf f0 f1 gg v cc T1 T e1 e S1 S zc mu B B0 L L0".split()
okm = True
for _ in range(12):
    X = {n: Fr(random.randint(0, 999), 997) for n in names}
    d, z, w, wf, f0, f1, gg, v, cc = (X[k] for k in "d z w wf f0 f1 gg v cc".split())
    T1, T, e1, e, S1, S, zc, mu, B, B0, L, L0 = (X[k] for k in "T1 T e1 e S1 S zc mu B B0 L L0".split())
    uf = gg + v + cc; w = wf + f0 + f1 + uf; tl = f0 + uf
    L1 = f1 + 4 * T1 + e1; Lp = v + 4 * T + e
    As = d + z - f1 - tl - L1 - Lp - 5 * cc
    Ms = 1 - 6 * tl - f1 - L1 - Lp - 5 * cc + T1 + 5 * v + T
    okm &= Ms == 1 - 6 * w + 6 * wf + 4 * (f1 + v) - 3 * (T1 + T) - e1 - e - 5 * cc
    ext = -f0 - L1 + T1 - S1 - uf - Lp + v - 4 * cc + T - S
    E0 = -w + wf - 3 * (T1 + T) - 4 * cc - S1 - S - e1 - e
    okm &= ext == E0
    DB0 = As + (L0 - 1) * (z - w - zc) - B0 * Ms
    eta0 = B - d - L * z
    Dref = (L0 - L) * z + (6 * B0 - L0) * w + B - B0
    okm &= DB0 - mu == (Dref - eta0 - 6 * B0 * wf - (1 + 4 * B0) * (f1 + v) - (4 - 3 * B0) * (T1 + T)
                        - (1 - B0) * (e1 + e + 5 * cc) - (L0 - 1) * zc - (mu - wf))
    for K in [Fr(1, 2), Fr(3, 4), X["B"]]:
        qd = f1 + v + cc
        okm &= E0 + K * (As - Ms) - (1 - K) * qd == (K * (d + z - 1) + (5 * K - 1) * tl - (2 - K) * f1 - (1 + 4 * K) * v
                                                    - (3 + K) * (T1 + T) - (5 - K) * cc - S1 - S - e1 - e)
    okm &= E0 + Fr(5, 6) * (As - Ms) == (Fr(5, 6) * (d + z - 1) + Fr(19, 6) * tl - f1 - Fr(25, 6) * v - Fr(23, 6) * (T1 + T)
                                        - 4 * cc - S1 - S - e1 - e)
    chi = Fr(random.randint(1, 99), 100)
    lhs = E0 + chi * (DB0 - mu) - (-w + chi * Dref)
    rhs = ((1 - 6 * B0 * chi) * wf - chi * (1 + 4 * B0) * (f1 + v) + (-3 - chi * (4 - 3 * B0)) * (T1 + T)
           + (-1 - chi * (1 - B0)) * (e1 + e) + (-4 - 5 * chi * (1 - B0)) * cc - (S1 + S) - chi * (L0 - 1) * zc
           - chi * (mu - wf) - chi * eta0)
    okm &= lhs == rhs
check("F6 marked identities: M_* expansion, E_0 = exterior, D_B0-mu (marked-defect), coefficient table, high-K/high-R planes", okm)
# F7 face-step and curve numerics
ok7 = True
for a in [Fr(1, 100), Fr(1, 200)]:
    for L in [Fr(2), Fr(201, 100), Fr(5, 2)]:
        ok7 &= Fr(1, 6) / (L - 1 - 3 * a) <= Fr(50, 291)
ok7 &= Fr(5, 6) - (4 + Fr(3, 100)) * Fr(50, 291) == Fr(41, 291)
ok7 &= Fr(5, 4) - 1 + (2 - 3 - Fr(9, 200)) * Fr(50, 291) == Fr(41, 582)
ok7 &= (Fr(1, 2) * (1 - Fr(5, 100)) + Fr(1, 100) >= Fr(485, 1000)) and (5 - Fr(3, 4) - (19 - Fr(3, 4)) / 100 == Fr(40675, 10**4))
for X_ in [Fr(k, 7) for k in range(8)]:
    for cc in [Fr(k, 9) for k in range(5)]:
        for l in [Fr(k, 11) for k in range(6)]:
            Phi = -Fr(9, 8) * X_ - Fr(29, 6) * cc - Fr(2, 3) * max(Fr(0), l - X_ - 5 * cc)
            ok7 &= Phi + 2 * l / 3 <= 0
for delta in [Fr(1, 2000), Fr(1, 4000), Fr(1, 10**6)]:
    B = 1 - delta; k = Fr(1, 2) + delta; a = Fr(1, 100)
    Ls = 6 * B * k / (1 + k)
    ok7 &= 6 * B * k - 2 * (1 + k) == delta - 6 * delta ** 2 and 2 < Ls < 3 and 3 <= 6 * B * k
    Ln = Fr(3)
    for _ in range(40):
        rho = a - k / Ln * (Ln - (Ln - a / 2 * ((1 + k) / k * Ln - 6 * B)) + (6 * B - Ln) * a)
        ok7 &= rho == a / 2 * (1 + k - 6 * B * k / Ln) and rho > 0 and 1 - k * (6 * B - Ln) / Ln > 0
        Lnew = Ln - a / 2 * ((1 + k) / k * Ln - 6 * B)
        ok7 &= Lnew - Ls == (1 - a * (1 + k) / (2 * k)) * (Ln - Ls) and Ls < Lnew < Ln <= 6 * B * k
        Ln = Lnew
    for kk in [k]:
        g_ = (2 * kk - 1) / (2 * (6 * kk - 1))
        ok7 &= Fr(1, 6) - 1 / (3 * (6 * kk - 1)) == g_ > 0
ok7 &= all(-Fr(2, 3) * v - max(Fr(0), Fr(1, 4) - v) <= -Fr(1, 6) for v in [Fr(k, 100) for k in range(101)])
check("F7 face-step constants 50/291, 41/291, 41/582, 0.485, 4.0675, Phi+2l/3<=0; curve rho_n, L_n->L_*, gamma_k; sup=-M/6", ok7)

# F8 returned-family identities (l.1420-1500) and Gauss-terminal payments (l.2190-2265)
ok8 = True
nm8 = "A c1 d1 R E m q K wa wao la g g2 t2 V cb2 p2 s0 BC Rl Ast cst dst u1 mst tau".split()
for _ in range(12):
    X = {n: Fr(random.randint(0, 999), 991) for n in nm8}
    A, c1, d1, R_, E_, m, q, K = (X[k] for k in "A c1 d1 R E m q K".split())
    wa, wao, la, g, g2, t2, V, cb2, p2, s0, BC = (X[k] for k in "wa wao la g g2 t2 V cb2 p2 s0 BC".split())
    M = m + q; qt = q + R_ + E_
    K0 = 2 * A - c1 - d1 + R_ + E_ - m
    a0 = A - c1 - wa; D0 = d1 + (K0 - K); J0 = D0 - c1 - 2 * wa + wao
    LamC = a0 + qt + wao + wa + BC - s0
    mp_ = 2 * a0 - K - g - g2 - V; qp = qt + wao + t2 + V; Mp = mp_ + qp
    lam0 = K + g - la - a0 + p2 - s0 + g2 - t2 - cb2
    ok8 &= Mp == M + J0 - g - g2 + t2
    ok8 &= LamC - lam0 == Mp + (cb2 - p2 + wa + BC + la)
    # Gauss-terminal (mom:gauss-target-data, gauss-linear/quadratic-payment)
    Ast, cs, ds, u1 = X["Ast"], X["cst"], X["dst"], X["u1"]
    mst = X["mst"]; qs = X["q"]; Mst = mst + qs
    al = Ast - cs; K0s = 2 * Ast - cs - ds + u1 - mst; Lams = al + qs + u1 + BC - s0
    Rl = X["Rl"]
    ok8 &= (Rl + K0s - s0) - Lams + s0 - s0 == Rl + Ast - Mst - ds - BC     # first term at K=K0
    ok8 &= Rl + 2 * al - 2 * s0 - Lams == Rl + Ast - Mst + mst - cs - u1 - BC - s0
    dh = min(ds, Ast); mf = mst - cs + dh - u1; K0h = K0s + ds - dh
    ok8 &= K0h == 2 * al - mf and max(K0h - s0, 2 * al - 2 * s0) == K0h - s0 + max(Fr(0), mf - s0)
# origin payment (mom:gauss-origin-payment) from the marked identities
for _ in range(12):
    X = {n: Fr(random.randint(0, 999), 997) for n in "d z wf f0 f1 gg v cc T1 T e1 e S1 S".split()}
    d, z, wf, f0, f1, gg, v, cc, T1, T, e1, e, S1, S = (X[k] for k in "d z wf f0 f1 gg v cc T1 T e1 e S1 S".split())
    uf = gg + v + cc; w = wf + f0 + f1 + uf; tl = f0 + uf; L1 = f1 + 4 * T1 + e1; Lp = v + 4 * T + e
    As = d + z - f1 - tl - L1 - Lp - 5 * cc; Ms = 1 - 6 * tl - f1 - L1 - Lp - 5 * cc + T1 + 5 * v + T
    Rloc = -f0 - L1 + T1 - S1 - uf - Lp + v - 4 * cc + T - S
    ok8 &= Rloc + As - Ms == (d + z) - 1 + 4 * tl - L1 - S1 - Lp - 4 * v - 4 * cc - S <= (d + z) - 1 + 4 * w
check("F8 returned family: M'=M+J_0-g-g_2+t_2, Lambda_C-lambda_0=M'+Delta' (mom:full-allowance); Gauss-terminal payments; origin payment", ok8)
# F9 Lemma mom:low-prefix: prefix blocks of the sorted list have slope <= total slope
ok9 = True
for _ in range(3000):
    J = random.randint(1, 12)
    zs = [random.uniform(0.001, 0.05) for _ in range(J)]; as_ = [random.uniform(-1, 1) for _ in range(J)]
    tot = sum(a * zz for a, zz in zip(as_, zs)) / sum(zs)
    order = sorted(range(J), key=lambda i: as_[i])
    cum_a = cum_z = 0.0
    for i in order:
        cum_a += as_[i] * zs[i]; cum_z += zs[i]
        ok9 &= cum_a / cum_z <= tot + 1e-12
check("F9 Lemma mom:low-prefix (sorted prefix slopes <= total slope): 3000 random factor lists", ok9)

# F10 width cap: the Lean Dirichlet faces carry a cap c < 1/6; with widths capped at 3/20 on
# the numerator band [0.30,0.35] the theta=0 saving is h g (3/40) + h theta_0 (m_0-1/2), and
# (that) - E(R_d) = g + delta x (3h/40 - T/4); so the cap is not load-bearing for B_2.
ok10 = True
for B, b in [(QZ(BRAT), QZ(bRAT)), (B2, QZ(b2))]:
    for _ in range(6):
        d = QZ(Fr(random.randint(220, 650), 1000)); x = QZ(Fr(random.randint(240, 500), 1000))
        m = model(B, b, d, x)
        h, ly, T, th0 = m["h"], m["ly"], m["T"], m["th0"]
        m0 = 1 - ly / h
        for mm in [QZ(Fr(30, 100)), QZ(Fr(33, 100)), QZ(Fr(35, 100))]:
            gcap = d * x / 2 * (Fr(3, 20) - (1 - 2 * mm) / P1)
            Acap = h * gcap + th0 * h * (m0 - mm)
            ok10 &= Acap - m["E"](m["Rd"]) == gQZ(B, b, d, x) + d * x * (3 * h / 40 - T / 4)
    G = geom(B, b)
    ok10 &= (3 * G["h"] / 40 - G["T"] / 4).sign() > 0
check("F10 a width cap 3/20 (< 1/6, as in the Lean Dirichlet faces) would not change B_2: capped band margin = g + dx(3h/40-T/4), 3h/40-T/4 > 0", ok10)

# F11 contour-lemma constants (l.3027-3125) and lower rows (sc:lower-rows)
s0_ = Fr(1, 7)
ok11 = (Fr(1, 2) - Fr(13, 20) / 2 == Fr(7, 40)) and s0_ < Fr(7, 40) and Fr(33, 40) < 1 - s0_
ok11 &= ((1 + Fr(11, 50)) / 2 == Fr(61, 100)) and ((1 + Fr(13, 20)) / 2 == Fr(33, 40)) and Fr(1, 2) <= Fr(61, 100) and Fr(33, 40) <= Fr(11, 12)
ok11 &= 1 - Fr(1, 2) - s0_ >= Fr(13, 40)                       # theta range of profiles-lemma covers delta/2
ok11 &= Fr(1, 250) - Fr(1, 40) * Fr(1, 5) == -Fr(1, 1000)
check("F11 contour constants: Re w in [7/40,1/2], 1/7 < 7/40, 33/40 < 6/7, a in [61/100,33/40] within NF3 Lemma 13.10's [1/2,11/12]; lower rows -1/1000", ok11)

# ----------------------------------------------------------------------------
print("== G. Counter-instance search: Gauss-mean Lemma mom:gauss-mean over Z (Jacobi symbols, all odd moduli)")
def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def gauss_table(u):
    w = [cmath.exp(2j * math.pi * k / u) for k in range(u)]
    chi = [jacobi(x, u) for x in range(u)]
    nz = [x for x in range(u) if chi[x]]
    return [sum(chi[x] * w[(hh * x) % u] for x in nz) / math.sqrt(u) for hh in range(u)]


# exact correlation facts used in the proof: F(u,v;0)=1_{u=v}phi(u), |F(u,v;j)| <= (u,v) 1_{(u,v)|j}
okc = True
GT = {}
for u in range(3, 46, 2):
    GT[u] = gauss_table(u)
phi = lambda n: sum(1 for k in range(1, n + 1) if math.gcd(k, n) == 1)
for u in range(3, 46, 2):
    for v in range(3, 46, 2):
        n = u * v
        for jj in range(0, n, max(1, n // 25)):
            Fv = sum(GT[u][hh % u] * GT[v][hh % v].conjugate() * cmath.exp(-2j * math.pi * jj * hh / n) for hh in range(n)) / math.sqrt(u * v)
            c = math.gcd(u, v)
            bound = c if jj % c == 0 else 0
            if jj == 0:
                okc &= abs(Fv - (phi(u) if u == v else 0)) < 1e-6
            okc &= abs(Fv) <= bound + 1e-6
check("G1 F(u,v;0)=1_{u=v}phi(u) and |F(u,v;j)| <= gcd 1_{gcd|j} for odd u,v<46 (incl. prime powers)", okc)
# Gauss-mean: ratio LHS/RHS for random and adversarial coefficients
worst = 0
Dlist = [40, 80, 120]
for Dm in Dlist:
    us = [u for u in range(Dm, 2 * Dm) if u % 2]
    tabs = {u: gauss_table(u) for u in us}
    for H in [Dm // 4, Dm, 4 * Dm, Dm * Dm // 4]:
        for trial in range(4):
            if trial == 0:
                Bc = {u: 1.0 for u in us}
            elif trial == 1:
                Bc = {u: random.choice([-1.0, 1.0]) for u in us}
            else:  # adversarial: align phases at a chosen frequency h0
                h0 = random.randint(1, H)
                Bc = {u: (tabs[u][h0 % u].conjugate() / abs(tabs[u][h0 % u]) if abs(tabs[u][h0 % u]) > 1e-9 else 0) for u in us}
            lhs = 0.0
            for hh in range(1, H + 1):
                s = sum(Bc[u] * tabs[u][hh % u] for u in us)
                lhs += abs(s) ** 2 / Dm
            rhs = H / Dm * sum(abs(Bc[u]) ** 2 for u in us) + sum(abs(Bc[u]) for u in us) ** 2
            worst = max(worst, lhs / rhs)
check("G2 Gauss-mean (mom:gauss-l1l2) over Z: LHS/RHS bounded (no power growth) for D=40..120, H=D/4..D^2/4, adversarial B",
      worst < 3, "max ratio %.3f" % worst)

# ----------------------------------------------------------------------------
el = time.time() - T0
mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
mem_mb = mem / 1e6 if sys.platform == "darwin" else mem / 1e3
nf = sum(1 for _, o in RESULTS if not o)
print("SUMMARY: %d checks, %d failed; runtime %.1f s; peak RSS %.0f MB; single thread" % (len(RESULTS), nf, el, mem_mb))
sys.exit(1 if nf else 0)
