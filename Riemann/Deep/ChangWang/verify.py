#!/usr/bin/env python3
"""Re-check the numerical and algebraic claims of Chang-Wang, arXiv 2611.00924v1.

Stdlib + mpmath only.  Exact identities use fractions.Fraction (random rational
points for rational-function identities in several variables, exact polynomial
arithmetic for one-variable identities).  Range claims use mpmath interval
arithmetic (mpmath.iv, outward rounding) with branch-and-bound subdivision, so a
PASS on a range claim is a proof up to the correctness of mpmath.iv.

Line numbers "CW l.N" refer to main.tex of the arXiv e-print 2611.00924v1;
"[H]" is Liu, Hexagon 2610.00655v1 (number_field_new_v1.tex).

Run:  python3 verify.py         (single-threaded, ~6 min, ~30 MB)
"""
import random
import sys
from fractions import Fraction as Fr

import mpmath
from mpmath import mp, iv, mpf

mp.dps = 50
iv.dps = 30
random.seed(20261105)

RESULTS = []


def ivfr(q):
    """Rigorous interval enclosure of a Fraction."""
    q = Fr(q)
    return iv.mpf(q.numerator) / q.denominator


def ivspan(lo, hi):
    return iv.mpf([ivfr(lo).a, ivfr(hi).b])


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("  [" + detail + "]") if detail else ""))
    sys.stdout.flush()


# ---------------------------------------------------------------- polynomials
def padd(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]


def pmul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def pscale(p, c):
    return [c * a for a in p]


def ptrim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def peval(p, x):
    s = 0
    for a in reversed(p):
        s = s * x + a
    return s


# p(r), D(r), K0(r), Num(r)  (coefficient lists, constant term first)
P_CUBIC = [61560, -11376, -376469, 218582]
D_R = [1800, -5011, 2286]
K0_R = [720, -1837, 600]
NUM_R = [2520, -7127, 3410]
# n_b, n_e, n_0 = alpha(r) * k + beta(r)
NB = ([-5400, 15033, -6858], [720, -1837, 600])
NE = ([-40500, 34995, -702], [5400, -3775, -408])
N0 = ([8100, -10887, 4806], [-1080, 1403, -696])


def nfun(pair, r, k):
    return peval(pair[0], r) * k + peval(pair[1], r)


def dnfun_dr(pair, r, k):
    al, be = pair
    dal = [i * al[i] for i in range(1, len(al))]
    dbe = [i * be[i] for i in range(1, len(be))]
    return peval(dal, r) * k + peval(dbe, r)


def section_A():
    print("\n== A. Endpoint arithmetic (CW l.37-58, l.439-452)")
    p = lambda r: peval(P_CUBIC, r)
    dp = lambda r: 3 * 218582 * r * r - 2 * 376469 * r - 11376
    a, b = Fr(9, 20), Fr(23, 50)
    # p' is convex (positive leading coeff) so max on [a,b] is at an endpoint
    check("A1 p(9/20)>0>p(23/50) and p'<0 on [9/20,23/50] (unique root)",
          p(a) > 0 > p(b) and dp(a) < 0 and dp(b) < 0)
    rs = mpmath.findroot(lambda r: peval(P_CUBIC, r), mpf("0.4505707467"))
    Dr = peval(D_R, rs)
    ls = peval(NUM_R, rs) / Dr
    Bs = mpf(11) / 12 - ls / 4
    ks = mpf(5) / 6 - ls / 2
    bs = -(ls * dnfun_dr(NE, rs, ks) + dnfun_dr(N0, rs, ks)) / dnfun_dr(NB, rs, ks)
    # [H] closed form for b_*
    bH = (ls ** 2 * (11665 - 468 * rs) + ls * (4528 * rs - 20554) + 5113 - 4412 * rs) / (
        ls * (4572 * rs - 5011) + 7127 - 6820 * rs)
    printed = {
        "r": ("0.45057074672406073136", rs), "l": ("0.16879343589254482845", ls),
        "b": ("0.09179235229895088958", bs), "k": ("0.74893661538706091910", ks),
        "B": ("0.87446830769353045955", Bs)}
    ok = all(abs(mpf(s) - v) < mpf(10) ** -19 for s, v in printed.values())
    check("A2 r_*, l_*, b_*, k_*, B_* agree with [H] printed 20 digits", ok,
          "B_*=" + mpmath.nstr(Bs, 25))
    check("A3 CW b_* (derivative definition) equals [H] closed form", abs(bs - bH) < mpf(10) ** -40)
    B = Bs - mpf(1) / 160000
    check("A4 B = B_* - 1/160000 = 0.87446205769353...", mpmath.nstr(B, 14).startswith("0.87446205769353"),
          mpmath.nstr(B, 25))
    check("A5 B(l_*+1/40000) = B_* - 1/160000 exactly (B(l)=11/12-l/4)",
          abs((mpf(11) / 12 - (ls + mpf(1) / 40000) / 4) - B) < mpf(10) ** -45)
    # Taylor expansion at .45057 (CW l.441-445)
    r0 = Fr(45057, 100000)
    c0 = p(r0)
    c1 = dp(r0)
    c2 = Fr(3 * 218582) * r0 - 376469
    printed_t = [Fr("0.162414095686326"), Fr("-217502.0989101246"), Fr("-81009.52478")]
    ok = (abs(c0 - printed_t[0]) < Fr(1, 10 ** 15) and abs(c1 - printed_t[1]) < Fr(1, 10 ** 10)
          and abs(c2 - printed_t[2]) < Fr(1, 10 ** 5))
    check("A6 Taylor coefficients of p at .45057", ok,
          "exact: %s, %s, %s" % (mpmath.nstr(mpf(c0.numerator) / c0.denominator, 16),
                                  mpmath.nstr(mpf(c1.numerator) / c1.denominator, 16), float(c2)))
    check("A7 p(.45057+74e-8)>0>p(.45057+75e-8)",
          p(r0 + Fr(74, 10 ** 8)) > 0 > p(r0 + Fr(75, 10 ** 8)))
    # rigorous enclosures on the root box (CW l.446-452)
    # hull of the interval images over 512 sub-boxes of the root box (removes overestimation)
    lo_r, hi_r = Fr(45057074, 10 ** 8), Fr(45057075, 10 ** 8)
    parts = []
    NSUB = 512
    for i in range(NSUB):
        R = ivspan(lo_r + (hi_r - lo_r) * i / NSUB, lo_r + (hi_r - lo_r) * (i + 1) / NSUB)
        lI = peval(NUM_R, R) / peval(D_R, R)
        kI = iv.mpf(5) / 6 - lI / 2
        dE, d0, dB = dnfun_dr(NE, R, kI), dnfun_dr(N0, R, kI), dnfun_dr(NB, R, kI)
        bI = -(lI * dE + d0) / dB
        parts.append((lI, bI, dE, d0, dB))
    hull = lambda j: iv.mpf([min(p[j].a for p in parts), max(p[j].b for p in parts)])
    lI, bI, dE, d0, dB = (hull(j) for j in range(5))
    check("A8 .16879<l_*<.16880 (interval)", lI.a > mpf("0.16879") and lI.b < mpf("0.16880"), str(lI))
    check("A9 .09175<b_*<.09185 (interval)", bI.a > mpf("0.09175") and bI.b < mpf("0.09185"), str(bI))
    check("A10 d_r n_e in (21592.47,21592.66), d_r n_0 in (-4134.33,-4134.28), d_r n_b in (5333.96,5334.03)",
          dE.a > mpf("21592.47") and dE.b < mpf("21592.66") and d0.a > mpf("-4134.33")
          and d0.b < mpf("-4134.28") and dB.a > mpf("5333.96") and dB.b < mpf("5334.03"),
          "%s %s %s" % (dE, d0, dB))
    # l_* numerator/denominator in t = r - .45 (CW l.450)
    t45 = Fr(45, 100)
    Dsh = [peval(D_R, t45), peval([D_R[1], 2 * D_R[2]], t45), D_R[2]]
    Nsh = [peval(NUM_R, t45), peval([NUM_R[1], 2 * NUM_R[2]], t45), NUM_R[2]]
    check("A11 shifted forms 7.965-2953.6t+2286t^2 and 3.375-4058t+3410t^2",
          Dsh == [Fr("7.965"), Fr("-2953.6"), 2286] and Nsh == [Fr("3.375"), -4058, 3410])
    # polynomial identities (CW l.484-489)
    lhs1 = ptrim(padd(pmul(NB[0], K0_R), pscale(pmul(D_R, NB[1]), 3)))
    check("A12 n_b(r,k(r)) = 0 identically, k(r)=K0/(3D)", lhs1 == [0])
    lhs2 = ptrim(padd(pmul(NUM_R, padd(pmul(NE[0], K0_R), pscale(pmul(D_R, NE[1]), 3))),
                      pmul(D_R, padd(pmul(N0[0], K0_R), pscale(pmul(D_R, N0[1]), 3)))))
    rhs2 = ptrim(pscale(pmul(pmul([0, 1], pmul([-5, 6], [-5, 6])), P_CUBIC), -1944))
    check("A13 3D^2{l(r)n_e+n_0} = -1944 r (6r-5)^2 p(r) identically", lhs2 == rhs2)
    lhs3 = ptrim(padd(pscale(K0_R, 2), padd(pscale(D_R, -5), pscale(NUM_R, 3))))
    check("A14 k(r) = 5/6 - l(r)/2 identically (so k(r_*) = k_*)", lhs3 == [0])
    AQ = (-6858 * ks + 600) * bs + (-702 * ks - 408) * ls + 4806 * ks - 696
    check("A15 A_Q > 2200", AQ > 2200, mpmath.nstr(AQ, 10))
    return dict(r=rs, l=ls, B=Bs, k=ks, b=bs, AQ=AQ)


# ------------------------------------------------------------ scalar system
def unit(d):
    if isinstance(d, (Fr, int)):
        return Fr(1)
    if hasattr(d, "delta"):  # mpmath interval
        return iv.mpf(1)
    return mpf(1)


def scal(d, x, b, l, k, one=None):
    """Return dict of the CW Section 6 scalar functions (works for Fr, mpf, iv)."""
    if one is None:
        one = unit(d)
    c = one / (3 * k)
    A = one * 5 / 6 - d
    D = 3 - (1 + 2 * c) * x
    P = 2 - 2 * (1 + c) * x + 2 * c * x * x
    J = A * D + d * P
    R = 1 - d + A * d * P / (2 * J)
    h = (1 + b + 3 * l) / 2
    ly = (1 + b - l) / 2
    C0 = -one / 4 + b / 6 + 5 * l / 4
    L = (23 * l - 1 - b) / 18
    th0 = 2 * d * x / 9
    E = C0 + d / 2 + d * x * L - (3 * h / 4) * (1 - R)
    S = (3 * R + 1) / 4 + th0 / 2 - one * 17 / 50
    aB = d * (189 * x - 25) / 139
    Q = E - (7 * aB / 108) * (2 * ly - h)
    return dict(c=c, A=A, D=D, P=P, J=J, R=R, h=h, ly=ly, C0=C0, L=L, th0=th0, E=E, S=S, aB=aB, Q=Q)


def rnd(lo, hi):
    return Fr(random.randint(0, 10 ** 6), 10 ** 6) * (Fr(hi) - Fr(lo)) + Fr(lo)


def section_B_identities():
    print("\n== B. Exact identities of the scalar comparison (CW l.410-495)")
    ok1 = ok2 = ok3 = ok4 = ok5 = True
    for _ in range(40):
        d, b, l, k = rnd("0.1", "0.8"), rnd("0.05", "0.12"), rnd("0.16", "0.18"), rnd("0.67", "0.75")
        s = scal(d, Fr(1, 2), b, l, k)
        nb, ne, n0 = nfun(NB, d, k), nfun(NE, d, k), nfun(N0, d, k)
        den = 432 * (54 * k * d - 75 * k - 6 * d + 10)
        ok1 &= (s["Q"] == (b * nb + l * ne + n0) / den)
        x = rnd("0", "0.5")
        s = scal(d, x, b, l, k)
        c, D, P = s["c"], s["D"], s["P"]
        Dp, Pp = -(1 + 2 * c), -2 * (1 + c) + 4 * c * x
        T = P * Dp - Pp * D
        ok2 &= (T == 4 + 2 * c - 12 * c * x + 2 * c * (1 + 2 * c) * x * x)
        # (T/D^2)' = (T' D - 2 T D')/D^3
        Tp = -12 * c + 4 * c * (1 + 2 * c) * x
        ok3 &= ((Tp * D - 2 * T * Dp) == 8 * (1 - c) ** 2)
        # R_x = -A^2 d T/(2 J^2)  ->  Q_x/d formula
        A, J = s["A"], s["J"]
        Jp = A * Dp + d * Pp
        Rx = A * d * (Pp * J - P * Jp) / (2 * J * J)
        ok4 &= (Rx == -A * A * d * T / (2 * J * J))
        Qx = d * s["L"] + (3 * s["h"] / 4) * Rx - (7 * d * 189 / (108 * 139)) * (2 * s["ly"] - s["h"])
        ok4 &= (Qx == d * (s["L"] - Fr(49, 556) * (2 * s["ly"] - s["h"])
                           - (3 * s["h"] / 8) * T / (D + (d / A) * P) ** 2))
        # dR/dk = d x A^2 (1-x)^2 / (3 k^2 J^2): check by exact derivative via c
        # R depends on k only through c = 1/(3k); dc/dk = -1/(3k^2)
        Dc, Pc = -2 * x, -2 * x + 2 * x * x
        Jc = A * Dc + d * Pc
        Rc = A * d * (Pc * J - P * Jc) / (2 * J * J)
        ok5 &= (Rc * (-Fr(1) / (3 * k * k)) == d * x * A * A * (1 - x) ** 2 / (3 * k * k * J * J))
    check("B1 Q(d,1/2) = (b n_b + l n_e + n_0)/(432(54kd-75k-6d+10)) (40 random rational points)", ok1)
    check("B2 T = P D' - P' D = 4+2c-12cx+2c(1+2c)x^2", ok2)
    check("B3 (T/D^2)' = 8(1-c)^2/D^3", ok3)
    check("B4 Q_x/d = L - (49/556)(2l_y-h) - (3h/8) T/[D+(d/A)P]^2", ok4)
    check("B5 d_k R = d x A^2 (1-x)^2/(3 k^2 J^2)", ok5)
    # E_b, E_l formulas and Q_l formula (affine in b, l)
    ok6 = True
    for _ in range(20):
        d, x, b, l, k = rnd("0.1", "0.8"), rnd("0", "0.5"), rnd("0.05", "0.12"), rnd("0.16", "0.18"), rnd("0.67", "0.75")
        s0, sb, sl = scal(d, x, b, l, k), scal(d, x, b + 1, l, k), scal(d, x, b, l + 1, k)
        R = s0["R"]
        ok6 &= (sb["E"] - s0["E"] == Fr(1, 6) - d * x / 18 - 3 * (1 - R) / 8)
        ok6 &= (sl["E"] - s0["E"] == Fr(5, 4) + Fr(23, 18) * d * x - 9 * (1 - R) / 8)
        ok6 &= (sl["Q"] - s0["Q"] == Fr(5, 4) + Fr(23, 18) * d * x - 9 * (1 - R) / 8 + 35 * s0["aB"] / 216)
    check("B6 E_b, E_l, Q_l formulas (E, Q affine in b, l)", ok6)
    # [H] sharp count identity behind R_{3/4} <= 1 - 2d/3 - 22dx/81
    ok7 = True
    for _ in range(20):
        d, x = rnd("0", "5/6"), rnd("0", "0.5")
        s = scal(d, x, Fr(0), Fr(0), Fr(3, 4))
        A, J, P = s["A"], s["J"], s["P"]
        y = Fr(1, 2) - x
        ok7 &= (2 * J - 3 * A * P - Fr(44, 27) * x * J
                == Fr(250, 729) * x * x + d / 243 * (352 * y ** 3 + 948 * y ** 2 + 984 * y + 199))
        ok7 &= ((1 - s["R"]) - 2 * d / 3 == d / (6 * J) * (2 * J - 3 * A * P))
    check("B7 [H] Lemma jsb:count identity (gives R_3/4 <= 1-2d/3-22dx/81)", ok7)
    # E(d,1/2), E(d,2/5) closed forms at the reference (b,l,k)=(.0918,.1688,3/4) (CW l.497-509)
    ok8 = True
    for _ in range(10):
        d = rnd("1/8", "5/6")
        e1 = scal(d, Fr(1, 2), Fr("0.0918"), Fr("0.1688"), Fr(3, 4))["E"]
        e2 = scal(d, Fr(2, 5), Fr("0.0918"), Fr("0.1688"), Fr(3, 4))["E"]
        ok8 &= (e1 == (-2659452 * d * d + 2425492 * d - 526140) / (120000 * (185 - 138 * d)))
        ok8 &= (e2 == (-201477144 * d * d + 167193610 * d - 35905500) / (600000 * (2525 - 1698 * d)))
    check("B8 closed forms of E(d,1/2), E(d,2/5) at b=.0918, l=.1688, k=3/4", ok8)
    s = scal(Fr(1, 2), Fr(1, 2), Fr("0.0918"), Fr("0.1688"), Fr(3, 4))
    check("B9 h_0=.7991, C_00=-.0237, L_0'=4651/30000",
          s["h"] == Fr("0.7991") and s["C0"] == Fr("-0.0237") and s["L"] == Fr(4651, 30000))
    # numerator expansions in y (CW l.512-514, l.527-531)
    ok10 = True
    for _ in range(10):
        y = rnd("-1", "1")
        d = y + Fr(2, 5)
        lhs = (-201477144 * d * d + 167193610 * d - 35905500) + Fr(600000, 1000) * (2525 - 1698 * d)
        ok10 &= (lhs == -201477144 * y * y + Fr(24965474, 5) * y - Fr(3922976, 25))
        d = y + Fr(23, 50)
        lhs = (-2659452 * d * d + 2425492 * d - 526140) - Fr(19 * 120000, 10000) * (185 - 138 * d)
        ok10 &= (lhs == -2659452 * y * y + Fr(256504, 25) * y - Fr(537677, 625))
    check("B10 numerator expansions of E(.,2/5)+1/1000 and E(.,1/2)-19/10000", ok10)
    # the two quadratic majorants are negative for all y
    f1 = lambda y: -156000 + 5000000 * abs(y) - 200000000 * y * y
    f2 = lambda y: -860 + 11000 * abs(y) - 2600000 * y * y
    check("B11 -156000+5e6|y|-2e8y^2<0 and -860+11000|y|-2.6e6y^2<0 for all y (discriminants)",
          5000000 ** 2 - 4 * 200000000 * 156000 < 0 and 11000 ** 2 - 4 * 2600000 * 860 < 0)
    # the dominations used in B11
    # E_dd formula (CW l.494-495) by numerical differentiation at the reference parameters
    ok13 = True
    for dd, xx in ((mpf("0.3"), mpf("0.45")), (mpf("0.5"), mpf("0.2")), (mpf("0.7"), mpf("0.5"))):
        f = lambda t: scal(t, xx, mpf("0.0918"), mpf("0.1688"), mpf(3) / 4)["E"]
        s = scal(dd, xx, mpf("0.0918"), mpf("0.1688"), mpf(3) / 4)
        formula = -(3 * s["h"] / 4) * (mpf(5) / 6) ** 2 * s["D"] * s["P"] ** 2 / s["J"] ** 3
        ok13 &= abs(mpmath.diff(f, dd, 2) - formula) < mpf(10) ** -30
    check("B13 E_dd = -(3h/4)(5/6)^2 D P^2/J^3 (numerical second derivative)", ok13)
    # values and slopes used for 'concavity and the derivative signs' (CW l.505-509)
    e1 = lambda d: (-2659452 * d * d + 2425492 * d - 526140) / (120000 * (185 - 138 * d))
    de1 = lambda d: ((-2 * 2659452 * d + 2425492) * (185 - 138 * d) + 138 * (-2659452 * d * d + 2425492 * d - 526140)) / (
        120000 * (185 - 138 * d) ** 2)
    check("B14 E(7/20,1/2)<-9/50000 with E'>0 there; E(3/5,1/2)<-1/500 with E'<0 there",
          e1(Fr(7, 20)) < -Fr(9, 50000) and de1(Fr(7, 20)) > 0 and e1(Fr(3, 5)) < -Fr(1, 500) and de1(Fr(3, 5)) < 0,
          "E(7/20)=%.3e E(3/5)=%.3e" % (float(e1(Fr(7, 20))), float(e1(Fr(3, 5)))))
    check("B12 coefficient dominations 24965474/5<5e6, 3922976/25>156000, 256504/25<11000, 537677/625>860",
          Fr(24965474, 5) < 5000000 and Fr(3922976, 25) > 156000 and Fr(256504, 25) < 11000
          and Fr(537677, 625) > 860)


# ------------------------------------------------------- interval machinery
def ivbox(lo, hi):
    return iv.mpf([lo, hi])


def bnb_max_below(f, box, bound, max_boxes=200000, min_width=None):
    """Prove sup f < bound on box (list of (lo,hi)) by subdivision.
    f takes a list of iv intervals and returns an iv interval.  Returns (ok, boxes, worst)."""
    eps = mpf(10) ** -28
    stack = [[(mpf(lo) - eps, mpf(hi) + eps) for lo, hi in box]]
    n = 0
    worst = None
    while stack:
        bx = stack.pop()
        n += 1
        if n > max_boxes:
            return False, n, ("budget", bx)
        val = f([ivbox(lo, hi) for lo, hi in bx])
        if val.b < bound:
            if worst is None or val.b > worst:
                worst = val.b
            continue
        # point check at the midpoint: a true counterexample?
        mid = [(mpf(lo) + mpf(hi)) / 2 for lo, hi in bx]
        vm = f([ivbox(m, m) for m in mid])
        if vm.a >= bound:
            return False, n, ("counterexample", mid, vm)
        widths = [mpf(hi) - mpf(lo) for lo, hi in bx]
        i = max(range(len(bx)), key=lambda j: widths[j])
        lo, hi = bx[i]
        m = (mpf(lo) + mpf(hi)) / 2
        b1 = list(bx)
        b2 = list(bx)
        b1[i] = (lo, m)
        b2[i] = (m, hi)
        stack.append(b1)
        stack.append(b2)
    return True, n, worst


def bnb_min_above(f, box, bound, **kw):
    ok, n, w = bnb_max_below(lambda X: -f(X), box, -bound, **kw)
    return ok, n, (-w if (w is not None and not isinstance(w, tuple)) else w)


def section_B_ranges(P):
    print("\n== B'. Range claims of the scalar comparison, interval branch-and-bound")
    ks = P["k"]
    kI = iv.mpf([ks - mpf(10) ** -30, ks + mpf(10) ** -30])
    bs, ls = P["b"], P["l"]
    # (eq:derivatives) CW l.433-437, k in [2/3,3/4]
    def Eb(X):
        d, x, k = X
        s = scal(d, x, iv.mpf(0), iv.mpf(0), k, one=iv.mpf(1))
        return iv.mpf(1) / 6 - d * x / 18 - 3 * (1 - s["R"]) / 8
    def El(X):
        d, x, k = X
        s = scal(d, x, iv.mpf(0), iv.mpf(0), k, one=iv.mpf(1))
        return iv.mpf(5) / 4 + iv.mpf(23) / 18 * d * x - 9 * (1 - s["R"]) / 8
    box = [(mpf(1) / 8, mpf(5) / 6), (mpf(0), mpf(1) / 2), (mpf(2) / 3, mpf(3) / 4)]
    o1 = bnb_max_below(Eb, box, mpf(1) / 5)
    o2 = bnb_min_above(Eb, box, -mpf(1) / 5)
    o3 = bnb_max_below(El, box, mpf(5) / 4 + mpf(10) ** -25)
    o4 = bnb_min_above(El, box, mpf(5) / 16 - mpf(10) ** -25)
    check("B'1 |E_b|<1/5 and 5/16<=E_l<=5/4 on [1/8,5/6]x[0,1/2], k in [2/3,3/4]",
          o1[0] and o2[0] and o3[0] and o4[0],
          "maxEb=%s minEb=%s maxEl=%s minEl=%s" % tuple(mpmath.nstr(o[2], 6) for o in (o1, o2, o3, o4)))
    # Lemma contact pieces at k_*: P/D>=3/8, T/D^2<2/3 for x in [2/5,1/2]
    def PD(X):
        (x,) = X
        c = 1 / (3 * kI)
        D = 3 - (1 + 2 * c) * x
        Pp = 2 - 2 * (1 + c) * x + 2 * c * x * x
        return Pp / D
    def TD2(X):
        (x,) = X
        c = 1 / (3 * kI)
        D = 3 - (1 + 2 * c) * x
        T = 4 + 2 * c - 12 * c * x + 2 * c * (1 + 2 * c) * x * x
        return T / (D * D)
    o5 = bnb_min_above(PD, [(mpf(2) / 5, mpf(1) / 2)], mpf(3) / 8)
    o6 = bnb_max_below(TD2, [(mpf(2) / 5, mpf(1) / 2)], mpf(2) / 3)
    check("B'2 P_c/D_c>=3/8 and T/D_c^2<2/3 on x in [2/5,1/2] (k=k_*)", o5[0] and o6[0],
          "min P/D=%s max T/D^2=%s" % (mpmath.nstr(o5[2], 8), mpmath.nstr(o6[2], 8)))
    check("B'3 (2/3)/(1+(21/29)(3/8))^2 = 107648/261075 < 5/12",
          Fr(2, 3) / (1 + Fr(21, 29) * Fr(3, 8)) ** 2 == Fr(107648, 261075) and Fr(107648, 261075) < Fr(5, 12))
    s = scal(bs, mpf(1) / 2, bs, ls, ks)  # only geometry used below
    check("B'4 at (b_*,l_*): L>3/20, 0<2l_y-h<1/8, h<4/5",
          s["L"] > mpf(3) / 20 and 0 < 2 * s["ly"] - s["h"] < mpf(1) / 8 and s["h"] < mpf(4) / 5,
          "L=%s 2ly-h=%s h=%s" % (mpmath.nstr(s["L"], 8), mpmath.nstr(2 * s["ly"] - s["h"], 8), mpmath.nstr(s["h"], 8)))
    def Qx_over_d(X):
        d, x = X
        s = scal(d, x, iv.mpf(bs), iv.mpf(ls), kI, one=iv.mpf(1))
        T = 4 + 2 * s["c"] - 12 * s["c"] * x + 2 * s["c"] * (1 + 2 * s["c"]) * x * x
        return s["L"] - iv.mpf(49) / 556 * (2 * s["ly"] - s["h"]) - (3 * s["h"] / 8) * T / (s["D"] + (d / s["A"]) * s["P"]) ** 2
    Rbox = [(mpf(7) / 20, mpf(3) / 5), (mpf(2) / 5, mpf(1) / 2)]
    o7 = bnb_min_above(Qx_over_d, Rbox, mpf(0))
    check("B'5 Q_x>0 on R at (b_*,l_*,k_*)", o7[0], "min Q_x/d=%s" % mpmath.nstr(o7[2], 8))
    def den(X):
        (d,) = X
        return 54 * kI * d - 75 * kI - 6 * d + 10
    o8 = bnb_max_below(den, [(mpf(7) / 20, mpf(3) / 5)], mpf(0))
    check("B'6 denominator 54k d-75k-6d+10<0 on [7/20,3/5] => Q(d,1/2)=A_Q(d-r_*)^2/den<=0", o8[0],
          "max=%s" % mpmath.nstr(o8[2], 8))
    def Sf(X):
        d, x = X
        return scal(d, x, iv.mpf(bs), iv.mpf(ls), kI, one=iv.mpf(1))["S"]
    o9 = bnb_min_above(Sf, Rbox, mpf(71) / 300 - mpf(10) ** -25)
    check("B'7 S>=71/300 on R", o9[0], "min S=%s" % mpmath.nstr(o9[2], 8))
    # direct check of Q(l_*)<=0 away from the contact (independent of the double-root argument)
    def Qf(lv):
        def f(X):
            d, x = X
            return scal(d, x, iv.mpf(bs), iv.mpf(lv), kI, one=iv.mpf(1))["Q"]
        return f
    far = [[(mpf(7) / 20, mpf("0.44")), (mpf(2) / 5, mpf(1) / 2)],
           [(mpf("0.46"), mpf(3) / 5), (mpf(2) / 5, mpf(1) / 2)],
           [(mpf("0.44"), mpf("0.46")), (mpf(2) / 5, mpf("0.49"))]]
    oks = [bnb_max_below(Qf(ls), bx, mpf(0)) for bx in far]
    check("B'8 Q(l_*)<0 on R minus a 0.01-neighbourhood of (r_*,1/2) (direct)", all(o[0] for o in oks),
          "max=%s" % mpmath.nstr(max(o[2] for o in oks), 6))
    # complement / full box / short sums at the actual parameters (CW l.494-560)
    corners = [(mpf("0.09175"), mpf("0.16883")), (mpf("0.09185"), mpf("0.16883"))]
    def Ef(bv, lv, kk):
        def f(X):
            d, x = X
            return scal(d, x, iv.mpf(bv), iv.mpf(lv), kk, one=iv.mpf(1))["E"]
        return f
    C1 = [(mpf(1) / 8, mpf(5) / 6), (mpf(0), mpf(2) / 5)]
    C2 = [(mpf(1) / 8, mpf(7) / 20), (mpf(2) / 5, mpf(1) / 2)]
    C3 = [(mpf(3) / 5, mpf(5) / 6), (mpf(2) / 5, mpf(1) / 2)]
    G = [(mpf(1) / 8, mpf(5) / 6), (mpf(0), mpf(1) / 2)]
    # reference (b_0,l_c,k=3/4): E_ref < -9/50000 on the three rectangles and < 19/10000 on G
    k34 = iv.mpf(3) / 4
    ref = [bnb_max_below(Ef(mpf("0.0918"), mpf("0.1688"), k34), bx, -mpf(9) / 50000) for bx in (C1, C2, C3)]
    check("B'9 reference E<-9/50000 on C1,C2,C3 (b=.0918,l=.1688,k=3/4)", all(o[0] for o in ref),
          "max=%s" % mpmath.nstr(max(o[2] for o in ref), 8))
    o10 = bnb_max_below(Ef(mpf("0.0918"), mpf("0.1688"), k34), G, mpf(19) / 10000)
    check("B'10 reference E<19/10000 on the full box", o10[0], "max=%s" % mpmath.nstr(o10[2], 8))
    # actual parameters: E affine in (b,l) with E_l>0, so the corners (b,.16883) dominate l<=.16883
    act = [bnb_max_below(Ef(bv, lv, kI), bx, -mpf(53) / 400000) for bv, lv in corners for bx in (C1, C2, C3)]
    check("B'11 E<-53/400000 on C1,C2,C3 for b in [.09175,.09185], l in [l_*,.16883], k=k_*",
          all(o[0] for o in act), "max=%s" % mpmath.nstr(max(o[2] for o in act), 8))
    actG = [bnb_max_below(Ef(bv, lv, kI), G, mpf(779) / 400000) for bv, lv in corners]
    check("B'12 E<779/400000 on the full box for the same parameter range", all(o[0] for o in actG),
          "max=%s" % mpmath.nstr(max(o[2] for o in actG), 8))
    def short(bv, lv):
        def f(X):
            d, x = X
            s = scal(d, x, iv.mpf(bv), iv.mpf(lv), kI, one=iv.mpf(1))
            Ds = (s["R"] - 1) / 4 + iv.mpf(3) / 80 + 17 * s["th0"] / 40
            return s["E"] + s["h"] * Ds
        return f
    sh = [bnb_max_below(short(bv, lv), Rbox, -mpf(1) / 100) for bv, lv in corners]
    check("B'13 short sums: E+h*Delta_s<-1/100 on R (same parameter range)", all(o[0] for o in sh),
          "max=%s" % mpmath.nstr(max(o[2] for o in sh), 8))
    def Ds_bound(X):
        d, x = X
        return iv.mpf(3) / 80 - d / 6 + iv.mpf(43) / 1620 * d * x
    o14 = bnb_max_below(Ds_bound, Rbox, -mpf(1049) / 64800 + mpf(10) ** -25)
    check("B'14 3/80-d/6+43dx/1620<=-1049/64800 on R", o14[0], "max=%s" % mpmath.nstr(o14[2], 10))
    # low-delta estimate and the x-monotonicity bound (CW l.496-501)
    s = scal(Fr(23, 100), Fr(1, 2), Fr("0.0918"), Fr("0.1688"), Fr(3, 4))
    val = s["C0"] + Fr(23, 100) * (1 - s["h"]) / 2
    check("B'15 C_00+(23/100)(1-h_0)/2 = -1193/2000000 and L_0'-11h_0/54<0",
          val == Fr(-1193, 2000000) and s["L"] - 11 * s["h"] / 54 < 0)
    check("B'16 856*181^2/7663^2<1/2 and 6697=37*181, 69/181=(23/100)/(5/6-23/100)",
          Fr(856 * 181 ** 2, 7663 ** 2) < Fr(1, 2) and 6697 == 37 * 181
          and Fr(23, 100) / (Fr(5, 6) - Fr(23, 100)) == Fr(69, 181))
    def Tratio(X):
        d, x = X
        s = scal(d, x, iv.mpf(0), iv.mpf(0), k34, one=iv.mpf(1))
        T = 4 + 2 * s["c"] - 12 * s["c"] * x + 2 * s["c"] * (1 + 2 * s["c"]) * x * x
        return T / (s["D"] + (d / s["A"]) * s["P"]) ** 2
    o17 = bnb_max_below(Tratio, [(mpf(23) / 100, mpf(5) / 6 - mpf(10) ** -20), (mpf(0), mpf(1) / 2)], mpf(1) / 2)
    check("B'17 T/[D+(d/A)P]^2<1/2 for d>=23/100, k=3/4 (E increasing in x)", o17[0],
          "max=%s" % mpmath.nstr(o17[2], 8))
    # Q_l bound (CW l.581-585)
    def Ql(X):
        d, x = X
        s = scal(d, x, iv.mpf(bs), iv.mpf(ls), kI, one=iv.mpf(1))
        return iv.mpf(5) / 4 + iv.mpf(23) / 18 * d * x - 9 * (1 - s["R"]) / 8 + 35 * s["aB"] / 216
    o18 = bnb_max_below(Ql, Rbox, mpf(593) / 500)
    check("B'18 Q_l<593/500 on R", o18[0], "max=%s" % mpmath.nstr(o18[2], 10))
    check("B'19 5/4-79(7/20)/432 < 593/500", Fr(5, 4) - Fr(79 * 7, 20 * 432) < Fr(593, 500))
    return o18[2]


def section_C():
    print("\n== C. Exhaustion norm and first-mean gain (CW l.198-405)")
    K = Fr(2247, 500)
    kap = Fr(749, 1000)
    cK = (K - 1) / K
    check("C1 K=6*kappa, c_K=1747/2247, kappa/K=1/6", K == 6 * kap and cK == Fr(1747, 2247) and kap / K == Fr(1, 6))
    ok = True
    for _ in range(30):
        v = {n: rnd("0", "0.05") for n in "wf f0 f1 gg vv cc T1 T e1 e S1 S z zC muextra".split()}
        m = rnd("0.41", "0.43")
        u = v["gg"] + v["vv"] + v["cc"]
        w = v["wf"] + v["f0"] + v["f1"] + u
        t = v["f0"] + u
        L1 = v["f1"] + 4 * v["T1"] + v["e1"]
        L = v["vv"] + 4 * v["T"] + v["e"]
        At = 2 * m + v["z"] + w
        Ast = At - v["f1"] - t - L1 - L - 5 * v["cc"]
        Mst = 1 - 6 * t - v["f1"] - L1 - L - 5 * v["cc"] + v["T1"] + 5 * v["vv"] + v["T"]
        q = v["f1"] + v["vv"] + v["cc"]
        Rloc = -v["f0"] - L1 + v["T1"] - v["S1"] - u - L + v["vv"] - 4 * v["cc"] + v["T"] - v["S"]
        etaK = 1 - 2 * m - K * v["z"]
        Dfree = Ast + (K - 1) * (v["z"] - v["zC"]) - Mst
        ok &= Dfree == 6 * w - etaK - (K - 1) * v["zC"] - 5 * v["wf"] - 5 * v["f1"] - v["T1"] - 5 * v["vv"] - v["T"]
        ok &= Rloc == -v["f0"] - v["f1"] - u - 3 * (v["T1"] + v["T"]) - 4 * v["cc"] - v["S1"] - v["S"] - v["e1"] - v["e"]
        mu = v["wf"] + v["muextra"]
        D0 = Dfree - mu
        ok &= Rloc + D0 / 6 == (-5 * v["f1"] / 6 - 19 * (v["T1"] + v["T"]) / 6 - 5 * v["vv"] / 6 - 4 * v["cc"]
                                - v["S1"] - v["S"] - v["e1"] - v["e"] - etaK / 6 - (K - 1) * v["zC"] / 6
                                + (v["wf"] - mu) / 6)
        # [H] small-formal identity (jtf:terminal) and its <= 3A_t-7/2+18w
        mst = Mst - q
        lhs = Rloc + 3 * Ast - Fr(7, 2) * mst - q
        rhs = (3 * At - Fr(7, 2) + 17 * v["f0"] + Fr(5, 2) * v["f1"] + 17 * v["gg"] + Fr(5, 2) * v["vv"]
               + 18 * v["cc"] - Fr(9, 2) * (v["T1"] + v["T"]) - v["S1"] - v["S"] - (v["e1"] + v["e"]) / 2)
        ok &= lhs == rhs and rhs <= 3 * At - Fr(7, 2) + 18 * w
        # M* and m* lower bounds when L1+L+5c<=2m
        if L1 + L + 5 * v["cc"] <= 2 * m:
            ok &= Mst >= 1 - 2 * m - 6 * w and mst >= 1 - 2 * m - 7 * w
    check("C2 D_free identity (l.261-263), R_loc simplification, R_loc+D_0/6 identity (l.313-319), [H] small-formal identity, M*/m* bounds", ok)
    check("C3 small-formal: -1/2+3d0/K+7/200 <= -93/200+90/749 < -1/3 (d0<=9/50)",
          -Fr(1, 2) + 3 * Fr(9, 50) / K + Fr(7, 200) == -Fr(93, 200) + Fr(90, 749)
          and -Fr(93, 200) + Fr(90, 749) < -Fr(1, 3))
    check("C4 origin/exceptional: c_K>5/7 (K>7/2)", cK > Fr(5, 7) and K > Fr(7, 2))
    # margins (l.374-377)
    W = lambda m: cK * (1 - 2 * m) / 6 + Fr(1, 1000)
    zK = lambda m: (1 - 2 * m) / K
    f1 = lambda m: 1 - 2 * m - 7 * W(m)
    f2 = lambda m: 1 - (Fr(5, 2) - 4 * m + K * W(m))
    m1, m2 = Fr(41, 100), Fr(43, 100)
    check("C5 margins: min 1-2m-7W_new = 40523/6741000 > 3/500; min 1-(5/2-4m+K W_new) = 15343/500000 > 3/100",
          min(f1(m1), f1(m2)) == Fr(40523, 6741000) > Fr(3, 500)
          and min(f2(m1), f2(m2)) == Fr(15343, 500000) > Fr(3, 100))
    # covering with intervals of width <1/1000, W=W_new(m_-): 1-2r_o-7W_new(m_-) > 1/250
    cov = min(1 - 2 * (mm + Fr(1, 1000)) - 7 * W(mm) for mm in (m1, m2 - Fr(1, 1000)))
    check("C6 covering: 1-2r_o-7W_new(m_-) > 1/250 for r_o<=m_-+1/1000 (linear in m_-)",
          cov > Fr(1, 250), "min=%s" % float(cov))
    check("C7 W_new<=1/40, z_K<41/1000 on [.41,.43]", W(m1) < Fr(1, 40) and zK(m1) < Fr(41, 1000))
    L0 = Fr(21, 100)
    a = lambda d, x, z: d * (x * L0 - z / 2) / (L0 - z)
    # a increases in d and x, decreases in z: minimum at d=7/20, x=2/5, z=z_K(.41)
    amin = a(Fr(7, 20), Fr(2, 5), zK(m1))
    check("C8 a(z_K)>13/100 on R x [.41,.43] (monotone; min at d=7/20,x=2/5,m=.41)", amin > Fr(13, 100),
          "min a=%s" % float(amin))
    check("C9 21/10000-9z/200>51/200000 for z<=41/1000", Fr(21, 10000) - 9 * Fr(41, 1000) / 200 == Fr(51, 200000))
    check("C10 [H] a_box = a(z=1/18) with L0=21/100: d(189x-25)/139",
          all(a(d, x, Fr(1, 18)) == d * (189 * x - 25) / 139 for d, x in ((Fr(1, 2), Fr(9, 20)), (Fr(2, 5), Fr(1, 2)))))
    okF = True
    for _ in range(30):
        d, x, v, s = rnd("7/20", "3/5"), rnd("2/5", "1/2"), rnd("7/50", "9/50"), rnd("0.2", "0.23")
        g = d * x
        F = lambda s: (g * s * v + a(d, x, s * v) * (1 - s) * v / 6) / 2
        # exact derivative of F via the quotient rule
        Nn = (1 - s) * (x * L0 - s * v / 2)
        Dn = L0 - s * v
        dN = -(x * L0 - s * v / 2) - (1 - s) * v / 2
        dF = (g * v + d * v / 6 * (dN * Dn + v * Nn) / Dn ** 2) / 2
        okF &= dF == d * v / 2 * (x - Fr(1, 12) + (Fr(1, 2) - x) * L0 * (L0 - v) / (6 * (L0 - s * v) ** 2))
        okF &= F(1 / K) >= F(Fr(2, 9))
    check("C11 F'(s) formula (l.387-390) and G_K>=G_0 (s=1/K>2/9)", okF and 1 / K > Fr(2, 9))
    # G_new < 1/100: g,a <= 3/10; z_K, W_new maximal at m=.41
    check("C12 G_new<1/100 (bound (3/10)(z_K+W_new)/2 at m=.41)", Fr(3, 10) * (zK(m1) + W(m1)) / 2 < Fr(1, 100))
    check("C13 inventory 2xL0-z_K>127/1000>W", 2 * Fr(2, 5) * L0 - zK(m1) > Fr(127, 1000) > W(m1))
    nu = Fr(13, 200000) - Fr(1, 10 ** 6) - Fr(6, 10 ** 6) - Fr(2, 10 ** 6)
    check("C14 net gain 13/200000-1e-6-6e-6-2e-6 > 1/20000", nu > Fr(1, 20000), "=%s" % float(nu))
    okR = True
    for _ in range(30):
        m = rnd("0.41", "0.43")
        z = rnd("0", "1") * (1 - 2 * m) / K
        w = rnd("0", "1") * W(m)
        okR &= (Fr(5, 2) - 4 * m + K * w) - (Fr(3, 2) - 2 * m + K * (z + w)) == 1 - 2 * m - K * z >= 0
    check("C15 restoration: 3/2-2m+K(z+w) <= 5/2-4m+Kw when z<=(1-2m)/K (l.342-343)", okR)
    # the [H] domain for comparison: W_K vs [H]'s 7(1-2m)/54
    check("C16 W_K=c_K v/6 < 7v/54 (CW's base width is below [H]'s; new gain is the +1/1000)",
          cK / 6 < Fr(7, 54))
    # cap kappa=749/1000 lies below the 7497/10000 of NF3 Remark 14.1(1); its two uses still hold
    check("C17 NF3 Rem.14.1(1) uses at kappa=749/1000: (21/100)(6k-1/8)-(1-1/12)>0 and 6k*l>3/4 (l>=.1687)",
          Fr(21, 100) * (6 * kap - Fr(1, 8)) - (1 - Fr(1, 12)) > 0 and 6 * kap * Fr(1687, 10000) > Fr(3, 4),
          "supply margin=%.2e (2263/3024=%.5f is Liu's threshold)" % (
              float(Fr(21, 100) * (6 * kap - Fr(1, 8)) - (1 - Fr(1, 12))), 2263 / 3024))


def section_D(P, Qlmax):
    print("\n== D. Global comparison (CW l.565-630)")
    bs, ls = P["b"], P["l"]
    l = ls + mpf(1) / 40000
    h = (1 + bs + 3 * l) / 2
    ly = (1 + bs - l) / 2
    lam = ly / h
    m0 = 1 - lam
    check("D1 lambda in (.577,.578), m_0 in (.422,.423), h>.79 at l=l_*+1/40000",
          mpf("0.577") < lam < mpf("0.578") and mpf("0.422") < m0 < mpf("0.423") and h > mpf("0.79"),
          "lambda=%s h=%s" % (mpmath.nstr(lam, 10), mpmath.nstr(h, 10)))
    check("D2 Q(l) <= Q(l_*) + Q_l/40000 < 593/20000000 (Q affine in l, Q(l_*)<=0, Q_l<593/500)",
          Qlmax / 40000 < mpf(593) / 20000000, "max increase=%s" % mpmath.nstr(Qlmax / 40000, 8))
    check("D3 593/20000000-.79/20000 = -9.85e-6", Fr(593, 20000000) - Fr(79, 100) / 20000 == Fr(-985, 10 ** 8))
    lamI = iv.mpf([lam - mpf(10) ** -30, lam + mpf(10) ** -30])
    m0I = 1 - lamI
    def envelope_excess(X):
        d, x = X
        th0 = 2 * d * x / 9
        a = d * (189 * x - 25) / 139
        GN = 7 * a * (2 * lamI - 1) / 108
        GL = th0 * (iv.mpf(5) / 8 - lamI)          # m=3/8 (and the central range m<3/8)
        e41 = (m0I - iv.mpf("0.41")) * (7 * a / 54 + th0)
        e43 = (iv.mpf("0.43") - m0I) * (d / 2 - th0 - 7 * a / 54)
        return iv.mpf([min(v.a for v in (GL - GN, e41, e43)), min(v.b for v in (GL - GN, e41, e43))])
    Rbox = [(mpf(7) / 20, mpf(3) / 5), (mpf(2) / 5, mpf(1) / 2)]
    o = bnb_min_above(envelope_excess, Rbox, mpf(7) / 150000)
    check("D4 [H] envelope F(m): excess over Gamma_N outside [.41,.43] >= 7/150000 on R", o[0],
          "min=%s" % mpmath.nstr(o[2], 8))
    def slopes(X):
        d, x = X
        th0 = 2 * d * x / 9
        a = d * (189 * x - 25) / 139
        return iv.mpf([min((th0 - iv.mpf(7) / 225).a, (d / 2 - th0 - 7 * a / 54 - iv.mpf(49) / 432).a,
                           (4 * a / 9 - th0).a),
                       max((th0 - iv.mpf(7) / 225).b, (d / 2 - th0 - 7 * a / 54 - iv.mpf(49) / 432).b,
                           (4 * a / 9 - th0).b)])
    o2 = bnb_min_above(slopes, Rbox, -mpf(10) ** -25)
    check("D5 slopes: theta_0>=7/225, d/2-theta_0-7a/54>=49/432, 4a/9-theta_0>0 on R", o2[0])
    check("D6 lambda<11/19 (so Gamma=Gamma_N, [H] Lemma jca:taper)", lam < mpf(11) / 19)
    check("D7 .79*7/150000-593/20000000-1e-7 = 427/60000000 > 7e-6",
          Fr(79, 100) * Fr(7, 150000) - Fr(593, 20000000) - Fr(1, 10 ** 7) == Fr(427, 60000000) > Fr(7, 10 ** 6))
    check("D8 779/400000-71/24000 = -1213/1200000", Fr(779, 400000) - Fr(71, 24000) == -Fr(1213, 1200000))
    # coarse-box conditions (l.589-594), b in (.0917,.0919), l in (.1687,.1689)
    bI = iv.mpf(["0.0917", "0.0919"])
    lI = iv.mpf(["0.1687", "0.1689"])
    hI = (1 + bI + 3 * lI) / 2
    lyI = (1 + bI - lI) / 2
    lxI = (1 - bI - lI) / 2
    C0I = -iv.mpf(1) / 4 + bI / 6 + 5 * lI / 4
    conds = [(lxI - lI).a > 0, (lyI - lI - 11 * bI / 6).a > 0, (2 * lyI - hI).a > 0, lyI.b < mpf(1) / 2,
             (lI * (iv.mpf(1) / 2 + 2 * lI / hI) - 11 * hI / 54).b < -mpf(1) / 200,
             (C0I + (1 - hI) / 16).b < -mpf(1) / 250,
             (C0I + (iv.mpf(1) / 2 + 3 * lI / 2) / 50).b < -mpf(1) / 125,
             (13 * hI / 75 - lyI / 2 + iv.mpf(1) / 50).b < -mpf(1) / 20,
             (lI / (hI + iv.mpf(1) / 2000)).a > mpf(21) / 100, (lI / hI).b < mpf(2) / 9,
             (11 * hI / 54 > 0)]
    check("D9 coarse-box conditions of l.589-594 plus l/(h+1/2000)>21/100, l/h<2/9", all(conds),
          "%s" % [int(c) for c in conds])
    check("D10 B>437/500 (H_eta nonzero there)", P["B"] - mpf(1) / 160000 > mpf(437) / 500)
    # [H] eq. (nf:outer-reserve): outer range m>=9/20 has excess >= 7 delta/2160 (lambda>=57/100)
    okO = True
    for _ in range(20):
        x, lm = rnd("2/5", "1/2"), rnd("57/100", "59/100")
        y, z = Fr(1, 2) - x, lm - Fr(57, 100)
        lhs = (Fr(1, 2) - 2 * x / 9) * (lm - Fr(11, 20)) - 7 * (189 * x - 25) * (2 * lm - 1) / (108 * 139) - Fr(7, 2160)
        okO &= lhs == Fr(997, 2502) * y * z + Fr(4199, 250200) * y + Fr(35, 108) * z
    check("D11 [H] outer-reserve identity (excess 7d/2160 for m>9/20) and lambda>=57/100 here", okO and lam > mpf(57) / 100)
    # [H] Lemma jca:taper identity for Gamma_L-Gamma_N (lambda<=11/19)
    okT = True
    for _ in range(20):
        x, lm = rnd("2/5", "1/2"), rnd("57/100", "11/19")
        y, z = Fr(1, 2) - x, Fr(11, 19) - lm
        GL = 2 * x / 9 * (Fr(5, 8) - lm)
        GN = 7 * (189 * x - 25) * (2 * lm - 1) / (139 * 108)
        okT &= GL - GN == Fr(175, 47538) * y + z * (Fr(19, 108) - Fr(997, 2502) * y)
    check("D12 [H] identity (Gamma_L-Gamma_N)/d = 175y/47538 + z(19/108-997y/2502)", okT)
    # CW's main new step at the contact class, explicit: Q_l at (r_*,1/2)
    s = scal(P["r"], mpf(1) / 2, bs, ls, P["k"])
    Ql_c = mpf(5) / 4 + mpf(23) / 18 * P["r"] / 2 - 9 * (1 - s["R"]) / 8 + 35 * s["aB"] / 216
    print("     info: Q_l at the contact class = %s; exponent at contact after the gain = %s"
          % (mpmath.nstr(Ql_c, 10), mpmath.nstr(Ql_c / 40000 - h / 20000, 6)))
    # informational (non-rigorous) true extrema on grids
    mp.dps = 20
    best = None
    for i in range(101):
        d = mpf(7) / 20 + (mpf(3) / 5 - mpf(7) / 20) * i / 100
        for j in range(21):
            x = mpf(2) / 5 + (mpf(1) / 10) * j / 20
            q = scal(d, x, bs, l, P["k"])["Q"]
            if best is None or q > best[0]:
                best = (q, d, x)
    qc = scal(P["r"], mpf(1) / 2, bs, l, P["k"])["Q"]
    print("     info: max Q(l_*+1/40000) on R ~ %s (grid; at contact %s) vs h*nu = %s"
          % (mpmath.nstr(best[0], 6), mpmath.nstr(qc, 6), mpmath.nstr(h / 20000, 6)))
    ex = []
    for d in (mpf(7) / 20, mpf(3) / 5):
        for x in (mpf(2) / 5, mpf(1) / 2):
            th0 = 2 * d * x / 9
            a = d * (189 * x - 25) / 139
            GN = 7 * a * (2 * lam - 1) / 108
            ex.append(min(th0 * (mpf(5) / 8 - lam) - GN, (m0 - mpf("0.41")) * (7 * a / 54 + th0),
                          (mpf("0.43") - m0) * (d / 2 - th0 - 7 * a / 54)))
    print("     info: true minimum envelope excess outside the band ~ %s (claimed >= %s)"
          % (mpmath.nstr(min(ex), 6), mpmath.nstr(mpf(7) / 150000, 6)))
    mp.dps = 50


def main():
    P = section_A()
    section_B_identities()
    Qlmax = section_B_ranges(P)
    section_C()
    section_D(P, Qlmax)
    nf = sum(1 for _, ok in RESULTS if not ok)
    print("\nSUMMARY: %d checks, %d failed" % (len(RESULTS), nf))
    return 1 if nf else 0


if __name__ == "__main__":
    sys.exit(main())

