#!/usr/bin/env python3
"""Re-check the numerical and algebraic claims of

  B. Liu, "An improved zero-free half-plane for Hecke L-functions over number
  fields and Quasi-Riemann Hypothesis", Hexagon 2610.00655v3 (hexagonmath.org,
  3 Nov 2026; source number_field_new_v1-4.tex), boundary B_tw = 0.874058503...

Stdlib + mpmath only.  "l.N" is line N of number_field_new_v1-4.tex (v3).

  * Exact identities use fractions.Fraction: rational identities in several
    variables are checked at random rational points, one-variable polynomial
    identities and root counts (Sturm) exactly.
  * Range claims use mpmath.iv (outward rounding) with forward-mode interval
    derivatives and branch-and-bound subdivision, so a PASS on a range claim is
    a proof up to the correctness of mpmath.iv.  Endpoint parameters that are
    algebraic numbers enter as intervals of radius 1e-50 around 90-digit values.
  * Contact (double-root) identities at the algebraic endpoints are checked
    numerically at 90 digits (the paper reduces them modulo the minimal
    polynomial; we do not redo that exact reduction).

Run:  python3 verify.py      (single-threaded, about 6 min, about 30 MB; 106 checks)
"""
import random
import sys
from fractions import Fraction as Fr

import mpmath
from mpmath import mp, iv, mpf

mp.dps = 90
iv.dps = 40
random.seed(20261105)
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("  [" + detail + "]") if detail else ""))
    sys.stdout.flush()


def info(name, detail):
    print("info " + name + "  [" + detail + "]")
    sys.stdout.flush()


def ivfr(q):
    q = Fr(q)
    return iv.mpf(q.numerator) / q.denominator


def ivnear(x, rad=mpf(10) ** -50):
    return iv.mpf([x - rad, x + rad])


def mpf2fr(x):
    s, man, ex, _ = x._mpf_
    v = Fr(int(man)) * (Fr(2) ** ex if ex >= 0 else Fr(1, 2 ** (-ex)))
    return -v if s else v


def rnd(lo, hi, den=10 ** 6):
    lo, hi = Fr(lo), Fr(hi)
    return lo + (hi - lo) * Fr(random.randint(1, den - 1), den)


# ------------------------------------------------------------ polynomials (low -> high)
def ptrim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def padd(p, q):
    n = max(len(p), len(q))
    return ptrim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def pscale(p, c):
    return ptrim([c * a for a in p])


def psub(p, q):
    return padd(p, pscale(q, -1))


def pmul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return ptrim(out)


def pderiv(p):
    return ptrim([i * p[i] for i in range(1, len(p))]) if len(p) > 1 else [0]


def peval(p, x):
    acc = 0
    for a in reversed(p):
        acc = acc * x + a
    return acc


def pdivmod(a, b):
    a = [Fr(c) for c in a]
    b = [Fr(c) for c in b]
    q = [Fr(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        c = a[-1] / b[-1]
        k = len(a) - len(b)
        q[k] = c
        a = ptrim(psub(a, [0] * k + pscale(b, c)))
        if len(a) == 1 and a[0] == 0:
            break
    return ptrim(q), a


def sturm_count(p, lo, hi):
    """Number of distinct real roots of p in (lo, hi]; p(lo), p(hi) != 0."""
    seq = [ptrim([Fr(c) for c in p]), pderiv([Fr(c) for c in p])]
    while len(seq[-1]) > 1 or seq[-1][0] != 0:
        _, r = pdivmod(seq[-2], seq[-1])
        if len(r) == 1 and r[0] == 0:
            break
        seq.append(pscale(r, -1))

    def var(x):
        vals = [peval(s, x) for s in seq]
        vals = [v for v in vals if v != 0]
        return sum(1 for i in range(len(vals) - 1) if (vals[i] > 0) != (vals[i + 1] > 0))

    return var(Fr(lo)) - var(Fr(hi))


def poly_of(f, deg, pts=None):
    """Interpolate the polynomial of degree <= deg through exact values f(x) (low->high)."""
    xs = pts or [Fr(1, 10) + Fr(i, 13) for i in range(deg + 1)]
    ys = [f(x) for x in xs]
    # Newton divided differences
    n = len(xs)
    coef = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (xs[i] - xs[i - j])
    p = [coef[-1]]
    for k in range(n - 2, -1, -1):
        p = padd(pmul(p, [-xs[k], 1]), [coef[k]])
    return p


def prop(p, q):
    """True if polynomials p, q are proportional (nonzero)."""
    p, q = ptrim(p), ptrim(q)
    if len(p) != len(q):
        return False
    i = max(range(len(q)), key=lambda k: abs(q[k]))
    c = Fr(p[i]) / Fr(q[i])
    return all(Fr(p[k]) == c * q[k] for k in range(len(p)))


# =============================================================== A. Part III algebra
def amax(p):
    return (4 * p - 6) / (p * (p + 24))


def rmax(p):
    return 6 * (p - 6) * (2 * p - 3) / (p * p - 193 * p + 198)


def Eal(p, al, rho):
    return (7 - 5 * p * al) * rho / 18 - p * al * (6 - p - p * rho) / (6 * (p - 1))


PLO, PHI = Fr(112, 25), Fr(9, 2)


def sectionA():
    print("== A. Part III: return estimate and returned-sum bound (l.8907-9459, 11400-11610)")
    ps = [rnd(PLO, PHI) for _ in range(6)]
    # A1 lambda_max (l.8925-8928)
    check("A1 lambda_max = 1 + p alpha_max = (5p+18)/(p+24) (l.8927)",
          all(1 + p * amax(p) == (5 * p + 18) / (p + 24) for p in ps))
    # A2 d rho_max / dp (l.8931)
    N, D = [108, -90, 12], [198, -193, 1]
    num = psub(pmul(pderiv(N), D), pmul(N, pderiv(D)))
    check("A2 d rho_max/dp = -42(53p^2-108p-72)/(p^2-193p+198)^2 (l.8931)",
          num == pscale([-72, -108, 53], -42))
    # A3 range of rho_max (l.8933-8936)
    ok = (rmax(PHI) == Fr(24, 289) and rmax(PLO) == Fr(5662, 67351) and Fr(5662, 67351) < Fr(1, 10)
          and sturm_count([-72, -108, 53], PLO, PHI) == 0 and peval([-72, -108, 53], PLO) > 0
          and sturm_count(D, PLO, PHI) == 0 and peval(D, PLO) != 0)
    check("A3 rho_max(9/2)=24/289, rho_max(112/25)=5662/67351<1/10, rho_max decreasing on [112/25,9/2] (l.8933)", ok)
    # A4 E_{alpha_max}(p, rho_max(p)) = 0 (l.10777)
    check("A4 E_{alpha_max}(p,rho_max(p)) = 0 (l.639, 10777)", all(Eal(p, amax(p), rmax(p)) == 0 for p in ps))
    # A5 derivatives of E (Lemma tw:lem:supremum, l.9283-9288); E is affine in rho and in alpha
    ok = True
    for p in ps:
        al, rho = rnd(0, amax(p)), rnd(0, 1)
        dr = Eal(p, al, rho + 1) - Eal(p, al, rho)
        da = Eal(p, al + 1, rho) - Eal(p, al, rho)
        ok &= dr == (7 - 5 * p * al) / 18 + p * p * al / (6 * (p - 1))
        ok &= da == -p * (18 - 3 * p + (2 * p - 5) * rho) / (18 * (p - 1))
        ok &= dr > 0 and da < 0
    check("A5 dE/drho=(7-5p a)/18+p^2 a/(6(p-1))>0, dE/dalpha=-p{18-3p+(2p-5)rho}/(18(p-1))<0 (l.9283-9288)", ok)
    # A6 dE_{amax}/drho (l.10780)
    check("A6 dE_{alpha_max}/drho = -(p^2-193p+198)/(18(p-1)(p+24)) > 0 (l.10780)",
          all(Eal(p, amax(p), 1) - Eal(p, amax(p), 0) == -(p * p - 193 * p + 198) / (18 * (p - 1) * (p + 24))
              and -(p * p - 193 * p + 198) > 0 for p in ps))
    # A7 alpha-choice margin (tw:eq:alpha-choice, l.9178; identity l.10806-10809)
    ok = True
    for p in ps:
        vr = 1 - 1 / p
        lhs = p * amax(p) / (6 * vr) - (7 - 5 * p * amax(p)) / 6000 - Fr(1, 12)
        ok &= lhs == (3513 * p * p - 17711 * p + 12198) / (6000 * (p - 1) * (p + 24))
    q = [12198, -17711, 3513]
    ok &= sturm_count(q, Fr(422, 100), PHI) == 0 and peval(q, PHI) > 0 and sturm_count(q, Fr(422, 100), 100) == 0
    check("A7 p a_max/(6 vr)-(7-5p a_max)/6000-1/12 = (3513p^2-17711p+12198)/(6000(p-1)(p+24)) > 0, roots < 4.22 (l.10806)", ok)
    # A8 alpha_max < vr/5, p > lambda_max, lambda_max < 6 (l.10783-10785) reduce to p^2+3p+6>0, p^2+19p-18>0, 5p+18<6p+144
    ok = all(amax(p) < (1 - 1 / p) / 5 and p > (5 * p + 18) / (p + 24) and (5 * p + 18) / (p + 24) < 6 for p in ps + [PLO, PHI])
    ok &= sturm_count([6, 3, 1], 0, 100) == 0 and sturm_count([-18, 19, 1], 1, 100) == 0
    check("A8 alpha_max<vr/5, p>lambda_max, lambda_max/6<1 on [112/25,9/2] (l.10783)", ok)
    # A9 numerical consequences of rho <= 5662/67351 (l.9077-9079, 9136-9141, 9349-9352, 10787-10790)
    rM = Fr(5662, 67351)
    ok = (Fr(503, 1500) * rM < Fr(3, 100) and 1 - rM / 4 > Fr(97, 100) and -1 + Fr(1381, 1500) * rM < Fr(-9, 10)
          and Fr(-97, 100) + (12 + 3 * rM) / 25 < Fr(-47, 100) and Fr(-93, 200) + (6 + 3 * rM) / 25 < Fr(-21, 100)
          and Fr(-47, 100) < Fr(-1, 5) and Fr(-21, 100) < Fr(-1, 5)
          and Fr(1, 3) + Fr(2, 1000) == Fr(503, 1500) and Fr(1, 4) + Fr(1006, 1500) == Fr(1381, 1500)
          and 2 - Fr(25, 6) - Fr(1, 6) == Fr(-7, 3) and 2 - 5 == -3 and 2 - Fr(7, 2) * 4 - 1 == -13)
    check("A9 ratios 503/1500, 97/100, 1381/1500; raw terminal ceilings < -47/100, -21/100; shifts -7/3, -3, -13 (l.9077-9142, 9349)", ok)
    # A10 strict margin identities (l.8992-9006)
    ok = True
    for _ in range(6):
        p = rnd(PLO, PHI); vr = 1 - 1 / p; rho = rnd(0, rM)
        dm, w = rnd(0, 1), rnd(0, Fr(1, 25))
        z = dm / p
        At = (1 - dm) + z + w
        G0 = dm - 6 * w
        kp = (6 - p - p * rho) / (p - 1)
        ok &= At - 1 + 4 * w == -vr * dm + 5 * w
        ok &= G0 == (vr * dm - (5 - rho) * w) / vr + kp * w
        ok &= (5 - rho) / vr - 6 == kp
    check("A10 A_t-1+4w=-vr d_m+5w at z=d_m/p; G_0=(vr d_m-(5-rho)w)/vr+k_p w; (5-rho)/vr-6=k_p (l.9002-9006, 10792)", ok)

    # A11 active estimate (tw:eq:active), exact identity incl. discarded terms (l.9236-9262)
    ok = True
    for _ in range(6):
        p = rnd(PLO, PHI); al = rnd(0, amax(p)); lam = 1 + p * al
        wf, f1, v, T1, T, cc, S1, S, e1, e = [rnd(0, Fr(1, 10)) for _ in range(10)]
        mu = wf + rnd(0, Fr(1, 10)); zC, z2, a, G0, w, ez, s = [rnd(0, Fr(1, 5)) for _ in range(7)]
        Dp = rnd(0, 1)
        Chom = -w + wf - 3 * (T1 + T) - 4 * cc - S1 - S - e1 - e
        pz = G0 + 6 * w - ez
        numer = pz - lam * (G0 + 5 * wf + 5 * f1 + 5 * v + T1 + T + mu + 5 * s) - (p - lam) * (zC + z2) + lam * a
        lhs = Chom + 2 * s + numer / 6 - Dp
        disc = (-p * al * wf - lam / 6 * (mu - wf) - 5 * lam / 6 * (f1 + v) - (3 + lam / 6) * (T1 + T)
                - 4 * cc - S1 - S - e1 - e - (p - lam) / 6 * (zC + z2))
        rhs = -p * al / 6 * G0 + (7 - 5 * p * al) / 6 * s - ez / 6 + (lam * a / 6 - Dp) + disc
        ok &= lhs == rhs
        ok &= p - lam > 0 and lam < 6
        # coefficient of w with G_0 = k_p w + nu/vr and s = rho w/3 + theta (tw:eq:Eendpoint)
        rho = rnd(0, rmax(p)); vr = 1 - 1 / p; kp = (6 - p - p * rho) / (p - 1); nu, th = rnd(0, 1), rnd(0, 1)
        G0b, sb = kp * w + nu / vr, rho * w / 3 + th
        ok &= (-p * al / 6 * G0b + (7 - 5 * p * al) / 6 * sb
               == Eal(p, al, rho) * w - p * al * nu / (6 * vr) + (7 - 5 * p * al) / 6 * th)
    check("A11 (tw:eq:active): exact identity with the displayed discarded terms; coefficient of w = E_alpha (l.9240-9262)", ok)
    # A12 untwisted homogeneous identity (l.9365-9373, 7390-7397) and R_loc = C_hom
    ok = True
    for _ in range(6):
        p = rnd(PLO, PHI)
        wf, f0, f1, g, v, cc, T1, T, S1, S, e1, e, ez, zC = [rnd(0, Fr(1, 10)) for _ in range(14)]
        mu = wf + rnd(0, Fr(1, 10))
        uf = g + v + cc; w = wf + f0 + f1 + uf
        L1, L = f1 + 4 * T1 + e1, v + 4 * T + e
        Rloc = -f0 - L1 + T1 - S1 - uf - L + v - 4 * cc + T - S
        Chom = -w + wf - 3 * (T1 + T) - 4 * cc - S1 - S - e1 - e
        Df = 6 * w - ez - (p - 1) * zC - 5 * wf - 5 * f1 - T1 - 5 * v - T
        ok &= Rloc == Chom
        ok &= Chom + (Df - mu) / 6 == ((wf - mu) / 6 - Fr(5, 6) * (f1 + v) - Fr(19, 6) * (T1 + T)
                                       - (p - 1) / 6 * zC - ez / 6 - 4 * cc - S1 - S - e1 - e)
    check("A12 R_loc = C_hom and C_hom+(D_f-mu)/6 identity (l.7363-7397, 9365-9373, 10984-10992)", ok)
    # A13 returned width/gap (tw:eq:gstar, tw:eq:greturn; l.9183-9210)
    ok = True
    for _ in range(6):
        wf, f0, f1, g, v, cc, T1, T, e1, e = [rnd(0, Fr(1, 20)) for _ in range(10)]
        m = rnd(Fr(3, 8), Fr(1, 2)); uf = g + v + cc; w = wf + f0 + f1 + uf; tloc = f0 + uf
        L1, L = f1 + 4 * T1 + e1, v + 4 * T + e
        Mst = 1 - 6 * tloc - f1 - L1 - L - 5 * cc + T1 + 5 * v + T
        Tst = 2 * m - L1 - L - 5 * cc
        G0 = 1 - 2 * m - 6 * w
        ok &= Mst - Tst == G0 + 6 * wf + 5 * f1 + 5 * v + T1 + T
        s, zC, z2, a, ao, D0, cC, c2, g2, t2, gaux, mu = [rnd(0, Fr(1, 5)) for _ in range(12)]
        J = D0 - cC - 2 * a + ao; Jp = max(J, 0); bC = wf + zC
        Gp = (Mst - Tst) + 5 * s + J - gaux - g2 + t2 + (cC - bC) + a + (c2 - z2)
        rhs = G0 + 5 * wf + 5 * f1 + 5 * v + T1 + T + mu + 5 * s - zC - a - z2
        ok &= D0 + ao - Jp == min(D0 + ao, cC + 2 * a)
        ok &= Gp - rhs == t2 + (c2 - g2) + (min(D0 + ao, cC + 2 * a) - mu) - (gaux - Jp)
    check("A13 G_*=M_*-T_*=G_0+6w_f+5f_1+5v+T_1+T; G' - (tw:eq:greturn) = t_2+(c_2-g_2)+(min(D_0+a_o,c_C+2a)-mu)-(g_aux-J_+) (l.9186-9210)", ok)
    # A14 moving-radical and terminal identities (App tw:app:moving/additive, l.10966-11010, 11333-11360; Part II l.7425-7445)
    ok = True
    for _ in range(6):
        At, q0, f0, f1, g, v, cc, T1, T, S1, S, e1, e, c, d, u, BC, tau = [rnd(0, Fr(1, 10)) for _ in range(18)]
        uf = g + v + cc; tloc = f0 + uf
        L1, L = f1 + 4 * T1 + e1, v + 4 * T + e
        A = At - f1 - tloc - L1 - L - 5 * cc
        M = 1 - 6 * tloc - f1 - L1 - L - 5 * cc + T1 + 5 * v + T
        q = q0 + f1 + v + cc; mst = M - q
        Rloc = -f0 - L1 + T1 - S1 - uf - L + v - 4 * cc + T - S
        ok &= Rloc + 3 * A - M == 3 * At - 4 + 20 * tloc + 4 * f1 + 4 * cc - 3 * T1 - 16 * v - 3 * T + 3 * mst + 3 * q0 - S1 - S
        q00 = f1 + v + cc
        ok &= (Rloc + 3 * A - Fr(7, 2) * mst - q) - (Rloc + 3 * A - Fr(7, 2) * (M - q00) - q00) == Fr(5, 2) * q0
        ok &= Rloc + A - M == At - 1 + 4 * tloc - 4 * v - L1 - L - 4 * cc - S1 - S
        ok &= Rloc + A - Fr(5, 6) * M == (At - Fr(5, 6) + 3 * tloc - f1 / 6 - Fr(7, 6) * (L1 + L) - Fr(29, 6) * cc
                                         + (T1 + T) / 6 - Fr(19, 6) * v - S1 - S)
        ok &= (T1 + T) / 6 - Fr(7, 6) * (L1 + L) <= 0
        at = A - c; mf = mst - c + d - u; K0 = 2 * at - mf
        ok &= K0 == 2 * A - c - d + u - mst
        ok &= K0 - at - q - u - BC == A - M - d - BC
        ok &= at - q - u - BC == A - M + (mst - c - u) - BC
    ok &= (Fr(-1) + 12 * Fr(1, 25) + Fr(3, 100) == Fr(-49, 100) and Fr(-93, 200) + 6 * Fr(1, 25) == Fr(-9, 40)
           and Fr(-1) + 12 * Fr(1, 25) + Fr(3, 100) + 3 * Fr(1, 200) == Fr(-19, 40)
           and Fr(-1, 2) + 6 * Fr(1, 25) + Fr(7, 200) + Fr(5, 2) * Fr(1, 200) == Fr(-17, 80))
    check("A14 R_loc+3A-M (+3q_0), +5q_0/2, R_loc+A-M<=A_t-1+4t_loc, raw R_loc+A-5M/6, additive-terminal identities; -49/100,-9/40,-19/40,-17/80", ok)
    # A15 Corollary tw:cor:return-bound (l.11497-11568)
    ok = True
    for _ in range(8):
        p = rnd(PLO, PHI); vr = 1 - 1 / p; al = rnd(0, amax(p)); gam = rnd(0, Fr(1, 2))
        W = al * gam; Amax = (1 - gam) + gam / p + W
        tlow = Fr(3, 2) * (Amax + 3 * W - Fr(5, 6))
        mminus = (1 - gam) / 2
        thigh = (1 - 2 + 4 * mminus - p * W) / 2   # from 2+2t-4m_-+pW<1
        Da = 4 + 2 * p * al - 6 * (vr - 4 * al)
        ok &= Amax == 1 - (vr - al) * gam and Amax + 4 * W == 1 - (vr - 5 * al) * gam
        ok &= tlow == Fr(3, 2) * (Fr(1, 6) - (vr - 4 * al) * gam) and thigh == (1 - (2 + p * al) * gam) / 2
        ok &= 4 * (thigh - tlow) == 1 - Da * gam
        ok &= Fr(1, 2) - gam - thigh == p * al * gam / 2
        ok &= 2 + p * al - 6 * (vr - 4 * al) == (p + 24) * (al - amax(p))
        ok &= amax(p) < vr / 5
        am = amax(p)
        ok &= 4 + 2 * p * am == 12 * (vr - 4 * am) and 1 / (2 + p * am) == 1 / (6 * (vr - 4 * am))
        g0 = 1 / (2 + p * am)
        ok &= Fr(3, 2) * (Fr(1, 6) - (vr - 4 * am) * g0) == 0 and (1 - (2 + p * am) * g0) / 2 == 0
    # coverage at alpha = 0.99 alpha_max, midpoint of (gamma_r, gamma_c), p = 112/25 and 9/2 (l.10797-10801)
    for p in (PLO, PHI):
        vr = 1 - 1 / p; al = Fr(99, 100) * amax(p)
        Da = 4 + 2 * p * al - 6 * (vr - 4 * al)
        gc = min(1 / (2 + p * al), 1 / Da) if Da > 0 else 1 / (2 + p * al)
        gr = 1 / (6 * (vr - 4 * al)); gm = (gr + gc) / 2
        tlow = Fr(3, 2) * (Fr(1, 6) - (vr - 4 * al) * gm); thigh = (1 - (2 + p * al) * gm) / 2
        ok &= gr < gc and max(tlow, 0) < min(thigh, Fr(1, 4)) and (vr - 4 * al) * gm > Fr(1, 6) and (vr - 5 * al) * gm > 0
    check("A15 Cor return-bound: A_max, t_low, t_high, 4(t_high-t_low)=1-D_a gamma, 2+pa-6(vr-4a)=(p+24)(a-a_max), coverage <=> a<a_max; margins at 0.99 a_max", ok)
    # A16 Lemma tw:lem:alpha-extremal (l.11586-11620)
    ok = True
    for _ in range(6):
        p = rnd(PLO, PHI); vr = 1 - 1 / p; a0 = amax(p) + rnd(0, Fr(1, 100))
        u, vv = 2 + p * a0, vr - 4 * a0; gam = 1 / u
        ok &= u - 6 * vv == (p + 24) * (a0 - amax(p)) and u - 6 * vv > 0
        ok &= Fr(3, 2) * (Fr(1, 6) - vv * gam) == (u - 6 * vv) / (4 * u) and (1 - u * gam) / 2 == 0
        ok &= (1 - gam) + (1 / p + a0) * gam + 3 * a0 * gam == 1 - vv * gam and vv * gam < Fr(1, 6) and gam < Fr(1, 2)
        t = rnd(0, Fr(1, 4))
        ok &= t + 1 - ((1 - gam) - t) + p * (1 / p + a0) * gam == 1 + 2 * t
    check("A16 alpha-extremal: u-6v=(p+24)(a0-a_max), t_high=0, t_low=(u-6v)/(4u), raw condition fails, reflection = 1+2t", ok)
    # A17 Lemma tw:lem:supremum configuration (l.9300-9318)
    ok = True
    for _ in range(4):
        p = rnd(PLO, PHI); vr = 1 - 1 / p; rho = rnd(rmax(p), Fr(1, 10)); w = Fr(1, 50); Dl = rnd(0, Fr(1, 10 ** 4))
        dm = ((5 - rho) * w + Dl) / vr; m = (1 - dm) / 2; z = dm / p; At = 2 * m + z + w
        ok &= At + 3 * w == 1 - (1 - rho) * w - Dl and At + 3 * w < 1
        ok &= 2 + 2 * Fr(1, 4) - 4 * m + p * w == Fr(1, 2) + 2 * dm + p * w and Fr(1, 2) + 2 * dm + p * w < 1
        ok &= (vr * (1 - 2 * m) - (5 - rho) * w) == Dl and Fr(1, 4) <= m
        ok &= Eal(p, amax(p), rho) > 0
    check("A17 supremum configuration: A_t+3w=1-(1-rho)w-Delta<1, reflection 1/2+2d_m+pw<1, origin margin Delta, E>0 for rho>rho_max", ok)
    # A18 proportional margins at p=9/2 (l.9845-9855)
    p = PHI; am = amax(p)
    ok = am == Fr(16, 171) and p * am == Fr(8, 19)
    for rho in [rnd(0, Fr(24, 289)) for _ in range(4)]:
        up = (p * am * (6 - p - p * rho) / (6 * (p - 1))) / ((7 - 5 * p * am) / 6)
        ok &= up == 8 * (1 - 3 * rho) / 217 and up - rho / 3 == (24 - 289 * rho) / 651
    check("A18 proportional margins at p=9/2: alpha_max=16/171, interval (rho/3, 8(1-3rho)/217) of length (24-289rho)/651", ok)


# =============================================================== B. endpoints
FTW = [-33539275, 46823455, 1111054532, -1437647898, 481077477, -19294821, 184194]
BTW = [36759635, -221603468, 471613016, -400709406, 82242441, 30485106, 1473552]
FEX = [-50, 68, 1657, -2064, 621]
BEX = [37099, -231405, 522846, -507168, 178848]
BPR = [77719, -525330, 1274868, -1308150, 481464]
QNAT = [-334400, -506420, 4106021, -3239236, 615108]
LNAT = [-200, 410, 5771, -7370, 2229]
BNAT = [399857, -2497497, 5643378, -5469552, 1925856]
RSTAR = [61560, -11376, -376469, 218582]
BFLAT = [26231, -167541, 384534, -376596, 133488]


def Fel(C, l):
    return 1600 * C**2 * l**2 - 320 * C**2 * l - 228 * C * l**2 + 184 * C * l - 20 * C - 103 * l**2 - 14 * l + 5


def Ctil(p):
    return -(p * p - 200 * p - 684) / (28 * p * (p + 126))


def Crho(p, rho):
    return (p + 4 - rho) / (4 * p * (5 - rho))


def Cex(p):
    return (p + 4) / (20 * p)


def RC(dl, C):
    exact = isinstance(dl, Fr) and isinstance(C, Fr)
    A = (Fr(5, 6) if exact else mpf(5) / 6) - dl
    return 1 - dl + A * dl * (1 - 4 * C) / (2 * (A * ((Fr(5, 2) if exact else mpf(5) / 2) - 8 * C) + dl * (1 - 4 * C)))


def Qred(dl, C, b, l):
    """(en:origin-Q), l.8402-8406."""
    exact = all(isinstance(t, (Fr, int)) for t in (dl, C, b, l))
    if not exact:
        dl, C, b, l = mpf(dl), mpf(C), mpf(b), mpf(l)
    h = (1 + b + 3 * l) / 2
    Ty = (1 + b - 5 * l) / 2
    q = Fr(1, 4) if exact else mpf(1) / 4
    return -q + b / 6 + 5 * l / 4 + dl / 2 + dl * l / 2 - C * dl * Ty - 3 * h / 4 * (1 - RC(dl, C))


def endpoint(fpoly, lo, hi, Cfun):
    f = lambda l: sum(c * l**i for i, c in enumerate(fpoly))
    l = mpmath.findroot(f, (mpf(lo), mpf(hi)), solver="anderson")
    l = mpmath.findroot(f, l)
    B = mpf(11) / 12 - l / 4
    kap = 2 * B - 1
    p = 6 * kap
    C = Cfun(p)
    r = (5 - 9 * l) / (6 * (1 + l + 8 * C * l))
    d0 = mpmath.diff(lambda d: Qred(d, C, 0, l), r)
    d1 = mpmath.diff(lambda d: Qred(d, C, 1, l), r)
    b = -d0 / (d1 - d0)
    return dict(l=l, B=B, kap=kap, p=p, C=C, r=r, b=b)


def digits_ok(x, s):
    """x agrees with the printed truncated decimal s (all printed digits correct)."""
    t = mpf(s)
    k = len(s.split(".")[1])
    return t - mpf(10) ** (-k) < x < t + 2 * mpf(10) ** (-k)


def sectionB():
    print("== B. Algebraic endpoints B_tw (Part III), B_ex (Part II), B', B_nat, B_* (l.116-218, 5363-5398, 8379-8430, 9501-9580, 10820-10850)")
    ps = [rnd(PLO, PHI) for _ in range(4)]
    check("B1 C~(p) = C_{rho_max(p)}(p) = -(p^2-200p-684)/(28p(p+126)) (l.9523, 10838)",
          all(Crho(p, rmax(p)) == Ctil(p) for p in ps))
    gt = poly_of(lambda l: Fel(Ctil(5 - 3 * l), l) * (28 * (5 - 3 * l) * (5 - 3 * l + 126)) ** 2, 6)
    lx = rnd(0, 1)
    ok = prop(gt, FTW) and peval(gt, lx) == Fel(Ctil(5 - 3 * lx), lx) * (28 * (5 - 3 * lx) * (131 - 3 * lx)) ** 2
    ge = poly_of(lambda l: Fel(Cex(5 - 3 * l), l) * (20 * (5 - 3 * l)) ** 2, 4)
    ok2 = prop(ge, FEX) and peval(ge, lx) == Fel(Cex(5 - 3 * lx), lx) * (20 * (5 - 3 * lx)) ** 2
    check("B2 F(C~(5-3l),l) ~ f_tw(l) (l.9564-9566) and F(C(5-3l),l) ~ f(l) of B_ex (l.8384)", ok and ok2)
    bt = poly_of(lambda B: peval(FTW, Fr(11, 3) - 4 * B), 6)
    be = poly_of(lambda B: peval(FEX, Fr(11, 3) - 4 * B), 4)
    check("B3 l=11/3-4B maps f_tw to the B_tw sextic (tw:eq:Bpoly) and f to the B_ex quartic", prop(bt, BTW) and prop(be, BEX))
    lm = Fr(170432651499230646679121715587598002408611527744903883509836, 10**60)
    ok = (sturm_count(FTW, Fr(17043265, 10**8), Fr(17043266, 10**8)) == 1
          and peval(FTW, lm) * peval(FTW, lm + Fr(1, 10**60)) < 0
          and sturm_count(BTW, Fr(87405850, 10**8), Fr(87405851, 10**8)) == 1
          and sturm_count(FEX, Fr(1703752759993550355, 10**19), Fr(1703752759993550356, 10**19)) == 1
          and sturm_count(BEX, Fr(87407284, 10**8), Fr(87407285, 10**8)) == 1
          and sturm_count(BPR, Fr(87422125, 10**8), Fr(87422126, 10**8)) == 1
          and sturm_count(QNAT, Fr(9, 20), Fr(23, 50)) == 1 and sturm_count(RSTAR, Fr(9, 20), Fr(23, 50)) == 1)
    check("B4 Sturm: unique roots in the stated isolating intervals (l_tw incl. the 1e-60 interval l.10827; B_tw, l_ex, B_ex, B', r_nat, r_*)", ok)
    tw = endpoint(FTW, "0.17043265", "0.17043266", lambda p: Crho(p, rmax(p)))
    ex = endpoint(FEX, "0.1703752759993550355", "0.1703752759993550356", Cex)
    for nm, E in (("tw", tw), ("ex", ex)):
        C, b, l, r = E["C"], E["b"], E["l"], E["r"]
        q0 = Qred(r, C, b, l)
        qd = mpmath.diff(lambda d: Qred(d, C, b, l), r)
        qb = Qred(r, C, 1, l) - Qred(r, C, 0, l)
        rc = RC(r, C) - (5 + 12 * C * r) / 9
        eps = mpf(10) ** -70
        check("B5 contact at B_%s: Q=dQ/ddelta=dQ/db=0 and R_C(r)=(5+12Cr)/9 to 1e-70 (l.9545, 8455-8459)" % nm,
              max(abs(q0), abs(qd), abs(qb), abs(rc)) < eps, "max residual %s" % mpmath.nstr(max(abs(q0), abs(qd), abs(qb), abs(rc)), 3))
    printed_tw = [("l_tw", tw["l"], "0.170432651499230646679121715587"), ("r_tw", tw["r"], "0.444284496630054769285578529082"),
                  ("b_tw", tw["b"], "0.095796105959869730628105439808"), ("kappa_tw", tw["kap"], "0.748117007583718009993772475539"),
                  ("p_tw", tw["p"], "4.48870204550230805"), ("rho_max(p_tw)", rmax(tw["p"]), "0.08362272823194174"),
                  ("alpha_max(p_tw)", amax(tw["p"]), "0.0934865580979076"), ("C_tw", tw["C"], "0.0952172959982581"),
                  ("B_tw", tw["B"], "0.874058503791859004996886237769767"),
                  ("B_tw(60)", tw["B"], "0.87405850379185900499688623776976716606"),
                  ("r_tw(38)", tw["r"], "0.44428449663005476928557852908290979870"),
                  ("b_tw(38)", tw["b"], "0.09579610595986973062810543980817815339")]
    bad = [n for n, x, s in printed_tw if not digits_ok(x, s)]
    check("B6 Part III endpoint numbers (l.9532-9541, 10829-10834) agree with all printed digits", not bad, ",".join(bad))
    printed_ex = [("B_ex", ex["B"], "0.8740728476668279077786693740"), ("r_ex", ex["r"], "0.4446941731145563527466866926"),
                  ("b_ex", ex["b"], "0.0956750076453764202754832248"), ("C_ex", ex["C"], "0.09455460151844812985")]
    bad = [n for n, x, s in printed_ex if not digits_ok(x, s)]
    check("B7 Part II endpoint numbers (l.8414-8419) agree with all printed digits", not bad, ",".join(bad))
    # other endpoints
    root = lambda P, lo, hi: mpmath.findroot(lambda x: sum(c * x**i for i, c in enumerate(P)), (mpf(lo), mpf(hi)), solver="anderson")
    Bp = root(BPR, "0.87422125", "0.87422126")
    rn = root(QNAT, "0.45", "0.46")
    ln = (-mpf(252833046) / 2942707225 * rn**3 + mpf(2932032572) / 2942707225 * rn**2
          - mpf(159506318657) / 76510387850 * rn + mpf(6990506117) / 7651038785)
    Bn = mpf(11) / 12 - ln / 4
    rs = root(RSTAR, "0.45", "0.46")
    ls = (3410 * rs**2 - 7127 * rs + 2520) / (2286 * rs**2 - 5011 * rs + 1800)
    Bs = mpf(11) / 12 - ls / 4
    Bf = root(BFLAT, "0.8748", "0.8749")
    ok = (digits_ok(Bp, "0.87422125049290208328") and digits_ok(rn, "0.45053392487884345250")
          and digits_ok(ln, "0.16879668066330557444") and digits_ok(Bn, "0.87446749650084027305")
          and digits_ok(5 * mpf(1) / 6 - ln / 2, "0.74893499300168054611")
          and abs(sum(c * ln**i for i, c in enumerate(LNAT))) < mpf(10) ** -60
          and abs(sum(c * Bn**i for i, c in enumerate(BNAT))) < mpf(10) ** -55
          and digits_ok(Bs, "0.87446830769353045955") and digits_ok((Bs - Bn) * 10**7, "8.1119269018")
          and digits_ok(Bf, "0.874886647039964"))
    check("B8 B'=0.87422125049290208328, B_nat (and its two minimal polynomials), B_*, B_*-B_nat, B_flat (l.177-218, 5363-5398, 280)", ok)
    # iteration numbers
    B0 = mpf("0.87407285")
    B1 = tw["B"] + (B0 - tw["B"]) / 50
    B1e = ex["B"] + (mpf(5463) / 6250 - ex["B"]) / 50
    k0 = 2 * B0 - 1
    ok = (digits_ok(B1, "0.8740587907160218248") and digits_ok(B1e, "0.8740729907134913496")
          and ex["B"] < B0 and k0 - tw["kap"] < mpf(3) * 10**-5 and 0 < k0 - tw["kap"]
          and mpf("0.083") < rmax(6 * k0) - (k0 - tw["kap"]) / 2 and rmax(6 * tw["kap"]) < mpf("0.0841")
          and Fr(7481457, 10**7) == 2 * Fr(87407285, 10**8) - 1 and tw["B"] > mpf(437) / 500)
    check("B9 B_1=0.87405879071602182 (l.9740), Part II B_1=0.87407299071349135 (l.8644), B_ex<B_0, kappa_0-kappa_tw<3e-5, 0.083<rho_n<0.0841 (l.9696-9700), B_tw>437/500", ok,
          "B_tw-B_ex=%s, B_ex margin to B_0=%s" % (mpmath.nstr(ex["B"] - tw["B"], 6), mpmath.nstr(B0 - ex["B"], 3)))
    # alpha = 1/25 weaker boundary (l.8978-8980)
    def rho_alpha(p, al):
        return (p * al * (6 - p) / (6 * (p - 1))) / ((7 - 5 * p * al) / 18 + p * p * al / (6 * (p - 1)))
    ok = all(rho_alpha(p, amax(p)) == rmax(p) for p in ps)
    l25 = mpmath.findroot(lambda l: Fel(Crho(5 - 3 * l, rho_alpha(5 - 3 * l, mpf(1) / 25)), l), mpf("0.1704"))
    ok &= digits_ok(mpf(11) / 12 - l25 / 4, "0.8740670193055")
    check("B10 with alpha=1/25 the same construction gives 0.8740670193055 (l.8978-8980)", ok,
          "B=%s" % mpmath.nstr(mpf(11) / 12 - l25 / 4, 16))
    return tw, ex


def sectionOpt(tw):
    print("== O. Optimality and the 13/15 limit (l.9743-9893, 10888-10911)")
    C, r, B = tw["C"], tw["r"], tw["B"]
    Dhi = 16 * C * r + 2 * r + 3
    lam = [(44 * C * r + 7 * r + 4) / (3 * Dhi), (4 * C * r - r + 2) / (3 * Dhi), 1 / Dhi]
    ok = all(x > 0 for x in lam) and abs(sum(lam) - 1) < mpf(10)**-80 and abs((44 * C * r + 7 * r + 7) / (3 * Dhi) - B) < mpf(10)**-70
    for _ in range(3):
        S, l = mpf(random.random()), mpf(random.random())
        Ldir, LK = mpf(5) / 6 + S / 12 - l / 6, mpf(1) / 3 + 7 * S / 12 + l / 3
        Hhi = 1 + (2 * C + mpf(1) / 2) * r - (mpf(1) / 2 + 2 * C * r) * S + (2 * C + mpf(1) / 2) * r * l
        ok &= abs(lam[0] * Ldir + lam[1] * LK + lam[2] * Hhi - B) < mpf(10)**-70
    md = 2 * (1 - 12 * C * r) / (9 * r * (1 - 4 * C))
    D0 = 19 - (12 * C + 9) * r
    ok &= mpf("0.39776") < md < mpf("0.39777") and digits_ok(13 / D0, "0.896935") and 3 * (2 - (4 * C + 3) * r) / D0 > 0
    ok &= (18 * C * r + 7 * r + 4) / (2 * D0) > 0 and -mpf(5) / 12 + r * (mpf(3) / 2 - 4 * C) > 0 and r * (1 - 4 * C) > 0
    check("O1 dual identity (tw:eq:dual) = B_tw for all Sigma,l; weights >0; m_d in (0.39776,0.39777); 13/D_hi0=0.896935>B_tw", ok)
    ok = True
    for _ in range(5):
        S, b, l = rnd(0, 1), rnd(0, 1), rnd(0, 1)
        lx = (S - b) / 2; h = 1 - lx + l; ly = (S + b) / 2
        Ldir, LK, J = Fr(5, 6) + S / 12 - l / 6, Fr(1, 3) + 7 * S / 12 + l / 3, S / 2 - l - Fr(4, 3) * b
        ok &= (1 - S + l) - Fr(2, 5) * h == Fr(341, 60) - Fr(35, 6) * Ldir - Fr(2, 3) * LK + Fr(3, 20) * J
        ok &= J == ly - l - Fr(11, 6) * b
        LJ = Fr(23, 24) - S / 6 + Fr(5, 24) * l
        ok &= Fr(17, 30) * Ldir + Fr(1, 30) * LK + Fr(2, 5) * LJ == Fr(13, 15)
        p, rho = rnd(PLO, PHI), rnd(0, Fr(1, 10))
        N = [684, 200, -1]; Dd = [0, 3528, 28]
        dC = peval(psub(pmul(pderiv(N), Dd), pmul(N, pderiv(Dd))), p) / peval(Dd, p) ** 2
        ok &= dC == -(163 * p * p + 684 * p + 43092) / (14 * p * p * (p + 126) ** 2)
        ok &= ((-(5 - rho) + (p + 4 - rho)) / (4 * p * (5 - rho) ** 2)) == (p - 1) / (4 * p * (5 - rho) ** 2)
        Cc = rnd(0, 1)
        ok &= Fel(Fr(1, 4), l) == -12 * l * (5 * l - 1) and Fel(Cc, Fr(1, 5)) == Fr(48, 25) * (4 * Cc - 1)
    ok &= Fr(341, 60) - Fr(13, 2) * Fr(8741, 10000) == Fr(101, 60000) and 684**2 - 4 * 163 * 43092 < 0
    ok &= sturm_count([-103, -228, 1600], Fr(9, 100), Fr(1, 4)) == 0 and peval([-103, -228, 1600], Fr(1, 4)) < 0
    check("O2 separation identity (en:geometry-separation), 101/60000, C~'(p), dC_rho/drho, F(1/4,l), F(C,1/5), a_C<0, 13/15 combination (l.8728-8740, 9773-9781, 9858-9890)", ok)


# =============================================================== C. comparison model, interval branch and bound
class AD:
    """Forward-mode derivative with interval (or exact) components."""
    __slots__ = ("v", "d")

    def __init__(self, v, d):
        self.v, self.d = v, d

    @staticmethod
    def c(x):
        return x if isinstance(x, AD) else AD(x, 0)

    def __add__(s, o):
        o = AD.c(o); return AD(s.v + o.v, s.d + o.d)
    __radd__ = __add__

    def __sub__(s, o):
        o = AD.c(o); return AD(s.v - o.v, s.d - o.d)

    def __rsub__(s, o):
        o = AD.c(o); return AD(o.v - s.v, o.d - s.d)

    def __mul__(s, o):
        o = AD.c(o); return AD(s.v * o.v, s.d * o.v + s.v * o.d)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = AD.c(o); q = s.v / o.v; return AD(q, (s.d - q * o.d) / o.v)

    def __rtruediv__(s, o):
        return AD.c(o).__truediv__(s)

    def __neg__(s):
        return AD(-s.v, -s.d)


def K(q, kind):
    """Constant q in the arithmetic of 'kind' ('iv' or 'fr')."""
    return ivfr(q) if kind == "iv" else Fr(q)


def comp(dl, x, kap, l, b, rho, br, kind="iv", want="Q"):
    """Part II/III comparison model (l.7851-7905, 7920-7960, 8029-8060, 9474-9486):
    E(R_i) - G_N for detector branch i in {o, j}; origin width vr(1-2m)/(5-rho)."""
    if kind == "iv" and any(isinstance(t, AD) for t in (dl, x, kap, l, b, rho)):
        dl, x, kap, l, b, rho = (AD.c(t) for t in (dl, x, kap, l, b, rho))
        k = lambda q: AD(ivfr(q), 0)
    else:
        k = lambda q: K(q, kind)
    p = 6 * kap
    vr = 1 - 1 / p
    c = 5 - rho
    ly = (1 + b - l) / 2
    h = (1 + b + 3 * l) / 2
    C0 = -k(Fr(1, 4)) + b / 6 + 5 * l / 4
    Ty = 2 * ly - h
    if br == "o":
        W0, Wp = vr / c, -2 * vr / c
    elif br == "j":
        W0, Wp = (k(Fr(1, 2)) + (1 - 3 / p)) / (p + 12), -2 * (1 - 3 / p) / (p + 12)
    th0 = dl * x / p
    if br in "oj":
        al = 2 * x / p + 2 * x * W0
        be = 2 - 4 * x / p + 2 * x * Wp
        A = k(Fr(5, 6)) - dl
        H = A / (1 - x) + dl
        mi = (k(Fr(5, 6)) + A * x / (1 - x) - H * al) / (A + H * be)
        Fi = al + be * mi
        R = 1 - dl * Fi
        Ti = mi + (Fi - x) / (1 - x)
    else:  # br == "s": earlier count R_kappa(delta,x) of Prop jdc:count (l.3995-4002)
        cc = 1 / (3 * kap)
        Dk = 3 - (1 + 2 * cc) * x
        Pk = (2 - 2 * cc * x) * (1 - x)
        J = (k(Fr(5, 6)) - dl) * Dk + dl * Pk
        R = 1 - dl + (k(Fr(5, 6)) - dl) * dl * Pk / (2 * J)
        mi = Ti = Fi = None
    E = C0 + dl / 2 + dl * x * (l + (h / 2 - ly) / p) - 3 * h * (1 - R) / 4
    zL = (1 - 2 * k(Fr(73, 200))) / p
    L0 = k(Fr(21, 100))
    a = dl * (L0 * x - zL / 2) / (L0 - zL)
    GN = a * vr * Ty / (2 * c)
    GL = th0 * ((1 - k(Fr(73, 200))) * h - ly)
    if want == "Q":
        return E - GN
    return dict(E=E, GN=GN, GL=GL, R=R, m=mi, T=Ti, F=Fi, p=p, a=a, th0=th0, h=h, ly=ly, vr=vr, c=c)


def bnb(fun, box, accept, init=8, max_boxes=60000):
    """Branch and bound over a box of (lo, hi) Fraction/mpf pairs; fun(list of iv) -> iv.
    Returns (ok, number of boxes, extreme accepted endpoint)."""
    lo0 = [ivfr(a).a if isinstance(a, Fr) else mpf(a) for a, _ in box]
    hi0 = [ivfr(b).b if isinstance(b, Fr) else mpf(b) for _, b in box]
    w0 = [hi0[i] - lo0[i] for i in range(len(box))]
    stack = []
    n0 = [init if (w0[i] > 0 and i < 2) else 1 for i in range(len(box))]
    import itertools
    for idx in itertools.product(*[range(n) for n in n0]):
        stack.append([(lo0[i] + w0[i] * idx[i] / n0[i], lo0[i] + w0[i] * (idx[i] + 1) / n0[i]) for i in range(len(box))])
    n, ext = 0, None
    while stack:
        bx = stack.pop()
        n += 1
        if n > max_boxes:
            return False, n, ext
        val = fun([iv.mpf([a, b]) for a, b in bx])
        good, e = accept(val)
        if good:
            e = lo(e) if hasattr(e, "_mpi_") else e
            ext = e if ext is None else (min(ext, e))
            continue
        j = max(range(len(bx)), key=lambda i: (bx[i][1] - bx[i][0]) / w0[i] if w0[i] > 0 else -1)
        if bx[j][1] - bx[j][0] < mpf(10) ** -12 * (w0[j] if w0[j] > 0 else 1):
            return False, n, ext
        mid = (bx[j][0] + bx[j][1]) / 2
        b1, b2 = list(bx), list(bx)
        b1[j] = (bx[j][0], mid); b2[j] = (mid, bx[j][1])
        stack += [b1, b2]
    return True, n, ext


def below(th):
    t = ivfr(th).a if isinstance(th, (Fr, int)) else mpf(th)
    return lambda v: (v.b < t, t - v.b)


def above(th):
    t = ivfr(th).b if isinstance(th, (Fr, int)) else mpf(th)
    return lambda v: (v.a > t, v.a - t)


DOM = [(Fr(8, 25), Fr(3, 5)), (Fr(9, 25), Fr(1, 2))]
BRIDGE = [(Fr(3, 10), Fr(8, 25)), (Fr(9, 25), Fr(1, 2))]
CRECTS = [[(Fr(1, 8), Fr(5, 6)), (Fr(0), Fr(9, 25))], [(Fr(1, 8), Fr(3, 10)), (Fr(9, 25), Fr(1, 2))],
          [(Fr(3, 5), Fr(5, 6)), (Fr(9, 25), Fr(1, 2))]]


def run(name, fun, box, acc, **kw):
    ok, n, ext = bnb(fun, box, acc, **kw)
    check(name, ok, "%d boxes, worst certified margin %s" % (n, mpmath.nstr(ext, 3) if ext is not None else "-"))


def quad_structure(C, b, l, r):
    """Q_{C,b,l}(delta) * 24(48C delta-80C-18 delta+25) is quadratic in delta: return its coefficients."""
    den = lambda d: 24 * (48 * C * d - 80 * C - 18 * d + 25)
    g = lambda d: Qred(d, C, b, l) * den(d)
    y0, y1, y2 = g(mpf(0)), g(mpf(1)), g(mpf(2))
    A2 = (y2 - 2 * y1 + y0) / 2
    A1 = y1 - y0 - A2
    return A2, A1, y0, g(r + mpf("0.37")) - (A2 * (r + mpf("0.37"))**2 + A1 * (r + mpf("0.37")) + y0)


def lo(x):
    return mp.make_mpf(x._mpi_[0])


def hi(x):
    return mp.make_mpf(x._mpi_[1])


def lift(t, *refs):
    return AD.c(t) if any(isinstance(r, AD) for r in refs) else t


def cform(f, B):
    """Mean-value (first-order Taylor) enclosure of f over the box B, intersected with
    the natural interval extension; derivatives by forward-mode interval AD."""
    nat = f(B)
    C = [b.mid for b in B]
    enc = f(C)
    for i in range(len(B)):
        if hi(B[i]) == lo(B[i]):
            continue
        args = list(B)
        args[i] = AD(B[i], 1)
        g = f(args)
        g = g.d if isinstance(g, AD) else 0
        enc = enc + g * (B[i] - C[i])
    return iv.mpf([max(lo(nat), lo(enc)), min(hi(nat), hi(enc))])


def sectionC(tw, ex):
    print("== C. Comparison inequalities by interval branch and bound (Lemma tw:lem:family l.9592-9675; Lemma en:endpoint-signs l.8436-8506)")
    ok = True
    for _ in range(4):
        dl, kap, l, b, rho = rnd(Fr(8, 25), Fr(3, 5)), rnd(Fr(3, 4) - Fr(1, 100), Fr(3, 4)), rnd(Fr(1, 6), Fr(1, 5)), rnd(0, Fr(1, 10)), rnd(0, Fr(1, 10))
        ok &= comp(dl, Fr(1, 2), kap, l, b, rho, "o", kind="fr") == Qred(dl, Crho(6 * kap, rho), b, l)
    check("C1 model: Q_o(delta,1/2) = Q_{C_rho(6 kappa),b,l}(delta) exactly (l.8442-8446, 9599-9601)", ok)
    for nm, E, kap0, ltag in (("tw", tw, Fr(7481457, 10**7), "Part III"), ("ex", ex, Fr(4676, 6250), "Part II")):
        rho_of = (lambda kap: rmax(6 * kap)) if nm == "tw" else (lambda kap: 0)
        lI, bI = ivnear(E["l"]), ivnear(E["b"])
        kI = 5 * ivfr(Fr(1, 6)) - lI / 2
        rhoI = rho_of(kI)
        A2, A1, A0, resid = quad_structure(E["C"], E["b"], E["l"], E["r"])
        disc = A1 * A1 - 4 * A2 * A0
        rt = -A1 / (2 * A2)
        ok = A2 < 0 and abs(disc) < mpf(10) ** -60 * A1 * A1 and abs(rt - E["r"]) < mpf(10) ** -30 and abs(resid) < mpf(10) ** -70
        ok &= all(48 * E["C"] * d - 80 * E["C"] - 18 * d + 25 > 0 for d in (mpf(8) / 25, mpf(3) / 5))
        check("C2 %s endpoint: Q_o(delta,1/2)=A(delta-r)^2/den with A<0, den>0 on [8/25,3/5] (90 digits)" % nm, ok,
              "A=%s, disc/A1^2=%s" % (mpmath.nstr(A2, 5), mpmath.nstr(disc / (A1 * A1), 3)))
        half = ivfr(Fr(1, 2))
        run("C3 %s endpoint: Q_j(delta,1/2) < 0 on [8/25,3/5]" % nm,
            lambda B: cform(lambda A: comp(A[0], half, kI, lI, bI, rhoI, "j"), B), [DOM[0]], below(0))
        for br in "oj":
            run("C4 %s endpoint: dQ_%s/dx > 1/1000 on D=[8/25,3/5]x[9/25,1/2] (l.9617, 8472)" % (nm, br),
                lambda B, br=br: comp(B[0], AD(B[1], 1), kI, lI, bI, rhoI, br).d, DOM, above(Fr(1, 1000)))
        kfam = [(lo(kI), ivfr(kap0).b)]

        def famd(B, br):
            kap = AD(B[2], 1)
            rho = rho_of(kap)
            return comp(B[0], B[1], kap, -(kap - kI) / 25 + lI, bI, rho, br).d
        for br in "oj":
            run("C5 %s family: dQ_%s/dkappa < -1/100 on D x [kappa_end,kappa_0] along l(kappa) (l.9618-9620, 8473-8475)" % (nm, br),
                lambda B, br=br: famd(B, br), DOM + kfam, below(Fr(-1, 100)), init=6)
        if nm == "tw":
            for br in "oj":
                run("C6 tw: dQ_%s/drho > -1/100 on D x P_tw, rho in [0.083,0.0841] free (tw:eq:rho-derivative)" % br,
                    lambda B, br=br: comp(B[0], B[1], B[2], lI - (B[2] - kI) / 25, bI, AD(B[3], 1), br).d,
                    DOM + kfam + [(Fr(83, 1000), Fr(841, 10000))], above(Fr(-1, 100)), init=6)

        def fam(A, br, key):
            d = comp(A[0], A[1], A[2], -(A[2] - kI) / 25 + lI, bI, rho_of(A[2]), br, want="all")
            return key(d, A)
        tests = [("E(R_%s)-G_L < -1e-7 on D", "oj", DOM, lambda d, A: d["E"] - d["GL"], below(mpf(-1) / 10**7)),
                 ("m_%s-19/50 > 0 on D", "oj", DOM, lambda d, A: d["m"] - ivfr(Fr(19, 50)), above(0)),
                 ("9/20-m_%s > 0 on D", "oj", DOM, lambda d, A: -d["m"] + ivfr(Fr(9, 20)), above(0)),
                 ("T_%s-m-3/5 > 0 on D", "oj", DOM, lambda d, A: d["T"] - d["m"] - ivfr(Fr(3, 5)), above(0)),
                 ("1-(T_%s-m) > 0 on D", "oj", DOM, lambda d, A: 1 - d["T"] + d["m"], above(0)),
                 ("T_%s-1 > 0 on D", "oj", DOM, lambda d, A: d["T"] - 1, above(0)),
                 ("3/2-T_%s > 0 on D", "oj", DOM, lambda d, A: -d["T"] + ivfr(Fr(3, 2)), above(0)),
                 ("9/10+x/(5p)-F_%s(m) > 0 on D", "oj", DOM, lambda d, A: -d["F"] + lift(A[1], d["p"]) / (5 * d["p"]) + ivfr(Fr(9, 10)), above(0)),
                 ("bridge E(R_%s)-G_N < -1e-7", "s", BRIDGE, lambda d, A: d["E"] - d["GN"], below(mpf(-1) / 10**7)),
                 ("bridge E(R_%s)-G_L < -1e-7", "s", BRIDGE, lambda d, A: d["E"] - d["GL"], below(mpf(-1) / 10**7))]
        for i, rect in enumerate(CRECTS):
            tests.append(("complementary E(R_%%s) < -1e-5 on C_%d" % (i + 1), "s", rect, lambda d, A: d["E"], below(mpf(-1) / 10**5)))
        for desc, brs, dom, key, acc in tests:
            for br in brs:
                run("C7 %s family (%s, Table): " % (nm, ltag) + (desc % br),
                    lambda B, br=br, key=key: cform(lambda A: fam(A, br, key), B), dom + kfam, acc, init=6)

        def geo(A, k):
            m, kap = A[0], A[1]
            p = 6 * kap; vr = 1 - 1 / p; c = 5 - rho_of(kap)
            Wo = vr * (1 - 2 * m) / c
            Wj = (ivfr(Fr(1, 2)) + (1 - 3 / p) * (1 - 2 * m)) / (p + 12)
            W = iv.mpf([min(lo(Wo), lo(Wj)), min(hi(Wo), hi(Wj))])   # encloses W_d = min(W_o, W_j)
            tl = 6 * W - 3 * vr * (1 - 2 * m) / 2 + ivfr(Fr(1, 4))
            th = 2 * m - ivfr(Fr(1, 2)) - p * W / 2
            return [tl, ivfr(Fr(1, 4)) - tl, th, tl - 5 * W, W, ivfr(Fr(1, 25)) - W, ivfr(Fr(21, 100)) - (1 - 2 * m) / p - W][k]
        names = ["t_low>0", "t_low<1/4", "t_high>0", "t_low>5W", "W>0", "W<1/25", "L_0-z_0-W>0"]
        okg, nb, worst = True, 0, None
        for k in range(7):
            o, n, e = bnb(lambda B, k=k: geo(B, k), [(Fr(19, 50), Fr(9, 20))] + kfam, above(0), init=32)
            okg &= o; nb += n
            if not o:
                print("   failed component", names[k])
            worst = e if worst is None else min(worst, e)
        check("C8 %s family: detector geometry %s at W=W_d=min(W_o,W_j) on [19/50,9/20] (l.7905-7925, 9640-9650)" % (nm, ", ".join(names)), okg, "%d boxes, worst margin %s" % (nb, mpmath.nstr(worst, 3)))



def sectionC9(tw, ex):
    print("== C9. Remaining entries of Table en:iteration-table on the families (l.8539-8575, with c=5-rho for Part III)")
    for nm, E, kap0 in (("tw", tw, Fr(7481457, 10**7)), ("ex", ex, Fr(4676, 6250))):
        rho_of = (lambda kap: rmax(6 * kap)) if nm == "tw" else (lambda kap: 0)
        lI, bI = ivnear(E["l"]), ivnear(E["b"])
        kI = 5 * ivfr(Fr(1, 6)) - lI / 2
        K0 = iv.mpf([lo(kI), ivfr(kap0).b])
        kap = K0; l = -(K0 - kI) / 25 + lI; b = bI; rho = rho_of(K0)
        p = 6 * kap; vr = 1 - 1 / p; c = 5 - rho
        lx = (1 - b - l) / 2; ly = (1 + b - l) / 2; h = (1 + b + 3 * l) / 2; C0 = -ivfr(Fr(1, 4)) + b / 6 + 5 * l / 4; Ty = 2 * ly - h
        tn = ivfr(Fr(23, 100)); mth = 1 - ly / h
        sc = {"p l-3/4": p * l - ivfr(Fr(3, 4)), "l/(h+1e-9)-21/100": l / (h + ivfr(Fr(1, 10**9))) - ivfr(Fr(21, 100)),
              "2/9-l/(h-1/80)": ivfr(Fr(2, 9)) - l / (h - ivfr(Fr(1, 80))), "l_x-l": lx - l, "l_y-l-11b/6": ly - l - 11 * b / 6,
              "1/2-l_y": ivfr(Fr(1, 2)) - ly, "T_y": Ty,
              "m_th-(p-1+c(1+2t_n))/(2(p-1)+4c)": mth - (p - 1 + c * (1 + 2 * tn)) / (2 * (p - 1) + 4 * c),
              "59/2000-(vr/c)(2l_y/h-1)": ivfr(Fr(59, 2000)) - vr / c * (2 * ly / h - 1),
              "vr(1-4/c)/10-(1/6-2t_n/3)": vr * (1 - 4 / c) / 10 - (ivfr(Fr(1, 6)) - 2 * tn / 3), "9/20-m_th": ivfr(Fr(9, 20)) - mth,
              "-Q_sm": -(l * (ivfr(Fr(1, 2)) + 2 * l / h) - 11 * h / 54), "-1/250-C0-(1-h)/16": -ivfr(Fr(1, 250)) - C0 - (1 - h) / 16,
              "-1/250-C0-(1/2+3l/2)/50": -ivfr(Fr(1, 250)) - C0 - (ivfr(Fr(1, 2)) + 3 * l / 2) / 50,
              "-1/20-(13h/75-l_y/2+1/50)": -ivfr(Fr(1, 20)) - (13 * h / 75 - ly / 2 + ivfr(Fr(1, 50))),
              "B-437/500": ivfr(Fr(11, 12)) - l / 4 - ivfr(Fr(437, 500))}
        bad = [k for k, v in sc.items() if not lo(v) > 0]
        worst = min(sc.items(), key=lambda kv: lo(kv[1]))
        check("C9 %s family: 16 scalar entries (p l>3/4, l/h bounds, geometry, origin branch at m_theta, 59/2000 cap, exceptional for t_n, Q_sm, small/lowest amplitudes, small parameters, B>437/500)" % nm,
              not bad, "smallest: %s = %s" % (worst[0], mpmath.nstr(lo(worst[1]), 3)) + ((" FAILED: " + ",".join(bad)) if bad else ""))
        kfam = [(lo(kI), ivfr(kap0).b)]

        def F2(A, br, key):
            isad = any(isinstance(t, AD) for t in A)
            k = (lambda q: AD(ivfr(q), 0)) if isad else ivfr
            dl, x, kp = (AD.c(t) for t in A) if isad else A
            d = comp(dl, x, kp, -(kp - kI) / 25 + lI, bI, rho_of(kp), br, want="all")
            R, th0 = d["R"], d["th0"]
            S = (3 * R + 1) / 4 + th0 / 2 - k(Fr(17, 50))
            Dsh = R + (-th0 + k(Fr(1, 2))) * k(Fr(3, 40)) - ((3 * R + 1) / 4 - th0 / 2)
            return key(d, R, S, Dsh, dl, x, d["p"], k)
        lst = [("E(R)<1/450", below(Fr(1, 450)), lambda d, R, S, Ds, dl, x, p, k: d["E"]),
               ("S>1/5", above(Fr(1, 5)), lambda d, R, S, Ds, dl, x, p, k: S),
               ("S<1", below(1), lambda d, R, S, Ds, dl, x, p, k: S),
               ("E+h D_short<-1e-5", below(mpf(-1) / 10**5), lambda d, R, S, Ds, dl, x, p, k: d["E"] + d["h"] * Ds),
               ("S+D_short>0", above(0), lambda d, R, S, Ds, dl, x, p, k: S + Ds),
               ("S+D_short<2", below(2), lambda d, R, S, Ds, dl, x, p, k: S + Ds),
               ("(p-2 delta x)L_0-R>0", above(0), lambda d, R, S, Ds, dl, x, p, k: (p - 2 * dl * x) * k(Fr(21, 100)) - R),
               ("G_L-G_N>0", above(0), lambda d, R, S, Ds, dl, x, p, k: d["GL"] - d["GN"]),
               ("a>0", above(0), lambda d, R, S, Ds, dl, x, p, k: d["a"]),
               ("2a/p-theta_0>0", above(0), lambda d, R, S, Ds, dl, x, p, k: 2 * d["a"] / p - d["th0"]),
               ("delta/2-theta_0-vr a/c>0", above(0), lambda d, R, S, Ds, dl, x, p, k: dl / 2 - d["th0"] - d["vr"] * d["a"] / d["c"]),
               ("a(59/2000)h/2-G_N>0", above(0), lambda d, R, S, Ds, dl, x, p, k: d["a"] * k(Fr(59, 2000)) * d["h"] / 2 - d["GN"]),
               ("(delta/2-theta_0)(l_y-11h/20)-G_N>0", above(0), lambda d, R, S, Ds, dl, x, p, k: (dl / 2 - d["th0"]) * (d["ly"] - 11 * d["h"] / 20) - d["GN"]),
               ("2xL_0-z_L>31/1000", above(Fr(31, 1000)), lambda d, R, S, Ds, dl, x, p, k: 2 * x * k(Fr(21, 100)) - (1 - 2 * k(Fr(73, 200))) / p)]
        okall, nb, fails = True, 0, []
        for desc, acc, key in lst:
            for dom, brs in ((DOM, "oj"), (BRIDGE, "s")):
                for br in brs:
                    o, n, _ = bnb(lambda B, br=br, key=key: cform(lambda A: F2(A, br, key), B), dom + kfam, acc, init=6)
                    okall &= o; nb += n
                    if not o:
                        fails.append(desc + "/" + br)
        check("C9 %s family: 14 two-dimensional Table entries on main (both detector branches) and bridge rectangles" % nm, okall,
              "%d boxes" % nb + ((" FAILED: " + ",".join(fails)) if fails else ""))
        full = [(Fr(1, 8), Fr(5, 6)), (Fr(0), Fr(1, 2))] + kfam
        o1, n1, _ = bnb(lambda B: cform(lambda A: F2(A, "s", lambda d, R, S, Ds, dl, x, p, k: S), B), full, above(0), init=8)
        o2, n2, _ = bnb(lambda B: cform(lambda A: F2(A, "s", lambda d, R, S, Ds, dl, x, p, k: 1 - S), B), full, above(0), init=8)
        o3, n3, _ = bnb(lambda B: cform(lambda A: F2(A, "s", lambda d, R, S, Ds, dl, x, p, k: (p - 2 * dl * x) * k(Fr(21, 100)) - R), B), full, above(0), init=8)
        check("C9 %s family: full class rectangle [1/8,5/6]x[0,1/2]: 0<S(R_s)<1 and (p-2 delta x)L_0-R_s>0" % nm, o1 and o2 and o3, "%d boxes" % (n1 + n2 + n3))


def main():
    import time
    t0 = time.time()
    sectionA()
    tw, ex = sectionB()
    sectionOpt(tw)
    sectionC(tw, ex)
    sectionC9(tw, ex)
    # Part II rational starting comparison B = 5463/6250 (Prop en:scalar-certificate, l.8162-8180)
    print("== D. Part II rational start 0.87408 (l.7851-7861, 8162-8180)")
    kap, l, b = Fr(37446831, 50000000), Fr(1597, 9375), Fr(1913, 20000)
    ok = Fr(11, 12) - l / 4 == Fr(5463, 6250) and (1 + kap) / 2 == Fr(87446831, 10**8) and (1 + b - l) / 2 == Fr(277591, 600000) and (1 + b + 3 * l) / 2 == Fr(160669, 200000)
    check("D1 rational parameters: B=11/12-l/4=0.87408, (1+kappa)/2=0.87446831>B_*, l_y, h (l.7851-7867)", ok)
    kI, lI, bI = ivfr(kap), ivfr(l), ivfr(b)
    for br in "oj":
        for key, nm in ((lambda d: d["E"] - d["GN"], "G_N"), (lambda d: d["E"] - d["GL"], "G_L")):
            run("D2 E(R_%s)-%s < -1e-6 on D" % (br, nm),
                lambda B, br=br, key=key: cform(lambda A: key(comp(A[0], A[1], kI, lI, bI, 0, br, want="all")), B),
                DOM, below(mpf(-1) / 10**6))
    for key, nm in ((lambda d: d["E"] - d["GN"], "G_N"), (lambda d: d["E"] - d["GL"], "G_L")):
        run("D3 bridge: E(R_s)-%s < -1e-6" % nm,
            lambda B, key=key: cform(lambda A: key(comp(A[0], A[1], kI, lI, bI, 0, "s", want="all")), B),
            BRIDGE, below(mpf(-1) / 10**6))
    for i, rect in enumerate(CRECTS):
        run("D4 complementary E(R_s) < -1e-5 on C_%d" % (i + 1),
            lambda B: cform(lambda A: comp(A[0], A[1], kI, lI, bI, 0, "s", want="all")["E"], B), rect, below(mpf(-1) / 10**5))
    nf = sum(1 for _, okk in RESULTS if not okk)
    print("\n%d checks, %d failed, %.0f s" % (len(RESULTS), nf, time.time() - t0))
    return nf


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
