#!/usr/bin/env python3
# Explicit constant in Holland's hyperbolicity wedge (arXiv 2608.08682), certificate.
# Run: uv run --no-project --with python-flint --with sympy python3 certify.py [Lstar] [K6] [NEPS]
# All inequalities are evaluated in Arb ball arithmetic (python-flint). "chk(cond, msg)" aborts on failure.
import sys
from flint import arb, acb, ctx
import sympy as sp
ctx.prec = 256

def I(lo, hi): return arb(lo).union(arb(hi))
def up(x): return arb(x.upper())
def lo_(x): return arb(x.lower())
def chk(c, msg):
    if not c:
        print("FAILED:", msg); sys.exit(1)
from fractions import Fraction
from decimal import Decimal, Context, ROUND_CEILING, ROUND_FLOOR
def _exact(a):                      # exact rational value of an exact (zero-radius) arb
    m, e = a.man_exp(); return Fraction(int(m)) * Fraction(2) ** int(e)
TINY = arb(10) ** -300
def fup(x, sig=5):                  # decimal string >= every point of the ball x (outward, upward)
    if x.lower() >= 0 and x.upper() < TINY: return "1e-300"   # astronomically small terms (e.g. e^{-1e41}): certified < 1e-300
    v = _exact(x.upper()); c = Context(prec=sig, rounding=ROUND_CEILING)
    d = c.divide(Decimal(v.numerator), Decimal(v.denominator)); assert Fraction(d) >= v; return str(d)
def flo(x, sig=5):                  # decimal string <= every point of the ball x (outward, downward)
    assert abs(x).upper() > TINY or x.lower() >= 0
    v = _exact(x.lower()); c = Context(prec=sig, rounding=ROUND_FLOOR)
    d = c.divide(Decimal(v.numerator), Decimal(v.denominator)); assert Fraction(d) <= v; return str(d)
DOC = []                            # documented constants, all rounded outward
pi = arb.pi()
Lstar = arb(sys.argv[1]) if len(sys.argv) > 1 else arb(100)
K6 = arb(sys.argv[2]) if len(sys.argv) > 2 else arb('9.5')
NEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 400
LAM = arb(sys.argv[4]) if len(sys.argv) > 4 and not sys.argv[4].startswith('-') else arb(2)   # Cauchy radius lambda*r (Holland: 2)
VERB = '-v' in sys.argv
d0 = arb(16140857488599731839205676)   # our all-n degree record (COROLLARIES): J^{d,n} hyperbolic for d <= d0, all n
report = {}

# ---------------- Section 0: regime ----------------
# u_0(x): the root of Z(u) = u (pi e^{2u} - 9/4) = x (Z increasing on u>0). L := 2 u_0(n), eps := 1/L.
a_s = Lstar / 2
Nstar = a_s * (pi * (2 * a_s).exp() - arb(9) / 4)        # n >= N*  <=>  u_0(n) >= a_s  <=>  eps <= 1/Lstar
nu = up(1 / Nstar)                                       # 1/n <= nu
Lnu = up(Nstar.log() / Nstar)                            # L/n <= log(n)/n <= log(N*)/N*  (n >= e^{2u_0} >= N*)
L5nu = up(Nstar.log() ** 5 / Nstar)                      # L^5/n bound
epsmax = 1 / Lstar
report['Nstar'] = Nstar

# ---------------- Section A: explicit complex saddle on D(x0) = {|z-x0| <= x0/10} ----------------
w1 = up(9 * (1 + arb('0.06') / a_s) / (arb('3.6') * (pi * (2 * a_s).exp() - arb(9) / 4)))
beta = up((-(arb('0.9').log()) - (1 - arb('0.06') / a_s).log() + 4 * w1) / 2)   # |u0(z)-u0(x0)| <= beta
chk(beta < arb('0.06'), "A1 disk invariance")
chk(1 / (2 * (a_s - arb('0.06')) * (1 - w1)) < 1, "A1 contraction")
alo = lo_(a_s - beta)                                     # Re u0(z) >= alo, |u0| >= alo
cos2b = lo_((2 * beta).cos())
chk(arb('0.1').asin() + 2 * (2 * beta / alo).atan() < pi / 2, "A2 Re(z/u^2)>0 on the ray")
phi2 = up(1 + (1 + arb(9) / (4 * pi * (2 * alo).exp())) / (2 * alo))      # |f''(u0)| <= 4 pi e^{2a} phi2
s1 = arb('0.3')
c3 = up((4 * (alo + beta) * (pi + arb(9) / 4 * (-2 * alo).exp()) / (alo - s1) ** 3 + 8 * pi * (2 * s1).exp()) / 6)
chk(c3 * s1 <= pi * cos2b, "A2 cubic remainder absorbed (C3 s1 <= kappa0/4)")
rel_core = up(c3 * (2 * phi2).sqrt() / (pi ** 2 * cos2b ** 2) * (-alo).exp())
E2a = (2 * alo).exp()
rel_gt = up(2 * (-2 * pi * cos2b * s1 ** 2 * E2a).exp() * (2 * phi2).sqrt() / (4 * pi * cos2b * s1 * alo.exp()))
E2r = (2 * alo - 2 * s1).exp()
rel_rt = up(2 * (-2 * pi * cos2b * s1 ** 2 * E2r).exp() * (2 * phi2).sqrt() * alo.exp() / (4 * pi * cos2b * E2r * s1))
rel_am = up(3 / (2 * pi) * (-alo).exp() * ((phi2 * (2 * s1).exp() / cos2b).sqrt() + rel_rt))
delta_main = up(rel_core + rel_gt + rel_rt + rel_am)
brk = up(arb('0.9') * (arb(1) / 4 + beta ** 2 / alo ** 2).log() + arb('0.6') * beta / alo + arb('1.1') / alo)
chk(brk < 0, "A2 pieces exponent negative")
rel_pc = up(arb('5.2').exp() * Nstar * (brk * Nstar).exp())
Ctheta = arb(0)
for k in range(2, 12):
    Ctheta += (2 * pi ** 2 * k ** 4 + 3 * pi * k ** 2) * (-pi * (k * k - 4)).exp()
Ctheta = up(Ctheta + (-pi * 100).exp())                    # tail k>=12 is < e^{-100 pi}
q2max = up(4 * (1 + arb('1e-30')) / (2 - (1 + arb('1e-30')) / alo))
gnet = lo_(arb('0.9') * arb(2).log() / (1 + beta / a_s) - q2max / (360 * (1 - beta / a_s)))
chk(gnet > 0, "A2 theta tail gain beats deficit")
rel_k = up(Ctheta * arb(2) ** arb('-4.5') / (2 * pi ** 2) * ((-gnet * (pi * (2 * a_s).exp() - arb(9) / 4)).exp()
            * ((2 * beta + 1).exp() * phi2.sqrt() + 1)))
deltaA = up(delta_main + rel_pc + rel_k)
chk(deltaA < arb('0.01'), "A2 relative error small")
etaG = up((-(arb('0.9').log()) - (1 - beta / a_s).log() + 2 * (beta / (2 * alo * a_s) + 2 * w1)) / 2 + 4 * deltaA)
report.update(beta=beta, deltaA=deltaA, etaG=etaG)
DOC.append("A: w1 <= %s, beta <= %s, alo >= %s, cos(2beta) >= %s, phi2 <= %s, c3 <= %s (limit (4pi/3)e^0.6 = %s..)" % (
    fup(w1, 3), fup(beta, 6), flo(alo, 6), flo(cos2b, 6), fup(phi2, 6), fup(c3, 6), flo(4 * pi / 3 * (2 * s1).exp(), 6)))
DOC.append("A: rel_core <= %s, rel_gt <= %s, rel_rt <= %s, rel_am <= %s, rel_pc <= %s, rel_k <= %s, deltaA <= %s" % tuple(
    fup(v, 3) for v in (rel_core, rel_gt, rel_rt, rel_am, rel_pc, rel_k, deltaA)))
DOC.append("A: b_rk <= %s, gnet >= %s, C_theta <= %s, q2max <= %s, etaG <= %s" % (fup(brk, 5), flo(gnet, 5), fup(Ctheta, 6), fup(q2max, 6), fup(etaG, 5)))
if VERB: print("A: beta", beta, "deltaA", deltaA, "etaG", etaG, "c3", c3, "brk", brk, "gnet", gnet)

# ---------------- Section B: Psi^{(k)} = Q_k(t,s) / (z^{k-1} * 2u),  t = 1/u0(z), s = 9/(4 pi e^{2u0}) -------------
u, E, t, s = sp.symbols('u E t s', positive=True)
Z = u * (E - sp.Rational(9, 4))
def Du(e): return sp.diff(e, u) + sp.diff(e, E) * 2 * E
Zp = Du(Z)
Psi = {2: 2 / (u * Zp)}
for k in range(2, 5):
    Psi[k + 1] = sp.together(Du(Psi[k]) / Zp)
QP = {}
for k in range(2, 6):
    q = sp.cancel(sp.together((Psi[k] * Z ** (k - 1) * 2 * u).subs({u: 1 / t, E: sp.Rational(9, 4) / s})))
    num, den = sp.fraction(q)
    QP[k] = (sp.Poly(num, t, s).terms(), sp.Poly(den, t, s).terms())
def evalpoly(terms, tv, sv, one):
    acc = one * 0
    for (i, j), c in terms:
        c = sp.Rational(c)
        acc += (one * int(c.p) / int(c.q)) * tv ** i * sv ** j
    return acc
def Qk(k, tv, sv, one):
    nt, dt = QP[k]
    return evalpoly(nt, tv, sv, one) / evalpoly(dt, tv, sv, one)
# complex polydisc bound for Q5 (used on n + Omega_F): |t| <= 1/(alo - 0.01), |s| <= 9/(4 pi) e^{-2 alo + 0.02}
tb = acb(0) + acb(arb(0, up(1 / (alo - arb('0.01')))), arb(0, up(1 / (alo - arb('0.01')))))
sb = acb(0) + acb(arb(0, up(arb(9) / (4 * pi) * (-2 * alo + arb('0.02')).exp())), arb(0, up(arb(9) / (4 * pi) * (-2 * alo + arb('0.02')).exp())))
q5 = up(abs(Qk(5, tb, sb, acb(1))))
report['q5'] = q5
DOC.append("B: |Q5| <= %s on |t| <= 1/(alo - 0.01), |s| <= (9/4pi) e^{-2 alo + 0.02}" % fup(q5, 5))
if VERB: print("B: |Q5| <=", q5)

# ---------------- Sections C-G over eps-subintervals ----------------
F = lambda k: arb(1) if k <= 1 else arb(k) * F(k - 1)
gk = {k: up(F(k) * etaG * (10 / (1 - 50 * nu)) ** k * Lnu) for k in (2, 3, 4)}
omega = up((d0 / Nstar).sqrt())   # sqrt(d/n) <= sqrt(d0/N*)  (wedge + K >= K_N)
worst = {}
def upd(name, val, mode):
    v = val.upper() if mode == 'max' else val.lower()
    if name not in worst: worst[name] = v
    else: worst[name] = max(worst[name], v) if mode == 'max' else min(worst[name], v)

for i in range(NEPS):
    el, eh = epsmax * i / NEPS, epsmax * (i + 1) / NEPS
    e = I(el, eh)
    tR = I(lo_(2 * el * (1 - 5 * nu)), up(2 * eh))
    sR = I(0, up(arb(9) / (4 * pi) * (-1 / eh).exp()))
    one = arb(1)
    fac = {k: I(lo_((1 + 5 * nu) ** (-(k - 1)) * (1 - 5 * nu)), 1) for k in (2, 3, 4)}
    X0 = Qk(2, tR, sR, one) * fac[2] + I(-gk[2], gk[2])     # n L delta_0   (delta_j = second difference of h at n+j)
    X1 = X0                                                  # same enclosure for n L delta_1
    Y = Qk(3, tR, sR, one) * fac[3] + I(-gk[3], gk[3])      # n^2 L (delta_1 - delta_0) = n^2 L Delta^3 h(n)
    W = Qk(4, tR, sR, one) * fac[4] + I(-gk[4], gk[4])      # n^3 L h''''(n+z), 0<=z<=4
    lam0 = I(1 - nu, 1) - e * X0                             # n log q0
    lam1 = I(1 - 2 * nu, 1) - e * X1                         # n log q1
    chk(lam0 > 0 and lam1 > 0, "q0>1, q1>1")
    ex = I(1, 1 + 2 * nu)
    nq0m1 = lam0 * ex; nq1m1 = lam1 * ex; q0 = ex
    # D_n = (delta1-delta0) + 2/(b+2) delta0 + (1/2) F''(xi) delta0^2,  F(x)=log(2-e^{-x}),  (1/2)F'' in [-1, -1+3/n]
    n2LDn = Y + 2 * I(1 - 3 * nu, 1) * X0 - I(1 - 3 * nu, 1) * e * X0 ** 2
    chk(n2LDn > 0, "Laguerre defect D_n > 0 (Lemma third-defect)")
    nq1Lm1 = nq0m1 / q0
    n2Lq1Lmq1 = I(1, 1 + 2 * nu) * n2LDn * I(1 - nu, 1)
    nLTmq0 = nq0m1 * n2Lq1Lmq1 / (nq1m1 * nq1Lm1)
    chk(nLTmq0 > 0, "T > q0")
    aJ = 2 * q0 / nLTmq0                                     # A_J / (n L)
    YB = I(1, 1 + 2 * nu) * X0 * I(1 - nu, 1) - q0 / aJ      # n L (1/b - 1/B_J)
    chk(YB > 0, "B_J > b")
    betaJ = 1 / (I(1 - nu, 1) - e * YB)                      # B_J / n
    bJ = I(1, 1 + nu) * betaJ * YB                           # (B_J - b) L / n
    # n^3 L [psi'''(B_J+z) - psi'''(b+z)] in [-(1+4nu) M, -M],  M = 2 bJ p q (p^2+pq+q^2), p = n/(b+z), q = n/(B_J+z)
    pp = I(1 - 5 * nu, 1)
    qq = 1 / (1 / pp + bJ * e)
    M = 2 * bJ * pp * qq * (pp ** 2 + pp * qq + qq ** 2)
    E4 = W - M * I(1, 1 + 4 * nu) - I(0, 2 * e ** 2 / aJ ** 3 * (1 + 3 * nu))
    tQ = -E4                                                 # n^3 L Q_n
    chk(tQ > 0, "Jacobi defect Q_n > 0 (Lemma signed-fourth-defect)")
    # ---- four-ratio matching via Brouwer on a box ----
    box = (I('1.10', '1.30'), I('0.25', '0.42'), I('0.015', '0.07'))
    def S(bx):
        a, b, c = bx
        lam = I(1 - nu, 1)
        D0 = (1 / aJ - 1) * (1 - 1 / a) + (1 / aJ) * ((lam - 1) - (lam - 1) / a)
        tB = betaJ + I(0, 1) * b * e
        tC = betaJ / 2 + I(0, 1) * c * e
        nn = I(0, nu)
        P = lambda tau, m: tau * (tau + nn) * (tau + 2 * nn) * (tau + 3 * nn) if m == 3 else (tau * (tau + nn) * (tau + 2 * nn) if m == 2 else tau * (tau + nn))
        D0 = D0 + b * (1 - 1 / P(tB, 1)) + c * (1 / P(tC, 1) - 4)
        amin = arb(min(a.lower(), 1))
        am1 = (a - 1).abs_upper()
        D1 = I(-1, 1) * 2 * am1 * e / (amin ** 3 * aJ ** 2) + b * (2 / P(tB, 2) - 2) + c * (16 - 2 / P(tC, 2))
        D2 = I(-1, 1) * 6 * am1 * e ** 2 / (amin ** 4 * aJ ** 3) + b * (6 - 6 / P(tB, 3)) + c * (6 / P(tC, 3) - 96) + (2 - tQ)
        bn = arb(1) / 3 - D1 - D2 / 6
        cn = arb(1) / 24 - D1 / 16 - D2 / 48
        an = 1 / (arb(5) / 6 + arb(3) / 4 * D1 + D2 / 12 + D0)
        S.last = (D0, D1, D2)
        return (an, bn, cn)
    Sb = S(box)
    if not all(Sb[k].lower() > box[k].lower() and Sb[k].upper() < box[k].upper() for k in range(3)):
        print("diag eps", e, "aJ", aJ, "betaJ", betaJ, "bJ", bJ, "tQ", tQ, "X0", X0, "Y", Y, "W", W)
        print("  D0,D1,D2:", S.last, " S(box):", Sb)
    for k in range(3):
        chk(Sb[k].lower() > box[k].lower() and Sb[k].upper() < box[k].upper(), "Brouwer box invariance comp %d (eps=%s)" % (k, eh))
    encl = Sb                       # the fixed point x_n lies in S(box); also x_n = S(x_n) lies in S(encl) for any encl containing it
    for it in range(6):
        nxt = S(encl)
        encl = tuple(nxt[k].intersection(encl[k]) for k in range(3))
    Sb = encl
    xa, xb, xc = Sb
    # ---- ranges (free-positive-range) ----
    Bn = betaJ + xb * e                                      # B/n
    Dn_ = betaJ / 2                                          # D/n
    chk(xa * aJ >= 8 * e * (betaJ + xb * e), "A >= 8B")
    chk(Bn >= 32 * omega ** 2 and Dn_ >= 32 * omega ** 2, "B >= 32d, D >= 32d")
    chk(xc > 4 * up(Nstar.log()) * omega ** 2, "C - D >= 4d  (C-D = xc n/L, d L/n <= log(N*) d0/N*)")
    chk(xc * e <= betaJ / 8, "C - D <= D/4")
    # ---- root localization (Lemma free-roots), sharpened Gershgorin constants ----
    sdB = omega / Bn.sqrt()                                  # sqrt(d/B)
    dB = omega ** 2 / Bn; dD = omega ** 2 / Dn_
    dH = 2 * omega ** 2 * e / (xa * aJ)                      # d/H <= 2d/A for p1 (H = A-d-1 >= A/2)
    VH = 2 * Bn * e / (xa * aJ)                              # B/H <= 2B/A
    cdiag1 = 2 * (1 + 2 * dH) ** 2 + VH                      # |diag - B| <= cdiag1 * d
    CJ1 = 2 * (1 + 2 * dH) * (1 + dB).sqrt() + cdiag1 * sdB  # roots(p1) in B +- CJ1 sqrt(Bd)
    cdiag2 = 3 + 2 * dD
    xcs = (xc * e).nonnegative_part().sqrt()                 # sqrt((C-D)/n)
    sGd_D = xcs * omega / Dn_                                # sqrt(G d)/D
    eta = cdiag2 * dD + 2 * (1 + 2 * dD) * (1 + dD) * sGd_D  # roots(p2) in [1-eta, 1+eta]
    chk(eta < arb('0.5'), "eta < 1/2")
    # |y - B| <= CJ1 sqrt(Bd) + eta (B + CJ1 sqrt(Bd));  eta sqrt(B/d) computed without dividing by omega:
    eta_sBd = cdiag2 * omega * Bn.sqrt() / Dn_ + 2 * (1 + 2 * dD) * (1 + dD) * xcs * Bn.sqrt() / Dn_
    Cloc = CJ1 + eta * CJ1 + eta_sBd
    chk(Cloc * sdB <= arb('0.75'), "roots in [B/4, 2B]")
    # ---- derivative ratios (Lemma free-ratios) ----
    yA = 2 * Bn * e / (xa * aJ)                              # y/A <= 2B/A
    dA = omega ** 2 * e / (xa * aJ)                          # d/A
    am = (Cloc + sdB * (1 + yA)) / (K6 * (1 - yA))
    bm = (1 + dA) * (1 + Cloc * sdB) / ((1 - yA) * K6 ** 2)
    rD = K6 * Bn.sqrt() * omega / (betaJ / 2)                # r/D
    rA = K6 * Bn.sqrt() * omega * e / (xa * aJ)              # r/A
    CV = 1 + rA + 3 * dA + 2 * sdB / K6 + 10 * dA * sdB / K6 + dB / (4 * K6 ** 2) + dA * dB / (4 * K6 ** 2)
    epsp = xc * e / (betaJ / 2 + xc * e)
    th = (Bn / (betaJ / 2)) * (1 + Cloc * sdB) * epsp * CV / ((1 - yA) * (1 - rD))
    tot = am + bm + th
    chk(tot < 1, "a+b+theta < 1")
    # ---- residual (Lemma free-defect) ----
    R0n = omega ** 2 + LAM * K6 * Bn.sqrt() * omega             # (d + lambda r)/n
    chk(R0n <= arb(1) / 20, "R0 <= n/20")
    rho1 = 2 * R0n + 4 * w1                                     # |u0(n+z) - u0(n)| <= rho1 on n + Omega_F
    # signed bound for n^4 L E_F^{(5)}(z), z in Omega_F (complex balls)
    zb = acb(arb(0, up(R0n + nu)), arb(0, up(R0n + nu)))      # (1/2 + z)/n and z/n lie in this ball
    one_c = acb(1)
    tc = acb(arb(2 * eh).union(arb(lo_(2 * el * (1 - 5 * nu)))), 0) + acb(arb(0, up(2 * rho1 * 4 / Lstar ** 2)), arb(0, up(2 * rho1 * 4 / Lstar ** 2)))
    sc = acb(arb(0, up(arb(9) / (4 * pi) * (-1 / eh + 2 * rho1).exp())), arb(0, up(arb(9) / (4 * pi) * (-1 / eh + 2 * rho1).exp())))
    Q5c = Qk(5, tc, sc, one_c)
    a0u = one_c + acb(arb(0, up(rho1 / (alo - rho1))), arb(0, up(rho1 / (alo - rho1))))      # alpha0/u0(n+z)
    nz = 1 / (one_c + zb)                                                                     # n/(n+z)
    term_h = Q5c * nz ** 4 * a0u
    P1 = 1 / (one_c + zb); P2 = 1 / (1 / P1 + acb(e * (bJ + xb)))
    term_B = 6 * acb(bJ + xb) * P1 * P2 * (P1 ** 3 + P1 ** 2 * P2 + P1 * P2 ** 2 + P2 ** 3)
    P3 = 1 / (acb(Dn_) + zb); P4 = 1 / (1 / P3 + acb(xc * e))
    term_C = -6 * acb(xc) * P3 * P4 * (P3 ** 3 + P3 ** 2 * P4 + P3 * P4 ** 2 + P4 ** 3)
    errs = up(arb(120) * etaG * arb(20) ** 5 * Lnu)                                          # G^{(5)}
    errs += up((bJ + xb) * 240 * nu / (1 - R0n) ** 6 + xc * 240 * nu / (Dn_ - R0n) ** 6)      # psi^{(5)} = 24 w^-5 + O(240 Re(w)^-6)
    errs += up(6 * e ** 3 / (xa * aJ) ** 4 * 2)                                               # psi^{(4)}(A+z)
    tot5 = term_h + term_B + term_C
    ce5 = arb(abs(tot5).upper()) + errs
    Cfin = ce5 / 120 * (LAM * K6) ** 5 * Bn ** arb('2.5') * (1 + sdB / (LAM * K6)) ** 5
    for nm, v, md in (('X0', X0, 'max'), ('tQ', tQ, 'min'), ('n2LDn', n2LDn, 'min'), ('aJ', aJ, 'max'), ('xc', xc, 'min'), ('tot', tot, 'max'),
                      ('Cloc', Cloc, 'max'), ('ce5', ce5, 'max'), ('Cfin', Cfin, 'max'), ('CJ1', CJ1, 'max'), ('CV', CV, 'max'),
                      ('eta', eta, 'max'), ('am', am, 'max'), ('bm', bm, 'max'), ('th', th, 'max'), ('YB', YB, 'min'), ('lam0', lam0, 'min')):
        upd(nm, v, md)
    if VERB and i == NEPS - 1:
        print("last eps-cell:", "X0", X0, "Y", Y, "W", W, "n2LDn", n2LDn, "aJ", aJ, "YB", YB, "betaJ", betaJ, "bJ", bJ, "tQ", tQ)
        print("  solution box:", xa, xb, xc, " Cloc", Cloc, " a+b+th", tot, " ce5", ce5, " Cfin", Cfin)

CfinMax = arb(worst['Cfin'])
thetaE = (1 + LAM ** 4 * (LAM - 1)).log() - arb('0.001')   # sup|E_F| <= thetaE  =>  sup|c-1| < lambda^4 (lambda-1)
lamL = up(((pi * Lstar / 2).log() + 2 * nu) / Lstar)
Kres = up((CfinMax / thetaE) ** 2 * (1 + lamL) ** 2)
KN = up(Nstar ** 3 * (Nstar + 2).log() ** 2 / d0 ** 5)          # with D_all = d0
d1 = arb(9) * arb(10) ** 24                                     # any all-n range D_all >= 9e24 (T^2, T = 3.000175e12)
KN1 = up(Nstar ** 3 * (Nstar + 2).log() ** 2 / d1 ** 5)         # omega^2 <= d1/N* <= d0/N* (used above) when K >= KN1
KSELF = up(Nstar ** 3 * (Nstar + 2).log() ** 2)                 # d >= 1 needs no all-n range; then d/n <= 1/N*
K = max(Kres.upper(), KN1.upper())
chk(arb(7000) >= K, "K <= 7000")
print("Lstar", Lstar.mid(), "K6", K6.mid(), "NEPS", NEPS)
print("N* =", Nstar.mid(), " nu =", nu.mid())
print("deltaA <=", deltaA.upper(), " etaG <=", etaG.upper(), " |Q5| <=", q5.upper())
for k_ in worst: print("  worst", k_, worst[k_])
print("Kres =", Kres.upper(), " KN(d0) =", KN.upper(), " KN(9e24) =", KN1.upper())
print("K = max(Kres, KN(9e24)) =", K, " <= 7000: PASS")
print("KSELF = max(Kres, N*^3 log^2(N*+2)) =", max(Kres.upper(), KSELF.upper()))
W = {k_: arb(v) for k_, v in worst.items()}
DOC.append("cells: tQ >= %s, n2LDn >= %s, aJ <= %s, xc >= %s, CJ1 <= %s, eta <= %s, Cloc <= %s, CV <= %s" % (
    flo(W['tQ'], 5), flo(W['n2LDn'], 5), fup(W['aJ'], 5), flo(W['xc'], 4), fup(W['CJ1'], 9), fup(W['eta'], 3), fup(W['Cloc'], 5), fup(W['CV'], 9)))
DOC.append("cells: a_m <= %s, b_m <= %s, theta <= %s, a+b+theta <= %s, ce5 <= %s, Cfin <= %s" % (
    fup(W['am'], 5), fup(W['bm'], 5), fup(W['th'], 5), fup(W['tot'], 5), fup(W['ce5'], 5), fup(W['Cfin'], 5)))
DOC.append("final: thetaE >= %s (log(1+lam^4(lam-1)) - 0.001), 1+lamL <= %s, Kres <= %s, KN(d0) <= %s, KN(9e24) <= %s, K <= %s, KSELF <= %s, N* >= %s" % (
    flo(thetaE, 6), fup(1 + lamL, 6), fup(Kres, 6), fup(KN, 4), fup(KN1, 5), fup(arb(K), 6), fup(max(Kres, KSELF), 4), flo(Nstar, 5)))
print("DOCUMENTED CONSTANTS (outward-rounded; DERIVATION.md/REPORT.md quote these):")
for line in DOC: print("  " + line)
