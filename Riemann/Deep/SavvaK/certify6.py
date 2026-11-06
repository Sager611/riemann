#!/usr/bin/env python3
# Explicit constant in Savva's two-thirds hyperbolicity wedge (Zenodo 22293642, J. Savva, 3 Sep 2026): certificate.
# Run: uv run --no-project --with python-flint --with sympy python3 certify6.py [Lstar] [Kr] [NEPS] [lambda] [-v]
# Architecture: Savva's six-match comparison (four Delta^2-log quotients + two normalizations), finite-free Jacobi
# factorization, MSS product localization, Savva's global-maximum critical radius, six-node Hermite-Genocchi residual,
# multiplier stability with Cauchy radius lambda*r.  Analytic input: OUR explicit complex saddle bound (HollandK Part A),
# not Savva's sector theorem.  All inequalities in Arb ball arithmetic (python-flint); chk() aborts on failure.
import sys
from math import comb
from flint import arb, acb, ctx
import sympy as sp
ctx.prec = 320

def I(lo, hi): return arb(lo).union(arb(hi))
def up(x): return arb(x.upper())
def lo_(x): return arb(x.lower())
def chk(c, msg):
    if not c:
        print("FAILED:", msg); sys.exit(1)
from fractions import Fraction
from decimal import Decimal, Context, ROUND_CEILING, ROUND_FLOOR
def _exact(a):
    m, e = a.man_exp(); return Fraction(int(m)) * Fraction(2) ** int(e)
def _huge(a):
    m, e = a.man_exp(); return abs(int(e)) > 4000
def _pow10(x, upper):                # outward power-of-ten bound for extreme magnitudes (positive values)
    import math
    v = arb(x.upper()) if upper else arb(x.lower())
    if not (v > 0): return str(v)
    lg = v.log() / arb(10).log()
    k = math.floor(float(lg.upper())) + 1 if upper else math.floor(float(lg.lower())) - 1
    return "1e%d" % k
def fup(x, sig=5):
    if _huge(x.upper()): return _pow10(x, True)
    v = _exact(x.upper()); c = Context(prec=sig, rounding=ROUND_CEILING)
    d = c.divide(Decimal(v.numerator), Decimal(v.denominator)); assert Fraction(d) >= v; return str(d)
def flo(x, sig=5):
    if _huge(x.lower()): return _pow10(x, False)
    v = _exact(x.lower()); c = Context(prec=sig, rounding=ROUND_FLOOR)
    d = c.divide(Decimal(v.numerator), Decimal(v.denominator)); assert Fraction(d) <= v; return str(d)
DOC = []
pi = arb.pi()
args = [a for a in sys.argv[1:] if not a.startswith('-')]
Lstar = arb(args[0]) if len(args) > 0 else arb(80)
KR = arb(args[1]) if len(args) > 1 else arb('2.4')        # critical radius r = KR sqrt(B d)
NEPS = int(args[2]) if len(args) > 2 else 400
LAM = arb(args[3]) if len(args) > 3 else arb('1.33')       # Cauchy radius lambda*r in the stability lemma
VERB = '-v' in sys.argv
d0 = arb(16140857488599731839205676)   # our all-n degree record (TREE N6)
D9 = arb(9) * arb(10) ** 24            # Chasse + RH to T = 3.000175e12: all-n range D_all >= 9.0e24
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
t, s = sp.symbols('t s', positive=True)
# u = u0(z), t = 1/u, s = 9/(4 pi e^{2u});  d/du t = -t^2, d/du s = -2s;  Z(u) = u (pi e^{2u} - 9/4), d/dz = (1/Z'(u)) d/du
Zp = sp.Rational(9, 4) * (t + 2 - t * s) / (t * s)          # Z'(u)
Zs = sp.Rational(9, 4) * (1 - s) / (t * s)                  # Z(u)
def Dz(Fx): return sp.cancel((sp.diff(Fx, t) * (-t ** 2) + sp.diff(Fx, s) * (-2 * s)) / Zp)
Psi = {2: sp.cancel(2 * t / Zp)}                            # Psi'' = 2/(u Z'(u))
for k in range(2, 6):
    Psi[k + 1] = Dz(Psi[k])
QP = {}
for k in range(2, 7):
    q = sp.cancel(Psi[k] * Zs ** (k - 1) * 2 / t)           # Q_k = Psi^(k) z^(k-1) 2u
    num, den = sp.fraction(q)
    QP[k] = (sp.Poly(num, t, s).terms(), sp.Poly(den, t, s).terms())
    dq = sp.cancel(sp.diff(q, t))
    dnum, dden = sp.fraction(dq)
    QP[('d', k)] = (sp.Poly(dnum, t, s).terms(), sp.Poly(dden, t, s).terms())
def evalpoly(terms, tv, sv, one):
    acc = one * 0
    for (i, j), c in terms:
        c = sp.Rational(c)
        acc += (one * int(c.p) / int(c.q)) * tv ** i * sv ** j
    return acc
def Qk(k, tv, sv, one):
    nt, dt = QP[k]
    return evalpoly(nt, tv, sv, one) / evalpoly(dt, tv, sv, one)
def Qk_mv(k, tv, sv, one):          # mean-value form in t:  Q(t) in Q(t_m) + Q'(tR)(tR - t_m)
    tm = arb(tv.mid())
    nt, dt = QP[('d', k)]
    dQ = evalpoly(nt, tv, sv, one) / evalpoly(dt, tv, sv, one)
    return Qk(k, tm, sv, one) + dQ * (tv - tm)

# ---------------- Section C: per-eps-cell pipeline (six-match branch, localization, radius, residual) ----------------
def fact(k): return 1 if k <= 1 else k * fact(k - 1)
gk = {k: up(arb(fact(k)) * etaG * (10 / (1 - 50 * nu)) ** k * Lnu) for k in (2, 3, 4, 5)}   # Cauchy, |w-n| <= 5
SELF = '-self' in sys.argv
if SELF:   # self-contained mode: K >= N*^2 log(N*+2) and n >= N*  =>  (d/n)^3 <= log(N*+2)/(K N*) <= 1/N*^3
    omega2 = up(1 / Nstar)
else:      # all-n-range mode: K >= 1 and n >= N*  =>  (d/n)^3 <= log(N*+2)/N*
    omega2 = up(((Nstar + 2).log() / Nstar) ** (arb(1) / 3))
omega = up(omega2.sqrt())                                  # sqrt(d/n) <= omega
MSER = 24
def Rk(k, T, U):
    """R_k(T,U) = int_0^1 (T+sU)^(-k-1) ds = [T^-k - (T+U)^-k]/(kU), with d/dT, d/dU; series in U/T + rigorous tail."""
    Tmin = lo_(T)
    chk(Tmin > 0, "Rk: T > 0")
    qa = up(abs(U) / Tmin)
    rho = up(arb(k + MSER + 2) / (MSER + 2) * qa * arb('1.1'))
    chk(rho < arb('0.5'), "Rk: series ratio < 1/2")
    Ti = 1 / T
    R = arb(0); RT = arb(0); RU = arb(0)
    Up = [arb(1)]
    for i in range(1, MSER + 1): Up.append(Up[-1] * U)
    for i in range(MSER):
        c = arb((-1) ** i * comb(k + i, i)) / (i + 1)
        tp = Ti ** (k + 1 + i)
        R += c * tp * Up[i]
        RT += c * (-(k + 1 + i)) * tp * Ti * Up[i]
        if i >= 1: RU += c * i * tp * Up[i - 1]
    M = MSER
    tauR = arb(comb(k + M, M)) / (M + 1) * Tmin ** (-(k + 1)) * qa ** M
    tauT = tauR * (k + 1 + M) / Tmin
    tauU = arb(comb(k + M, M)) * M / (M + 1) * Tmin ** (-(k + 2)) * qa ** (M - 1)
    return R + I(-2 * tauR, 2 * tauR), RT + I(-2 * tauT, 2 * tauT), RU + I(-2 * tauU, 2 * tauU)

# Matching system (k = 1..4):  a^k eps^(k-1)/k + w R_k(t, w eps) + delta R_k(1+nu/2, delta eps) = Y_(k-1)  (+- 6k L/n)
# unknowns y = (a, t, w, delta):  A = n/(a eps), C = t n, B = (t + w eps) n, D = n + 1/2 + delta eps n.
def Gfun(y, e, Yd, nuI):
    a, t, w, de = y
    T1 = 1 + nuI / 2
    out = []
    for k in range(1, 5):
        Rb = Rk(k, t, w * e)[0]
        Rd = Rk(k, T1, de * e)[0]
        ak = a ** k * (e ** (k - 1) if k > 1 else arb(1)) / k
        out.append(ak + w * Rb + de * Rd - Yd[k - 1] + I(-6 * k * Lnu, 6 * k * Lnu))
    return out
def DGfun(Y, e, nuI):
    a, t, w, de = Y
    T1 = 1 + nuI / 2
    rows = []
    for k in range(1, 5):
        Rb, RbT, RbU = Rk(k, t, w * e)
        Rd, RdT, RdU = Rk(k, T1, de * e)
        da = (a * e) ** (k - 1) if k > 1 else arb(1)
        rows.append([da, w * RbT, Rb + w * e * RbU, Rd + de * e * RdU])
    return rows
def mat_inv(Mx):
    n = len(Mx)
    A = [[arb(Mx[i][j]) for j in range(n)] + [arb(1 if i == j else 0) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(float(A[r][c].mid())))
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c:
                f = A[r][c]
                A[r] = [A[r][j] - f * A[c][j] for j in range(2 * n)]
    return [[arb(A[i][n + j].mid()) for j in range(n)] for i in range(n)]
def newton(y, em, Ym):
    P = None
    for it in range(60):
        g = [arb(x.mid()) for x in Gfun(y, em, Ym, arb(0))]
        J = [[arb(x.mid()) for x in row] for row in DGfun(y, em, arb(0))]
        P = mat_inv(J)
        dy = [sum((P[i][j] * g[j] for j in range(4)), arb(0)) for i in range(4)]
        y = [arb((y[i] - dy[i]).mid()) for i in range(4)]
        if max(abs(float(dy[i].mid())) for i in range(4)) < 1e-70: break
    return y, P
def krawczyk(ym, P, e, Yd, nuI, r):
    Ybox = [ym[i] + I(-r[i], r[i]) for i in range(4)]
    Gm = Gfun(ym, e, Yd, nuI)
    Jb = DGfun(Ybox, e, nuI)
    Kx = []
    for i in range(4):
        v = ym[i] - sum((P[i][j] * Gm[j] for j in range(4)), arb(0))
        for j in range(4):
            cij = (arb(1) if i == j else arb(0)) - sum((P[i][l] * Jb[l][j] for l in range(4)), arb(0))
            v += cij * I(-r[j], r[j])
        Kx.append(v)
    inside = all(Kx[i].lower() > Ybox[i].lower() and Kx[i].upper() < Ybox[i].upper() for i in range(4))
    return inside, Kx, Ybox
def CJ(uinv, Vp, nuI):
    """Gershgorin for the transported Jacobi matrix of q_{U,V}(Y) = 2F1(-d,U;V;Y/U), U >= V + d:
       every root Y has |Y - V| <= CJ sqrt(V d).  uinv = n/U, Vp = V/n, d/n in [0, omega2]."""
    dd = I(0, omega2)
    s1 = 1 - (dd + nuI) * uinv                       # S/U, S = U - d - 1
    chk(s1 > 0, "S > 0")
    diag = omega * (Vp * s1 * uinv + 2 * (s1 + (dd + nuI) * uinv) * abs(1 - 2 * Vp * uinv)) / (s1 ** 2 * Vp.sqrt())
    offsq = (1 + dd / Vp) * (1 - Vp * uinv) * (s1 + dd * uinv) / (s1 ** 2 * (s1 ** 2 - nuI ** 2 * uinv ** 2))
    chk(1 - Vp * uinv > 0, "U > V")
    return up(diag + 2 * up(offsq).sqrt())

import time
T0 = time.time()
if VERB: print('setup done', flush=True)
worst = {}
def upd(name, val, mode):
    v = val.upper() if mode == 'max' else val.lower()
    if name not in worst: worst[name] = v
    else: worst[name] = max(worst[name], v) if mode == 'max' else min(worst[name], v)
yprev = [arb(1) / 3, arb(2), arb(16) / 3, arb(1) / 3]       # Savva's limit y* (alpha = 3)
CELLROWS = ["cell,eps_lo,eps_hi,C6_up,Bn_up,Cfin_up,radius_budget_up,Cloc_up"]   # per-cell witnesses (outward-rounded)
nuI = I(0, nu)
samples = []
def cell(el, eh, i):
    global yprev
    e = I(el, eh)
    one = arb(1)
    tR = I(lo_(2 * el * (1 - 6 * nu)), up(2 * eh))
    sR = I(0, up(arb(9) / (4 * pi) * (-1 / eh).exp()))
    fac = {k: I(lo_((1 + 6 * nu) ** (-(k - 1)) * (1 - 6 * nu)), 1) for k in (2, 3, 4, 5)}
    X = [Qk_mv(j + 2, tR, sR, one) * fac[j + 2] + I(-gk[j + 2], gk[j + 2]) for j in range(4)]   # n^(j+1) L Delta^(j+2) h(n)
    Yd = [X[0], -X[1] / 2, X[2] / 6, -X[3] / 24]
    em = arb(e.mid()); Ym = [arb(x.mid()) for x in Yd]
    ym, P = newton(yprev, em, Ym)
    yprev = ym
    r = [arb('1e-12')] * 4
    ins = False
    for attempt in range(12):
        ins, Kx, Yb = krawczyk(ym, P, e, Yd, nuI, r)
        if ins: break
        r = [up(arb('1.5') * max(r[k], arb((Kx[k] - ym[k]).abs_upper()))) for k in range(4)]
        if max(float(x) for x in r) > 0.2: break
    if not ins: return False
    # Brouwer: T(y) = y - P G_exact(y) maps Y into K(Y) subset int(Y), because G_exact = G_approx + rho(y) with
    # |rho_k(y)| <= (k/eps) sum_U (n/U)^(k+1)/n <= 5k L/n <= 6k L/n for every y in Y (needs n/U <= 1 for all five U):
    chk(Yb[3] > 0 and Yb[1] > 1 and Yb[0] > 0 and Yb[0] * e < 1, "box: delta > 0, t > 1, 0 < a eps < 1 (so every n/U <= 1)")
    br = [Kx[k].intersection(Yb[k]) for k in range(4)]
    a, t, w, de = br
    chk(a > 0 and w > 0 and t > 1, "a = 1/alpha > 0, w > 0, t > 1")
    rA = a * e                                   # n/A
    b = t + w * e                                # B/n
    c = 1 + nuI / 2 + de * e                     # D/n
    chk(1 - rA * (b + omega2 + nu) > 0, "A - B - d > 0")
    chk(t - c - omega2 > 0, "C - D - d > 0")
    # ---- localization: Jacobi Gershgorin for both factors, MSS product in both orientations ----
    CJ1 = CJ(rA, b, nuI)                         # Q^(1)_{A,B}: roots in B +- CJ1 sqrt(Bd)
    CJ2 = CJ(1 / t, c, nuI)                      # q_{C,D}: roots in D +- CJ2 sqrt(Dd); Q^(D)_{C,D}: [1-eta, 1+eta]
    eta = up(CJ2 * omega / c.sqrt())
    chk(eta < arb('0.5'), "eta < 1/2 (positive lower endpoint, second factor)")
    chk(CJ1 * omega / b.sqrt() < arb('0.5'), "B - CJ1 sqrt(Bd) > 0 (positive lower endpoint, first factor)")
    Cloc = up(CJ1 * (1 + eta) + CJ2 * (b / c).sqrt())
    # ---- critical radius (Savva's four-term recurrence + global maximum), r = KR sqrt(Bd) ----
    sb = b.sqrt()
    yh = b + I(-1, 1) * Cloc * sb * omega         # y/n at every root / critical point of p_F
    chk(yh > 0, "y > 0")
    ky = c * yh * rA / t                          # D y/(A C)
    P3 = 1 - ky
    p2 = b + c + I(0, 2 * omega2 + nu) - c * yh / t - ky * (t + I(-omega2, 3 * omega2 + 3 * nu))   # P_2/n
    chk(p2 > 0, "P2 > 0")
    P1n = c * Cloc + (omega * (b + c) + omega ** 3 + c * yh * omega * (rA + 1 / t + omega2 * rA / t)
                      + (c * yh / t + ky * (t + 2 * omega2 + nu)) * omega) / sb                  # |P_1|/(n sqrt(Bd))
    P0n = c * yh * (1 + omega2 * rA) * (1 + omega2 / t)                                         # P_0/(n^2 d)
    totr = abs(P3) * KR * sb * omega / p2 + abs(P1n) / (KR * p2) + abs(P0n) / (KR ** 2 * b * p2)
    chk(totr < 1, "critical-radius budget < 1")
    # ---- six-node residual: n^5 L E_F^(6)(z), |z| <= d + lambda r ----
    R0n = up(omega2 + LAM * KR * sb * omega)
    chk(R0n <= arb(1) / 20, "R0 <= n/20")
    rho1 = up(2 * R0n + 4 * w1)
    zb = acb(arb(0, R0n), arb(0, R0n))
    one_c = acb(1)
    rad_t = up(2 * rho1 * 4 / Lstar ** 2)
    tc = acb(tR, 0) + acb(arb(0, rad_t), arb(0, rad_t))
    srad = up(arb(9) / (4 * pi) * (-1 / eh + 2 * rho1).exp())
    sc = acb(arb(0, srad), arb(0, srad))
    Q6c = Qk(6, tc, sc, one_c)
    a0u = one_c + acb(arb(0, up(rho1 / (alo - rho1))), arb(0, up(rho1 / (alo - rho1))))
    term_h = Q6c * (1 / (one_c + zb)) ** 5 * a0u
    bz = acb(b) + zb; tz = acb(t) + zb
    term_B = -24 * acb(w) * sum((bz ** (-(j + 1)) * tz ** (-(5 - j)) for j in range(5)), acb(0))
    cz = acb(c) + zb; hz = acb(1 + nuI / 2) + zb
    term_D = -24 * acb(de) * sum((cz ** (-(j + 1)) * hz ** (-(5 - j)) for j in range(5)), acb(0))
    term_A = -24 * acb(a ** 5 * e ** 4) * (one_c + acb(rA) * zb) ** (-5)
    errs = up(arb(720) * etaG * arb(20) ** 6 * Lnu)                                              # G^(6), Cauchy
    errs += up(240 * Lnu * (1 / (1 - R0n) ** 6 + 1 / (lo_(c) - R0n) ** 6 + 2 / (lo_(t) - R0n) ** 6 + 1))  # psi^(5) - 24 w^-5
    tot6 = term_h + term_B + term_D + term_A
    C6 = up(abs(tot6)) + errs
    Cfin = up(C6 * (LAM * KR * sb + omega) ** 6 / 720)
    CELLROWS.append("%d,%s,%s,%s,%s,%s,%s,%s" % (i, flo(arb(el), 8), fup(arb(eh), 8), fup(C6, 6), fup(b, 6), fup(Cfin, 6), fup(arb(totr.upper()), 6), fup(Cloc, 6)))
    for nm, v, md in (('C6', C6, 'max'), ('Cfin', Cfin, 'max'), ('totr', totr, 'max'), ('Cloc', Cloc, 'max'), ('CJ1', CJ1, 'max'),
                      ('CJ2', CJ2, 'max'), ('eta', eta, 'max'), ('p2', p2, 'min'), ('P1n', abs(P1n), 'max'), ('P0n', P0n, 'max'),
                      ('a', a, 'min'), ('t_lo', t, 'min'), ('t_hi', t, 'max'), ('w_lo', w, 'min'), ('w_hi', w, 'max'),
                      ('de_lo', de, 'min'), ('de_hi', de, 'max'), ('b_hi', b, 'max'), ('R0n', R0n, 'max'), ('tmc', t - c, 'min')):
        upd(nm, v, md)
    if i in (0, NEPS // 2, NEPS - 1):
        samples.append((i, el, eh, [x.mid() for x in Yd], [x for x in br], tot6, C6, totr, Cloc))
    return True
ncells = 0; nsplit = 0; maxdepth = 0
for i in range(NEPS):
    stack = [(epsmax * i / NEPS, epsmax * (i + 1) / NEPS, 0)]
    while stack:
        el, eh, dep = stack.pop(0)
        if cell(el, eh, i if dep == 0 else -1):
            ncells += 1; maxdepth = max(maxdepth, dep)
        else:
            chk(dep < 10, "Krawczyk: unique matching solution (cell %d, depth %d)" % (i, dep))
            mid = (el + eh) / 2; nsplit += 1
            stack = [(el, mid, dep + 1), (mid, eh, dep + 1)] + stack


thetaE = (1 + LAM ** 5 * (LAM - 1)).log() - arb('0.001')    # |E_F| <= thetaE  =>  sup|c_F - 1| < lambda^5 (lambda - 1)
lamL = up(((pi * Lstar / 2).log() + 2 * nu) / Lstar)          # log(n+2) <= (1 + lamL) L
CfinMax = arb(worst['Cfin'])
Kres = up(CfinMax / thetaE * (1 + lamL))
KN9 = up(Nstar ** 2 * (Nstar + 2).log() / D9 ** 3)
KNd0 = up(Nstar ** 2 * (Nstar + 2).log() / d0 ** 3)
KSELF = up(Nstar ** 2 * (Nstar + 2).log())
K = max(Kres.upper(), KN9.upper(), arb(1).upper())
print("Lstar", Lstar.mid(), "KR", KR.mid(), "NEPS", NEPS, "lambda", LAM.mid(), "prec", ctx.prec, "cells", ncells, "splits", nsplit, "maxdepth", maxdepth)
print("N* =", Nstar.mid(), " nu =", nu.mid(), " omega2 (d/n max) =", omega2.mid())
print("deltaA <=", deltaA.upper(), " etaG <=", etaG.upper())
for k_ in worst: print("  worst", k_, worst[k_])
for smp in samples:
    print("  sample cell", smp[0], "eps in [%s, %s]" % (smp[1].mid(), smp[2].mid()))
    print("     Y =", smp[3]); print("     (a,t,w,delta) =", smp[4]); print("     n^5 L E6 =", smp[5], " C6", smp[6], " totr", smp[7], " Cloc", smp[8])
print("thetaE =", thetaE.lower(), " 1+lamL =", (1 + lamL).upper())
print("Kres =", Kres.upper(), " KN(9e24) =", KN9.upper(), " KN(d0) =", KNd0.upper())
print("K = max(Kres, KN(9e24), 1) =", K)
if not SELF:
    chk(arb(24) >= K, "K <= 24")
    print("K <= 24: PASS  (n^2 log(n+2) >= 24 d^3 => J^{d,n} hyperbolic, every d >= 1, with the all-n range d <= 9.0e24)")
else:
    KS = max(Kres.upper(), KSELF.upper())
    chk(arb('2.76e25') >= KS, "KSELF <= 2.76e25")
    print("self-contained: KSELF <= 2.76e25: PASS (no all-n range used)")
print("KSELF = max(Kres, N*^2 log(N*+2)) =", max(Kres.upper(), KSELF.upper()))
W = {k_: arb(v) for k_, v in worst.items()}
DOC.append("cells: a >= %s, t in [%s, %s], w in [%s, %s], delta in [%s, %s], t - D/n >= %s" % (
    flo(W['a'], 5), flo(W['t_lo'], 6), fup(W['t_hi'], 6), flo(W['w_lo'], 5), fup(W['w_hi'], 5), flo(W['de_lo'], 5), fup(W['de_hi'], 5), flo(W['tmc'], 5)))
DOC.append("cells: CJ1 <= %s, CJ2 <= %s, eta <= %s, Cloc <= %s, P2/n >= %s, |P1|/(n sqrt(Bd)) <= %s, P0/(n^2 d) <= %s, radius budget <= %s" % (
    fup(W['CJ1'], 6), fup(W['CJ2'], 6), fup(W['eta'], 3), fup(W['Cloc'], 6), flo(W['p2'], 5), fup(W['P1n'], 5), fup(W['P0n'], 5), fup(W['totr'], 5)))
DOC.append("cells: C6 = sup n^5 L |E_F^(6)| <= %s, Cfin <= %s, R0/n <= %s" % (fup(W['C6'], 5), fup(W['Cfin'], 5), fup(W['R0n'], 3)))
DOC.append("final: thetaE >= %s, 1+lamL <= %s, Kres <= %s, KN(9e24) <= %s, KN(d0) <= %s, K <= %s, KSELF <= %s, N* >= %s" % (
    flo(thetaE, 6), fup(1 + lamL, 6), fup(Kres, 6), fup(KN9, 5), fup(KNd0, 5), fup(arb(K), 6), fup(max(Kres, KSELF), 4), flo(Nstar, 5)))
cellfile = 'cells_self.csv' if SELF else 'cells_main.csv'
with open(cellfile, 'w') as fh: fh.write("\n".join(CELLROWS) + "\n")
print("per-cell witnesses written to", cellfile, "(%d rows; cell = -1 marks bisected sub-cells)" % (len(CELLROWS) - 1))
print("DOCUMENTED CONSTANTS (outward-rounded):")
for line in DOC: print("  " + line)
