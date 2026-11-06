#!/usr/bin/env python3
# cert.py -- Ingham3: a prime between n^k and (n+1)^k for every n with n^k >= x0(k), k = 3..7, and for every n >= 1 when k = 7.
# One command (from the repo root):
#   uv run --no-project --with python-flint --with mpmath python3 Riemann/Deep/Ingham3/cert.py
# Every inequality is checked in Arb ball arithmetic (python-flint); decimal inputs are entered as exact decimals
# and every derived constant keeps its ball, so all roundings are safe. See DERIVATION.md for the mathematics.
import sys, time, json, math
from flint import arb, acb, ctx
ctx.prec = 128
T_START = time.time()
QUICK = 'quick' in sys.argv
LOG = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); LOG.append(s)
def up(x):   # an exact arb >= every point of the ball x
    return arb(x.upper())
def lo(x):
    return arb(x.lower())
def certain_lt(a, b):
    return bool(a < b)
PI = arb.pi()
ONE = arb(1)
# ------------------------------------------------------------------ inputs (DERIVATION.md sec. 1)
B2 = (arb(777) - arb(1081).sqrt())/852                 # Liu's contact B_2
TSTAR = B2 + arb('2.2e-254')                            # Bebel (Lean): zeta(s) != 0 for Re s > T*
B = arb('0.8733819666132076')                           # B >= T*  (checked below)
T0 = arb('3e12')                                        # Platt-Trudgian: RH up to height 3*10^12
def Nplus(T):   # Hasanalizade-Shen-Wong Cor. 1.2: N(T) <= T/(2pi) log(T/(2 pi e)) + 0.1038 log T + 0.2573 loglog T + 9.3675, T >= e
    L = T.log(); return T/(2*PI)*(T/(2*PI*arb(1).exp())).log() + arb('0.1038')*L + arb('0.2573')*L.log() + arb('9.3675')
# ------------------------------------------------------------------ Kadiri-Lumley-Ng (arXiv 2101.12263, JMAA 465 (2018)), Lemma 4.14
H0 = arb('3.0610046e10'); HH = H0 - 1
b1, b2, b3, b4 = arb('0.62'), arb('1.048'), arb('0.605'), arb('0.529')
a1, a2 = arb('0.63'), arb('2.851')     # |zeta(1/2+it)| <= a1 t^{1/6} log t (t >= 3) holds with 0.618 (Hiary-Patel-Yang 2024)
m0 = (1 + arb(2)/3*(arb(6)/5).sqrt()).sqrt()
EG = arb.const_euler()
LH0 = H0.log(); llH0 = LH0.log()
def zeta_r(s): return acb(s).zeta().real
def Ifun(Ae, n, al):     # I(A,n) = int_0^oo x^A e^{-2 al x^2} (log x)^n dx  (KLN (4.17) with beta = 2)
    s = (Ae + 1)/2
    g = s.gamma(); ps = s.digamma(); p1 = acb(s).polygamma(1).real
    l = (2*al).log()
    core = [g, g*ps - l*g, g*(ps*ps + p1) - 2*l*g*ps + l*l*g][n]
    return (2*al)**(-s)*core/arb(2)**(n + 1)
def C1k(k): return 6/PI**2 + b2/(k*H0).log()
def C2k(k):
    L = (k*H0).log(); return PI*m0*b1/L + 6*m0/(PI*k*H0) + PI*m0*b2/(k*H0*L)
a3 = -6*a1/arb(1).exp() + a2
def C4(k):
    C3 = a3**2*C1k(k)*(k*H0).log(); return C1k(k)*a1**2*(1 + 1/C3.sqrt())**2
def Jparts(k, al):
    r = C2k(k)/C1k(k)*k; b = arb(2); t3, t6 = arb(1)/3, arb(1)/6
    N0 = Ifun(b+t3,0,al) + r*Ifun(b-2*t3,0,al)
    N1 = 2*Ifun(b+t3,1,al) + 2*r*Ifun(b-2*t3,1,al)
    N2 = Ifun(b+t3,2,al) + r*Ifun(b-2*t3,2,al)
    N3 = 2*a2*(Ifun(b+t6,0,al) + r*Ifun(b-5*t6,0,al))/a1
    N4 = 2*a2*(Ifun(b+t6,1,al) + r*Ifun(b-5*t6,1,al))/a1
    N5 = a2**2*(Ifun(b,0,al) + r*Ifun(b-1,0,al))/a1**2
    return N0, N1, N2, N3, N4, N5
def Jfun(k, T, al):
    N0, N1, N2, N3, N4, N5 = Jparts(k, al); LT = T.log(); T6 = T**(arb(1)/6)
    return N0 + N1/LT + N2/LT**2 + N3/(T6*LT) + N4/(T6*LT**2) + N5/(T6**2*LT**2)
def om1(s, T, al): return (al*(s/T)**2).exp()
def Ufun(al, k, T): return 4*al*2*C4(k)*om1(arb(1)/2, T, al)**2*Jfun(k, T, al)
def C5(k, de):
    L = (k*H0).log(); return PI*m0*b4/(2*de)*(1 + 2*de/L)**2*(2*de*EG/L).exp()
def C6(k, de):
    L = (k*H0).log(); return b4/(5*de*de.exp())*(1 + de/L)**2 + b3*(-2*de).exp()/L**2
def Kfun(k, de, T, al): return (C5(k, de) + C6(k, de)*PI*m0/(k*T))*Ifun(arb(1),0,al) + C6(k, de)/k*Ifun(arb(2),0,al)
def Vfun(al, k, de, T): return 8*al*om1(1 + de/(k*T).log(), T, al)**2*Kfun(k, de, T, al)
def Vinf(al, k, de): return 8*al*(C5(k, de)*Ifun(arb(1),0,al) + C6(k, de)/k*Ifun(arb(2),0,al))
def Mkd(k, de): return (LH0/((k*H0).log() + 2*de)).max(ONE)
def b12(H): return 1/(2*(1 - 1/H)**2)
def logC1(al, d, de, k, H, s):
    LkH = (k*H0).log()
    e1 = arb(8)/3*de*(2*s - 1)*Mkd(k, de) + 4*de*(2*s - 1)*llH0/(LkH + 2*de)
    eU = 2*(1 - s) + 2*d/LH0 + 2*de*(2*s - 1)/(LkH + 2*de)
    e3 = 2*d*(2*llH0 - LkH.log())/LH0 + arb(8)*d/3 + 2*al
    return b12(H).log() + e1 + eU*Ufun(al, k, H0).log() + (2*s - 1)*Vfun(al, k, de, H0).log() + e3
def b5(eta): return zeta_r(1 + eta)**4/zeta_r(2 + 2*eta)**2 + 2*zeta_r(1 + eta)**2/zeta_r(2 + 2*eta)
def b6(X, eta):
    L = X.log(); return (1 + eta)*L/(eta*X**eta)*(1 + 1/(eta*L) + EG/L + 7*eta/(12*(1 + eta)*X*L))
def b7(k, eta):
    z = zeta_r(1 + eta); return (1 + 2/(arb('3.006')*z*(k*H0)**(1 + eta)))*(arb('3.006')*z)**2
def b8(eta, H): return ((2 + eta)**2/H**2 + ((1 + 2*eta)/H + 1)**2).sqrt()
def b9(eta, H, k):
    l2 = arb(2).log()
    return (PI*b7(k, eta).log()/l2 + PI*b5(eta).log()/l2 - 2*PI*(1 - b6(arb(10)**9, eta)**2).log()/l2 + PI
            + 2*(1 + 2*eta)/l2*(b8(eta, H)/(2*PI)).log())
def C7(eta, H, k): return (2*(1 + 2*eta) + 2*PI*(1 + eta))/arb(2).log() + b9(eta, H, k)/LH0
def b10(k, mu):
    b = b6(k*H0, mu - 1); return -(1 - b**2).log()/b**2
def b11(X, tau):
    L = X.log(); return 1 + 3/((tau - 1)*L) + 6/((tau - 1)**2*L**2) + 6/((tau - 1)**3*L**3)
def C8(k, mu):
    X = k*H0; L = X.log()
    return b10(k, mu)*L**2/X**(2*mu - 2)*(4*mu*b11(X, 2*mu)/(k*(2*mu - 1)) + 2*PI*m0*(2*mu - 1)*b11(X, 2*mu - 1)/(mu - 1))
ETA = arb('0.2561882')        # any eta in [eta0, 1/2) with mu <= 1 + eta is admissible; this one minimises C7
# KLN Table 1 rows: sigma0, k, mu, alpha, delta, d, A (table), B (table)
KLN_ROWS = [('0.60','0.5','1.251','0.288','0.3140','0.341','2.177','5.663'),('0.65','0.6','1.249','0.256','0.3070','0.340','2.963','5.249'),
('0.70','0.8','1.247','0.222','0.3040','0.339','3.983','4.824'),('0.75','1.0','1.245','0.189','0.3030','0.338','5.277','4.403'),
('0.80','1.0','1.245','0.160','0.3030','0.337','6.918','3.997'),('0.85','1.0','1.245','0.133','0.3030','0.336','8.975','3.588'),
('0.86','1.0','1.245','0.127','0.3030','0.335','9.441','3.514'),('0.87','1.0','1.245','0.122','0.3030','0.335','9.926','3.430'),
('0.88','1.0','1.245','0.116','0.3030','0.335','10.431','3.346'),('0.89','1.0','1.245','0.111','0.3030','0.335','10.955','3.262'),
('0.90','1.0','1.245','0.105','0.3030','0.334','11.499','3.186'),('0.91','1.0','1.245','0.100','0.3030','0.334','12.063','3.102'),
('0.92','1.0','1.245','0.095','0.3030','0.334','12.646','3.017'),('0.93','1.0','1.245','0.089','0.3030','0.333','13.250','2.941'),
('0.94','1.0','1.245','0.084','0.3030','0.333','13.872','2.856'),('0.95','1.0','1.245','0.079','0.3030','0.333','14.513','2.772'),
('0.96','1.0','1.245','0.074','0.3030','0.332','15.173','2.694'),('0.97','1.0','1.245','0.069','0.3030','0.332','15.850','2.609'),
('0.98','1.0','1.245','0.064','0.3030','0.331','16.544','2.532'),('0.99','1.0','1.245','0.060','0.3030','0.331','17.253','2.446')]
class KP:
    """One admissible KLN parameter set; for every sigma > 1/2 + d/log H0 and T >= H0:
       N(sigma,T) <= A(sigma) (log kT)^{2 sigma} (log T)^{5-4 sigma} T^{8(1-sigma)/3} + Bc(sigma) (log T)^2,
       log A(sigma) = c0 + c1 sigma (affine), Bc(sigma) = q0 - q1 sigma (affine, q1 >= 0)."""
    def __init__(self, row, de_override=None):
        s0, k, mu, al, de, d, At, Bt = [arb(v) for v in row]
        if de_override is not None: de = arb(de_override)
        self.name = row[0]; self.k, self.mu, self.al, self.de, self.d = k, mu, al, de, d
        l0, l1 = logC1(al, d, de, k, HH, arb(0)), logC1(al, d, de, k, HH, arb(1))
        self.c0 = l0 - (2*PI*d).log(); self.c1 = l1 - l0
        c7 = C7(ETA, HH, k); c8 = C8(k, mu)
        self.q0 = (c7*(mu + d/LH0) + c8)/(2*PI*d); self.q1 = c7/(2*PI*d)
        self.smin = arb(1)/2 + d/LH0
        # admissibility (KLN Lemma 4.14 hypotheses)
        N0, N1, N2, N3, N4, N5 = Jparts(k, al)
        self.ok = all([Ufun(al, k, H0) > 1, N1 > 0, N4 > 0, Vinf(al, k, de) > 1, de > 0, d > 0, al > 0,
                       k >= arb(10)**9/H0, k <= 1, b6(arb(10)**9, ETA) < 1, ETA < arb(1)/2, b6(arb(10)**9, mu - 1) < 1, mu <= 1 + ETA])
        self.Atab, self.Btab, self.s0 = At, Bt, s0
    def logA(self, s): return self.c0 + self.c1*s
    def Bc(self, s): return self.q0 - self.q1*s
def build_kln(de_override=None):
    P = [KP(r, de_override) for r in KLN_ROWS]
    if de_override is not None:
        for p in P:
            say('KLN row %s with delta=%s: admissible=%s  A=%s  B=%s' % (p.name, de_override, p.ok, p.logA(p.s0).exp().str(6, radius=False), p.Bc(p.s0).str(6, radius=False)))
            assert p.ok
        return P
    for p in P:
        Am, Bm = p.logA(p.s0).exp(), p.Bc(p.s0)
        good = p.ok and (abs(Am - p.Atab) < arb('0.002')) and (abs(Bm - p.Btab) < arb('0.002'))   # we use A_mine, B_mine (our own evaluation of KLN's formulas)
        say('KLN row %s: admissible=%s  A_mine=%s (table %s)  B_mine=%s (table %s)  reproduces=%s' % (
            p.name, p.ok, Am.str(6, radius=False), p.Atab.str(5, radius=False), Bm.str(6, radius=False), p.Btab.str(4, radius=False), good))
        assert good
    for p in P: assert p.c1 > 0          # A(sigma) increasing in sigma (used for Amax in the tail)
    return P
# ------------------------------------------------------------------ kernel (DERIVATION.md sec. 3)
S1 = arb('0.2173')     # >= sup_{y >= pi} |sin y|/y = 1/sqrt(1 + y*^2) = 0.2172336..., tan y* = y*, y* = 4.4934
def check_S1():
    # |sin y|/y <= S1 on [pi, 2 pi] (interval subdivision) and 1/y <= 1/(2 pi) < S1 beyond
    n = 8000
    for i in range(n):
        y = lo(PI + PI*i/n).union(up(PI + PI*(i + 1)/n))     # a ball containing the whole subinterval
        v = abs(y.sin())/y
        assert v < S1, (i, v)
    assert 1/(2*PI) < S1
    say('kernel: sup_{y>=pi} |sin y|/y <= 0.2173 checked on 8000 subintervals of [pi, 2pi]')
def S_up(y):
    """upper bound of S(y) = sup_{y'>=y} |sin y'|/y' (non-increasing in y)."""
    if y.upper() <= PI.lower():
        v = (y.sin()/y) if y > 0 else ONE
        return up(v).max(S1)
    if y.lower() >= PI.upper():
        return up(1/y).min(S1)
    return up(y.sin()/y).max(S1)
class Kernel:
    """Box spline Phi = 1_[0,1/m] * ... * 1_[0,1/m] (m factors). For Re z = delta*beta in [0, delta_max],
       |Phi~(z)|/int Phi <= G(Im z), G(l) = (S(l/2m) + 2 delta_max min(1/m, 1/l))^m, non-increasing.
       Mom(p; l0) = int_{l0}^oo max(1,l)^p d(-G(l)) = max(1,l0)^p G(l0) + p int_{max(1,l0)}^oo G(l) l^{p-1} dl."""
    def __init__(self, m, dmax, pgrid, r=arb(1) + arb(1)/500, Lam_fac=10000, lmax_store=30):
        self.m, self.dmax, self.pgrid = m, dmax, pgrid
        mm = arb(m); self.mm = mm
        logr = r.log(); self.logr = logr
        N = int(math.ceil(math.log(2*m*Lam_fac)/float(logr.mid()))) + 1
        self.N = N
        lam = [(logr*i).exp() for i in range(N + 1)]
        G = [self.G(l) for l in lam]
        self.G0 = (1 + 2*dmax/mm)**m
        self.lam, self.Gv = lam, G
        Lam = lam[N]
        istore = int(math.ceil(math.log(lmax_store)/float(logr.mid())))
        self.istore = istore
        self.T = {}
        for p in pgrid:
            # tail: G(l) <= ((2m + 2 dmax)/l)^m for l >= 2 m pi; int_Lam^oo l^{p-m-1} dl * p = p Lam^{p-m}/(m-p)
            assert Lam > 2*mm*PI and p < m
            acc = p*((2*mm + 2*dmax)**m)*Lam**(p - m)/(m - p)
            pw = [((logr*i)*p).exp() for i in range(N + 1)]
            store = [None]*(istore + 1)
            for i in range(N - 1, -1, -1):
                acc = acc + G[i]*(pw[i + 1] - pw[i])          # p * int_{lam_i}^{lam_i+1} G l^{p-1} dl <= G_i (lam_{i+1}^p - lam_i^p)
                if i <= istore: store[i] = up(acc)
            self.T[str(p)] = store
    def G(self, l):
        mm = self.mm
        corr = 2*self.dmax*(ONE/mm).min(ONE/l) if l > 0 else 2*self.dmax/mm
        return up((S_up(l/(2*mm)) + corr)**self.m)
    def mom(self, p, l0):
        """upper bound for Mom(p; l0); p must be on the grid."""
        mm = self.mm
        if l0.lower() >= 2*mm*PI and (l0.lower() >= self.lam[self.istore].upper()):
            # G(l) <= ((2m + 2 dmax)/l)^m for l >= 2 m pi:  Mom <= (2m+2dmax)^m l0^{p-m} (1 + p/(m-p))
            return up((2*mm + 2*self.dmax)**self.m*lo(l0)**(p - self.m)*mm/(mm - p))
        T = self.T[str(p)]
        if l0.upper() <= 1:
            return self.G0 + T[0]
        i = int(math.floor(float((l0.log()/self.logr).lower())))
        assert 0 <= i <= self.istore
        return up(l0**p*self.Gv[i] + T[i])
def p_up(p, pgrid):
    for q in pgrid:
        if q >= p: return q
    raise ValueError('p beyond grid: %s' % p)

# ------------------------------------------------------------------ the bound for one cell of x (DERIVATION.md sec. 4)
LN2 = arb(2).log()
class Setup:
    def __init__(self, P, kernels, pgrid, sig0, dsig):
        self.P, self.kernels, self.pgrid = P, kernels, pgrid
        self.dsig = dsig
        self.sig0 = sig0
        # sigma cells [s_j, s_{j+1}] covering [sig0, B]; each with the KLN parameter set minimising A at both ends
        cells = []
        s = sig0
        while s < B:
            t = (s + dsig).min(B) if (s + dsig) < B else B
            best = None
            for p in P:
                if not (s > p.smin): continue
                v = p.logA(s) + p.logA(t)
                if best is None or v.mid() < best[0].mid(): best = (v, p)
            cells.append((s, t, best[1]))
            if t is B: break
            s = t
        self.cells = cells
        self.Amax = max((c[2].logA(c[1]).exp().upper() for c in cells), key=lambda z: float(arb(z).mid()))
        self.Bmax = max((c[2].Bc(c[0]).upper() for c in cells), key=lambda z: float(arb(z).mid()))
        self.Amax, self.Bmax = arb(self.Amax), arb(self.Bmax)
NT0 = up(Nplus(T0))
def cell_bound(S, k, La, Lb, m, detail=False):
    """Upper bound R for 1 - (main term)^{-1} * (prime sum) over all x = n^k in [e^La, e^Lb] (DERIVATION sec. 4)."""
    kk = arb(k); mm = arb(m)
    La, Lb = arb(La), arb(Lb)
    na, nb = (La/kk).exp(), (Lb/kk).exp()
    da = up(kk*(1/na).log1p()); db = lo(kk*(1/nb).log1p())
    K = None
    for KK in S.kernels[m]:                      # kernels sorted by dmax; take the first with da <= dmax
        if da < KK.dmax: K = KK; break
    assert K is not None
    Db = -db.log(); lDb = Db.log(); lam0 = db*T0
    assert Db > 1
    # every zero (on or off the line), weight x^{-1/2}, and the off-line excess for beta in (1/2, sig0]  (DERIVATION sec. 4.2)
    assert (1/db) >= 14                         # t1 = 1/delta >= 14 > every height below the first zero
    pS = p_up(1 + 1/Db, S.pgrid)
    Sall = Db/db*K.mom(pS, arb(0))
    tot = ((S.sig0 - 1)*La).exp()*Sall
    Zlow = arb(0)
    # sigma in [sig0, B]: KLN.  Case (i) lam0 >= 1: expand log t around log T0; case (ii) lam0 < 1: around D = log(1/delta)
    ell0 = T0.log()
    if lam0.lower() >= 1:
        Lam = ell0.log() - lam0.log()/ell0
        pden = ell0
        momB = (ell0**2)*lam0**(-2/ell0)*K.mom(p_up(2/ell0, S.pgrid), lam0)
    else:
        assert lam0.upper() < 1
        Lam = lDb
        pden = Db
        momB = Db**2*K.mom(p_up(2/Db, S.pgrid), lam0)
    worst = None
    for (s, t, P) in S.cells:
        p = p_up(arb(8)/3*(1 - s) + (5 - 2*s)/pden, S.pgrid)
        mom = K.mom(p, lam0)
        E0 = -La + P.c0 + arb(8)/3*Db + 5*Lam
        E1 = La + P.c1 - arb(8)/3*Db - 2*Lam
        if abs(E1) > arb('1e-6'):
            integ = ((E0 + E1*t).exp() - (E0 + E1*s).exp())/E1
        else:
            integ = (t - s)*(E0 + E1*s).exp().max((E0 + E1*t).exp())*arb('1.000001')
        main = La*mom*integ
        Bt = P.Bc(s)*momB*(((t - 1)*La).exp() - ((s - 1)*La).exp())
        tot = tot + main + Bt
        if detail:
            c = main + Bt
            if worst is None or c.mid() > worst[0].mid(): worst = (c, s)
    Zhigh = 2*tot
    # prime powers and trivial zeros
    J = int(math.floor(float(((Lb + da)/LN2).upper())))
    PP = mm*(Lb + da)*((da/2).exp()*(-La/2).exp()/4 + J*da.exp()*(-2*La/3).exp()/9 + arb(J).log()/(db*La.exp()))
    Et = 2*(-3*La).exp()
    R = up(Zlow + Zhigh + PP + Et)
    if detail:
        return R, dict(Zlow=Zlow, Zhigh=Zhigh, PP=PP, lam0=lam0, worst_sigma=worst[1], worst=worst[0], Sall_part=2*((S.sig0 - 1)*La).exp()*Sall)
    return R
def tail_bound(S, k, Ls, m):
    """R_an(L) at L = Ls: an upper bound for R(x), valid for every x with log x = L >= Ls, and non-increasing in L (DERIVATION sec. 5)."""
    K = S.kernels[m][0]; kk = arb(k); mm = arb(m); L = arb(Ls)
    ns = (L/kk).exp()
    cst = 1/(1 - 1/(2*ns))                      # 1/delta <= cst * n/k for n >= ns
    assert up(kk*T0*cst/ns) < 1, 'need lambda0 <= 1 on the tail'
    Dlow = L/kk - kk.log()                      # D = log(1/delta) >= L/k - log k
    assert Dlow > 1
    pmax = p_up(arb(8)/3*(1 - S.sig0) + 4/Dlow, S.pgrid)
    Mst = K.mom(pmax, arb('0.5'))               # Mom(p; l0) for l0 <= 1, any p <= pmax
    p1 = p_up(1 + 1/Dlow, S.pgrid); M1 = K.mom(p1, arb(0))
    one_m = 1 - B
    kap = 1 - arb(8)/(3*kk)
    Xi = L*kap + arb(8)/3*(kk/cst).log() - 2*(L/kk).log()
    assert Xi > 0
    # monotonicity conditions for L >= Ls (DERIVATION sec. 5)
    mono = [one_m*(kap - 2/L) > 4/L, (1 - S.sig0 - 1/kk) > 2/L, (1 - 1/kk) > 3/L]
    assert all(mono), mono
    dmax = K.dmax
    main = L*S.Amax*Mst*(L/kk)**3*(-one_m*Xi).exp()/Xi
    Bt = S.Bmax*Mst*(L/kk)**2*(-one_m*L).exp()
    Sallp = (-(1 - S.sig0)*L).exp()*(cst/kk)*(L/kk).exp()*(L/kk)*M1
    Zhigh = 2*(main + Bt + Sallp)
    Zlow = arb(0)
    assert up(kk/ns) < dmax
    J = (L + 1)/LN2
    PP = mm*(L + 1)*((-L/2).exp()*arb(1).exp()/4 + J*arb(1).exp()*(-2*L/3).exp()/9 + J.log()*(cst/kk)*(L/kk - L).exp())
    Et = 2*(-3*L).exp()
    return up(Zlow + Zhigh + PP + Et)

# ------------------------------------------------------------------ k = 7, n <= 100: a direct check (DERIVATION sec. 6)
def direct_small_k7(nmax):
    """for 1 <= n < nmax: an explicit prime in (n^7, (n+1)^7) (deterministic Miller-Rabin; n^7 < 2^64)."""
    import sympy
    out = []
    for n in range(1, nmax):
        assert (n + 1)**7 < 2**64
        p = sympy.nextprime(n**7)
        assert p < (n + 1)**7 and sympy.isprime(p)
        out.append((n, int(p)))
    return out
# ------------------------------------------------------------------ driver
def scan(S, k, m, Lstart, Lend, wfun):
    """certify every cell [L, L + w(L)] from Lstart to Lend; returns (ok, cells, worst R, its L)."""
    L = Lstart; worst = (-1.0, None); n = 0
    while L < Lend:
        Lb = min(L + wfun(L), Lend)
        R = cell_bound(S, k, L, Lb, m)
        if not (R < 1): return False, n, (float(R.mid()), L)
        if float(R.mid()) > worst[0]: worst = (float(R.mid()), L)
        L = Lb; n += 1
    return True, n, worst
def find_L0(S, k, m, lo_L, hi_L, w):
    """bisection for the smallest L (to 0.01) with cell_bound(L, L + w) < 1 (lo_L fails, hi_L passes)."""
    assert cell_bound(S, k, hi_L, hi_L + w, m) < 1 and not (cell_bound(S, k, lo_L, lo_L + w, m) < 1)
    a, b = lo_L, hi_L
    while b - a > 0.005:
        c = (a + b)/2
        if cell_bound(S, k, c, c + w, m) < 1: b = c
        else: a = c
    return math.ceil(b*100)/100
PLAN = {   # k: (bracket for L0, cell width, tail start L*)
    3: ((1800, 1950), 0.25, 2400),
    4: ((440, 500), 0.1, 640),
    5: ((270, 300), 0.05, 400),
    6: ((205, 230), 0.05, 300),
}
def run_k(S, k, m=2, bracket=None, Lstar=None):
    t = time.time()
    (lo_L, hi_L), w, Ls = PLAN[k]
    if bracket: lo_L, hi_L = bracket
    if Lstar: Ls = Lstar
    L0 = find_L0(S, k, m, lo_L, hi_L, w)
    ok, n, worst = scan(S, k, m, L0, Ls, lambda L: w)
    Rt = tail_bound(S, k, Ls, m)
    ok = ok and bool(Rt < 1)
    say('k=%d: x0 = e^%.2f (n0 = e^%.4f = %s); %d cells on [%.2f, %d] all R < 1 (max %.4f at L=%.2f); tail from L*=%d: R <= %.3g; %s  %.0fs' % (
        k, L0, L0/k, '%.4e' % math.exp(L0/k) if L0/k < 700 else '10^%.2f' % (L0/k/math.log(10)), n, L0, Ls, worst[0], worst[1] if worst[1] else -1, Ls, float(Rt.mid()), 'CERTIFIED' if ok else 'FAILED', time.time() - t))
    assert ok
    return dict(k=k, L0=L0, logn0=L0/k, cells=n, worstR=worst[0], worstL=worst[1], Lstar=Ls, tailR=float(Rt.mid()))
def run_k7(S, m=2, nsmall=100, Lstar=240):
    t = time.time()
    k = 7
    L1 = 7*math.log(nsmall)
    ok1, n1, w1 = scan(S, 7, m, L1, 180, lambda L: 0.25)
    ok2, n2, w2 = scan(S, 7, m, 180, Lstar, lambda L: 0.02)
    Rt = tail_bound(S, 7, Lstar, m)
    small = direct_small_k7(nsmall + 1)
    ok = ok1 and ok2 and bool(Rt < 1)
    say('k=7: every n >= %d certified: %d + %d cells on [%.3f, %d] (max R %.4f at L=%.2f), tail from L*=%d: R <= %.3g; n <= %d checked directly (largest prime found %d); %s  %.0fs' % (
        nsmall, n1, n2, L1, Lstar, max(w1[0], w2[0]), (w1 if w1[0] > w2[0] else w2)[1], Lstar, float(Rt.mid()), nsmall, small[-1][1], 'CERTIFIED' if ok else 'FAILED', time.time() - t))
    assert ok
    return dict(k=7, nsmall=nsmall, cells=n1 + n2, worstR=max(w1[0], w2[0]), Lstar=Lstar, tailR=float(Rt.mid()))
def main():
    say('Ingham3 certificate, python-flint', __import__('flint').__version__, 'prec', ctx.prec, 'quick' if QUICK else 'full')
    assert B > TSTAR
    say('half-plane: B2 =', B2.str(20), '; T* = B2 + 2.2e-254 < B =', B.str(17, radius=False))
    check_S1()
    # N(t) <= t log t for t >= 14 (used in sec. 4.2): T/(2pi) log(T/(2 pi e)) <= T log T/(2 pi) and the rest <= (1 - 1/(2 pi)) T log T at T = 14, growing slower
    t14 = arb(14); assert arb('0.1038')*t14.log() + arb('0.2573')*t14.log().log() + arb('9.3675') < (1 - 1/(2*PI))*t14*t14.log()
    P = build_kln()
    pgrid = [arb(i)/200 for i in range(1, 301)]          # 0.005 .. 1.5
    kernels = {}
    for m in (2,):
        kernels[m] = []
        for dmax in (arb('1e-9'), arb('0.025'), arb('0.15')):
            t = time.time(); kernels[m].append(Kernel(m, dmax, pgrid))
            say('kernel m=%d dmax=%s: %d lambda-points, %.1fs' % (m, dmax.str(3, radius=False), kernels[m][-1].N, time.time() - t))
    S = Setup(P, kernels, pgrid, arb('0.515'), arb('0.0025'))
    say('sigma cells:', len(S.cells), ' Amax', S.Amax.str(6, radius=False), ' Bmax', S.Bmax.str(6, radius=False))
    return S
def setup_variant(S, de):
    """same kernels, KLN parameter sets with delta = de (KLN Thm 1.1 as printed asks delta >= 1)."""
    P = build_kln(de)
    S2 = Setup(P, S.kernels, S.pgrid, S.sig0, arb('0.0025'))
    say('variant delta=%s: sigma cells %d, Amax %s' % (de, len(S2.cells), S2.Amax.str(6, radius=False)))
    return S2
# ------------------------------------------------------------------ sanity check (not part of the proof): Ford's zero-free region instead of B
def vk_sanity(S, k=3, c_vk=57.54):
    """Replace the half-plane by beta < 1 - nu(gamma), nu(t) = 1/(c (log t)^{2/3} (log log t)^{1/3}) (Ford 2002, Dudek's input), keep
       RH(T0), KLN and the same kernel; return log log x0 (floats; a sanity check only)."""
    import mpmath as mp
    mp.mp.dps = 30
    K = S.kernels[2][0]
    Prow = [p for p in S.P if p.name == '0.99'][0]         # one admissible parameter set, valid for every sigma > 0.5141
    c0, c1 = float(Prow.c0.mid()), float(Prow.c1.mid())
    l0 = float(T0.log().mid())
    def nu(y): return 1.0/(c_vk*y**(2/3.0)*math.log(y)**(1/3.0))     # nu as a function of y = log t
    def nu_inv_logt(u):      # log t with nu(t) = u, or l0 if nu(T0) <= u (then T0 is the binding bound)
        if u >= nu(l0): return l0
        a, b = l0, l0
        while nu(b) > u: b *= 2
        for _ in range(200):
            c = (a + b)/2
            if nu(c) > u: a = c
            else: b = c
        return b
    pg = [float(p.mid()) for p in S.pgrid]
    mom0 = [math.log(float(K.mom(q, arb(0)).mid())) for q in S.pgrid]     # log Mom(p; 0) >= log Mom(p; l0) for every l0
    def mom(p, l0log):       # log of an upper bound for Mom(p; e^{l0log})
        i = [j for j, q in enumerate(pg) if q >= p][0]
        if l0log > math.log(4*math.pi):
            mm = 2; return min(mom0[i], math.log((2*mm)**mm*mm/(mm - pg[i])) + (pg[i] - mm)*l0log)
        return mom0[i]
    US = [0.485*math.exp(-i*0.005) for i in range(14000)]
    LT = [max(l0, nu_inv_logt(u)) for u in US]
    def g(u, ltstar, L, D):  # log of L x^{-u} F(1-u)  (F from KLN at sigma = 1-u, zeros only above max(T0, t_sigma))
        s = 1 - u
        lam0log = ltstar - D
        a = 8*u/3; e = 5 - 2*s
        if lam0log >= 0:
            pp = a + e/ltstar
            lF = c0 + c1*s + a*ltstar + e*math.log(ltstar) - pp*lam0log + mom(pp, lam0log)
        else:
            pp = a + e/D
            lF = c0 + c1*s + a*D + e*math.log(D) + mom(pp, -1)
        return math.log(L) - u*L + lF
    def logR(Lx):
        L = float(Lx); D = L/k - math.log(k)
        # integrate exp(g) over u in (0, 1 - 0.515] on a fine log grid (trapezoid in log u; a sanity check, not a bound)
        us = US
        vals = [g(u, lt, L, D) for (u, lt) in zip(US, LT)]
        mx = max(vals)
        tot = 0.0
        for j in range(len(us) - 1):
            tot += 0.5*(math.exp(vals[j] - mx) + math.exp(vals[j + 1] - mx))*(us[j] - us[j + 1])
        return math.log(2) + mx + math.log(tot)
    if 'probe' in sys.argv:
        for lL in (25, 28, 30, 32, 33, 34, 35, 36, 38, 40):
            print('VK probe log log x =', lL, ' log R =', logR(math.exp(lL)), flush=True)
    lo, hi = 25.0, 40.0
    assert logR(math.exp(lo)) > 0 and logR(math.exp(hi)) < 0
    while hi - lo > 0.01:
        mid = (lo + hi)/2
        if logR(math.exp(mid)) < 0: hi = mid
        else: lo = mid
    return hi
if __name__ == '__main__':
    S = main()
    res = {}
    res[7] = run_k7(S)
    for k in (6, 5, 4, 3):
        res[k] = run_k(S, k)
    say('--- robustness: KLN parameter sets with delta = 1 (Theorem 1.1 as printed) ---')
    S1v = setup_variant(S, '1')
    rv = {7: run_k7(S1v)}
    for k in (6, 5, 4, 3):
        L0m = res[k]['L0']
        rv[k] = run_k(S1v, k, bracket=(L0m, L0m + 150), Lstar=None if k != 3 else 2600)
    say('--- sanity check (floats): Ford zero-free region (c = 57.54) in place of B, k = 3 ---')
    lL = vk_sanity(S)
    say('VK: the bound first drops below 1 at log x = e^%.2f, i.e. n0 = exp(e^%.2f)  (Dudek 2016: exp(e^33.217); Cully-Hugill 2023: exp(e^32.892), thesis exp(e^32.537))' % (lL, lL - math.log(3)))
    json.dump(dict(main=res, delta1=rv, vk_loglogx0=lL), open(__file__.replace('cert.py', 'logs/cert_result.json'), 'w'), indent=1)
    say('total %.0fs' % (time.time() - T_START))
    open(__file__.replace('cert.py', 'logs/cert.log'), 'w').write('\n'.join(LOG) + '\n')
