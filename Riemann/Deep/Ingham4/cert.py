#!/usr/bin/env python3
# cert.py -- Ingham4: a prime in (x, x + x^theta] for every real x >= 1 (theta from KLN only), and the conditional gains
# from the Chourasiya-Simonic density (arXiv 2507.15184, unrefereed; b4 raised to Ramare-Zuniga Alterman's published 0.5501).
# One command (from the repo root):
#   uv run --no-project --with python-flint --with mpmath --with sympy python3 Riemann/Deep/Ingham4/cert.py
# Built on Riemann/Deep/Ingham3/cert.py: inputs, KLN machinery (sec. 1), kernel and moments (sec. 3) are copied unchanged;
# new: delta(x) for h = x^theta (cells in log x), the theta tail lemma, small x (sieve + Oliveira e Silva-Herzog-Pardi),
# the CS density as an alternative bound in every sigma-cell (conditional runs only). See DERIVATION.md.
import sys, time, json, math
from flint import arb, acb, ctx
ctx.prec = 128
T_START = time.time()
LOG = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); LOG.append(s)
def up(x):   # an exact arb >= every point of the ball x
    return arb(x.upper())
def lo(x):
    return arb(x.lower())
def exact_max(vals):
    """the maximum of exact numbers (radius-0 balls), computed exactly (no float keys)."""
    m = vals[0]
    for v in vals[1:]:
        assert v.rad() == 0 and m.rad() == 0
        m = m.max(v)
    return m
def fup(x):
    """a Python float >= every point of the ball x (the stored moment tables use floats: half the memory of arbs)."""
    u = x.upper(); f = float(u)
    while arb(f) < u: f = math.nextafter(f, math.inf)
    return f
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
                if i <= istore: store[i] = fup(acc)        # float rounded up
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
            return self.G0 + arb(T[0])
        assert l0.lower() >= 1
        i = int(math.floor(float((l0.log()/self.logr).lower())))
        while i > 0 and not (self.lam[i].upper() <= l0.lower()): i -= 1     # exact: lam_i <= l0 (float rounding cannot pick a point above l0)
        assert 0 <= i <= self.istore and self.lam[i].upper() <= l0.lower()
        return up(l0**p*self.Gv[i] + arb(T[i]))
def p_up(p, pgrid):
    for q in pgrid:
        if q >= p: return q
    raise ValueError('p beyond grid: %s' % p)


# ------------------------------------------------------------------ (new) trapezoid weight: Phi = 1_[0,a] * 1_[0,b], a + b = 1 (DERIVATION sec. 3')
class TrapKernel:
    """Phi = 1_[0,a] * 1_[0,b] (0 < a <= b = 1 - a): Phi >= 0 Lipschitz on [0,1], Phi(0) = Phi(1) = 0, int Phi = ab, sup Phi = a,
       sup/int = 1/b <= 2 (so the PP bound with m = 2 holds). For Re z in [0, dmax] (dmax <= 1):
       |Phi~(z)|/int Phi <= G(l), l = |Im z|, G(l) = (S(a l/2) + 2 dmax min(a, 1/l)) (S(b l/2) + 2 dmax min(b, 1/l)), non-increasing;
       G(l) <= C2/l^2 for l >= 2 pi/a, C2 = (2/a + 2 dmax)(2/b + 2 dmax).  Same moment machinery as Kernel (a = 1/2 is the triangle)."""
    def __init__(self, a, dmax, pgrid, r=arb(1) + arb(1)/500, Lam_fac=10000, lmax_store=None):
        self.a = arb(a); self.b = 1 - self.a; self.dmax, self.pgrid = dmax, pgrid
        if lmax_store is None: lmax_store = max(30.0, float((2*PI/self.a).upper()) + 1)   # stored range must reach the tail threshold
        assert self.a > 0 and self.a <= self.b and dmax <= 1
        logr = r.log(); self.logr = logr
        self.lthr = 2*PI/self.a
        N = int(math.ceil(math.log(float((self.lthr*Lam_fac).mid()))/float(logr.mid()))) + 1
        self.N = N
        lam = [(logr*i).exp() for i in range(N + 1)]
        G = [self.G(l) for l in lam]
        self.G0 = up((1 + 2*dmax*self.a)*(1 + 2*dmax*self.b))
        self.C2 = up((2/self.a + 2*dmax)*(2/self.b + 2*dmax))
        self.lam, self.Gv = lam, G
        Lam = lam[N]
        istore = int(math.ceil(math.log(lmax_store)/float(logr.mid())))
        self.istore = istore
        self.T = {}
        for p in pgrid:
            assert Lam > self.lthr and p < 2
            acc = p*self.C2*Lam**(p - 2)/(2 - p)
            pw = [((logr*i)*p).exp() for i in range(N + 1)]
            store = [None]*(istore + 1)
            for i in range(N - 1, -1, -1):
                acc = acc + G[i]*(pw[i + 1] - pw[i])
                if i <= istore: store[i] = fup(acc)        # float rounded up
            self.T[str(p)] = store
    def G(self, l):
        a, b, dm = self.a, self.b, self.dmax
        c1 = 2*dm*a.min(ONE/l) if l > 0 else 2*dm*a
        c2 = 2*dm*b.min(ONE/l) if l > 0 else 2*dm*b
        return up((S_up(a*l/2) + c1)*(S_up(b*l/2) + c2))
    def mom(self, p, l0):
        if l0.lower() >= self.lthr.upper() and (l0.lower() >= self.lam[self.istore].upper()):
            return up(self.C2*lo(l0)**(p - 2)*2/(2 - p))      # l0^p G(l0) + p int_l0^oo G l^{p-1} <= C2 l0^{p-2} (1 + p/(2-p))
        T = self.T[str(p)]
        if l0.upper() <= 1:
            return self.G0 + arb(T[0])                            # G(l0) <= G(0) for every l0 in [0, 1]
        assert l0.lower() >= 1
        i = int(math.floor(float((l0.log()/self.logr).lower())))
        while i > 0 and not (self.lam[i].upper() <= l0.lower()): i -= 1     # exact: lam_i <= l0 (float rounding cannot pick a point above l0)
        assert 0 <= i <= self.istore and self.lam[i].upper() <= l0.lower()
        return up(l0**p*self.Gv[i] + arb(T[i]))
def build_kln_H0(H0new, de_override=None):
    """KLN Lemma 4.14 rows (Table 1 parameters, H = H0 - 1) evaluated with RH(H0) at H0 = H0new (KLN set H0 = 3.0610046e10 'unless
       specified', TeX l. 622; Platt-Trudgian give RH(3e12)); admissibility re-checked; Table 1 is not reproduced (different H0)."""
    global H0, HH, LH0, llH0
    old = (H0, HH, LH0, llH0)
    H0 = arb(H0new); HH = H0 - 1; LH0 = H0.log(); llH0 = LH0.log()
    try:
        P0 = [KP(r, de_override) for r in KLN_ROWS]
        bad = [p.name for p in P0 if not p.ok]
        if bad:
            N0, N1, N2, N3, N4, N5 = Jparts(P0[0].k, P0[0].al)
            say('KLN rows at H0 = %s: not admissible (dropped): %s' % (H0new, bad))
        P = [p for p in P0 if p.ok]
        for p in P:
            assert p.c1 > 0
        say('KLN rows at H0 = %s: A(sigma0) = %s' % (H0new, ' '.join('%s:%s' % (p.name, p.logA(p.s0).exp().str(5, radius=False)) for p in P)))
    finally:
        H0, HH, LH0, llH0 = old
    return P

# ------------------------------------------------------------------ (new) Chourasiya-Simonic, arXiv 2507.15184 v2, Corollary 1.1 + Table 1 (CONDITIONAL)
# N(sigma,T) <= B1 T^{3(1-s)/(2-s)} (log T)^{(7-5s)/(2-s)} + B2 log^2 T + B3 log T for T >= 3e12 and s in [s1, s2] (N counts beta >= sigma).
# Their Lemma 4.2 uses b4 = 0.54 (private communication); the published Ramare-Zuniga Alterman bound (arXiv 2603.25961 v3, Thm 1.2) is
# S_eps(X) <= 0.5501 for 1e12 <= X < 1e16, eps <= 1/50 (0.5103, 0.4905 beyond). K3 is linear in b4 and enters B1 only through
# K3^{(2 s2 - 1)/(2 - s2)} (K3 > 1), so B1' = B1 (0.5501/0.54)^{(2 s2-1)/(2-s2)} (DERIVATION sec. 7).
CS_ROWS = [('0.50000','0.53125','5.069','9.461','167.8','0.012'), ('0.53125','0.56250','6.046','9.168','162.6','0.025'),
('0.56250','0.59375','7.077','8.875','157.4','0.038'), ('0.59375','0.62500','8.185','8.582','152.2','0.052'),
('0.62500','0.65625','9.387','8.290','147.0','0.066'), ('0.65625','0.68750','10.70','7.997','141.9','0.081'),
('0.68750','0.71875','12.11','7.704','136.7','0.097'), ('0.71875','0.75000','13.66','7.411','131.5','0.114'),
('0.75000','0.78125','15.35','7.118','126.3','0.132'), ('0.78125','0.81250','17.21','6.826','121.1','0.150'),
('0.81250','0.84375','19.31','6.533','115.9','0.170'), ('0.84375','0.87500','21.77','6.240','110.7','0.192'),
('0.87500','0.90625','24.83','5.745','101.9','0.215')]
B4_CS, B4_RZA = arb('0.54'), arb('0.5501')
class CSRow:
    def __init__(self, row):
        s1, s2, b1, b2, b3, de = [arb(v) for v in row]
        self.s1, self.s2, self.name = s1, s2, row[0] + '-' + row[1]
        # their Lemma 4.2 needs zeta(1+e) S_e(X) <= b4 e^{gamma e}/e for e = delta/log X and 2 delta/log X, X = T >= T0:
        self.eps = 2*de/T0.log()                          # largest e used
        assert self.eps < arb(1)/50                       # RZA Thm 1.2 (1e12 <= X < 1e16, e <= 1/50: 0.5501; larger X: smaller)
        # zeta(1+e) <= 1 + 1/e and e^{gamma e} >= 1 + gamma e  =>  zeta(1+e) <= (1 + (1-gamma) e) e^{gamma e}/e
        self.b4 = up(B4_RZA*(1 + (1 - arb.const_euler())*self.eps))
        ex = (2*s2 - 1)/(2 - s2)
        self.B1 = up(b1*(self.b4/B4_CS)**ex)              # raised B1
        self.B2 = up(b2 + b3/T0.log())                    # B2 log^2 t + B3 log t <= (B2 + B3/log T0) log^2 t for t >= T0
def a_cs(s): return 3 - 3/(2 - s)        # = 3(1-s)/(2-s)   (concave, decreasing)
def e_cs(s): return 5 - 3/(2 - s)        # = (7-5s)/(2-s)   (concave, decreasing)
def cs_for_cell(s, t, rows):
    """CS data on the sigma-cell [s,t]: max B1', max B2' over overlapping rows; tangent lines at the midpoint m
       (a(sig) <= al0 + al1 sig and e(sig) <= ep0 + ep1 sig for every sig < 2, by concavity)."""
    fs, ft = float(s.mid()), float(t.mid())
    ov = [r for r in rows if float(r.s1.mid()) <= ft + 1e-9 and float(r.s2.mid()) >= fs - 1e-9]   # rows meeting [s,t] (extra rows only raise the max)
    if not ov: return None
    r_lo = min(ov, key=lambda r: float(r.s1.mid())); r_hi = max(ov, key=lambda r: float(r.s2.mid()))
    if not (s >= r_lo.s1 and t <= r_hi.s2): return None      # [s,t] is covered by the (contiguous) rows in ov
    B1 = exact_max([r.B1 for r in ov]); B2 = exact_max([r.B2 for r in ov])
    m = arb((fs + ft)/2)                 # an exact point
    sl = -3/(2 - m)**2
    al1 = sl; al0 = a_cs(m) - sl*m; ep1 = sl; ep0 = e_cs(m) - sl*m
    return dict(B1=B1, lB1=B1.log(), B2=B2, al0=al0, al1=al1, ep0=ep0, ep1=ep1, a_s=a_cs(s), e_s=e_cs(s), rows=[r.name for r in ov])

# ------------------------------------------------------------------ the sigma-cells (Ingham3 sec. 4.4; CS added as an alternative)
LN2 = arb(2).log()
class Setup:
    def __init__(self, P, kernels, pgrid, sig0, dsig, cs_rows=None):
        self.P, self.kernels, self.pgrid = P, kernels, pgrid
        self.dsig = dsig
        self.sig0 = sig0
        self.use_cs = cs_rows is not None
        cells = []
        s = sig0
        while s < B:
            t = (s + dsig).min(B) if (s + dsig) < B else B
            best = None
            for p in P:
                if not (s > p.smin): continue
                v = p.logA(s) + p.logA(t)
                if best is None or v.mid() < best[0].mid(): best = (v, p)
            cs = cs_for_cell(s, t, cs_rows) if self.use_cs else None
            if self.use_cs: assert cs is not None, (s, t)
            cells.append((s, t, best[1], cs))
            if t is B: break
            s = t
        self.cells = cells
        self.Amax = exact_max([arb(c[2].logA(c[1]).exp().upper()) for c in cells])
        self.Bmax = exact_max([arb(c[2].Bc(c[0]).upper()) for c in cells])
        self.ncs = 0
# ------------------------------------------------------------------ the bound for one cell of x (Ingham3 sec. 4-5; delta now from a function of L)
def delta_pow(k):
    kk = arb(k)
    return lambda L: kk*(1/(L/kk).exp()).log1p()            # x = n^k: delta = k log(1 + 1/n), n = e^{L/k}
def delta_theta(th):
    om = 1 - arb(th)
    return lambda L: (-om*L).exp().log1p()                   # h = x^theta: delta = log(1 + x^{theta-1})
def cell_bound(S, dfun, La, Lb, m, detail=False):
    """Upper bound R over all x in [e^La, e^Lb] (delta = dfun(log x) is decreasing in log x)."""
    mm = arb(m)
    La, Lb = arb(La), arb(Lb)
    assert La*(1 - B) > 1                       # log x * x^{sigma-1} (and x^{t-1} - x^{s-1}) decrease in x on the cell for every sigma <= B
    da = up(dfun(La)); db = lo(dfun(Lb))
    K = None
    for KK in S.kernels[m]:
        if da < KK.dmax: K = KK; break
    assert K is not None
    Db = -db.log(); lDb = Db.log(); lam0 = db*T0
    assert Db > 1
    assert (1/db) >= 14
    pS = p_up(1 + 1/Db, S.pgrid)
    Sall = Db/db*K.mom(pS, arb(0))
    tot = ((S.sig0 - 1)*La).exp()*Sall
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
    ncs = 0; parts = []
    for (s, t, P, cs) in S.cells:
        p = p_up(arb(8)/3*(1 - s) + (5 - 2*s)/pden, S.pgrid)
        mom = K.mom(p, lam0)
        E0 = -La + P.c0 + arb(8)/3*Db + 5*Lam
        E1 = La + P.c1 - arb(8)/3*Db - 2*Lam
        main = La*mom*expint(E0, E1, s, t)
        Bt = P.Bc(s)*momB*(((t - 1)*La).exp() - ((s - 1)*La).exp())
        c = up(main + Bt)
        if cs is not None:
            pc = p_up(cs['a_s'] + cs['e_s']/pden, S.pgrid)
            momc = K.mom(pc, lam0)
            E0c = -La + cs['lB1'] + cs['al0']*Db + cs['ep0']*Lam
            E1c = La + cs['al1']*Db + cs['ep1']*Lam
            mainc = La*momc*expint(E0c, E1c, s, t)
            Btc = cs['B2']*momB*(((t - 1)*La).exp() - ((s - 1)*La).exp())
            cc = up(mainc + Btc)
            if cc < c: c = cc; ncs += 1
        tot = tot + c
        if detail: parts.append((float(s.mid()), float(c.mid())))
    Zhigh = 2*tot
    J = int(math.floor(float(((Lb + da)/LN2).upper())))
    PP = mm*(Lb + da)*((da/2).exp()*(-La/2).exp()/4 + J*da.exp()*(-2*La/3).exp()/9 + arb(J).log()/(db*La.exp()))
    Et = 2*(-3*La).exp()
    R = up(Zhigh + PP + Et)
    if detail:
        return R, dict(Zhigh=float(Zhigh.mid()), PP=float(PP.mid()), lam0=float(lam0.mid()), ncs=ncs, Sall_part=float((2*((S.sig0 - 1)*La).exp()*Sall).mid()), parts=parts)
    return R
def expint(E0, E1, s, t):
    """int_s^t exp(E0 + E1 sig) dsig (exact formula, or a crude bound when E1 ~ 0)."""
    if abs(E1) > arb('1e-6'):
        return ((E0 + E1*t).exp() - (E0 + E1*s).exp())/E1
    return (t - s)*(E0 + E1*s).exp().max((E0 + E1*t).exp())*arb('1.000001')
def tail_pow(S, k, Ls, m):
    """Ingham3's R_an(L*) for x = n^k (unchanged; KLN only)."""
    K = S.kernels[m][0]; kk = arb(k); mm = arb(m); L = arb(Ls)
    ns = (L/kk).exp()
    cst = 1/(1 - 1/(2*ns))
    assert up(kk*T0*cst/ns) < 1, 'need lambda0 <= 1 on the tail'
    Dlow = L/kk - kk.log()
    assert Dlow > 1
    pmax = p_up(arb(8)/3*(1 - S.sig0) + 4/Dlow, S.pgrid)
    Mst = K.mom(pmax, arb('0.5'))           # the l0 <= 1 branch: G(0) + p int_1^oo G l^{p-1} dl >= Mom(p; l0) for EVERY l0 in [0, 1]
    p1 = p_up(1 + 1/Dlow, S.pgrid); M1 = K.mom(p1, arb(0))
    one_m = 1 - B
    kap = 1 - arb(8)/(3*kk)
    Xi = L*kap + arb(8)/3*(kk/cst).log() - 2*(L/kk).log()
    assert Xi > 0
    mono = [one_m*(kap - 2/L) > 4/L, (1 - S.sig0 - 1/kk) > 2/L, (1 - 1/kk) > 3/L]
    assert all(mono), mono
    main = L*S.Amax*Mst*(L/kk)**3*(-one_m*Xi).exp()/Xi
    Bt = S.Bmax*Mst*(L/kk)**2*(-one_m*L).exp()
    Sallp = (-(1 - S.sig0)*L).exp()*(cst/kk)*(L/kk).exp()*(L/kk)*M1
    Zhigh = 2*(main + Bt + Sallp)
    assert up(kk/ns) < K.dmax
    J = (L + 1)/LN2
    PP = mm*(L + 1)*((-L/2).exp()*arb(1).exp()/4 + J*arb(1).exp()*(-2*L/3).exp()/9 + J.log()*(cst/kk)*(L/kk - L).exp())
    Et = 2*(-3*L).exp()
    return up(Zhigh + PP + Et)
def tail_theta(S, th, Ls, m):
    """R_an(L) for h = x^theta at L = Ls: valid for every x with log x = L >= Ls and non-increasing in L (DERIVATION sec. 5)."""
    K = S.kernels[m][0]; mm = arb(m); L = arb(Ls); thb = arb(th); om = 1 - thb
    ys = (-om*L).exp()                          # y = x^{theta-1} <= ys
    cst = 1/(1 - ys/2)                          # log(1+y) >= y(1 - y/2): 1/delta <= cst/y
    assert up(T0*ys) < 1, 'need lambda0 = delta T0 <= y T0 < 1 on the tail'
    Dlow = om*L                                 # D = log(1/delta) >= log(1/y) = (1-theta) L
    Dmax = om*L + cst.log()                     # D <= log(cst/y)
    assert Dlow > 1
    pmax = p_up(arb(8)/3*(1 - S.sig0) + 4/Dlow, S.pgrid)
    Mst = K.mom(pmax, arb('0.5'))           # the l0 <= 1 branch: G(0) + p int_1^oo G l^{p-1} dl >= Mom(p; l0) for EVERY l0 in [0, 1]
    p1 = p_up(1 + 1/Dlow, S.pgrid); M1 = K.mom(p1, arb(0))
    one_m = 1 - B
    kap = 1 - arb(8)*om/3
    Xi = L - arb(8)/3*Dmax - 2*Dmax.log()
    assert Xi > 0
    mono = [one_m*(kap - 2/L) > 4/L, (thb - S.sig0) > 2/L, thb > 3/L, kap - 2/L > 0]
    assert all(mono), mono
    main = L*S.Amax*Mst*Dmax**3*(-one_m*Xi).exp()/Xi
    Bt = S.Bmax*Mst*Dmax**2*(-one_m*L).exp()
    Sallp = (-(1 - S.sig0)*L).exp()*cst*(om*L).exp()*Dmax*M1
    Zhigh = 2*(main + Bt + Sallp)
    assert up(ys) < K.dmax
    J = (L + 1)/LN2
    PP = mm*(L + 1)*((-L/2).exp()*arb(1).exp()/4 + J*arb(1).exp()*(-2*L/3).exp()/9 + J.log()*cst*(-thb*L).exp())
    Et = 2*(-3*L).exp()
    return up(Zhigh + PP + Et)

# ------------------------------------------------------------------ small x for h = x^theta (DERIVATION sec. 6)
GAP_OSHP = 1476            # every gap between consecutive primes below 4e18 is <= 1476 (Oliveira e Silva-Herzog-Pardi, Math. Comp. 83 (2014))
X_OSHP = 4*10**18
SIEVE_N = 10**7
def gap_records(N):
    """maximal-gap records (g, p) for primes p < N (sieve of Eratosthenes), checked against the known initial records."""
    sv = bytearray([1])*(N + 1); sv[0] = sv[1] = 0
    for i in range(2, int(N**0.5) + 1):
        if sv[i]: sv[i*i::i] = bytearray(len(range(i*i, N + 1, i)))
    recs = []; prev = 2; gmax = 0; i = 3
    while True:
        j = sv.find(1, i)
        if j < 0: break
        g = j - prev
        if g > gmax: gmax = g; recs.append((g, prev))
        prev = j; i = j + 1
    return recs, prev
def small_x_theta(th_num, th_den, L1):
    """for 1 <= x < e^L1: a prime in (x, x + x^theta], theta = th_num/th_den (exact integer arithmetic)."""
    recs, plast = gap_records(SIEVE_N)
    known = [(1,2),(2,3),(4,7),(6,23),(8,89),(14,113),(18,523),(20,887),(22,1129),(34,1327),(36,9551),(44,15683),(52,19609),
             (72,31397),(86,155921),(96,360653),(112,370261),(114,492113),(118,1349533),(132,1357201),(148,2010733),(154,4652353)]
    assert recs[:len(known)] == known, recs[:25]
    # (i) 1 <= x < 2: x + x^theta >= 2.  (ii) p_n <= x < p_{n+1} <= SIEVE_N: need g_n <= p_n^theta; enough for the records:
    for (g, p) in recs:
        assert g**th_den <= p**th_num, (g, p)
    # (iii) plast <= x < e^L1 <= 4e18 (plast = the largest prime below SIEVE_N, whose gap is not in recs):
    #       p_n <= x < p_{n+1} with p_n < 4e18, so p_{n+1} - x <= p_{n+1} - p_n <= 1476 <= plast^theta <= x^theta
    assert GAP_OSHP**th_den <= plast**th_num
    assert up(arb(L1).exp()) < X_OSHP
    worst = max(recs, key=lambda r: math.log(r[0])/math.log(r[1]) if r[1] > 2 else 0)
    say('small x (theta = %d/%d): %d gap records below %d (largest %s), all g <= p^theta (worst log g/log p = %.4f at p = %d); 1476 <= %d^theta; e^%s < 4e18' % (
        th_num, th_den, len(recs), SIEVE_N, recs[-1], math.log(worst[0])/math.log(worst[1]), worst[1], plast, L1))
    return dict(records=len(recs), last=recs[-1], worst_ratio=math.log(worst[0])/math.log(worst[1]))

# ------------------------------------------------------------------ drivers
def scan(S, dfun, m, Lstart, Lend, wfun, stop=True):
    L = Lstart; worst = (-1.0, None); n = 0; first_fail = None
    while L < Lend - 1e-12:
        Lb = min(L + wfun(L), Lend)
        R = cell_bound(S, dfun, L, Lb, m)
        if not (R < 1):
            if stop: return False, n, (float(R.mid()), L)
            if first_fail is None: first_fail = L
        if float(R.mid()) > worst[0]: worst = (float(R.mid()), L)
        L = Lb; n += 1
    return first_fail is None, n, worst
FINE = [0.02]
def theta_width(L):
    return 0.25 if L < 150 else (FINE[0] if 165 <= L < 215 else 0.05)
def direct_small_k(k, nmax):
    """for 1 <= n < nmax: an explicit prime in (n^k, (n+1)^k) (sympy; deterministic below 2^64)."""
    import sympy
    out = []
    for n in range(1, nmax):
        assert (n + 1)**k < 2**64
        p = sympy.nextprime(n**k)
        assert p < (n + 1)**k and sympy.isprime(p)
        out.append((n, int(p)))
    return out
def run_k_alln(S, k, nsmall=100, Lstar=300, w=0.05, m=2, tag=''):
    """every n >= 1: cells from log x = k log(nsmall) to L*, tail from L*, direct check for n < nsmall."""
    t = time.time()
    L1 = k*math.log(nsmall)
    ok, n, worst = scan(S, delta_pow(k), m, L1, Lstar, lambda L: w)
    Rt = tail_pow(S, k, Lstar, m)
    small = direct_small_k(k, nsmall + 1)        # n <= nsmall directly (covers the joint at x = nsmall^k)
    ok = ok and bool(Rt < 1)
    say('k=%d%s: every n >= %d: %d cells on [%.3f, %d] (max R %.4f at L=%.2f), tail from L*=%d: R <= %.3g; 1 <= n <= %d checked directly (largest prime %d); %s  %.0fs' % (
        k, tag, nsmall, n, L1, Lstar, worst[0], worst[1], Lstar, float(Rt.mid()), nsmall, small[-1][1], 'CERTIFIED' if ok else 'FAILED', time.time() - t))
    return dict(k=k, nsmall=nsmall, ok=ok, cells=n, worstR=worst[0], worstL=worst[1], Lstar=Lstar, tailR=float(Rt.mid()))
def run_theta(S, th_str, L1=42, Lstar=300, m=2, tag=''):
    t = time.time()
    num = int(round(float(th_str)*1000)); den = 1000
    assert abs(num/den - float(th_str)) < 1e-12
    dfun = delta_theta(arb(num)/den)
    ok, n, worst = scan(S, dfun, m, L1, Lstar, theta_width)
    if not ok:
        say('theta=%s%s: FAILS at the cell starting at log x = %.2f (R = %.4f)  %.0fs' % (th_str, tag, worst[1], worst[0], time.time() - t))
        return dict(theta=th_str, ok=False, failL=worst[1], failR=worst[0])
    Rt = tail_theta(S, arb(num)/den, Lstar, m)
    ok = bool(Rt < 1)
    small = small_x_theta(num, den, L1)
    say('theta=%s%s: %d cells on [%d, %d] all R < 1 (max %.4f at log x = %.2f); tail from L*=%d: R <= %.3g; small x below e^%d: gap records + OeS-H-P; %s  %.0fs' % (
        th_str, tag, n, L1, Lstar, worst[0], worst[1], Lstar, float(Rt.mid()), L1, 'CERTIFIED' if ok else 'FAILED', time.time() - t))
    return dict(theta=th_str, ok=ok, cells=n, worstR=worst[0], worstL=worst[1], Lstar=Lstar, tailR=float(Rt.mid()), small=small)
def fine_peak(S, th, a, b, w=0.004, tag=''):
    """the largest cell bound on [a, b] with cells of width w (no stop): shows that theta fails by a margin, not by cell width."""
    t = time.time(); num = int(round(float(th)*1000)); dfun = delta_theta(arb(num)/1000)
    ok, n, worst = scan(S, dfun, 2, a, b, lambda L: w, stop=False)
    say('theta=%s%s, cells of width %g on [%g, %g]: max R %.4f at log x = %.3f (%s)  %.0fs' % (th, tag, w, a, b, worst[0], worst[1], 'all < 1' if ok else 'exceeds 1', time.time() - t))
    return dict(theta=th, width=w, maxR=worst[0], at=worst[1], ok=ok)
def probe_theta(S, th, Ls):
    """max over the given cell starts of the cell bound (width 0.02) -- a quick locator, not a certificate."""
    dfun = delta_theta(arb(th))
    best = (-1, None)
    for L in Ls:
        R = float(cell_bound(S, dfun, L, L + 0.02, 2).mid())
        if R > best[0]: best = (R, L)
    return best
def find_L0(S, dfun, m, lo_L, hi_L, w):
    assert cell_bound(S, dfun, hi_L, hi_L + w, m) < 1 and not (cell_bound(S, dfun, lo_L, lo_L + w, m) < 1)
    a, b = lo_L, hi_L
    while b - a > 0.005:
        c = (a + b)/2
        if cell_bound(S, dfun, c, c + w, m) < 1: b = c
        else: a = c
    return math.ceil(b*100)/100
PLAN = {3: ((1700, 1950), 0.25, 2400), 4: ((380, 500), 0.1, 640), 5: ((230, 300), 0.05, 400), 6: ((150, 230), 0.05, 300)}
def run_k(S, k, m=2, bracket=None, Lstar=None, tag=''):
    t = time.time()
    (lo_L, hi_L), w, Ls = PLAN[k]
    if bracket: lo_L, hi_L = bracket
    if Lstar: Ls = Lstar
    dfun = delta_pow(k)
    L0 = find_L0(S, dfun, m, lo_L, hi_L, w)
    ok, n, worst = scan(S, dfun, m, L0, Ls, lambda L: w)
    Rt = tail_pow(S, k, Ls, m)
    ok = ok and bool(Rt < 1)
    say('k=%d%s: x0 = e^%.2f (n0 = e^%.4f = %s); %d cells on [%.2f, %d] all R < 1 (max %.4f at L=%.2f); tail from L*=%d: R <= %.3g; %s  %.0fs' % (
        k, tag, L0, L0/k, '%.4e' % math.exp(L0/k) if L0/k < 700 else '10^%.2f' % (L0/k/math.log(10)), n, L0, Ls, worst[0], worst[1] if worst[1] else -1, Ls, float(Rt.mid()), 'CERTIFIED' if ok else 'FAILED', time.time() - t))
    assert ok
    return dict(k=k, L0=L0, logn0=L0/k, cells=n, worstR=worst[0], worstL=worst[1], Lstar=Ls, tailR=float(Rt.mid()))
A_TRAP = '0.1'                 # the trapezoid's short side a (pipeline A1)
H0_RH = '3e12'                 # RH height for KLN in pipeline A1 (Platt-Trudgian verified RH to 3.000175e12)
def build(cs=False, de_override=None, pipe='A1'):
    """pipe 'A0' = Ingham3 unchanged (triangle, KLN at H0 = 3.0610046e10 with Table 1 reproduction);
       pipe 'A1' = trapezoid a = 1/10 and KLN Lemma 4.14 evaluated with RH(3e12)."""
    rows = [CSRow(r) for r in CS_ROWS]
    for r1, r2 in zip(rows, rows[1:]): assert r1.s2 == r2.s1          # contiguous rows [n/32, (n+1)/32] (exact dyadic decimals)
    if pipe == 'A0':
        P = build_kln(de_override); kern = KERNELS; lh = LH0
    else:
        P = build_kln_H0(H0_RH, de_override); kern = KERNELS_T; lh = arb(H0_RH).log()
    for p in P: assert p.de < lh*(lh.log() - 1)/2     # KLN TeX l. 1077: loglog T/(log kT + 2 delta) decreasing needs this
    S = Setup(P, kern, PGRID, arb('0.515'), arb('0.0025'), rows if cs else None)
    return S
def main_setup(kerns=('tri', 'trap'), dmaxes=('1e-9', '0.025', '0.15')):
    """build only the kernels a part needs (memory: each process stays well below 300 MB)."""
    global KERNELS, KERNELS_T, PGRID
    say('Ingham4 certificate, python-flint', __import__('flint').__version__, 'prec', ctx.prec)
    assert B > TSTAR
    say('half-plane: B2 =', B2.str(20), '; T* = B2 + 2.2e-254 < B =', B.str(17, radius=False))
    check_S1()
    t14 = arb(14); assert arb('0.1038')*t14.log() + arb('0.2573')*t14.log().log() + arb('9.3675') < (1 - 1/(2*PI))*t14*t14.log()
    PGRID = [arb(i)/200 for i in range(1, 301)]
    KERNELS = {2: []}
    for dmax in ([arb(v) for v in dmaxes] if 'tri' in kerns else []):
        t = time.time(); KERNELS[2].append(Kernel(2, dmax, PGRID))
        say('kernel m=2 dmax=%s: %d lambda-points, %.1fs' % (dmax.str(3, radius=False), KERNELS[2][-1].N, time.time() - t))
    KERNELS_T = {2: []}
    for dmax in ([arb(v) for v in dmaxes] if 'trap' in kerns else []):
        t = time.time(); KERNELS_T[2].append(TrapKernel(arb(A_TRAP), dmax, PGRID))
        say('trapezoid kernel a=%s dmax=%s: %d lambda-points, stored to %.1f, %.1fs' % (A_TRAP, dmax.str(3, radius=False), KERNELS_T[2][-1].N, float(KERNELS_T[2][-1].lam[KERNELS_T[2][-1].istore].mid()), time.time() - t))
KERNELS = None; KERNELS_T = None; PGRID = None

def rss_mb():
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20     # macOS: bytes
def run_part(name):
    """one part of the certificate in its own process; returns its results (logged), peak memory included."""
    res = {}
    if name == 'main':
        main_setup(('trap',), ('1e-9', '0.025'))
        build_kln()                                  # KLN implementation check: Table 1 reproduced at H0 = 3.0610046e10 (asserted)
        say('=== MAIN CLAIM (pipeline A1: trapezoid weight a = 1/10, KLN Lemma 4.14 with RH(3e12); refereed zero counts only) ===')
        S1 = build(False, pipe='A1')
        say('sigma cells:', len(S1.cells), ' Amax', S1.Amax.str(6, radius=False), ' Bmax', S1.Bmax.str(6, radius=False))
        res['main'] = run_theta(S1, '0.856'); assert res['main']['ok']
        res['main_minus'] = run_theta(S1, '0.855'); assert not res['main_minus']['ok']
        res['main_minus_fine'] = fine_peak(S1, '0.855', 184, 198)
    elif name == 'a0':
        main_setup(('tri',), ('1e-9', '0.025'))
        say('=== Ingham3 pipeline unchanged (A0: triangle, KLN at H0 = 3.0610046e10) ===')
        S0 = build(False, pipe='A0')
        res['A0'] = run_theta(S0, '0.857'); assert res['A0']['ok']
        res['A0_minus'] = run_theta(S0, '0.856'); assert not res['A0_minus']['ok']
        res['A0_minus_fine'] = fine_peak(S0, '0.856', 183, 197)
    elif name == 'delta1':
        main_setup(('trap',), ('1e-9', '0.025'))
        say('=== robustness: every KLN row with delta = 1 (KLN Thm 1.1 as printed), pipeline A1 ===')
        Sd = build(False, '1', pipe='A1')
        for th in ('0.860', '0.861', '0.862', '0.863'):
            rr = run_theta(Sd, th, tag=' (delta=1)')
            res['delta1_' + th] = rr
            if rr['ok']: break
    elif name == 'cs':
        main_setup(('trap',), ('1e-9', '0.025', '0.15'))
        say('=== CONDITIONAL: A1 + Chourasiya-Simonic 2507.15184 Cor. 1.1 (unrefereed), B1 raised for b4 = 0.5501 (RZA 2603.25961) ===')
        S1c = build(True, pipe='A1')
        for c_ in S1c.cells[::12] + [S1c.cells[-1]]:
            say('  sigma-cell [%.4f, %.4f]: KLN row %s, CS rows %s, B1 %s, B2 %s' % (float(c_[0].mid()), float(c_[1].mid()), c_[2].name, c_[3]['rows'], c_[3]['B1'].str(6, radius=False), c_[3]['B2'].str(6, radius=False)))
        res['cs'] = run_theta(S1c, '0.842', tag=' (A1+CS)'); assert res['cs']['ok']
        res['cs_minus'] = run_theta(S1c, '0.841', tag=' (A1+CS)'); assert not res['cs_minus']['ok']
        res['cs_k6_alln'] = run_k_alln(S1c, 6, tag=' (A1+CS)'); assert res['cs_k6_alln']['ok']
        for k, br in ((5, (200, 240)), (4, (330, 400)), (3, (1300, 1500))):
            res['cs_k%d' % k] = run_k(S1c, k, bracket=br, tag=' (A1+CS)')
    say('part %s: peak memory %.0f MB, %.0fs' % (name, rss_mb(), time.time() - T_START))
    res['_rss_MB_' + name] = rss_mb()
    return res
PARTS = ('main', 'a0', 'delta1', 'cs')
if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'full'
    if mode == 'part':
        import os
        name = sys.argv[2]; res = run_part(name)
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs')
        json.dump(dict(res=res, log=LOG), open(os.path.join(d, 'part_%s.json' % name), 'w'), indent=1)
        sys.exit(0)
    if mode == 'full':
        # one command: the parts run one after another, each in its own small process (< 300 MB); this parent stays idle and small
        import os, subprocess
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs'); allres = {}; alllog = []
        for name in PARTS:
            t = time.time()
            subprocess.run([sys.executable, os.path.abspath(__file__), 'part', name], check=True)
            pj = json.load(open(os.path.join(d, 'part_%s.json' % name)))
            allres.update(pj['res']); alllog += pj['log']
            os.remove(os.path.join(d, 'part_%s.json' % name))
            assert allres['_rss_MB_' + name] < 300, name
        say('all parts done: peak memory per part (MB) %s; total %.0fs' % (', '.join('%s %.0f' % (n, allres['_rss_MB_' + n]) for n in PARTS), time.time() - T_START))
        json.dump(allres, open(os.path.join(d, 'cert_result.json'), 'w'), indent=1)
        open(os.path.join(d, 'cert.log'), 'w').write('\n'.join(alllog + LOG) + '\n')
        sys.exit(0)
    main_setup()
    if mode == 'probe':
        S0 = build(False, pipe='A0'); S1c = build(True, pipe='A0')
        Ls = [float(v) for v in sys.argv[3].split(':')] if len(sys.argv) > 3 else [140 + 2*i for i in range(60)]
        if len(Ls) == 3: Ls = [Ls[0] + Ls[2]*i for i in range(int((Ls[1] - Ls[0])/Ls[2]) + 1)]
        for th in sys.argv[2].split(','):
            r0 = probe_theta(S0, th, Ls); r1 = probe_theta(S1c, th, Ls)
            say('probe theta=%s: KLN max R %.4f at L=%s;  KLN+CS max R %.4f at L=%s' % (th, r0[0], r0[1], r1[0], r1[1]))
        sys.exit(0)
    if mode == 'probek':
        S1c = build(True, pipe='A0')
        for k in [int(v) for v in sys.argv[2].split(',')]:
            a, b, w = [float(v) for v in sys.argv[3].split(':')]
            dfun = delta_pow(k); L = a; fails = []; mx = (0, None)
            while L < b:
                R = float(cell_bound(S1c, dfun, L, L + w, 2).mid())
                if R >= 1: fails.append(L)
                if R > mx[0]: mx = (R, L)
                L += w
            say('probek k=%d (KLN+CS) on [%g,%g] w=%g: max R %.4f at L=%s; failing cells: %d, first %s, last %s' % (k, a, b, w, mx[0], mx[1], len(fails), fails[0] if fails else None, fails[-1] if fails else None))
        sys.exit(0)
    if mode == 'probe3':
        S1c = build(True, pipe='A0'); dfun = delta_pow(3)
        for L in [1200, 1300, 1400, 1500, 1600, 1700, 1800]:
            say('probe3 L=%d: R=%.4f' % (L, float(cell_bound(S1c, dfun, L, L + 0.25, 2).mid())))
        Sd = build(False, '1', pipe='A0')
        for th in ('0.857', '0.858', '0.859', '0.860'):
            r = probe_theta(Sd, th, [170 + i for i in range(50)])
            say('probe delta=1 theta=%s: max R %.4f at L=%s' % (th, r[0], r[1]))
        sys.exit(0)
    if mode == 'fine':
        th = sys.argv[2]; FINE[0] = float(sys.argv[3]); a, b = float(sys.argv[4]), float(sys.argv[5])
        S0 = build(sys.argv[6] == 'cs' if len(sys.argv) > 6 else False, pipe=sys.argv[7] if len(sys.argv) > 7 else 'A1')
        num = int(round(float(th)*1000)); dfun = delta_theta(arb(num)/1000)
        ok, n, worst = scan(S0, dfun, 2, a, b, lambda L: FINE[0], stop=False)
        say('fine theta=%s width %g on [%g,%g]: ok=%s cells=%d max R %.5f at L=%.3f  %.0fs' % (th, FINE[0], a, b, ok, n, worst[0], worst[1], time.time() - T_START))
        sys.exit(0)
    if mode == 'probe2':
        P2 = build_kln_H0('3e12')
        for a in [v for v in sys.argv[2].split(',')]:
            kern = {2: [TrapKernel(arb(a), dm, PGRID) for dm in (arb('1e-9'), arb('0.025'), arb('0.15'))]}
            S2 = Setup(P2, kern, PGRID, arb('0.515'), arb('0.0025'))
            for th in sys.argv[3].split(','):
                r = probe_theta(S2, th, [float(v) for v in sys.argv[4].split(',')])
                say('probe2 a=%s theta=%s (KLN at H0=3e12, trapezoid): max R %.4f at L=%s' % (a, th, r[0], r[1]))
        sys.exit(0)
    if mode == 'probe2b':
        a = arb(sys.argv[2])
        kern = {2: [TrapKernel(a, dm, PGRID) for dm in (arb('1e-9'), arb('0.025'), arb('0.15'))]}
        Ls = [float(v) for v in sys.argv[4].split(',')]
        for (tag, de, cs) in [('A1', None, False), ('A1 delta=1', '1', False), ('A1+CS', None, True)]:
            P2 = build_kln_H0('3e12', de)
            S2 = Setup(P2, kern, PGRID, arb('0.515'), arb('0.0025'), [CSRow(r) for r in CS_ROWS] if cs else None)
            for th in sys.argv[3].split(';')[['A1', 'A1 delta=1', 'A1+CS'].index(tag)].split(','):
                r = probe_theta(S2, th, Ls)
                say('probe2b %s a=%s theta=%s: max R %.4f at L=%s' % (tag, sys.argv[2], th, r[0], r[1]))
        sys.exit(0)
