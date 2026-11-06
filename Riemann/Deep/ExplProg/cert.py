#!/usr/bin/env python3
"""
cert.py -- Arb certificate for explicit prime number theorems in progressions from the zero-free half-plane.
Run (one command, one thread, < 100 MB):
    uv run --no-project --with python-flint --with mpmath python3 cert.py all
Sections: core | linnik (Thms 2, 2') | nonres (Thm 3) | pnt (Thm 1) | lprime (Thm 4) | all.
Every inequality used as a proof step is a certified sign of an arb ball;
floats only choose parameters.  Notation and proofs: DERIVATION.md (section numbers in comments).
Inputs (DERIVATION.md sec. 0):
  H   L(s,chi) != 0 for Re s > B for zeta and every Dirichlet chi (Bebel, Lean-checked: Re s > B_2 + 2.2e-254 < B).
  BMOR  Bennett-Martin-O'Bryant-Rechnitzer, Math. Comp. 90 (2021), Thm 1.1 (zero count, conductor > 1, T >= 5/7).
  ROS  Rosser 1941, Thm 19: |N(T) - (T/2pi)log(T/2pi e) - 7/8| < 0.137 log T + 0.443 log log T + 1.588 (T >= 2), zeta.
  ZETA  sum_rho Re(1/rho) = 1 + gamma/2 - log(4 pi)/2 (Davenport ch. 12); first zeta zero 14.1347... > 14.
  RS62  Rosser-Schoenfeld 1962: theta(x) < 1.01624 x, psi(x) < 1.03883 x (x > 0);
        phi(n) > n/(e^gamma loglog n + 2.50637/loglog n) (n >= 3).
  DUS  Dusart (arXiv 1002.0442, Ramanujan J. 45 (2018)): sum_{p<=x} log p/p <= log x + E + 0.2/log x + 0.2/log^2 x (x >= 2974).
The sieve step of Theorem 2 for q < 100 is in sanity.py (exact integers).
"""
import sys, math, time
from fractions import Fraction as Fr
from flint import arb, ctx

ctx.prec = 200
PI = arb.pi(); EUL = arb.const_euler(); E1 = arb(1).exp()
LOG2 = arb(2).log()

def A(x):
    if isinstance(x, Fr):
        return arb(x.numerator) / x.denominator
    return arb(x)

def pos(x, msg):
    """certified x > 0, else abort"""
    if not (x > 0):
        raise SystemExit("CERT FAIL: %s  (%s)" % (msg, x.str(10)))

def up(x):  # float upper bound (for printing only)
    return float(x.upper())

# ---------------------------------------------------------------- sec 0: the half-plane constant
B2 = (777 - arb(1081).sqrt()) / 852
B = A(Fr(8733819666132076, 10 ** 16))
b = 1 - B
INV_b = 1 / b

def core():
    pos(B - B2 - arb(10) ** -17, "B > B_2 + 1e-17 > T_* = B_2 + 2.2e-254")
    pos(B2 - A(Fr(87338196661320752, 10 ** 17)) + arb(10) ** -17, "B_2 = 0.87338196661320752 (to 1e-17)")
    print("[core] B_2 = %s, B = %s (> B_2 + 1e-17), 1/(1-B) = %s" % (B2.str(20), B.str(17), INV_b.str(12)))
    print("[core] ZETA: sum_rho 1/gamma^2 <= %s" % SUM_INV_G2.str(8))
    # sec 2.1: BMOR right-hand side at conductor 9 is positive for t >= 5/7 (each term non-decreasing there)
    R9 = (S0 / PI) * (9 * S0 / (2 * PI * E1)).log() + A(Fr(1, 4)) + CB * (9 * (S0 + 2) / (2 * PI)).log() \
         + 2 * (1 + (9 * (S0 + 2) / (2 * PI)).log()).log() - A(Fr(1, 2))
    pos(R9, "R(5/7; 9) > 0")
    pos(S0 - 2 * PI / 9, "5/7 > 2 pi/9 (monotonicity of (t/pi) log(9t/2 pi e))")
    pos(9 * (S0 + 2) / (2 * PI) - 1, "ell(5/7; 9) > 0")
    print("[core] BMOR at conductor 9: R(5/7;9) = %s > 0" % R9.str(8))

# ---------------------------------------------------------------- sec 2: zero sums by partial summation
CB = A(Fr(22737, 100000))
UE = 2 * PI - 2                      # log((t+2)/2pi) > 0 iff t > UE
S0 = A(Fr(5, 7))

def Nplus_dir(t, L):
    """BMOR upper bound, decomposed (sec 2.1): N(t,chi) <= L a(t) + 2 log(1+L) + b(t) + e(t)/(1+L), t >= 5/7, L >= log 9"""
    lt = ((t + 2) / (2 * PI)).log()
    e = 2 * lt if t > UE else arb(0)
    return L * (t / PI + CB) + 2 * (1 + L).log() + (t / PI) * (t / (2 * PI * E1)).log() + A(Fr(1, 4)) + CB * lt - A(Fr(1, 2)) + e / (1 + L)

def envelope(n0, k1, k2):
    """G(t) = min(n0, k1/t, k2/t^2) as pieces [(lo, hi, c, m)] meaning c t^-m on [lo, hi] (hi = None for infinity)"""
    if k1 * k1 <= n0 * k2:           # three pieces (decided on the balls; if undecided, fall back)
        ta = k1 / n0; tb = k2 / k1
        return [(arb(0), ta, n0, 0), (ta, tb, k1, 1), (tb, None, k2, 2)]
    if k1 * k1 > n0 * k2:
        tc = (k2 / n0).sqrt()
        return [(arb(0), tc, n0, 0), (tc, None, k2, 2)]
    raise SystemExit("envelope undecided")

def Gval(pieces, t):
    for (lo, hi, c, m) in pieces:
        if hi is None or t <= hi:
            return c / t ** m if m else c
    raise SystemExit("Gval")

def F_dir(t, m, L):
    """antiderivative of rho_main(t) t^-m, rho_main = L/pi + log(t/2pi)/pi + CB/(t+2)"""
    if m == 0:
        return (L / PI) * t + (t * (t / (2 * PI)).log() - t) / PI + CB * (t + 2).log()
    if m == 1:
        return (L / PI) * t.log() + (t / (2 * PI)).log() ** 2 / (2 * PI) + CB / 2 * (t / (t + 2)).log()
    if m == 2:
        if t is None:
            return arb(0)
        return -(L / PI) / t - ((t / (2 * PI)).log() + 1) / (PI * t) + CB * (-1 / (2 * t) + ((t + 2) / t).log() / 4)

def F_e(t, m, L):
    """antiderivative of 2/((1+L)(t+2)) t^-m"""
    k = 2 / (1 + L)
    if m == 0:
        return k * (t + 2).log()
    if m == 1:
        return k / 2 * (t / (t + 2)).log()
    if m == 2:
        if t is None:
            return arb(0)
        return k * (-1 / (2 * t) + ((t + 2) / t).log() / 4)

def amax(x, y):
    return x if x > y else (y if y > x else arb.union(x, y).upper() + 0 * x)

def zsum_coeffs(n0, k1, k2):
    """sec 2.2: (za, zb, ze) with  sum_{rho in Z(chi)} G(|gamma_rho|) <= Z(L) := L za + 2 n0 log(1+L) + zb + ze/(1+L)
    for every primitive chi of conductor <= e^L, L >= log 9;  G = min(n0, k1/t, k2/t^2)."""
    pcs = envelope(n0, k1, k2)
    a0 = S0 / PI + CB
    b0 = (S0 / PI) * (S0 / (2 * PI * E1)).log() + A(Fr(1, 4)) + CB * ((S0 + 2) / (2 * PI)).log() - A(Fr(1, 2))
    za = a0 * n0; zb = b0 * n0; ze = arb(0)          # N+(5/7) G(0); e(5/7) = 0
    def Pa(t, m):   # antiderivative of t^-m / pi
        if t is None:
            return arb(0)
        return (t if m == 0 else (t.log() if m == 1 else -1 / t)) / PI
    def Pb(t, m):   # antiderivative of (log(t/2pi)/pi + CB/(t+2)) t^-m
        if t is None:
            return arb(0)
        if m == 0:
            return (t * (t / (2 * PI)).log() - t) / PI + CB * (t + 2).log()
        if m == 1:
            return (t / (2 * PI)).log() ** 2 / (2 * PI) + CB / 2 * (t / (t + 2)).log()
        return -((t / (2 * PI)).log() + 1) / (PI * t) + CB * (-1 / (2 * t) + ((t + 2) / t).log() / 4)
    def Pe(t, m):   # antiderivative of 2 t^-m/(t+2)
        if t is None:
            return arb(0)
        if m == 0:
            return 2 * (t + 2).log()
        if m == 1:
            return (t / (t + 2)).log()
        return 2 * (-1 / (2 * t) + ((t + 2) / t).log() / 4)
    for (lo, hi, c, m) in pcs:
        u = lo if lo > S0 else S0
        if hi is not None and not (hi > u):
            continue
        if hi is None and m != 2:
            raise SystemExit("unbounded piece")
        za += c * (Pa(hi, m) - Pa(u, m))
        zb += c * (Pb(hi, m) - Pb(u, m))
        ue = u if u > UE else UE
        if hi is None or hi > ue:
            ze += c * (Pe(hi, m) - Pe(ue, m))
    return za, zb, ze

def Zof(co, n0, L):
    za, zb, ze = co
    return L * za + 2 * n0 * (1 + L).log() + zb + ze / (1 + L)

def zsum_dir(L, n0, k1, k2):
    return Zof(zsum_coeffs(n0, k1, k2), n0, L)

# zeta (ROS): N_z(t) <= (t/pi) log(t/2pi e) + 4.926 + 0.274 log t + 0.886 log log t  (t >= 14; N_z = 0 below 14.13)
CR0, CR1, CR2 = A(Fr(4926, 1000)), A(Fr(274, 1000)), A(Fr(886, 1000))
def Nplus_z(t):
    return (t / PI) * (t / (2 * PI * E1)).log() + CR0 + CR1 * t.log() + CR2 * t.log().log()

def zsum_zeta(n0, k1, k2):
    """sum over zeta zeros of G(|gamma|) <= N+(14) G(14) + int_14^oo rho_z G"""
    pcs = envelope(n0, k1, k2)
    T = arb(14)
    tot = Nplus_z(T) * Gval(pcs, T)
    for (lo, hi, c, m) in pcs:
        u = lo if lo > T else T
        if hi is not None and not (hi > u):
            continue
        # rho_z = log(t/2pi)/pi + 0.274/t + 0.886/(t log t)
        def P(t):
            if m == 0:
                return (t * (t / (2 * PI)).log() - t) / PI + CR1 * t.log() + CR2 * t.log().log()
            if m == 1:
                return (t / (2 * PI)).log() ** 2 / (2 * PI) - CR1 / t
            if m == 2:
                return -((t / (2 * PI)).log() + 1) / (PI * t) - CR1 / (2 * t * t)
        if hi is None:
            main = -P(u)
            corr = (CR2 / u.log() / u if m == 1 else (CR2 / u.log() / (2 * u * u) if m == 2 else None))
        else:
            main = P(hi) - P(u)
            corr = (CR2 / u.log() * (1 / u - 1 / hi) if m == 1 else (CR2 / u.log() * (1 / (2 * u * u) - 1 / (2 * hi * hi)) if m == 2 else arb(0)))
        if corr is None:
            raise SystemExit("zeta piece m=0 unbounded")
        tot += c * (main + corr)
    return tot

# ZETA: sum_rho beta/|rho|^2 = 1 + gamma/2 - log(4pi)/2;  sum_rho 1/gamma^2 <= 2 (1 + B^2/14^2) * that
SUM_BETA = 1 + EUL / 2 - (4 * PI).log() / 2
SUM_INV_G2 = 2 * (1 + B * B / 196) * SUM_BETA

# ---------------------------------------------------------------- sec 3: the Linnik / non-residue kernel
# g(u) on [-ell, 0]: rises linearly from 0 at -ell to 1 at -c, falls to 0 at 0 (rise length a = ell - c);
# f(u) = g(u) e^{-Bu}, w(t) = f(log t), W(s) = int f(u) e^{su} du.
class Tri:
    def __init__(self, ell, a):
        self.ell, self.a = A(ell), A(a)
        self.c = self.ell - self.a
        ell, a, c = self.ell, self.a, self.c
        self.n0 = ell / 2                                   # ||g||_1
        self.n1 = arb(2)                                    # ||g'||_1
        self.n2 = 2 * (1 / a + 1 / c)                       # total variation of g'
        self.W1 = self.I(b)                                 # W(1) = int g e^{bu}
        self.W0 = self.I(-B)                                # W(0) = int f
        m = A(Fr(1, 2)) - B                                 # left kernel g_L = g e^{mu}
        self.m0 = self.I(m)
        J1 = (1 / a) * ((-m * c).exp() - (-m * ell).exp()) / m + (1 / c) * (1 - (-m * c).exp()) / m   # int |g'| e^{mu}
        J2 = (-m * ell).exp() / a + (1 / a + 1 / c) * (-m * c).exp() + 1 / c                            # int e^{mu} |dg'|
        self.m1 = J1 + abs(m) * self.m0
        self.m2 = J2 + 2 * abs(m) * J1 + m * m * self.m0
        pos(a - 1 / B, "peak beyond 1/B")
        self.winf = amax((B * c).exp(), (B * ell - 1).exp() / (a * B))   # ||w||_inf
    def I(self, m):
        ell, a, c = self.ell, self.a, self.c
        return ((-m * ell).exp() / a - (1 / a + 1 / c) * (-m * c).exp() + 1 / c) / (m * m)
    def ks(self, Y):
        """sec 3.2: sup over beta in [1/2,B] of x^{beta-B}|W(beta+i gamma)| <= min(n0, k1/g, k2/g^2); left analogues"""
        eY = E1 * Y
        k1 = self.n1 + self.n0 / eY
        k2 = self.n2 + 2 * self.n1 / eY + (2 / eY) ** 2 * self.n0
        k1L = self.m1 + self.m0 / eY
        k2L = self.m2 + 2 * self.m1 / eY + (2 / eY) ** 2 * self.m0
        return (self.n0, k1, k2), (self.m0, k1L, k2L)

def omega_parts(K, Y, L):
    """Omega(x;q) <= (x^B Z + x^(1/2) ZL)/2 and Omega_1 <= x^B z1 + x^(1/2) z1L, for x with log(alpha x) >= Y"""
    (n0, k1, k2), (m0, k1L, k2L) = K.ks(Y)
    Lq = L if L > arb(9).log() else arb(9).log()
    Z = zsum_dir(Lq, n0, k1, k2)
    ZL = zsum_dir(Lq, m0, k1L, k2L)
    z1 = SUM_INV_G2 * k2 / 2
    z1L = SUM_INV_G2 * k2L / 2
    return Z, ZL, z1, z1L

KER = Tri(7, Fr(49, 10))

# ---------------------------------------------------------------- sec 4: Linnik
def linnik_cell(K, c2, L0, L1):
    """for every q with e^L0 <= q <= e^L1 (q >= 3) and (a,q)=1 there is a prime p = a (q), p <= x := (c2 e^L0 L0)^(1/b).
    Returns certified margin (LHS - RHS)/LHS of sec 4 (1)."""
    X = (A(c2) * L0.exp() * L0).log() / b
    x = X.exp(); Y = X - K.ell
    pos(Y - 10, "Y >= 10")
    Z, ZL, z1, z1L = omega_parts(K, Y, L1)
    xB = (B * X).exp(); xh = (X / 2).exp()
    Q1 = L1.exp()
    om1 = xB * z1 + xh * z1L
    triv = K.W0 / (1 - (-2 * Y).exp())
    Dq = K.winf * ((L1 / LOG2) * K.ell + L1)
    pp = K.winf * A(Fr(101624, 100000)) * (xh + (X / 3).exp() * X / LOG2)
    rhs = om1 + Q1 * ((xB * Z + xh * ZL) / 2 + triv + Dq + pp)
    lhs = x * K.W1
    return (lhs - rhs) / lhs, Z

def Z_over_L_majorant(c, nn, Lm):
    """sec 4.3: sup_{L >= Lm} Z(L)/L.  Z(L)/L = za + h(L)/L + ze/(L(1+L)), h(L) = 2 nn log(1+L) + zb.
    h(L)/L is non-increasing on [Lm, oo) when 2 nn (log(1+Lm) - Lm/(1+Lm)) + zb >= 0 (its derivative has numerator
    2 nn (L/(1+L) - log(1+L)) - zb, and log(1+L) - L/(1+L) increases); otherwise zb is replaced by max(zb, 0)."""
    za, zb, ze = c
    cond = 2 * nn * ((1 + Lm).log() - Lm / (1 + Lm)) + zb
    zbb = zb if cond >= 0 else amax(zb, arb(0))
    return za + (2 * nn * (1 + Lm).log() + zbb) / Lm + amax(ze, arb(0)) / (Lm * (1 + Lm))

def small_terms_linnik(K, c2, Lm):
    """tail (sec 4.3): for all q >= e^Lm, with x = (c2 q log q)^(1/b): returns certified margin of
    c2 W1 > Z(L)/(2L) + [small]/L using monotone majorants at L = Lm."""
    X = (A(c2) * Lm.exp() * Lm).log() / b
    Y = X - K.ell
    pos(Y - 1, "alpha x > e"); pos(Lm - arb(9).log(), "log q >= log 9")
    (n0, k1, k2), (m0, k1L, k2L) = K.ks(Y)
    co = zsum_coeffs(n0, k1, k2); coL = zsum_coeffs(m0, k1L, k2L)
    def ZL_over_L(c, nn):
        return Z_over_L_majorant(c, nn, Lm)
    xB = (B * X).exp(); xhB = ((A(Fr(1, 2)) - B) * X).exp(); x3B = ((A(Fr(1, 3)) - B) * X).exp()
    z1 = SUM_INV_G2 * k2 / 2; z1L = SUM_INV_G2 * k2L / 2
    rhs = ZL_over_L(co, n0) / 2 + xhB * ZL_over_L(coL, m0) / 2
    rhs += (z1 + xhB * z1L) / (Lm.exp() * Lm)
    rhs += K.W0 / ((1 - (-2 * Y).exp()) * xB * Lm)
    rhs += K.winf * (K.ell / LOG2 + 1) / xB
    rhs += K.winf * A(Fr(101624, 100000)) * (xhB + x3B * X / LOG2) / Lm
    lhs = A(c2) * K.W1
    pos(X - 2, "X > 1/(B-1/3)")
    return (lhs - rhs) / lhs

def phi_lo_q(L):
    """phi(q) >= f(q) := q/(e^gamma log log q + 2.50637/log log q) for q >= 3 (RS62 Thm 15); f is increasing for q >= 16
    (f' has the sign of D - qD', D - qD' = e^gamma (log log q - 1/log q) + 2.50637 (1/log log q + 1/(log q (log log q)^2)) > 0)."""
    return L.exp() / (EUL.exp() * L.log() + PHI_C / L.log())

def small_terms_linnik_phi(K, c2, Lm):
    """Theorem 2': for all q >= e^Lm (>= 16), with x = (c2 phi(q) log q)^(1/b) >= x_lo := (c2 f(q) log q)^(1/b):
    c2 W1 > Z(L)/(2L) + small/L with every small term evaluated at x_lo(e^Lm) (non-increasing in q since f increases)."""
    pos(Lm - arb(16).log(), "q >= 16")
    fq = phi_lo_q(Lm)
    X = (A(c2) * fq * Lm).log() / b
    Y = X - K.ell
    pos(Y - 1, "alpha x > e"); pos(Lm - arb(9).log(), "log q >= log 9")
    (n0, k1, k2), (m0, k1L, k2L) = K.ks(Y)
    co = zsum_coeffs(n0, k1, k2); coL = zsum_coeffs(m0, k1L, k2L)
    xB = (B * X).exp(); xhB = ((A(Fr(1, 2)) - B) * X).exp(); x3B = ((A(Fr(1, 3)) - B) * X).exp()
    z1 = SUM_INV_G2 * k2 / 2; z1L = SUM_INV_G2 * k2L / 2
    rhs = Z_over_L_majorant(co, n0, Lm) / 2 + xhB * Z_over_L_majorant(coL, m0, Lm) / 2
    rhs += (z1 + xhB * z1L) / (fq * Lm)
    rhs += K.W0 / ((1 - (-2 * Y).exp()) * xB * Lm)
    rhs += K.winf * (K.ell / LOG2 + 1) / xB
    rhs += K.winf * A(Fr(101624, 100000)) * (xhB + x3B * X / LOG2) / Lm
    pos(X - 2, "X > 1/(B-1/3)")
    lhs = A(c2) * K.W1
    return (lhs - rhs) / lhs

def linnik_range(K, c2, Lstart, Lm, step=Fr(1, 200)):
    """cells [L0, L1] covering [Lstart, Lm], then the tail; returns worst cell margin"""
    L0 = A(Lstart); worst = None; n = 0
    st = A(step)
    while L0 < A(Lm):
        L1 = L0 + st * amax(arb(1), L0 / 10)
        mg, _ = linnik_cell(K, c2, L0, L1)
        pos(mg, "Linnik cell L0=%s c2=%s" % (L0.str(8), c2))
        worst = mg if worst is None or mg < worst else worst
        L0 = L1; n += 1
    mt = small_terms_linnik(K, c2, L0)
    pos(mt, "Linnik tail L>=%s c2=%s" % (L0.str(8), c2))
    return worst, mt, n, L0

def limit_constant(K):
    """c_inf = za/(2 W(1)) at Y = oo (k1 = ||g'||, k2 = TV(g')): for every eps > 0, p(a,q) <= ((c_inf + eps) q log q)^(1/b)
    and n_chi <= ((c_inf + eps) log q)^(1/b) for q >= q(eps)"""
    za, zb, ze = zsum_coeffs(K.n0, K.n1, K.n2)
    return za / (2 * K.W1)

def linnik():
    print("[linnik] limit constant (Y = oo): c_inf = %s" % limit_constant(KER).str(6))
    print("[linnik] kernel ell=7, rise 4.9: W(1) = %s, W(0) = %s, ||w||_inf <= %s, n2 = %s"
          % (KER.W1.str(10), KER.W0.str(10), KER.winf.str(8), KER.n2.str(8)))
    rows = [(Fr(1), 2), (Fr(79, 100), 2), (Fr(75, 100), 3), (Fr(72, 100), 4), (Fr(67, 100), 6), (Fr(62, 100), 10),
            (Fr(56, 100), 20), (Fr(53, 100), 30), (Fr(51, 100), 50), (Fr(49, 100), 100), (Fr(46, 100), 1000)]
    for (c2, e10) in rows:
        Ls = arb(10) ** e10
        Ls = Fr(int(float(Ls.log().mid()) * 10 ** 6), 10 ** 6)          # L_start <= log(10^e10)
        pos(arb(10 ** e10).log() - A(Ls), "L_start <= log 10^e")
        mt = small_terms_linnik(KER, c2, A(Ls))
        pos(mt, "Linnik q >= 10^%d c2=%s" % (e10, c2))
        mg, _ = linnik_cell(KER, c2, A(Ls), A(Ls))           # cross-check: direct evaluation at q = e^Ls
        print("[linnik] all q >= 10^%d: p(a,q) <= (%s q log q)^(1/(1-B)); margin %s (direct at q = 10^%d: %s)"
              % (e10, c2, mt.str(4), e10, mg.str(4)), flush=True)
        mp_ = small_terms_linnik_phi(KER, c2, A(Ls))
        pos(mp_, "Linnik phi-version q >= 10^%d c2=%s" % (e10, c2))
        print("[linnik]   and p(a,q) <= (%s phi(q) log q)^(1/(1-B)) for q >= 10^%d; margin %s" % (c2, e10, mp_.str(4)), flush=True)

# ---------------------------------------------------------------- sec 5: least character non-residue
def nonres_cell(K, c3, L0, L1):
    """every nonprincipal chi mod q, e^L0 <= q <= e^L1, has chi(n) not in {0,1} for some n <= x := (c3 L0)^(1/b)"""
    X = (A(c3) * L0).log() / b
    x = X.exp(); Y = X - K.ell
    pos(Y - 10, "Y >= 10 (nonres)")
    Z, ZL, z1, z1L = omega_parts(K, Y, L1)
    xB = (B * X).exp(); xh = (X / 2).exp()
    om1 = xB * z1 + xh * z1L
    triv = K.W0 / (1 - (-2 * Y).exp())
    Dq = K.winf * ((L1 / LOG2) * K.ell + L1)
    rhs = om1 + (xB * Z + xh * ZL) / 2 + 2 * triv + 2 * Dq
    lhs = x * K.W1
    return (lhs - rhs) / lhs

def nonres_tail(K, c3, Lm):
    X = (A(c3) * Lm).log() / b
    Y = X - K.ell
    pos(Y - 1, "alpha x > e"); pos(Lm - arb(9).log(), "log q >= log 9")
    (n0, k1, k2), (m0, k1L, k2L) = K.ks(Y)
    co = zsum_coeffs(n0, k1, k2); coL = zsum_coeffs(m0, k1L, k2L)
    def ZL_over_L(c, nn):
        return Z_over_L_majorant(c, nn, Lm)
    xB = (B * X).exp(); xhB = ((A(Fr(1, 2)) - B) * X).exp()
    z1 = SUM_INV_G2 * k2 / 2; z1L = SUM_INV_G2 * k2L / 2
    rhs = ZL_over_L(co, n0) / 2 + xhB * ZL_over_L(coL, m0) / 2 + (z1 + xhB * z1L) / Lm
    rhs += 2 * K.W0 / ((1 - (-2 * Y).exp()) * xB * Lm) + 2 * K.winf * (K.ell / LOG2 + 1) / xB
    lhs = A(c3) * K.W1
    return (lhs - rhs) / lhs

def nonres_range(K, c3, Lstart, Lm, step=Fr(1, 100)):
    L0 = A(Lstart); worst = None; n = 0; st = A(step)
    while L0 < A(Lm):
        L1 = L0 + st * amax(arb(1), L0 / 10)
        mg = nonres_cell(K, c3, L0, L1)
        pos(mg, "nonres cell L0=%s c3=%s" % (L0.str(8), c3))
        worst = mg if worst is None or mg < worst else worst
        L0 = L1; n += 1
    mt = nonres_tail(K, c3, L0)
    pos(mt, "nonres tail")
    return worst, mt, n, L0

def trivial_nonres(LT):
    """q - 1 <= (log q)^(1/b) for 3 <= q <= e^LT: on cells [l0,l1], (1/b) log l0 >= log(e^l1 - 1)"""
    l0 = arb(3).log(); n = 0; worst = None
    while l0 < LT:
        h = A(Fr(1, 1000)) if l0 > 20 else A(Fr(1, 100))
        l1 = l0 + h
        d = INV_b * l0.log() - ((l1).exp() - 1).log()
        pos(d, "trivial nonres at l0=%s" % l0.str(6))
        worst = d if worst is None or d < worst else worst
        l0 = l1; n += 1
    return worst, n

def nonres():
    LT = A(Fr(256, 10))
    w, n = trivial_nonres(LT)
    print("[nonres] 3 <= q <= e^25.6: n_chi <= q - 1 <= (log q)^(1/(1-B)) (%d cells, worst log-margin %s)" % (n, w.str(4)))
    rows = [(Fr(1), Fr(25)), (Fr(62, 100), Fr(25)), (Fr(60, 100), Fr(30)), (Fr(56, 100), Fr(50)), (Fr(52, 100), Fr(100)),
            (Fr(47, 100), Fr(1000)), (Fr(46, 100), Fr(10000))]
    for (c3, Ls) in rows:
        mt = nonres_tail(KER, c3, A(Ls))
        pos(mt, "nonres log q >= %s c3=%s" % (Ls, c3))
        mg = nonres_cell(KER, c3, A(Ls), A(Ls))
        print("[nonres] all q >= e^%s: n_chi <= (%s log q)^(1/(1-B)); margin %s (direct at q = e^%s: %s)"
              % (float(Ls), c3, mt.str(4), float(Ls), mg.str(4)), flush=True)

# ---------------------------------------------------------------- sec 6: uniform PNT in progressions
HALF = A(Fr(1, 2))
PSI_RS = A(Fr(103883, 100000))       # psi(y) < 1.03883 y  (RS62 Thm 12)
PHI_C = A(Fr(250637, 100000))        # RS62 Thm 15
def phi_lo(L0, L1):
    """lower bound for phi(q), e^L0 <= q <= e^L1, L0 > 1:  max(sqrt(q/2), q/(e^gamma loglog q + 2.50637/loglog q))"""
    v1 = (L0 / 2).exp() / arb(2).sqrt()
    v2 = L0.exp() / (EUL.exp() * L1.log() + PHI_C / L0.log())
    return v2 if v2 > v1 else v1          # a valid LOWER bound for max(v1, v2) in every case

def trap_consts(d, al):
    """sec 6.2: trapezoid kernels w_+/w_- (ramps of relative width d, left end alpha <= al):
    sup_{beta in [1/2,B]} x^{beta-B}|W| <= min(WB, A1/t, A2/t^2); left (beta in [b,1/2], x^{beta-1/2}): (Wh, A1L, A2L)"""
    WB = ((1 + d) ** (B + 1) - 1) / (d * B * (B + 1))
    A1 = (1 + d) ** B * (1 + al ** B)
    A2 = ((1 + d) ** (B + 1) + 1) * (1 + al ** B) / d
    Wh = ((1 + d) ** A(Fr(3, 2)) - 1) / (d * A(Fr(3, 4)))
    A1L = (1 + d) ** HALF * (1 + al ** HALF)
    A2L = ((1 + d) ** A(Fr(3, 2)) + 1) * (1 + al ** HALF) / d
    return (WB, A1, A2), (Wh, A1L, A2L)

def pnt_cell(kind, L0, L1, X0, X1, thetas=(Fr(1, 4), Fr(1), Fr(4), Fr(16), Fr(64))):
    """certified upper bound for |psi(x;q,a) - x/phi(q)| / (x^B log^2(qx)) over x in [e^X0, e^X1] and
    kind = 'q1' (q = 1), 'q2' (q = 2) or 'gen' (integers q in [e^L0, e^L1], q >= 3); requires e^X0 >= 64."""
    pos(X0 - arb(64).log(), "x >= 64")
    if kind == 'gen':
        phl = phi_lo(L0, L1); den = (L0 + X0) ** 2
    else:
        phl = arb(1); den = (X0 if kind == 'q1' else (LOG2 + X0)) ** 2
    xB0 = (B * X0).exp(); xb1 = (b * X1).exp(); xhB0 = ((HALF - B) * X0).exp()
    al = (-X0 / 2).exp()
    best = None
    for th in thetas:
        dfl = min(0.5, float((A(th) * phl / xb1).mid()))
        d = Fr(dfl).limit_denominator(10 ** 12)
        if d > Fr(1, 2): d = Fr(1, 2)
        if d <= 0: continue
        dd = A(d)
        pos(arb(1) / 2 - dd + 0 * dd + arb(0) if d < Fr(1, 2) else arb(1), "delta <= 1/2")
        (WB, A1, A2), (Wh, A1L, A2L) = trap_consts(dd, al)
        T1 = PSI_RS * xhB0
        T2 = dd * (xb1 + xhB0) / (2 * phl)
        Zz = zsum_zeta(WB, A1, A2); ZzL = zsum_zeta(Wh, A1L, A2L)
        T3 = (Zz + xhB0 * ZzL) / (2 * phl)
        lg = ((1 + dd) / (1 - dd)).log()
        T5 = (X1 / 2 + lg) / (1 - 1 / ((X0 / 2).exp() * (1 - dd)) ** 2) / xB0
        if kind == 'gen':
            Lq = L1 if L1 > arb(9).log() else arb(9).log()
            Z = zsum_dir(Lq, WB, A1, A2); ZL = zsum_dir(Lq, Wh, A1L, A2L)
            T4 = (Z + xhB0 * ZL) / 2
            T6 = ((L1 / LOG2) * (X1 / 2 + lg) + L1) / xB0
        elif kind == 'q2':
            T4 = arb(0); T6 = (X1 / 2 + lg + LOG2) / xB0
        else:
            T4 = arb(0); T6 = arb(0)
        R = (T1 + T2 + T3 + T4 + T5 + T6) / den
        if best is None or R.upper() < best.upper():
            best = R
    return best

def pnt_tail(XT, th=Fr(1)):
    """sec 6.4: certified bound for the ratio, valid for all X = log x >= XT and all 0 <= log q <= X"""
    XT = A(XT); bX = b * XT
    dT = A(th) * (-bX).exp()
    pos(A(Fr(1, 7)) - dT, "delta_T <= 1/7")
    N0 = (1 + dT) ** B / B          # W_B(delta) = ((1+d)^(B+1)-1)/(d B (B+1)) <= (1+d)^B/B (mean value theorem); no cancellation
    K1 = (1 + dT) ** B * (1 + (-B * XT / 2).exp())
    pos(1 / N0 - S0, "t_a >= 5/7")
    pos(2 / (B * dT) - K1 * K1, "three-piece envelope")
    tmax = K1 * B                                                 # t_a <= K1 B
    pos(-((tmax / (2 * PI)).log() / PI + CB / (tmax + 2)), "rho_0 < 0 on [5/7, t_a]")
    cth = ((2 + dT) / (2 * PI * A(th))).log()                     # lambda = log(t_b/2pi) <= bX + cth
    cp = amax(cth, arb(0))
    lamX = b + cp / XT                                            # lambda/X <= b + cth^+/X  (non-increasing)
    lam = lamX * XT
    pos(b * XT + ((2 / (2 * PI * A(th)))).log(), "lambda >= 0")
    a0 = S0 / PI + CB
    b0 = (S0 / PI) * (S0 / (2 * PI * E1)).log() + A(Fr(1, 4)) + CB * ((S0 + 2) / (2 * PI)).log() - HALF
    pos(-b0, "b0 < 0")
    al1 = a0 * N0 + (2 * K1 + K1 * ((2 * PI).log() + N0.log())) / PI - S0 / (PI * B)
    be0 = b0 / B + K1 * (CB / 2 * (1 + 2 * N0).log() + 1 / PI + CB / 2) + K1 * ((1 + 2 / UE).log() + dT / 2)
    T14 = arb(14)
    ga0 = K1 * (Nplus_z(T14) / 14 + CR1 / 14 + CR2 / (14 * T14.log()) + 1 / PI + CR1 / 2 / 14 + CR2 / (28 * T14.log()))
    epsX = ((HALF - B) * XT).exp()
    M = amax(amax(N0_left(dT) * B, (1 + dT) ** HALF * (1 + (-XT / 4).exp())), ((1 + dT) ** A(Fr(3, 2)) + 1) * (1 + (-XT / 4).exp()) / 2)
    fac = (1 + M * epsX) / 2
    L9 = arb(9).log()
    R = fac * ((al1 + K1 * lam / PI) / (4 * XT) + L9 * (al1 + K1 * lam / PI) / XT ** 2 + 2 * N0 * (1 + XT).log() / XT ** 2
               + (amax(be0, arb(0)) + ga0) / XT ** 2 + K1 * lamX ** 2 / PI + 2 * K1 * lamX / (PI * XT))
    xB = (B * XT).exp()
    R += (PSI_RS * epsX + A(th) * (1 + (-XT / 2).exp()) / 2) / XT ** 2
    R += (XT / 2 + ((1 + dT) / (1 - dT)).log()) / ((1 - 4 / XT.exp()) * xB * XT ** 2)
    R += ((XT / 2 + ((1 + dT) / (1 - dT)).log()) / LOG2 + 1) / (4 * XT * xB)
    return R

def N0_left(d):
    return 2 * (1 + d) ** HALF      # W_h(delta) = ((1+d)^(3/2)-1)/(3d/4) <= 2 (1+d)^(1/2) (mean value theorem)

def exact_small(c1, xmax=65):
    """1 <= q <= x <= xmax, (a,q) = 1: sup over x in [max(2,q), xmax] of |psi(x;q,a) - x/phi(q)|/(x^B log^2(qx)) (Arb)"""
    lam = [arb(0)] * (xmax + 2)
    for n in range(2, xmax + 2):
        m = n; p = None
        for d in range(2, n + 1):
            if m % d == 0:
                p = d; break
        while m % p == 0: m //= p
        if m == 1: lam[n] = arb(p).log()
    worst = arb(0); arg = None; worst64 = arb(0)
    for q in range(1, xmax + 1):
        ph = sum(1 for k in range(1, q + 1) if math.gcd(k, q) == 1)
        for a in range(1, q + 1):
            if math.gcd(a, q) != 1: continue
            psi = arb(0)
            for n in range(1, xmax + 1):
                if n % q == a % q: psi += lam[n]
                if n < max(2, q): continue
                pieces = 2000 if n == 2 else (200 if n < 8 else 4)
                for j in range(pieces):
                    s0 = arb(n) + arb(j) / pieces; s1 = arb(n) + arb(j + 1) / pieces
                    e = amax(abs(psi - s0 / ph), abs(psi - s1 / ph))
                    r = e / ((B * s0.log()).exp() * (arb(q) * s0).log() ** 2)
                    if r > worst: worst = r; arg = (q, a, n, j)
                    if n >= 64 and r > worst64: worst64 = r
    pos(A(c1) - worst, "exact small range")
    return worst, arg, worst64

def pnt_x_below_q(c1, x=None):
    """sec 6.1: x < q: |E| <= max(log x, sqrt(2x)) <= c1 x^B log^2(qx); certified at x = x0 (default 2), monotone beyond"""
    x = arb(2) if x is None else x
    v1 = 1 / (4 * (B * x.log()).exp() * x.log())                  # need log x <= 4 c1 x^B log^2 x
    v2 = arb(2).sqrt() * ((HALF - B) * x.log()).exp() / (4 * x.log() ** 2)
    pos(A(c1) - amax(v1, v2), "x < q case")
    return amax(v1, v2)

def pnt(c1=Fr(3, 2), XT=Fr(60)):
    w, arg, w64 = exact_small(c1)
    print("[pnt] exact, 2 <= q <= x <= 66 (x < 66): max ratio %s at (q,a,n,piece) = %s; on 64 <= x < 66: %s" % (w.str(6), arg, w64.str(4)), flush=True)
    v = pnt_x_below_q(c1)
    print("[pnt] x < q (any x >= 2): ratio <= %s" % v.str(4), flush=True)
    worst_by_X = []
    X0 = A(Fr(416, 100)); ncell = 0; worst = arb(0)          # e^4.16 = 64.07 < 66: exact range overlaps
    XTa = A(XT)
    while X0 < XTa:
        X1 = X0 + amax(A(Fr(1, 20)), X0 / 50)
        slab = arb(0)
        for kind in ('q1', 'q2'):
            R = pnt_cell(kind, None, None, X0, X1)
            pos(A(c1) - R, "pnt cell %s X0=%s" % (kind, X0.str(6)))
            slab = amax(slab, R); ncell += 1
        L0 = arb(3).log()
        while L0 < X1:
            L1 = L0 + amax(A(Fr(1, 10)), L0 / 20)
            R = pnt_cell('gen', L0, L1, X0, X1)
            pos(A(c1) - R, "pnt cell gen L0=%s X0=%s" % (L0.str(6), X0.str(6)))
            slab = amax(slab, R); ncell += 1
            L0 = L1
        worst_by_X.append((float(X0.mid()), float(X1.mid()), float(slab.upper())))
        worst = amax(worst, slab)
        X0 = X1
    Rt = pnt_tail(XT)
    pos(A(c1) - Rt, "pnt tail")
    print("[pnt] x >= 64, q <= x: %d cells to log x = %s, max ratio %s; tail log x >= %s: ratio <= %s"
          % (ncell, float(XTa.mid()), worst.str(4), float(XT), Rt.str(4)), flush=True)
    print("[pnt]   on 64 <= x < 66 also the exact range applies: %s" % w64.str(4))
    for X in [4.16, 10, 20, 30, 40, 50]:
        m = max([r for (x0, x1, r) in worst_by_X if x1 > X] + [float(Rt.upper())])
        vq = pnt_x_below_q(c1, A(Fr(int(X * 100), 100)).exp())
        print("[pnt]   sup ratio over log x >= %.2f: q <= x (cells and tail) %.4f; q > x %s" % (X, m, vq.str(3)))
    for XT2 in (Fr(100), Fr(300), Fr(1000), Fr(10 ** 4), Fr(10 ** 6)):
        print("[pnt]   tail log x >= %s (all q; q > x part negligible): ratio <= %s" % (float(XT2), pnt_tail(XT2).str(4)))
    return worst_by_X, Rt

# ---------------------------------------------------------------- sec 7: L'/L(1, chi)  (item iv)
DUS = A(Fr(2, 10))                   # Dusart (arXiv 1002.0442, Ramanujan J. 2018): sum_{p<=y} log p/p <= log y + E + 0.2/log y + 0.2/log^2 y, y >= 2974
def lp_kernels(h):
    """sec 7.1: one-sided ramp v (width h in n/x), V(s) = ((1+h)^(s+1) - 1)/(h s (s+1)).
    right zeros (beta in [1/2,B]): x^{beta-1}|V(rho-1)| <= x^{-b} min(1/b, 1/|g|, ((1+h)^B+1)/(h g^2));
    left zeros (beta in [b,1/2]):  x^{beta-1}|V(rho-1)| <= x^{-1/2} min(2, 1/|g|, ((1+h)^(1/2)+1)/(h g^2))"""
    hh = A(h)
    return (1 / b, arb(1), ((1 + hh) ** B + 1) / hh), (arb(2), arb(1), ((1 + hh) ** HALF + 1) / hh)

def lp_value(L0, L1, kappa, h):
    """certified upper bound for |L'/L(1,chi)|, chi primitive mod q, e^L0 <= q <= e^L1, using x = max(2974, (kappa L0)^(1/b))"""
    X = (A(kappa) * L0).log() / b
    if X < arb(2974).log():
        X = arb(2974).log()
    x = X.exp(); hh = A(h)
    G, GL = lp_kernels(h)
    Lq = L1 if L1 > arb(9).log() else arb(9).log()
    Z = zsum_dir(Lq, *G); ZL = zsum_dir(Lq, *GL)
    prime = X + (1 + hh).log() - EUL + DUS / X + DUS / X ** 2
    return prime + ((-b * X).exp() * Z + (-X / 2).exp() * ZL) / 2 + 1 / (x - 1)

def lp_tail(Lm, kappa, h):
    """all q >= e^Lm with x = (kappa log q)^(1/b): |L'/L(1,chi)| - (1/b) log log q <= returned constant"""
    Lm = A(Lm)
    X = (A(kappa) * Lm).log() / b
    pos(X - arb(2974).log(), "x >= 2974 in the tail")
    x = X.exp(); hh = A(h)
    G, GL = lp_kernels(h)
    pos(Lm - arb(9).log(), "log q >= log 9")
    co = zsum_coeffs(*G); coL = zsum_coeffs(*GL)
    mG = Z_over_L_majorant(co, G[0], Lm); mL = Z_over_L_majorant(coL, GL[0], Lm)
    pos(1 / (2 * b) - 1, "(kappa L)^(-1/(2b)) L decreasing")
    return INV_b * A(kappa).log() + (1 + hh).log() - EUL + DUS / X + DUS / X ** 2 + 1 / (x - 1) \
        + mG / (2 * A(kappa)) + (-X / 2).exp() * Lm * mL / 2

def lprime():
    """sec 7: |L'/L(1,chi)| <= (1/(1-B)) log log q + C for primitive chi mod q"""
    grid = [(Fr(k, 100), Fr(hh, 10)) for k in (20, 25, 30, 40, 60, 100, 200) for hh in (3, 5, 7, 10)]
    rows = []
    for Lm in (Fr(25), Fr(100), Fr(1000), Fr(10000), Fr(10 ** 6)):
        best = None
        for (kap, h) in grid:
            try:
                c = lp_tail(Lm, kap, h)
            except SystemExit:
                continue
            if best is None or c.upper() < best[0].upper():
                best = (c, kap, h)
        rows.append((Lm, best))
        print("[lprime] log q >= %s: |L'/L(1,chi)| <= (1/(1-B)) log log q + %s  (kappa = %s, h = %s)"
              % (float(Lm), best[0].str(5), best[1], best[2]), flush=True)
    # cells from log 3 to 25 with the per-cell best grid point: constant C' for all q >= 3
    Ctail = rows[0][1][0]
    L0 = arb(3).log(); worstC = Ctail; n = 0
    while L0 < A(25):
        L1 = L0 + amax(A(Fr(1, 20)), L0 / 100)
        bestv = None
        for (kap, h) in grid:
            v = lp_value(L0, L1, kap, h) - INV_b * L0.log()
            if bestv is None or v.upper() < bestv.upper():
                bestv = v
        worstC = amax(worstC, bestv); n += 1
        L0 = L1
    print("[lprime] all q >= 3 (%d cells to log q = 25, then the tail): |L'/L(1,chi)| <= (1/(1-B)) log log q + %s"
          % (n, worstC.str(5)), flush=True)
    return rows, worstC

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if mode in ("core", "all"):
        core()
    if mode in ("linnik", "all"):
        linnik()
    if mode in ("nonres", "all"):
        nonres()
    if mode in ("pnt", "all"):
        pnt()
    if mode in ("lprime", "all"):
        lprime()
    if mode == "probe":
        for L in [9.2103, 23.026, 50, 230.3, 1000]:
            for c2 in [Fr(73, 100), Fr(62, 100), Fr(55, 100), Fr(49, 100)]:
                mg, Z = linnik_cell(KER, c2, A(L), A(L) + A(Fr(1, 100)))
                print(L, c2, mg.str(5), Z.str(6))
    print("runtime %.1f s" % (time.time() - t0))
