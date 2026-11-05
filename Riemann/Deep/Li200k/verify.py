#!/usr/bin/env python3
"""
verify.py -- certificate for the explicit semi-local Li criterion (node N3).

Claim: lambda_n > 0 for all 1 <= n <= N(a), where a bounds |Re(rho) - 1/2| for the zeros above
T0 = 3,000,175,332,800 (Platt-Trudgian: RH holds for 0 < gamma <= T0):
  variant A: a = 1/2                                (no input beyond 0 < beta < 1)
  variant B: a(t) = 1/2 - 1/(5.558691 log t)          (Mossinghoff-Trudgian-Yang zero-free region)
  variant C: a = 3/8                                (zero-free half-plane Re s > 7/8)
Proof: REPORT.md.  Key inequality  lambda_n >= M_n - Q_n  with
  M_n = sum_{gamma>0} 2(1 - cos(n theta(gamma))),   theta(t) = 2 arctan(1/(2t))   (ordinates only)
  Q_n = sum_{gamma>=T0} q_n(gamma),  q_n(t) = 4 sinh^2(n a(t)/(2t^2)) + 4 n a(t)^2/t^3.
Inputs used numerically: Trudgian 2014 |N(t) - t/(2pi) log(t/(2pi e)) - 7/8| <= E(t),
  E(t) = 0.111 log t + 0.275 log log t + 2.450 + 0.2/t  (t >= e).
All arithmetic is ball arithmetic (python-flint / Arb); every comparison is a certified sign.
Needs: pip install python-flint (tested 0.9.0, FLINT 3.6).  Runtime about 1 min, peak memory about 40 MB, one thread.
Option --conditional adds two conditional half-plane variants (D, E).
"""
import sys, time
from flint import arb, acb, ctx

ctx.prec = 160
t_start = time.time()
PI = arb.pi()
T0 = arb(3000175332800)
LOG2PI = (2 * PI).log()

def r(p, q=1):
    return arb(p) / q

C111, C275, C245, C02 = r(111, 1000), r(275, 1000), r(245, 100), r(2, 10)

def E(t):            # Trudgian 2014 bound on |N(t) - Ntilde(t)|, t >= e
    return C111 * t.log() + C275 * t.log().log() + C245 + C02 / t

def Ntilde(t):
    return t / (2 * PI) * (t / (2 * PI)).log() - t / (2 * PI) + r(7, 8)

def certified_pos(x, what):
    if not (x > 0):
        raise SystemExit("FAIL: could not certify %s > 0 (value %s)" % (what, x.str(10, radius=True)))

# ---------------------------------------------------------------- constants K0, K1
# K0 = int_0^pi (1-cos u)/u^2 du ,  K1 = int_0^pi (1-cos u) log(u)/u^2 du   (power series, rigorous tail)
def K_constants(terms=80):
    K0 = arb(0); K1 = arb(0); lp = PI.log()
    fact = arb(1)
    for k in range(1, terms + 1):
        fact = fact * (2 * k - 1) * (2 * k)            # (2k)!
        s = (-1) ** (k + 1)
        pw = PI ** (2 * k - 1)
        K0 += s * pw / (fact * (2 * k - 1))
        K1 += s * pw / fact * (lp / (2 * k - 1) - arb(1) / (2 * k - 1) ** 2)
    # tail: sum_{k>terms} pi^(2k-1)/(2k)! * (1 + log pi + 1)  <=  2 * first omitted term bound
    K = terms + 1
    fK = fact * (2 * K - 1) * (2 * K)
    tail = 2 * PI ** (2 * K - 1) / fK * (2 + lp)
    K0 = K0 + arb(0, tail.upper()); K1 = K1 + arb(0, tail.upper())
    return K0, K1

K0, K1 = K_constants()
# cross-check K0 = Si(pi) - 2/pi
K0_check = PI.si() - 2 / PI
assert (K0 - K0_check).contains(0) or abs(float(K0 - K0_check)) < 1e-30, "K0 cross-check failed"

# ---------------------------------------------------------------- lower bound for M_n
def t_of(n):          # t_n = cot(pi/(2n))/2 : n*theta(t_n) = pi
    return (PI / (2 * arb(n))).cot() / 2

def main_term(n):     # (n/pi) [K0 (log(n/2pi) - pi^2/(4n^2)) - K1], increasing in n
    n = arb(n)
    return n / PI * (K0 * ((n / (2 * PI)).log() - PI ** 2 / (4 * n * n)) - K1)

def M_lo_interval(n1, n2):
    """certified lower bound for M_m valid for every integer m in [n1, n2] (n1 >= 89)."""
    t1, t2 = t_of(n1), t_of(n2)
    corr = (C111 + C275 / t1.log()) * arb(n1) ** 2 / (2 * t1 * t1)   # m^2/t_m^2 decreasing in m
    return main_term(n1) - 8 * E(t2) - corr

K1INF = PI / 2 * (1 - arb.const_euler())      # int_0^inf (1-cos u) log(u)/u^2 du (Mellin transform)

def M_lo_full(n1, n2):
    """lower bound for M_m, all integers m in [n1, n2] (n1 >= 10^4): uses every zero with gamma >= t_low = sqrt(n1),
    cut into pieces on which G = 2-2cos(m theta(t)) is monotone; each piece costs <= 8 E (Stieltjes with |N - Ntilde| <= E)."""
    tl = arb(n1).sqrt()
    th = 2 * (1 / (2 * tl)).atan()
    U1, U2 = arb(n1) * th, arb(n2) * th
    K0U = PI / 2 - 2 / U1                                # int_0^U (1-cos u)/u^2 du >= pi/2 - 2/U
    n = arb(n1)
    main = n / PI * (K0U * ((n / (2 * PI)).log() - 1 / (4 * tl * tl)) - K1INF)
    t1, t2 = t_of(n1), t_of(n2)
    pieces = U2 / PI + 3
    corr = (C111 + C275 / t1.log()) * arb(n1) ** 2 / (2 * t1 * t1)
    return main - 8 * E(t2) * pieces - corr

def M_lo(n1, n2):
    a = M_lo_interval(n1, n2)
    if n1 >= 10 ** 4:
        b = M_lo_full(n1, n2)
        if b > a:
            return b
    return a

# ---------------------------------------------------------------- upper bound for Q_n
A_HALF = r(1, 2); A_38 = r(3, 8); R0_MTY = r(5558691, 1000000)
def a_const(v):
    return (lambda t: v), v
def a_classical():
    return (lambda t: r(1, 2) - 1 / (R0_MTY * t.log())), r(1, 2)
VARIANTS = {
    'A (a=1/2)': a_const(A_HALF),
    'B (a=1/2-1/(5.558691 log t))': a_classical(),
    'C (a=3/8, 7/8 half-plane)': a_const(A_38),
}
if '--conditional' in sys.argv:      # side results: conditional zero-free half-planes (see REPORT.md)
    VARIANTS['D (a=0.374957, Liu, conditional)'] = a_const(r(374957, 10 ** 6))
    VARIANTS['E (a=0.37446205769353, Chang-Wang, cond.)'] = a_const(r(37446205769353, 10 ** 14))

def W(t):              # antiderivative of w(t) = (1/2pi) log(t/2pi) + 0.13/t  (>= Ntilde' + E' for t >= T0)
    return t / (2 * PI) * (t / (2 * PI)).log() - t / (2 * PI) + r(13, 100) * t.log()

def tail_bound(n, amax, tK):
    """int_{tK}^inf q(t) w(t) dt with q <= C n^2 amax^2/t^4 + 4 n amax^2/t^3, C=(sinh xK/xK)^2."""
    n = arb(n)
    xK = n * amax / (2 * tK * tK)
    C = (xK.sinh() / xK) ** 2 if xK > 0 else arb(1)
    L = (tK / (2 * PI)).log()
    i4 = (3 * L + 1) / (9 * tK ** 3) / (2 * PI) + r(13, 100) / (4 * tK ** 4)
    i3 = (2 * L + 1) / (4 * tK ** 2) / (2 * PI) + r(13, 100) / (3 * tK ** 3)
    return C * n * n * amax * amax * i4 + 4 * n * amax * amax * i3

def q(n, a, t):
    n = arb(n); at = a(t)
    x = n * at / (2 * t * t)
    return 4 * x.sinh() ** 2 + 4 * n * at * at / t ** 3

def Q_hi(n, var):
    a, amax = VARIANTS[var]
    X = arb(n) * amax / (2 * T0 * T0)
    bound = 2 * q(n, a, T0) * E(T0)
    if X < r(1, 1000):
        return bound + tail_bound(n, amax, T0)
    # partition [T0, 100 T0]: fine steps in x = X T0^2/t^2 (dx <= 0.004) while x > 1, then geometric 1%
    Xf = float(X)
    pts = [T0]
    xs = Xf
    while xs > 1.0:
        xs = max(xs - 0.004, 1.0)
        pts.append(T0 * arb(Xf / xs).sqrt())
    t = pts[-1]
    while t < 100 * T0:
        t = t * r(101, 100)
        pts.append(t)
    # partition points are balls; every term encloses its value for each choice of points inside the balls
    # (e.g. the increasing midpoints), so the sum still majorises the integral for an exact partition
    s = arb(0)
    Wprev = W(pts[0])
    for k in range(len(pts) - 1):
        Wnext = W(pts[k + 1])
        s += q(n, a, pts[k]) * (Wnext - Wprev)        # q decreasing, w increasing
        Wprev = Wnext
    return bound + s + tail_bound(n, amax, pts[-1])

# ---------------------------------------------------------------- part 1: small n via the zero near 14.1347
def hardy_Z(t):
    t = arb(t)
    s = acb(r(1, 2), t)
    theta = (acb(r(1, 4), t / 2).lgamma()).imag - t / 2 * PI.log()
    return ((acb(0, theta)).exp() * s.zeta()).real

def part_small():
    za, zb = hardy_Z(r(1413, 100)), hardy_Z(r(1414, 100))
    certified_pos(-(za * zb), "sign change of Z on [14.13,14.14]")
    th_lo = 2 * (1 / (2 * r(1414, 100))).atan()    # theta(gamma*) in [th_lo, th_hi]
    th_hi = 2 * (1 / (2 * r(1413, 100))).atan()
    Qmax = max((Q_hi(88, v) for v in VARIANTS), key=lambda v: float(v))
    worst = None
    for n in range(1, 89):
        lo, hi = n * th_lo, n * th_hi          # n*theta(gamma*) lies in [lo, hi] (0 < lo, hi < 2 pi)
        certified_pos(lo, "n theta > 0")
        certified_pos(2 * PI - hi, "n theta < 2 pi")
        for g in (2 - 2 * lo.cos(), 2 - 2 * hi.cos()):   # g = 2-2cos is unimodal on (0, 2pi): min at an end
            certified_pos(g - Qmax, "g(n theta(gamma*)) - Q_88 (n=%d)" % n)
            worst = g if worst is None or float(g) < float(worst) else worst
    print("part 1: zero gamma* in (14.13,14.14) on the line (Z sign change %s, %s);"
          % (za.str(5), zb.str(5)))
    print("        for 1 <= n <= 88: lambda_n >= %s - %s > 0" % (worst.str(6), Qmax.str(3)))

# ---------------------------------------------------------------- part 0: the per-zero inequality (random tests)
def part0(trials=3000, seed=1):
    """Off-line quartet rho, rho* = 1-conj(rho), conj(rho), conj(rho*) with rho = 1/2 + d + i y (0 < d <= a): its total
    contribution to lambda_n is 4 - 2(r^n + r^-n) cos(n phi) = 2 * lhs below.  BOTH rho and rho* have height y > 0, so
    each of them carries one term g(n theta(y)) in M_n and one term q_n(y) in Q_n (N(t) counts both).  Check, per zero,
    lhs = 2 - Re(z_rho^n + z_rho*^n) >= g(n theta(y)) - 4 sinh^2(n a/(2y^2)) - 4 n a^2/y^3."""
    import random
    rnd = random.Random(seed)
    worst = None
    for _ in range(trials):
        a = arb(rnd.choice([1, 2, 3, 4])) / 8
        d = a * arb(rnd.random())
        y = arb(10) ** arb(rnd.uniform(0.5, 4))
        n = int(float(y * y) * rnd.uniform(0.01, 120)) + 1
        rho = acb(r(1, 2) + d, y)
        z = 1 - 1 / rho
        zs = 1 - 1 / (1 - rho.conjugate())
        lhs = (2 - z ** n - zs ** n).real                     # per-zero share of the 4-zero group, (4 - 2Re(...))/2
        th = 2 * (1 / (2 * y)).atan()
        rhs = 2 - 2 * (n * th).cos() - 4 * (n * a / (2 * y * y)).sinh() ** 2 - 4 * n * a * a / y ** 3
        certified_pos(lhs - rhs + arb(10) ** -30 * (1 + abs(rhs)), "per-zero inequality")
        sl = float(lhs - rhs)
        worst = sl if worst is None or sl < worst else worst
    print("part 0: per-zero inequality holds in %d random cases (min slack %.3g)" % (trials, worst))

# ---------------------------------------------------------------- part 3: consistency with directly computed lambda_n
def part3(NMAX=600):
    """lambda_n for n <= NMAX by Johansson's formula (1.1n+80 bits): checks lambda_1..3 and lower bound <= lambda_n."""
    from flint import arb_series
    old = ctx.prec
    ctx.prec = int(1.1 * NMAX) + 80; ctx.cap = NMAX + 1
    x = arb_series([0, 1]); one = arb_series([1])
    f = (-arb_series.zeta(x)).log() + (one + x / 2).lgamma() + (one - x).log() - x * arb.pi().log() / 2
    co = f.coeffs()
    lam = [None]
    for n in range(1, NMAX + 1):
        s = arb(0); c = arb(1)
        for k in range(1, n + 1):
            if k > 1:
                c = c * (n - k + 1) / (k - 1)
            s += (-1) ** k * c * co[k]
        lam.append(n * s)
    ctx.prec = old
    assert abs(float(lam[1]) - 0.0230957089661210338) < 1e-15
    mx = 0.0
    for n in range(89, NMAX + 1):
        lb = M_lo_interval(n, n) - Q_hi(n, 'A (a=1/2)')
        certified_pos(lam[n] - lb, "lambda_n - lower bound (n=%d)" % n)
        mx = max(mx, float(lb / lam[n]))
    print("part 3: lambda_1 = %s; for 89 <= n <= %d our lower bound <= %.3f lambda_n; min lambda_n (n<=%d) = lambda_1"
          % (lam[1].str(12), NMAX, mx, NMAX))

# ---------------------------------------------------------------- part 2: certification loop
RATIOS = [2.0, 1.5, 1.2, 1.1, 1.05, 1.02, 1.01, 1.005, 1.002, 1.001, 1.0005, 1.0002, 1.0001,
          1.00005, 1.00002, 1.00001, 1.000005, 1.000002, 1.000001, 1.0000005, 1.0000002, 1.0000001]

def certify(var):
    n = 89; steps = 0
    while True:
        advanced = False
        for rr in RATIOS:
            m = int(n * rr)
            if m <= n:
                continue
            gap = M_lo(n, m) - Q_hi(m, var)
            if gap > 0:
                n = m; steps += 1; advanced = True
                break
        if not advanced:
            return n, steps

def asym(n):
    n = arb(n)
    return n / 2 * (n.log() - 1 - LOG2PI + arb.const_euler())

if __name__ == "__main__":
    print("inputs: T0 = 3000175332800, E(T0) = %s, K0 = %s, K1 = %s"
          % (E(T0).str(8), K0.str(15, radius=False), K1.str(15, radius=False)))
    part0()
    part_small()
    part3()
    # consistency: our lower bound must lie below the known value lambda_100000 (Johansson 2015, Li normalisation)
    lam_1e5 = arb("462580.782406902231409416038")
    lb = M_lo(100000, 100000) - Q_hi(100000, 'A (a=1/2)')
    certified_pos(lam_1e5 - lb, "known lambda_100000 above our lower bound")
    print("check : at n = 1e5 our lower bound %s <= lambda_n = 462580.78 (Johansson)" % lb.str(8))
    T0sq = T0 * T0
    for var in VARIANTS:
        t1 = time.time()
        N, steps = certify(var)
        print("part 2: %-30s lambda_n > 0 for all n <= N = %d  (N/T0^2 = %.6f, %d steps, %.1f s)"
              % (var, N, float(arb(N) / T0sq), steps, time.time() - t1))
        print("        at N: M_lo = %s, Q_hi = %s, asymptotic (n/2)(log n-1-log 2pi+gamma) = %s"
              % (M_lo(N, N).str(6), Q_hi(N, var).str(6), asym(N).str(6)))
    print("total runtime %.1f s" % (time.time() - t_start))
