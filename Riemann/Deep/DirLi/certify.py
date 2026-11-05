#!/usr/bin/env python3
"""
certify.py -- certificate for explicit Li-coefficient positivity ranges of Dirichlet L-functions (node N3).

Claim: for every primitive Dirichlet character chi mod q (3 <= q <= 400,000) and 1 <= n <= N(q),
Re lambda_n(chi) > 0, where lambda_n(chi) = sum_rho [1 - (1 - 1/rho)^n] over the nontrivial zeros of L(s, chi).
Variants: a = 1/2 (zeros above t0 only known to lie in 0 < beta < 1) and a = 3/8 (7/8 half-plane for every
Dirichlet L-function, OpenAI 2026, Lean port in Riemann/Deep/QRH; with the functional equation 1/8 <= beta <= 7/8).
Inputs: Platt (Math. Comp. 85 (2016), arXiv 1305.3087): GRH for |Im rho| <= t0(q);
        Bennett-Martin-O'Bryant-Rechnitzer (Math. Comp. 90 (2021), arXiv 2005.02989, Thm 1.1): zero counting.
Proof: PROOF.md.  Key inequality  Re lambda_n >= M_n/2 - Q_n/2,
  M_n = sum_rho G_n(|gamma|), G_n(t) = 2 - 2cos(n theta(t)), theta(t) = 2 atan(1/(2t))   (all zeros; ordinates only)
  Q_n = sum_{|gamma| > t0} q_n(|gamma|),  q_n(t) = 4 sinh^2(n a/(2t^2)) + 4 n a^2/t^3.
Ball arithmetic (python-flint / Arb); every comparison used as a proof step is a certified sign; floats are used only
to choose block end points (guidance), never as proof.
Run:  uv run --no-project --with python-flint --with numpy python3 certify.py <mode>   (one thread, ~40 MB)
  modes: part0 (pair inequality, random tests) | small (q = 3,4,5,7,8) | q1000 [lo hi] (each modulus) |
         blocks [ratio qmin qmax] (blocks of moduli, odd and 0 mod 4 separately)
Outputs of the runs: results_q1000.txt (q1000 3 1000), results_blocks.txt (blocks 1.02 1000 400000).
"""
import sys, time, math
from fractions import Fraction
import numpy as np
from flint import arb, ctx

ctx.prec = 128
PI = arb.pi()
EG = arb.const_euler()
E1 = arb(1).exp()

def R(p, q=1):
    return arb(p) / q

def A(fr):                       # exact Fraction -> arb
    return arb(fr.numerator) / fr.denominator

# ---------------------------------------------------------------- inputs
def t0_of(q):
    """Platt's height (exact rational): GRH for all primitive chi mod q <= 400,000 and |Im rho| <= t0(q)."""
    assert 3 <= q <= 400000
    if q % 2 == 0:
        return max(Fraction(10 ** 8, q), Fraction(75 * 10 ** 6, q) + 200)
    return max(Fraction(10 ** 8, q), Fraction(375 * 10 ** 5, q) + 200)

CA = R(22737, 100000)            # BMOR Thm 1.1: |N - Nt| <= CA*l + 2 log(1+l) - 1/2, l = log(q(T+2)/2pi) > 1.567
LTH = R(1567, 1000)
T57 = R(5, 7)

def ell(t, q):
    return (arb(q) * (t + 2) / (2 * PI)).log()

def E(t, q):
    l = ell(t, q)
    return CA * l + 2 * (1 + l).log() - R(1, 2)

def cE(t, q):                    # E'(s) <= cE(t,q)/s for s >= t
    return CA + 2 / (1 + ell(t, q))

def Nt0(t, q):                   # Nt without the parity constant -chi(-1)/4
    return t / PI * (arb(q) * t / (2 * PI * E1)).log()

def theta(t):
    return 2 * (1 / (2 * t)).atan()

def t_of(n):                     # n theta(t_n) = pi
    return (PI / (2 * arb(n))).cot() / 2

def G(n, t):
    return 2 - 2 * (arb(n) * theta(t)).cos()

# ---------------------------------------------------------------- constants
def K_constants(terms=80):
    """K0 = int_0^pi (1-cos u)/u^2 du, K1 = int_0^pi (1-cos u) log(u)/u^2 du (power series, rigorous tail)."""
    K0 = arb(0); K1 = arb(0); lp = PI.log(); fact = arb(1)
    for k in range(1, terms + 1):
        fact = fact * (2 * k - 1) * (2 * k)
        s = (-1) ** (k + 1)
        pw = PI ** (2 * k - 1)
        K0 += s * pw / (fact * (2 * k - 1))
        K1 += s * pw / fact * (lp / (2 * k - 1) - arb(1) / (2 * k - 1) ** 2)
    K = terms + 1
    tail = 2 * PI ** (2 * K - 1) / (fact * (2 * K - 1) * (2 * K)) * (2 + lp)
    return K0 + arb(0, tail.upper()), K1 + arb(0, tail.upper())

K0, K1 = K_constants()
assert abs(float(K0 - (PI.si() - 2 / PI))) < 1e-30
K1INF = PI / 2 * (1 - EG)        # int_0^oo (1-cos u) log(u)/u^2 du

# ---------------------------------------------------------------- parameter set (a q-block)
class P:
    """qlo <= q <= qhi; t0 <= t0(q) for every q of the block; a = strip half-width for zeros above t0."""
    def __init__(self, qlo, qhi, t0, a, label=''):
        self.qlo, self.qhi, self.t0f, self.a, self.label = qlo, qhi, t0, a, label
        self.t0 = A(t0)
        self.af = float(a)
        self.aa = A(a)

# ---------------------------------------------------------------- gain: lower bounds for M_m/2, all m in [n1, n2]
def gainA(n1, n2, p):
    """zeros with |gamma| > t_m (G_m decreasing 4 -> 0 there)."""
    if n1 < 3:
        return None
    t1, t2 = t_of(n1), t_of(n2)
    if not (t1 > T57 and ell(t1, p.qlo) > LTH and arb(p.qlo) * t1 > 2 * PI):
        return None
    n = arb(n1)
    br = K0 * ((arb(p.qlo) * n / (2 * PI)).log() - PI ** 2 / (4 * n * n)) - K1
    if not (br > 0):
        return None
    return n / PI * br - 4 * E(t2, p.qhi) - cE(t1, p.qlo) * n * n / (4 * t1 * t1)

def gainB(n1, n2, p, kap=1):
    """zeros with |gamma| > t_low = kap sqrt(n1), cut into pieces on which G_m is monotone."""
    tl = arb(n1).sqrt() * kap
    t1, t2 = t_of(n1), t_of(n2)
    if not (tl < t1 and tl > 1 and ell(tl, p.qlo) > LTH and arb(p.qlo) * tl > 2 * PI):
        return None
    th = theta(tl)
    U1, U2 = arb(n1) * th, arb(n2) * th
    if not (U1 > 4 / PI):
        return None
    n = arb(n1)
    br = (PI / 2 - 2 / U1) * ((arb(p.qlo) * n / (2 * PI)).log() - 1 / (4 * tl * tl)) - K1INF
    if not (br > 0):
        return None
    return n / PI * br - 4 * E(t2, p.qhi) * (U2 / PI + 2) - cE(t1, p.qlo) * n * n / (4 * t1 * t1)

def gain(n1, n2, p):
    best = None
    for g in (gainA(n1, n2, p), gainB(n1, n2, p) if n1 >= 10 ** 4 else None):
        if g is not None and (best is None or g > best):
            best = g
    return best

def ceil_lo(x):                  # N >= x, N integer  =>  N >= k
    lo = x.lower(); k = math.ceil(float(lo))
    while not (arb(k - 1) < lo):
        k -= 1
    return k

def floor_hi(x):                 # N <= x, N integer  =>  N <= k
    hi = x.upper(); k = math.floor(float(hi))
    while not (arb(k + 1) > hi):
        k += 1
    return k

def Nlo(s, p):                   # certified integer lower bound for N(s), every q of the block, both parities
    if s > T57 and ell(s, p.qlo) > LTH:
        return max(0, ceil_lo(Nt0(s, p.qlo) - R(1, 4) - E(s, p.qhi)))
    return 0

def Nhi(s, p):                   # certified integer upper bound for N(s), s >= 5/7
    if ell(s, p.qhi) <= LTH:
        return 0
    return max(0, floor_hi(Nt0(s, p.qhi) + R(1, 4) + E(s, p.qhi)))

def gainC(n, p, ratio=R(103, 100), span=2000):
    """M_n/2 >= (1/2) int_{s0}^oo (N(s) - N(s0)) |G_n'(s)| ds, s0 = max(t_n, 5/7): lower Stieltjes sum."""
    if n <= 2:
        s0 = T57                                   # t_1 = 0, t_2 = 1/2 < 5/7
    else:
        s0 = t_of(n)
        assert s0 > T57
    N0 = Nhi(s0, p)
    S = s0 * span
    tot = arb(0); s = s0; Gs = G(n, s)
    while s < S:
        s2 = s * ratio; G2 = G(n, s2)
        c = Nlo(s, p) - N0
        if c > 0:
            tot += c * (Gs - G2)                   # G decreasing on [t_n, oo)
        s, Gs = s2, G2
    c = Nlo(s, p) - N0
    if c > 0:
        tot += c * Gs
    return tot / 2

# ---------------------------------------------------------------- loss: upper bound for Q_n/2, every q of the block
def qn(n, a, t):
    x = arb(n) * a / (2 * t * t)
    return 4 * x.sinh() ** 2 + 4 * arb(n) * a * a / t ** 3

def tail(n, a, tK, qh, c):
    """int_{tK}^oo q_n w dt, w = log(qh t/2pi)/pi + c/t; q_n <= C n^2 a^2/t^4 + 4 n a^2/t^3, C = (sinh xK/xK)^2."""
    n = arb(n)
    xK = n * a / (2 * tK * tK)
    C = (xK.sinh() / xK) ** 2
    L = (arb(qh) * tK / (2 * PI)).log()
    i4 = (3 * L + 1) / (9 * tK ** 3) / PI + c / (4 * tK ** 4)
    i3 = (2 * L + 1) / (4 * tK ** 2) / PI + c / (3 * tK ** 3)
    return C * n * n * a * a * i4 + 4 * n * a * a * i3

def loss(n, p, dx=0.004, width=6.0):
    """Q_n/2 <= (1/2)[int_{t0}^oo q_n (Nt' + E') dt + 2 q_n(t0) E(t0)]  (Stieltjes; q_n decreasing)."""
    a, t0, qh = p.aa, p.t0, p.qhi
    c = cE(t0, p.qlo)
    bound = 2 * qn(n, a, t0) * E(t0, qh)
    X = arb(n) * a / (2 * t0 * t0)
    if X < R(1, 1000):
        return (bound + tail(n, a, t0, qh, c)) / 2
    def W(t):                                  # antiderivative of w = log(qh t/2pi)/pi + c/t (w > 0)
        return t / PI * ((arb(qh) * t / (2 * PI)).log() - 1) + c * t.log()
    Xf = float(X); xs = Xf; xmin = max(Xf - width, 1e-3)
    pts = [t0]
    while xs > xmin:                           # fine steps in x = X t0^2/t^2 near t0
        xs = max(xs - dx, xmin)
        pts.append(t0 * arb(Xf / xs).sqrt())   # balls; every term encloses its value for any points inside them
    s = arb(0); Wp = W(pts[0])
    for k in range(len(pts) - 1):
        Wn = W(pts[k + 1])
        s += qn(n, a, pts[k]) * (Wn - Wp)      # upper sum: q_n decreasing, w > 0
        Wp = Wn
    return (bound + s + tail(n, a, pts[-1], qh, c)) / 2

# ---------------------------------------------------------------- float mirrors (guidance only)
def f_ell(t, q): return math.log(q * (t + 2) / (2 * math.pi))
def f_E(t, q):
    l = f_ell(t, q); return 0.22737 * l + 2 * math.log1p(l) - 0.5
def f_cE(t, q): return 0.22737 + 2 / (1 + f_ell(t, q))
def f_t(n): return 0.5 / math.tan(math.pi / (2 * n))
F_K0, F_K1 = float(K0), float(K1)
F_K1INF = float(K1INF)

def f_gain(n1, n2, p):
    best = -1e300
    t1, t2 = f_t(n1), f_t(n2)
    if n1 >= 3:
        g = n1 / math.pi * (F_K0 * (math.log(p.qlo * n1 / (2 * math.pi)) - math.pi ** 2 / (4 * n1 * n1)) - F_K1) \
            - 4 * f_E(t2, p.qhi) - f_cE(t1, p.qlo) * n1 * n1 / (4 * t1 * t1)
        best = max(best, g)
    if n1 >= 10 ** 4:
        tl = math.sqrt(n1); th = 2 * math.atan(1 / (2 * tl)); U1, U2 = n1 * th, n2 * th
        g = n1 / math.pi * ((math.pi / 2 - 2 / U1) * (math.log(p.qlo * n1 / (2 * math.pi)) - 1 / (4 * tl * tl)) - F_K1INF) \
            - 4 * f_E(t2, p.qhi) * (U2 / math.pi + 2) - f_cE(t1, p.qlo) * n1 * n1 / (4 * t1 * t1)
        best = max(best, g)
    return best

def f_qn(n, a, t):
    x = n * a / (2 * t * t)
    if x > 340: return float('inf')
    sh = math.sinh(x)
    return 4 * sh * sh + 4 * n * a * a / t ** 3

def f_loss(n, p, dx=0.004, width=6.0):
    a, t0, qh = p.af, float(p.t0f), p.qhi
    n = float(n)
    c = f_cE(t0, p.qlo)
    bound = 2 * f_qn(n, a, t0) * f_E(t0, qh)
    X = n * a / (2 * t0 * t0)
    def tl(tK):
        xK = n * a / (2 * tK * tK)
        if xK > 340: return float('inf')
        C = (math.sinh(xK) / xK) ** 2
        L = math.log(qh * tK / (2 * math.pi))
        return C * n * n * a * a * ((3 * L + 1) / (9 * tK ** 3) / math.pi + c / (4 * tK ** 4)) \
            + 4 * n * a * a * ((2 * L + 1) / (4 * tK ** 2) / math.pi + c / (3 * tK ** 3))
    if X < 1e-3:
        return (bound + tl(t0)) / 2
    if X > 340:
        return float('inf')
    xmin = max(X - width, 1e-3)
    k = np.arange(0, int(math.ceil((X - xmin) / dx)) + 1, dtype=float)
    xs = np.maximum(X - dx * k, xmin)
    ts = t0 * np.sqrt(X / xs)
    Ws = ts / math.pi * (np.log(qh * ts / (2 * math.pi)) - 1) + c * np.log(ts)
    x = n * a / (2 * ts[:-1] ** 2)
    qv = 4 * np.sinh(x) ** 2 + 4 * n * a * a / ts[:-1] ** 3
    s = float(np.sum(qv * np.diff(Ws)))
    return (bound + s + tl(float(ts[-1]))) / 2

# ---------------------------------------------------------------- certification
def small_n(p, nmax=400):
    """regime C for n = 1, 2, ... until gainA certifies a single n; returns (last n done, worst margin)."""
    n = 0; worst = None
    while True:
        m = n + 1
        gA = gainA(m, m, p)
        if gA is not None and gA - loss(m, p) > 0 and m > 3:
            return n, worst
        d = gainC(m, p) - loss(m, p)
        if not (d > 0):
            raise SystemExit("FAIL regime C n=%d %s: %s" % (m, p.label, d.str(5)))
        worst = d if worst is None or float(d) < float(worst) else worst
        n = m
        if n >= nmax:
            raise SystemExit("regime C limit %s" % p.label)

def certify(p, cap=10 ** 4, rel_stop=1e-9):
    """chain of blocks [n1, m]: gain(n1, m) - loss(m) > 0 certified; floats only choose m."""
    n0, worstC = small_n(p)
    n = n0; steps = 0; arbcalls = 0
    while True:
        n1 = n + 1
        f = lambda m: f_gain(n1, m, p) - f_loss(m, p) * (1 + 1e-9) - 1e-9 * abs(f_gain(n1, m, p))
        if f(n1) <= 0 or steps >= 200:              # frontier reached (relative 1e-9 of the float bound)
            return n, n0, steps, worstC
        lo, hi = n1, cap * n1
        if f(hi) > 0:
            lo = hi
        else:
            while hi - lo > max(1, lo // 10 ** 12):  # bisection, geometric while hi/lo is large
                mid = (lo + hi) // 2 if hi < 4 * lo else int(math.sqrt(lo) * math.sqrt(hi))
                mid = min(max(mid, lo + 1), hi - 1)
                if f(mid) > 0: lo = mid
                else: hi = mid
        m = lo
        ok = False
        for _ in range(12):
            g = gain(n1, m, p); arbcalls += 1
            if g is not None and g - loss(m, p) > 0:
                ok = True; break
            if m == n1:
                break
            m = n1 + (m - n1) // 2
        if not ok:
            return n, n0, steps, worstC
        n = m; steps += 1
        if n1 > 10 ** 5 and m - n1 < rel_stop * n1:  # near the frontier: steps no longer move it
            return n, n0, steps, worstC

def run_q(q, a, label=None):
    t0 = t0_of(q)
    p = P(q, q, t0, a, label or ('q=%d a=%s' % (q, a)))
    t1 = time.time()
    N, n0, steps, worstC = certify(p)
    return N, float(N / t0 ** 2), n0, steps, time.time() - t1

A_HALF, A_38 = Fraction(1, 2), Fraction(3, 8)

def has_primitive(q):
    return q >= 3 and q % 4 != 2

# ---------------------------------------------------------------- part 0: the per-zero inequality (random tests)
def part0(trials=4000, seed=1):
    """pair rho = 1/2 + d + i y, rho* = 1 - conj(rho) (0 < d <= a < 1/2 or d < a = 1/2): per zero
    (1/2)[2 - Re(z^n + z*^n)] >= (1/2)[G_n(y) - q_n(y)];  also for y < 0 (same real part)."""
    import random
    from flint import acb
    rnd = random.Random(seed); worst = None
    for _ in range(trials):
        a = arb(rnd.choice([1, 2, 3, 4])) / 8
        d = a * arb(rnd.random())
        y = arb(10) ** arb(rnd.uniform(1.0, 4))
        if rnd.random() < 0.5:
            y = -y
        n = int(float(y * y) * rnd.uniform(0.01, 120)) + 1
        rho = acb(R(1, 2) + d, y)
        z = 1 - 1 / rho
        zs = 1 - 1 / (1 - rho.conjugate())
        lhs = (2 - z ** n - zs ** n).real
        ay = abs(y)
        rhs = G(n, ay) - 4 * (n * a / (2 * ay * ay)).sinh() ** 2 - 4 * n * a * a / ay ** 3
        if not (lhs - rhs + arb(10) ** -30 * (1 + abs(rhs)) > 0):
            raise SystemExit("FAIL per-zero inequality")
        sl = float(lhs - rhs); worst = sl if worst is None or sl < worst else worst
    print("part 0: per-zero inequality holds in %d random cases (|gamma| in [10, 10^4], min slack %.3g)" % (trials, worst))

def blocks(kind, qmin, qmax, ratio):
    """moduli with primitive characters in (qmin, qmax]: kind 'odd' or 'even' (q = 0 mod 4); blocks [q1, q2], q2 <= ratio q1."""
    out = []; q1 = qmin + 1
    while q1 <= qmax:
        q2 = min(qmax, int(q1 * ratio))
        mem = [m for m in range(q1, q2 + 1) if (m % 2 == 1 if kind == 'odd' else m % 4 == 0)]
        if mem:
            out.append((mem[0], mem[-1]))
        q1 = q2 + 1
    return out

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else 'small'
    print("K0 = %s, K1 = %s, K1INF = %s" % (K0.str(15), K1.str(15), K1INF.str(15)))
    T_START = time.time()
    if mode == 'part0':
        part0()
    if mode == 'small':
        for q in [3, 4, 5, 7, 8]:
            for a in (A_HALF, A_38):
                N, r, n0, steps, dt = run_q(q, a)
                print("q=%-6d a=%-4s t0=%-14s N=%d  N/t0^2=%.4f  (regime C n<=%d, %d steps, %.1fs)"
                      % (q, a, float(t0_of(q)), N, r, n0, steps, dt), flush=True)
    if mode == 'q1000':      # every modulus 3 <= q <= 1000 with primitive characters, one at a time
        lo = int(sys.argv[2]) if len(sys.argv) > 2 else 3
        hi = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
        for q in range(lo, hi + 1):
            if not has_primitive(q):
                continue
            Nh, rh, n0h, sth, dth = run_q(q, A_HALF)
            N3, r3, n03, st3, dt3 = run_q(q, A_38)
            print("q %d t0 %s N_half %d r_half %.4f N_38 %d r_38 %.4f  nC %d %d  steps %d %d"
                  % (q, t0_of(q), Nh, rh, N3, r3, n0h, n03, sth, st3), flush=True)
    if mode == 'blocks':     # blocks of moduli in (qmin, qmax], odd and q = 0 mod 4 separately
        ratio = float(sys.argv[2]) if len(sys.argv) > 2 else 1.02
        qmin = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
        qmax = int(sys.argv[4]) if len(sys.argv) > 4 else 400000
        for kind in ('odd', 'even'):
            for (q1, q2) in blocks(kind, qmin, qmax, ratio):
                t0 = t0_of(q2)                      # t0 decreasing in q within a parity class
                res = []
                for a in (A_HALF, A_38):
                    p = P(q1, q2, t0, a, '%s [%d,%d] a=%s' % (kind, q1, q2, a))
                    N, n0, steps, w = certify(p)
                    res.append(N)
                print("block %s %d %d t0min %s N_half %d N_38 %d  N_half/t0(q1)^2 %.4f N_38/t0(q1)^2 %.4f"
                      % (kind, q1, q2, t0, res[0], res[1], float(res[0] / t0_of(q1) ** 2), float(res[1] / t0_of(q1) ** 2)),
                      flush=True)
    print("total runtime %.1f s" % (time.time() - T_START))
