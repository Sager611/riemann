#!/usr/bin/env python3
"""Re-check the certificates of NOTE.md.  Standard library + mpmath only.

c_n(k,m) are the Taylor coefficients of 1/q_{k,m}(z),
q_{k,m}(z) = sum_j (-1)^j C(m,j) (k+m-j)!/(k+m)! z^j  (Pade denominator of e^z).

A  m = 3, every k: exact polynomial identity + sign argument behind the theorem.
B  first failures: c_n(kmin(m), m) < 0 exactly, for every odd m in 5..101, and
   c_n(kmin(m)-1, m) > 0 for all n <= 4(k+m)+200.
C  diagonal: c_n(m,m) > 0 for ALL n >= 0, odd m <= MV (default 101): exact signs
   for n < n1 and a rigorous interval-arithmetic pole bound for n >= n1.

Exact signs use an integer fixed-point recurrence with a proven error bound:
Y_n = round(sum_j -b_j s^j Y_{n-j} / b_0) is 2^P c_n s^n up to an error
|E_n| <= (1/2) sum_{i<n} |c_i s^i|  (errors propagate through the recurrence's
own impulse response, which is c_n s^n itself), so the sign is certain when
|Y_n| exceeds the running bound.  Usage: python3 verify.py [MV]
"""
import sys, time
from math import comb
from fractions import Fraction
from operator import mul
from mpmath import mp, iv, mpf, mpc

T0 = time.time()
FAIL = []

def check(cond, msg):
    print(('ok   ' if cond else 'FAIL ') + msg, flush=True)
    if not cond:
        FAIL.append(msg)

def bcoef(k, m):
    """integers b_j = (-1)^j C(m,j) (k+m-j)!/k!, so that q_{k,m}(z) = sum_j b_j z^j / b_0"""
    rising = [1]                                   # rising[i] = (k+1)(k+2)...(k+i)
    for i in range(1, m + 1):
        rising.append(rising[-1] * (k + i))
    return [(-1) ** j * comb(m, j) * rising[m - j] for j in range(m + 1)]

def scale(k, m):
    """largest integer s with q_{k,m}(s) > 0 (q decreases on R for odd m), at least 1"""
    b = bcoef(k, m)
    Q = lambda x: sum(bj * x ** j for j, bj in enumerate(b))
    lo, hi = 0, 1
    while Q(hi) > 0:
        lo, hi = hi, 2 * hi
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if Q(mid) > 0: lo = mid
        else: hi = mid
    return max(1, lo)

def signs(k, m, nmax, P=None):
    """certified signs (+1/-1, 0 = undecided) of c_n(k,m) for 0 <= n <= nmax"""
    s = scale(k, m)
    P = P or (64 + 2 * m)
    while True:
        b = bcoef(k, m); b0 = b[0]
        alpha = [-b[j] * s ** j for j in range(1, m + 1)]   # alpha[j-1] multiplies Y_{n-j}
        Y = [1 << P]; out = [1]
        S = 1 << P                     # S / 2^P >= sum_{i<n} |c_i s^i|
        sh = P + 1
        for n in range(1, nmax + 1):
            j = min(m, n)
            N = sum(map(mul, alpha[:j], Y[n - 1::-1][:j] if n - 1 - j < 0 else Y[n - 1:n - 1 - j:-1]))
            q, r = divmod(N, b0)
            if 2 * r >= b0: q += 1
            Y.append(q)
            B = -((-S) >> sh)          # ceil(S / 2^(P+1)) >= |error of Y_n|
            out.append(1 if q > B else (-1 if q < -B else 0))
            S += abs(q) + B
        if 0 not in out or P > 1 << 16: return out
        P *= 2

# ---------------------------------------------------------------- A: m = 3
def part_A():
    print('== A: m = 3 for every k')
    # polynomials in d as lists of Fractions (index = power)
    def pmul(a, b):
        out = [Fraction(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return out
    def padd(a, b):
        n = max(len(a), len(b)); a = a + [0] * (n - len(a)); b = b + [0] * (n - len(b))
        return [Fraction(x) + y for x, y in zip(a, b)]
    def psc(a, c): return [c * x for x in a]
    def trim(a):
        while a and a[-1] == 0: a = a[:-1]
        return a
    d = [0, 1]; one = [1]
    e = [2, -3]                                   # e = 2 - 3d,  N = d^3 / e
    d3 = pmul(pmul(d, d), d)
    Nd = padd(d3, pmul(d, e))                     # e*(N+d)
    A1 = padd(pmul(Nd, Nd), psc(pmul(pmul(d3, padd(one, psc(d, -1))), e), 3))  # e^2*((N+d)^2+3N(1-d))
    A2 = padd(psc(d3, 3), pmul(padd(padd(pmul(d, d), d), one), e))            # e*(3N+d^2+d+1)
    om = padd(one, psc(d, -1))
    lhs = padd(pmul(d3, pmul(A1, A2)), psc(pmul(d3, pmul(pmul(Nd, Nd), pmul(om, om))), -8))
    d5 = pmul(d3, pmul(d, d))
    rhs = pmul(pmul(d5, pmul(om, om)), pmul([-2, 3], [12, -30, 12, -3]))
    check(trim(lhs) == trim(rhs),
          'e^3*(d^3[(N+d)^2+3N(1-d)][3N+d^2+d+1] - 8N(N+d)^2(1-d)^2) = d^5(1-d)^2(3d-2)(-3d^3+12d^2-30d+12), N=d^3/(2-3d)')
    g = lambda x: -3 * x ** 3 + 12 * x ** 2 - 30 * x + 12
    check(24 ** 2 - 4 * 9 * 30 < 0 and g(Fraction(1, 2)) < 0,
          "g'(d) = -9d^2+24d-30 < 0 everywhere and g(1/2) = -3/8 < 0, so g < 0 on [1/2, inf)")
    # d(N) in (1/2, 2/3) for N >= 1:  F(d) = d^3+3Nd-2N increasing, F(1/2) = 1/8-N/2 < 0, F(2/3) = 8/27 > 0
    check(Fraction(1, 8) - Fraction(1, 2) < 0, 'F(1/2) < 0 for N >= 1, F(2/3) = 8/27 > 0  =>  1/2 < d < 2/3')
    # numerical cross-check of the closed forms against roots of the cubic
    mp.dps = 40
    worst = 0; allok = True
    for k in [0, 1, 2, 7, 30, 500, 10 ** 6]:
        N = mpf(k + 1)
        dd = mp.findroot(lambda x: x ** 3 + 3 * N * x - 2 * N, mpf(2) / 3 - mpf(2) / (27 * N + 1))
        ratio2 = 8 * N * (N + dd) ** 2 * (1 - dd) ** 2 / (dd ** 3 * ((N + dd) ** 2 + 3 * N * (1 - dd)) * (dd ** 2 + dd + 1 + 3 * N))
        b = bcoef(k, 3)
        rts = mp.polyroots([mpf(x) for x in b[::-1]], maxsteps=200, extraprec=200)
        dQ = lambda z: b[1] + 2 * b[2] * z + 3 * b[3] * z ** 2
        rts = sorted(rts, key=lambda z: abs(z))
        R = [-b[0] / (z * dQ(z)) for z in rts]
        direct = (2 * abs(R[1]) / abs(R[0])) ** 2
        worst = max(worst, abs(direct - ratio2) / ratio2)
        allok = allok and ratio2 < 1 and abs(rts[1]) > abs(rts[0]) and abs(mp.im(rts[0])) < 1e-30
    check(allok and worst < mpf(10) ** -30, 'closed form (2|R+|/R0)^2 < 1 matches roots of q_{k,3} (k = 0..1e6 sample), rel.err %.1e' % float(worst))
    bad = [(k, n) for k in range(0, 201) for n, sg in enumerate(signs(k, 3, 400)) if sg != 1]
    check(not bad, 'exact: c_n(k,3) > 0 for 0 <= k <= 200, 0 <= n <= 400')

# ------------------------------------------------- B: first failures, odd m <= 101
KMIN = {5: (66, 38), 7: (44, 37), 9: (42, 41), 11: (43, 46), 13: (44, 51), 15: (48, 57), 17: (50, 62),
        19: (54, 68), 21: (57, 74), 23: (61, 80), 25: (65, 86), 27: (69, 92), 29: (73, 98), 31: (77, 104),
        33: (81, 110), 35: (83, 115), 37: (87, 121), 39: (91, 127), 41: (95, 133), 43: (99, 139),
        45: (103, 145), 47: (107, 151), 49: (111, 157), 51: (115, 163), 53: (119, 169), 55: (123, 175),
        57: (127, 181), 59: (131, 187), 61: (135, 193), 63: (139, 199), 65: (143, 205), 67: (147, 211),
        69: (151, 217), 71: (155, 223), 73: (160, 230), 75: (164, 236), 77: (168, 242), 79: (172, 248),
        81: (176, 254), 83: (180, 260), 85: (184, 266), 87: (188, 272), 89: (192, 278), 91: (196, 284),
        93: (200, 290), 95: (202, 295), 97: (206, 301), 99: (210, 307), 101: (214, 313)}

def part_B():
    print('== B: first failures k_min(m), odd m = 5..101')
    for (k, m, n) in [(66, 5, 38), (44, 7, 37), (45, 9, 42)]:
        check(signs(k, m, n)[n] == -1, 'paper example c_%d(%d,%d) < 0' % (n, k, m))
    bad = []
    for m, (k, n) in sorted(KMIN.items()):
        sg = signs(k, m, n)
        if not (sg[n] == -1 and all(x == 1 for x in sg[:n])): bad.append(('neg', m, k, n))
        L = 4 * (k - 1 + m) + 200
        if any(x != 1 for x in signs(k - 1, m, L)): bad.append(('pos', m, k - 1))
    check(not bad, 'for all 49 odd m in 5..101: first negative of c_n(kmin,m) at the listed n; '
                   'c_n(kmin-1,m) > 0 for n <= 4(kmin-1+m)+200   %s' % bad[:3])

# ------------------------------------------------ C: diagonal, all n, odd m <= MV
def roots_aberth(b, guesses, dps):
    """polish approximations of all roots of sum b_j z^j (Aberth-Ehrlich, Gauss-Seidel)"""
    mp.dps = dps
    m = len(b) - 1
    cf = [mpf(x) for x in b]
    z = [mpc(g) for g in guesses]
    done = 0
    for it in range(200):
        mx = 0
        for j in range(m):
            zj = z[j]
            p = cf[m]; dp = mpf(0)
            for c in cf[m - 1::-1]:
                dp = dp * zj + p; p = p * zj + c
            ratio = p / dp
            sm = sum(1 / (zj - z[i]) for i in range(m) if i != j)
            w = ratio / (1 - ratio * sm)
            z[j] = zj - w
            mx = max(mx, abs(w) / abs(zj))
        if mx < mpf(10) ** (-(dps // 2)):
            done += 1
            if done >= 2: break          # one extra sweep past the target (cubic convergence)
    return z

def guess_from(prev, m):
    """initial guesses for degree m from the sorted roots of degree m-2 (continuation along the curve)"""
    if prev is None: return None
    mo = len(prev)
    up = sorted([z for z in prev if mp.im(z) > 0], key=lambda z: mp.arg(z))   # upper half, by angle
    sc = mpf(m) / mo
    pts = [mpc(mp.re(prev_r) * sc, 0) for prev_r in prev if abs(mp.im(prev_r)) < 1e-20][:1]
    t_old = [(i + 1) / mpf(len(up) + 1) for i in range(len(up))]
    nu = (m - 1) // 2
    out = list(pts)
    curve = [mpc(mp.re(pts[0]), 0)] + [u * sc for u in up]
    tt = [mpf(0)] + t_old
    for i in range(nu):
        t = (i + 1) / mpf(nu + 1)
        # piecewise-linear interpolation, linear extrapolation at the end
        h = max(1, min(len(tt) - 1, next((q for q in range(1, len(tt)) if tt[q] >= t), len(tt) - 1)))
        a, bb = tt[h - 1], tt[h]
        z = curve[h - 1] + (curve[h] - curve[h - 1]) * (t - a) / (bb - a)
        out += [z, mp.conj(z)]
    return out

def tail_cert(m, approx, dps):
    """rigorous: returns n1 such that c_n(m,m) > 0 for all n >= n1 (or None)"""
    b = bcoef(m, m)
    iv.dps = dps
    B = [iv.mpf(x) for x in b]
    zs = [iv.mpc(mp.re(z), mp.im(z)) for z in approx]
    up = lambda x: mpf(x.b)          # exact endpoints of an interval
    lo_ = lambda x: mpf(x.a)
    rho = []
    for j, zj in enumerate(zs):
        p = B[m]
        for c in B[m - 1::-1]: p = p * zj + c
        pr = B[m]
        for i, zi in enumerate(zs):
            if i != j: pr = pr * (zj - zi)
        rho.append(up(m * abs(p / pr)))   # root j lies in the disc |z - z_j| <= m|W_j| (Gerschgorin)
    for i in range(m):                     # disjoint discs => exactly one root in each
        for j in range(i + 1, m):
            if not lo_(abs(zs[i] - zs[j])) > up(iv.mpf(rho[i]) + iv.mpf(rho[j])): return None
    j1 = min(range(m), key=lambda j: abs(approx[j]))
    if mp.im(approx[j1]) != 0: return None  # disc symmetric about R => its root is real
    unit = iv.mpf([-1, 1])
    boxes = [iv.mpc(iv.mpf(mp.re(z)) + unit * rr, iv.mpf(mp.im(z)) + unit * rr) for z, rr in zip(approx, rho)]
    def zdQ(Z):
        dp = iv.mpf(0); p = B[m]
        for c in B[m - 1::-1]:
            dp = dp * Z + p; p = p * Z + c
        return Z * dp
    R1 = (-B[0] / zdQ(boxes[j1])).real
    if not lo_(R1) > 0: return None
    rhi = up(abs(boxes[j1]))
    terms = []
    for j in range(m):
        if j == j1: continue
        lo = lo_(abs(boxes[j]))
        if not lo > rhi: return None       # D1 (unique minimal modulus) not certified
        terms.append((iv.mpf(up(abs(-B[0] / zdQ(boxes[j])))) / iv.mpf(lo_(R1)), iv.mpf(rhi) / iv.mpf(lo)))
    def T(n):
        return up(sum(a * g ** n for a, g in terms))
    lo, hi = 0, 1
    while not T(hi) < 1:
        hi *= 2
        if hi > 10 ** 7: return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if T(mid) < 1: hi = mid
        else: lo = mid
    return hi

def part_C(MV):
    print('== C: diagonal k = m, all n, odd m <= %d' % MV)
    prev = None; bad = []
    for m in range(3, MV + 1, 2):
        b = bcoef(m, m)
        dps = 2 * m + 40
        if prev is None:
            mp.dps = dps
            g0 = mp.polyroots([mpf(x) for x in b[::-1]], maxsteps=300, extraprec=300)
        else:
            g0 = guess_from(prev, m)
        rts = roots_aberth(b, g0, dps)
        rts = [mpc(mp.re(z), 0) if abs(mp.im(z)) < mpf(10) ** (-dps // 2) else z for z in rts]
        n1 = tail_cert(m, rts, dps)
        prev = rts
        if n1 is None: bad.append((m, 'tail')); continue
        sg = signs(m, m, max(n1, 1))
        if any(x != 1 for x in sg): bad.append((m, 'sign'))
        print('     m=%3d  n1=%5d  (%.0fs)' % (m, n1, time.time() - T0), flush=True)
    check(not bad, 'c_n(m,m) > 0 for every n >= 0 and every odd m <= %d   %s' % (MV, bad[:3]))

if __name__ == '__main__':
    MV = int(sys.argv[1]) if len(sys.argv) > 1 else 101
    parts = sys.argv[2] if len(sys.argv) > 2 else 'ABC'
    if 'A' in parts: part_A()
    if 'B' in parts: part_B()
    if 'C' in parts: part_C(MV)
    print('ALL OK' if not FAIL else 'FAILURES: %d' % len(FAIL), ' (%.0f s)' % (time.time() - T0))
    sys.exit(1 if FAIL else 0)
