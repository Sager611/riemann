# Fixes after the reviewer (aca2af55): C3 with continued-fraction units (no skips) and independent class numbers by
# cycles of reduced indefinite forms; B2 and S1 by independent routes (Hurwitz Fourier series; functional equation from
# L(k)); A5 with direct lattice sums for G4 and G6 at a non-symmetric tau.
import sys, math, random
from fractions import Fraction
import mpmath as mp
from sympy import factorint, Symbol, Poly, bernoulli as bp
random.seed(3)
mp.mp.dps = 30
X = Symbol('x')
def rep(tag, ok, msg): print(f"{tag}: {'PASS' if ok else 'FAIL'} {msg}"); sys.stdout.flush()
def sqfree(m): return all(e == 1 for e in factorint(abs(m)).values())
def fund(d): return (d % 4 == 1 and sqfree(d)) or (d % 4 == 0 and (d // 4) % 4 in (2, 3) and sqfree(d // 4))
def kron(d, n):
    r = 1
    while n % 2 == 0:
        n //= 2
        if d % 2 == 0: return 0
        r *= 1 if d % 8 in (1, 7) else -1
    a = d % n; res = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: res = -res
        a %= n
    return r * res if n == 1 else 0
def fund_unit_cf(D):
    # returns (x, y, norm) with eps = (x + y sqrt D)/2 fundamental, via continued fractions
    if D % 4 == 1: P, Q, N = 1, 2, D
    else: P, Q, N = 0, 1, D // 4
    r = math.isqrt(N)
    pm2, pm1, qm2, qm1 = 0, 1, 1, 0
    for _ in range(100000):
        a = (P + r) // Q
        p, q = a * pm1 + pm2, a * qm1 + qm2
        if D % 4 == 1:
            nrm = p * p - p * q + q * q * (1 - D) // 4
            if nrm in (1, -1): return (2 * p - q, q, nrm)
        else:
            nrm = p * p - N * q * q
            if nrm in (1, -1): return (2 * p, q, nrm)
        pm2, pm1, qm2, qm1 = pm1, p, qm1, q
        P = a * Q - P; Q = (N - P * P) // Q
    raise RuntimeError("no unit")
def hplus_cycles(D):
    r = math.isqrt(D)
    red = set()
    for b in range(1, r + 1):
        if (b - D) % 2: continue
        t = (b * b - D) // 4  # a*c = t < 0
        for a0 in range(1, -t + 1):
            if t % a0: continue
            for a in (a0, -a0):
                c = t // a
                A = abs(a)
                if not (b * b < D and D < (2 * A + b) ** 2 and (2 * A - b <= 0 or (2 * A - b) ** 2 < D)): continue
                if math.gcd(math.gcd(A, b), abs(c)) != 1: continue
                red.add((a, b, c))
    seen = set(); cycles = 0
    for f in red:
        if f in seen: continue
        cycles += 1; g = f
        while g not in seen:
            seen.add(g)
            a, b, c = g; C = abs(c)
            b2 = r - ((r + b) % (2 * C))
            a2 = (b2 * b2 - D) // (4 * c)
            g = (c, b2, a2)
            if g not in red: raise RuntimeError(f"rho left reduced set at D={D}: {g}")
    return cycles
real = [d for d in range(5, 801) if fund(d)]
bad = 0; mxdev = 0; diffs = []
for d in real:
    x, y, nrm = fund_unit_cf(d)
    assert x * x - d * y * y == 4 * nrm
    eps = (mp.mpf(x) + y * mp.sqrt(d)) / 2
    L1 = -mp.fsum(kron(d, a) * mp.digamma(mp.mpf(a) / d) for a in range(1, d + 1)) / d
    hp = mp.sqrt(d) * L1 / (2 * mp.log(eps))
    hpl = hplus_cycles(d); hind = hpl if nrm == -1 else hpl // 2
    dev = abs(hp - hind); mxdev = max(mxdev, dev)
    if dev > 1e-12: bad += 1; diffs.append((d, mp.nstr(hp, 8), hind))
rep("C3' h R = sqrt(d) L(1,chi)/2 vs independent h (cycles of reduced forms, N(eps) from continued fractions)",
    bad == 0, f"all {len(real)} fundamental d in [5,800], no skips, max |h_formula - h_cycles| {mp.nstr(mxdev,3)}, mismatches {diffs[:5]}")
# B2': zeta(1-k, x) = -B_k(x)/k via Hurwitz's Fourier series zeta(1-k,x) = 2 Gamma(k)/(2 pi)^k sum_n cos(pi k/2 - 2 pi n x)/n^k  (k >= 2)
err = 0
for k in range(2, 13):
    for xv in [Fraction(1, 10), Fraction(37, 100), Fraction(1, 2), Fraction(93, 100), Fraction(1)]:
        xx = mp.mpf(xv.numerator) / xv.denominator
        ser = 2 * mp.gamma(k) / (2 * mp.pi) ** k * mp.nsum(lambda n: mp.cos(mp.pi * k / 2 - 2 * mp.pi * n * xx) / n ** k, [1, mp.inf])
        ex = -Fraction(Poly(bp(k, X), X).eval(xv)) / k
        err = max(err, abs(ser - mp.mpf(ex.numerator) / ex.denominator))
rep("B2' Hurwitz zeta at 1-k by Fourier series vs -B_k(x)/k (exact)", err < 1e-15, f"k = 2..12, 5 values of x, max err {mp.nstr(err,3)}")
# S1': zeta_F(1-k) for F = Q(sqrt q): zeta(1-k) and L(1-k,chi) from the functional equations and L(k) by direct summation
err = 0
for q in (5, 8, 12, 13, 17, 21, 24, 28, 29):
    for k in (2, 4, 6):
        Lk = mp.nsum(lambda n: kron(q, int(n)) / n ** k, [1, mp.inf], method='direct', steps=[200000]) if False else mp.fsum(kron(q, n) / mp.mpf(n) ** k for n in range(1, 200001))
        # tail bound for k >= 2: < 1/(k-1) * 200000^(1-k) ~ 5e-6 for k=2; use Euler-Maclaurin-free periodic acceleration instead:
        Lk = mp.fsum(kron(q, a) * mp.zeta(k, mp.mpf(a) / q) for a in range(1, q + 1)) / mp.mpf(q) ** k  # positive-s Hurwitz (convergent series)
        L1mk = (mp.mpf(q) / mp.pi) ** (k - mp.mpf(1) / 2) * mp.gamma(mp.mpf(k) / 2) / mp.gamma((1 - mp.mpf(k)) / 2) * Lk
        z1mk = 2 * (2 * mp.pi) ** (-k) * mp.cos(mp.pi * k / 2) * mp.gamma(k) * mp.zeta(k)
        Bk = sum(kron(q, a) * Fraction(q) ** (k - 1) * Fraction(Poly(bp(k, X), X).eval(Fraction(a, q))) for a in range(1, q + 1))
        exact = Fraction(bp(k)) / k * Bk / k  # zeta(1-k) L(1-k,chi) = (-B_k/k)(-B_{k,chi}/k)
        err = max(err, abs(z1mk * L1mk - mp.mpf(exact.numerator) / exact.denominator))
rep("S1' zeta_F(1-k) via functional equations from s = k vs Bernoulli rationals", err < 1e-20, f"9 fields x k = 2,4,6, max err {mp.nstr(err,3)}")
# A5': direct lattice sums G4, G6 at non-symmetric tau vs q-expansions, then Delta = (2 pi)^12 eta^24
def sig(k, n): return sum(dd ** k for dd in range(1, n + 1) if n % dd == 0)
tau = mp.mpc(0.31, 1.07); R = 150
G4d = mp.fsum(1 / (m * tau + n) ** 4 for m in range(-R, R + 1) for n in range(-R, R + 1) if (m, n) != (0, 0))
G6d = mp.fsum(1 / (m * tau + n) ** 6 for m in range(-R, R + 1) for n in range(-R, R + 1) if (m, n) != (0, 0))
q = mp.exp(2j * mp.pi * tau)
G4q = 2 * mp.zeta(4) * (1 + 240 * mp.fsum(sig(3, n) * q ** n for n in range(1, 60)))
G6q = 2 * mp.zeta(6) * (1 - 504 * mp.fsum(sig(5, n) * q ** n for n in range(1, 60)))
eta = mp.exp(2j * mp.pi * tau / 24) * mp.qp(q)
Dd = (60 * G4d) ** 3 - 27 * (140 * G6d) ** 2
e4, e6 = abs(G4d - G4q) / abs(G4q), abs(G6d - G6q) / abs(G6q)
eD = abs(Dd - (2 * mp.pi) ** 12 * eta ** 24) / abs((2 * mp.pi) ** 12 * eta ** 24)
rep("A5' lattice sums G4, G6 (box 150) at tau = 0.31+1.07i; Delta from lattice sums vs (2pi)^12 eta^24", e4 < 1e-4 and e6 < 1e-8 and eD < 1e-3,
    f"rel diffs G4 {mp.nstr(e4,3)}, G6 {mp.nstr(e6,3)}, Delta {mp.nstr(eD,3)} (box truncation)")
