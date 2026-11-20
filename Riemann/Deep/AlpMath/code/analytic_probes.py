# Probes A1-A8 (Chowla-Selberg / Kronecker limit folder), C1-C4 (class number formulas), B1-B4 (Bernoulli)
# statements quoted in PROBES.md; independent numerics with mpmath / sympy / exact rationals.
import sys, math, cmath, random
from fractions import Fraction
import mpmath as mp
from sympy import bernoulli as sb, Symbol, expand, Poly, factorint, isprime, primerange
random.seed(7)
mp.mp.dps = 30
out = []
def rep(tag, ok, msg):
    print(f"{tag}: {'PASS' if ok else 'FAIL'} {msg}"); sys.stdout.flush(); out.append((tag, ok))

def eta(tau):
    q = mp.exp(2j * mp.pi * tau)
    return mp.exp(2j * mp.pi * tau / 24) * mp.qp(q)  # qp(q) = prod (1-q^n)

# A1 CotangentSum: sum_n ((n+u)^2+a^2)^-1 = pi/a (1 + 2 Re sum_k exp(2 pi i (u + i a))^(k+1))
err = 0
for _ in range(20):
    a = mp.mpf(random.uniform(0.05, 3)); u = mp.mpf(random.uniform(-2, 2))
    lhs = mp.nsum(lambda n: 1 / ((n + u) ** 2 + a ** 2), [-mp.inf, mp.inf])
    z = mp.exp(2j * mp.pi * (u + 1j * a))
    rhs = mp.pi / a * (1 + 2 * mp.re(z / (1 - z)))
    err = max(err, abs(lhs - rhs) / abs(rhs))
rep("A1 CotangentSum", err < 1e-20, f"20 random (a,u), max rel err {mp.nstr(err,3)}")
# A2 PowerIntegral: int_R (t^2+a^2)^-s dt = sqrt(pi) Gamma(s-1/2)/Gamma(s) a^(1-2s)
err = 0
for _ in range(12):
    a = mp.mpf(random.uniform(0.2, 3)); s = mp.mpf(random.uniform(0.55, 4))
    lhs = 2 * mp.quad(lambda t: (t * t + a * a) ** (-s), [0, 1, mp.inf])
    rhs = mp.sqrt(mp.pi) * mp.gamma(s - 0.5) / mp.gamma(s) * a ** (1 - 2 * s)
    err = max(err, abs(lhs - rhs) / rhs)
rep("A2 PowerIntegral", err < 1e-12, f"12 random (a,s), max rel err {mp.nstr(err,3)}")
# A3 GammaZetaLimit: lim_{s->1+} 2 sqrt(pi)Gamma(s-1/2)/Gamma(s) y^(1-s) zeta(2s-1) - pi/(s-1) = 2 pi (gamma - log 2) - pi log y
mp.mp.dps = 60
err = 0
for y in [0.3, 1, 2.5, 7]:
    y = mp.mpf(y); s = 1 + mp.mpf(10) ** -25
    lhs = 2 * mp.sqrt(mp.pi) * mp.gamma(s - 0.5) / mp.gamma(s) * y ** (1 - s) * mp.zeta(2 * s - 1) - mp.pi / (s - 1)
    rhs = 2 * mp.pi * (mp.euler - mp.log(2)) - mp.pi * mp.log(y)
    err = max(err, abs(lhs - rhs))
rep("A3 GammaZetaLimit", err < 1e-15, f"y in {{0.3,1,2.5,7}}, s-1=1e-25, max abs err {mp.nstr(err,3)}")
mp.mp.dps = 30
# A4 LogEta: sum_m 1/m Re sum_k q^(m k) = -pi Im tau/12 - log|eta(tau)|
err = 0
for _ in range(10):
    tau = mp.mpc(random.uniform(-1, 1), random.uniform(0.3, 2))
    q = mp.exp(2j * mp.pi * tau)
    lhs = mp.nsum(lambda m: mp.re(q ** m / (1 - q ** m)) / m, [1, mp.inf])
    rhs = -mp.pi * mp.im(tau) / 12 - mp.log(abs(eta(tau)))
    err = max(err, abs(lhs - rhs))
rep("A4 LogEta", err < 1e-20, f"10 random tau, max abs err {mp.nstr(err,3)}")
# A5 DeltaEta: g2^3 - 27 g3^2 of Z+Z tau = (2 pi)^12 eta^24, g2 = 60 G4, g3 = 140 G6 (lattice sums via classical q-expansions)
def sig(k, n): return sum(d ** k for d in range(1, n + 1) if n % d == 0)
err = 0
for _ in range(8):
    tau = mp.mpc(random.uniform(-0.5, 0.5), random.uniform(0.8, 1.6))
    q = mp.exp(2j * mp.pi * tau)
    G4 = 2 * mp.zeta(4) * (1 + 240 * sum(sig(3, n) * q ** n for n in range(1, 80)))
    G6 = 2 * mp.zeta(6) * (1 - 504 * sum(sig(5, n) * q ** n for n in range(1, 80)))
    lhs = (60 * G4) ** 3 - 27 * (140 * G6) ** 2
    rhs = (2 * mp.pi) ** 12 * eta(tau) ** 24
    err = max(err, abs(lhs - rhs) / abs(rhs))
rep("A5 DeltaEta", err < 1e-15, f"8 random tau, max rel err {mp.nstr(err,3)}")
# lattice-sum cross-check of G4 at tau = i (independent of q-expansion)
tau = mp.mpc(0, 1); R = 60
G4d = sum(1 / (m * tau + n) ** 4 for m in range(-R, R + 1) for n in range(-R, R + 1) if (m, n) != (0, 0))
q = mp.exp(2j * mp.pi * tau); G4q = 2 * mp.zeta(4) * (1 + 240 * sum(sig(3, n) * q ** n for n in range(1, 60)))
rep("A5b G4 lattice sum vs q-expansion (tau=i)", abs(G4d - G4q) / abs(G4q) < 1e-4, f"rel diff {mp.nstr(abs(G4d-G4q)/abs(G4q),3)} (cutoff {R})")
# A6 KroneckerLimit via Fourier-Bessel expansion of E(tau,s) = sum' y^s/|m tau+n|^(2s)
def Efour(tau, s, N=60):
    x, y = mp.re(tau), mp.im(tau)
    t = 2 * mp.zeta(2 * s) * y ** s + 2 * mp.sqrt(mp.pi) * mp.gamma(s - 0.5) / mp.gamma(s) * mp.zeta(2 * s - 1) * y ** (1 - s)
    t += 8 * mp.pi ** s * mp.sqrt(y) / mp.gamma(s) * sum(mp.mpf(n) ** (s - 0.5) * sum(mp.mpf(d) ** (1 - 2 * s) for d in range(1, n + 1) if n % d == 0)
                                                         * mp.besselk(s - 0.5, 2 * mp.pi * n * y) * mp.cos(2 * mp.pi * n * x) for n in range(1, N))
    return t
# A8 KroneckerSplit/expansion sanity at s = 2 by direct summation
tau = mp.mpc(0.3, 1.1); s = 2; R = 300
direct = mp.fsum(mp.im(tau) ** s / abs(m * tau + n) ** (2 * s) for m in range(-R, R + 1) for n in range(-R, R + 1) if (m, n) != (0, 0))
split = 2 * mp.im(tau) ** s * mp.zeta(2 * s) + 2 * mp.nsum(lambda m: mp.im(tau) ** s * mp.nsum(lambda n: ((n + m * mp.re(tau)) ** 2 + (m * mp.im(tau)) ** 2) ** (-s), [-mp.inf, mp.inf]), [1, mp.inf])
four = Efour(tau, mp.mpf(s))
rep("A8 KroneckerSplit (s=2, tau=0.3+1.1i)", abs(split - four) < 1e-15 and abs(direct - split) < 1e-4, f"split-fourier {mp.nstr(abs(split-four),3)}, direct(R=300)-split {mp.nstr(abs(direct-split),3)}")
mp.mp.dps = 50
err = 0
for _ in range(6):
    tau = mp.mpc(random.uniform(-0.5, 0.5), random.uniform(0.9, 2))
    s = 1 + mp.mpf(10) ** -20
    lhs = Efour(tau, s) - mp.pi / (s - 1)
    rhs = 2 * mp.pi * (mp.euler - mp.log(2) - mp.log(mp.sqrt(mp.im(tau)) * abs(eta(tau)) ** 2))
    err = max(err, abs(lhs - rhs))
rep("A6 KroneckerLimit", err < 1e-12, f"6 random tau, s-1=1e-20, max abs err {mp.nstr(err,3)}")
mp.mp.dps = 30
# class numbers
def kron(d, n):
    if n == 0: return 1 if abs(d) == 1 else 0
    r = 1
    while n % 2 == 0:
        n //= 2
        if d % 2 == 0: return 0
        r *= 1 if d % 8 in (1, 7) else -1
    # Jacobi (d/n), n odd positive
    a = d % n; res = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: res = -res
        a %= n
    return r * res if n == 1 else 0
def sqfree(m): return all(e == 1 for e in factorint(abs(m)).values())
def fund(d):
    if d % 4 == 1: return sqfree(d)
    if d % 4 == 0: return (d // 4) % 4 in (2, 3) and sqfree(d // 4)
    return False
def hforms(d):  # d < 0, primitive reduced forms
    h = 0; a = 1
    while 3 * a * a <= -d:
        for b in range(-a + 1, a + 1):
            if (b * b - d) % (4 * a): continue
            c = (b * b - d) // (4 * a)
            if c < a: continue
            if c == a and b < 0: continue
            if math.gcd(math.gcd(a, abs(b)), c) != 1: continue
            h += 1
        a += 1
    return h
imag = [d for d in range(-3, -3001, -1) if fund(d)]
bad1 = bad2 = bad4 = 0; mx2 = 0
for d in imag:
    h = hforms(d); D = -d; w = 6 if d == -3 else 4 if d == -4 else 2
    S = sum(kron(d, a) * a for a in range(D))
    if Fraction(-w * S, 2 * D) != h: bad1 += 1
    L1 = -mp.fsum(kron(d, a) * mp.digamma(mp.mpf(a) / D) for a in range(1, D + 1)) / D
    e2 = abs(w * mp.sqrt(D) * L1 / (2 * mp.pi) - h); mx2 = max(mx2, e2)
    if e2 > 1e-12: bad2 += 1
    if not h <= math.sqrt(D) * (1 + math.log(math.sqrt(D))): bad4 += 1
rep("C1 classNumber_eq_sum_of_discr_neg", bad1 == 0, f"{len(imag)} fundamental d in [-3000,-3], h by reduced forms, failures {bad1}")
rep("C2 classNumber_eq_of_discr_neg", bad2 == 0, f"{len(imag)} fields, L(1,chi) by digamma, max |err| {mp.nstr(mx2,3)}")
rep("C4 classNumber_le_sqrt_discr_mul_pow_pred (n=2, imaginary)", bad4 == 0, f"{len(imag)} fields, failures {bad4}")
# C3 real quadratic: h R = sqrt(d) L(1,chi)/2 -> predicted h integral; known values
def fund_unit(d):
    # smallest x,y>0 with x^2 - d y^2 = +-4
    y = 1
    while True:
        for sgn in (-4, 4):
            x2 = d * y * y + sgn
            if x2 > 0:
                x = math.isqrt(x2)
                if x * x == x2: return (mp.mpf(x) + y * mp.sqrt(d)) / 2
        y += 1
        if y > 10 ** 6: return None
real = [d for d in range(5, 801) if fund(d)]
bad3 = 0; mx3 = 0; known = {5: 1, 40: 2, 60: 2, 316: 3, 328: 4}; kn_ok = True; hs = {}
for d in real:
    eps = fund_unit(d)
    if eps is None: continue
    L1 = -mp.fsum(kron(d, a) * mp.digamma(mp.mpf(a) / d) for a in range(1, d + 1)) / d
    hp = mp.sqrt(d) * L1 / (2 * mp.log(eps)); hr = int(mp.nint(hp)); hs[d] = hr
    e = abs(hp - hr); mx3 = max(mx3, e)
    if e > 1e-10 or hr < 1: bad3 += 1
for d, h in known.items():
    if hs.get(d) != h: kn_ok = False
rep("C3 classNumber_mul_regulator_eq_of_discr_pos", bad3 == 0 and kn_ok, f"{len(real)} fundamental d<=800: predicted h integral (max dev {mp.nstr(mx3,3)}), table {known} matches={kn_ok}")
# Bernoulli
B = {k: sb(k) for k in range(2, 203)}
bad = 0
for k in range(0, 201):
    v = B[k + 2]
    for p in primerange(2, 100):
        if (Fraction(int(v.p), int(v.q)) * p).denominator % p == 0: bad += 1
rep("B1 prime_mul_bernoulli_mem (p*B_{k+2} p-integral)", bad == 0, f"k<=200, p<100, failures {bad}")
err = 0
for k in range(1, 13):
    for x in [0.1, 0.37, 0.5, 0.93, 1.0]:
        lhs = mp.zeta(1 - k, x); rhs = -mp.bernpoly(k, x) / k
        err = max(err, abs(lhs - rhs))
rep("B2 hurwitzZeta_one_sub_eq_bernoulli", err < 1e-20, f"k<=12, 5 values of x in (0,1], max err {mp.nstr(err,3)}")
X = Symbol('x'); from sympy import bernoulli as bpoly
bad = sum(1 for k in range(0, 41) if expand(bpoly(k, 1 - X) - (-1) ** k * bpoly(k, X)) != 0)
rep("B3 bernoulli_poly_reflection", bad == 0, f"k<=40 exact, failures {bad}")
# B4 generalized Bernoulli change of level
def bcoef(N, k, a): return Fraction(N) ** (k - 1) * Fraction(Poly(bpoly(k, X), X).eval(Fraction(a, N)))
def genB(N, k, phi): return sum(phi(a) * complex(bcoef(N, k, a)) for a in range(1, N + 1))
def chars_mod(N):
    if N == 1: return [lambda a: 1.0]
    g = next(x for x in range(2, N) if math.gcd(x, N) == 1 and len({pow(x, e, N) for e in range(N)}) == sum(1 for t in range(1, N) if math.gcd(t, N) == 1))
    phiN = sum(1 for t in range(1, N) if math.gcd(t, N) == 1)
    dl = {pow(g, e, N): e for e in range(phiN)}
    return [(lambda j: (lambda a: cmath.exp(2j * math.pi * j * dl[a % N] / phiN) if math.gcd(a, N) == 1 else 0))(j) for j in range(phiN)]
err = 0; cnt = 0
for N in [1, 3, 4, 5, 7, 9, 11]:
    for chi in chars_mod(N):
        for m in [1, 2, 3, 4, 6]:
            M = N * m
            ind = lambda a, chi=chi, M=M: chi(a) if math.gcd(a, M) == 1 else 0
            for k in range(1, 6):
                lhs = genB(M, k, ind)
                prod = 1
                for p in factorint(M):
                    prod *= (1 - chi(p) * p ** (k - 1)) if N > 1 else (1 - p ** (k - 1))
                rhs = genB(N, k, chi) * prod
                err = max(err, abs(lhs - rhs)); cnt += 1
rep("B4 generalized_bernoulli_change_level", err < 1e-6, f"{cnt} cases (N<=11 cyclic, m<=6, k<=5), max |err| {err:.2e}")
# A7 ClassMean (Chowla-Selberg): mean over classes of log(4 pi Im tau |eta|^4) = gamma + log(2 pi/sqrt p) - L'(1)/L(1)
mp.mp.dps = 80
ps = [p for p in primerange(3, 270) if p % 4 == 3]
err = 0
for p in ps:
    d = -p; forms = []
    a = 1
    while 3 * a * a <= p:
        for b in range(-a + 1, a + 1):
            if (b * b + p) % (4 * a): continue
            c = (b * b + p) // (4 * a)
            if c < a or (c == a and b < 0) or math.gcd(math.gcd(a, abs(b)), c) != 1: continue
            forms.append((a, b, c))
        a += 1
    lp = [mp.log(4 * mp.pi * (mp.sqrt(p) / (2 * a)) * abs(eta(mp.mpc(mp.mpf(-b) / (2 * a), mp.sqrt(p) / (2 * a)))) ** 4) for (a, b, c) in forms]  # patched after the run: was a float -b/(2a); a7_rerun.py is authoritative
    mean = mp.fsum(lp) / len(lp)
    s = 1 + mp.mpf(10) ** -30
    chi = [kron(d, t) for t in range(1, p + 1)]
    L = p ** (-s) * mp.fsum(chi[t - 1] * mp.zeta(s, mp.mpf(t) / p) for t in range(1, p + 1))
    Ld = -mp.log(p) * L + p ** (-s) * mp.fsum(chi[t - 1] * mp.zeta(s, mp.mpf(t) / p, 1) for t in range(1, p + 1))
    rhs = mp.euler + mp.log(2 * mp.pi / mp.sqrt(p)) - Ld / L
    err = max(err, abs(mean - rhs))
rep("A7 ClassMean (Chowla-Selberg, prime discriminants)", err < 1e-20, f"{len(ps)} primes p=3 mod 4, p<270, max |err| {mp.nstr(err,3)}")
print("SUMMARY", sum(1 for _, ok in out if ok), "pass of", len(out))
