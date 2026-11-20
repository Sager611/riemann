# Probes Q1 (QuadraticClassFieldTheory ProductFormula, Q-instance = Hilbert reciprocity), W1 (WeylBranching via Schur functions),
# L1 (HeckeLEntireAndNonzeroAtOne, Dirichlet instance), S1 (SiegelKlingen rationality instances)
import sys, math, random, cmath, itertools
import mpmath as mp
from fractions import Fraction
from sympy import factorint, primerange, bernoulli
random.seed(11)
def rep(tag, ok, msg): print(f"{tag}: {'PASS' if ok else 'FAIL'} {msg}"); sys.stdout.flush()
# Q1 Hilbert symbols over Q
def leg(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1
def val(a, p):
    v = 0
    while a % p == 0: a //= p; v += 1
    return v, a
def hilb(a, b, p):
    if p == 0: return -1 if (a < 0 and b < 0) else 1
    va, ua = val(a, p); vb, ub = val(b, p)
    if p != 2:
        eps = (p - 1) // 2
        s = (-1) ** (va * vb * eps)
        return s * (leg(ua, p) ** vb) * (leg(ub, p) ** va)
    e = lambda u: ((u - 1) // 2) % 2
    w = lambda u: ((u * u - 1) // 8) % 2
    ex = e(ua) * e(ub) + va * w(ub) + vb * w(ua)
    return (-1) ** ex
bad = 0; n = 0
for _ in range(3000):
    a = random.choice([-1, 1]) * random.randint(1, 10 ** 5); b = random.choice([-1, 1]) * random.randint(1, 10 ** 5)
    ps = set(factorint(abs(a))) | set(factorint(abs(b))) | {2}
    prod = hilb(a, b, 0)
    for p in ps: prod *= hilb(a, b, p)
    n += 1
    if prod != 1: bad += 1
rep("Q1 product formula (Hilbert reciprocity) for L = Q(sqrt a), t = b", bad == 0, f"{n} random pairs |a|,|b| <= 1e5, failures {bad}")
# W1 Weyl branching U(n+1) -> U(n): s_lam(x1..xn, 1) = sum over interlacing mu of s_mu(x1..xn)
mp.mp.dps = 30
def schur(lam, xs):
    N = len(xs)
    num = mp.matrix(N, N); den = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            num[i, j] = xs[i] ** (lam[j] + N - 1 - j); den[i, j] = xs[i] ** (N - 1 - j)
    return mp.det(num) / mp.det(den)
bad = 0; cnt = 0; err = 0
for N in (2, 3, 4):
    for lam in itertools.product(range(4), repeat=N):
        if any(lam[i] < lam[i + 1] for i in range(N - 1)): continue
        xs = [mp.mpc(random.uniform(0.5, 1.5), random.uniform(-0.5, 0.5)) for _ in range(N - 1)]
        lhs = schur(lam, xs + [mp.mpf(1)])
        rhs = 0
        for mu in itertools.product(range(4), repeat=N - 1):
            if all(lam[i] >= mu[i] >= lam[i + 1] for i in range(N - 1)): rhs += schur(mu, xs)
        e = abs(lhs - rhs) / max(1, abs(lhs)); err = max(err, e); cnt += 1
        if e > 1e-15: bad += 1
rep("W1 WeylBranching (multiplicity one on interlacing weights)", bad == 0, f"{cnt} dominant weights, N = 2..4, entries <= 3, max rel err {mp.nstr(err, 3)}")
# L1 L(1, chi) != 0 for nontrivial Dirichlet chi mod q <= 60 (instance of HeckeLEntireAndNonzeroAtOne over Q)
def chars(q):
    units = [a for a in range(1, q + 1) if math.gcd(a, q) == 1]
    # generic construction: characters as homomorphisms found by brute force on a generating set (q <= 60 small)
    # use structure: enumerate all maps defined via discrete logs on a basis found greedily
    gens = []; sub = {1 % q}
    for a in units:
        if a % q not in sub:
            gens.append(a)
            new = set(sub); cur = set(sub)
            while True:
                cur = {(x * a) % q for x in cur}
                if cur <= new: break
                new |= cur
            sub = new
            # close under products
            changed = True
            while changed:
                changed = False
                for x in list(sub):
                    for y in list(sub):
                        z = (x * y) % q
                        if z not in sub: sub.add(z); changed = True
    # brute force: a character is determined by values on gens (roots of unity of order dividing ord(g)); check consistency
    ords = []
    for g in gens:
        o = 1; x = g % q
        while x != 1 % q: x = (x * g) % q; o += 1
        ords.append(o)
    res = []
    for exps in itertools.product(*[range(o) for o in ords]):
        val = {1 % q: 1}; ok = True; frontier = [1 % q]
        while frontier:
            x = frontier.pop()
            for g, o, ex in zip(gens, ords, exps):
                y = (x * g) % q; v = val[x] * cmath.exp(2j * math.pi * ex / o)
                if y in val:
                    if abs(val[y] - v) > 1e-9: ok = False; break
                else: val[y] = v; frontier.append(y)
            if not ok: break
        if ok and len(val) == len(units): res.append(val)
    # dedupe
    uniq = []
    for v in res:
        if not any(all(abs(v[a] - u[a]) < 1e-9 for a in v) for u in uniq): uniq.append(v)
    return uniq
mn = mp.inf; tot = 0; okc = True
for q in range(3, 61):
    cs = chars(q)
    phi = sum(1 for a in range(1, q + 1) if math.gcd(a, q) == 1)
    if len(cs) != phi: okc = False
    for v in cs:
        if all(abs(v[a] - 1) < 1e-9 for a in v): continue
        L1 = -mp.fsum(v.get(a % q, 0) * mp.digamma(mp.mpf(a) / q) for a in range(1, q + 1)) / q
        mn = min(mn, abs(L1)); tot += 1
rep("L1 L(1,chi) != 0 (nontrivial chi mod q <= 60)", okc and mn > 0.01, f"{tot} characters (count = phi(q) for all q: {okc}), min |L(1,chi)| = {mp.nstr(mn, 4)}")
# S1 SiegelKlingen instances: zeta_F(1-k) rational for F = Q(sqrt d), via zeta(1-k) L(1-k, chi_d) (Hurwitz, numerics) vs exact Bernoulli
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
from sympy import Symbol, Poly, bernoulli as bp
X = Symbol('x')
err = 0; vals = {}
for d in (5, 8, 12, 13, 17, 21, 24, 28, 29):
    for k in (2, 4, 6):
        Lnum = d ** (k - 1) * mp.fsum(kron(d, a) * mp.zeta(1 - k, mp.mpf(a) / d) for a in range(1, d + 1))
        zf = mp.zeta(1 - k) * Lnum
        Bk = sum(kron(d, a) * Fraction(d) ** (k - 1) * Fraction(Poly(bp(k, X), X).eval(Fraction(a, d))) for a in range(1, d + 1))
        exact = Fraction(-1, k) * Fraction(bernoulli(k)) * Fraction(-1, k) * Bk  # zeta(1-k) = -B_k/k, L(1-k,chi) = -B_{k,chi}/k
        err = max(err, abs(zf - mp.mpf(exact.numerator) / exact.denominator)); vals[(d, k)] = exact
ok = vals[(5, 2)] == Fraction(1, 30) and vals[(8, 2)] == Fraction(1, 12) and vals[(13, 2)] == Fraction(1, 6)
rep("S1 SiegelKlingen: zeta_F(1-k) rational, F = Q(sqrt d)", err < 1e-20 and ok, f"9 fields x k in (2,4,6): max |numeric - Bernoulli rational| {mp.nstr(err,3)}; zeta_F(-1) = 1/30, 1/12, 1/6 for d = 5, 8, 13: {ok}")
