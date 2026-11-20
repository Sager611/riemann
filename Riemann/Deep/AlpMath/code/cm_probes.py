# Probes H1-H4 (HilbertClassFieldJ / CM: ExistsEqOfNatCastMulMem.lean:66, ComplexMultiplicationJPremises.lean:38) and R1 (RadfordShufflePolynomial)
import sys, math, itertools
import mpmath as mp
from sympy import Symbol, Poly, factor_list, sqrt as ssqrt, isprime, primerange, factorint, discriminant, ZZ
X = Symbol('X')
def rep(tag, ok, msg): print(f"{tag}: {'PASS' if ok else 'FAIL'} {msg}"); sys.stdout.flush()
def sqfree(m): return all(e == 1 for e in factorint(abs(m)).values())
def fund(d): return (d % 4 == 1 and sqfree(d)) or (d % 4 == 0 and (d // 4) % 4 in (2, 3) and sqfree(d // 4))
def forms(d):
    D = -d; out = []; a = 1
    while 3 * a * a <= D:
        for b in range(-a + 1, a + 1):
            if (b * b - d) % (4 * a): continue
            c = (b * b - d) // (4 * a)
            if c < a or (c == a and b < 0) or math.gcd(math.gcd(a, abs(b)), c) != 1: continue
            out.append((a, b, c))
        a += 1
    return out
mp.mp.dps = 150
ds = [d for d in range(-3, -400, -1) if fund(d)]
bad_int = bad_dist = bad_irrQ = 0; maxh = 0; polys = {}
for d in ds:
    F = forms(d); h = len(F); maxh = max(maxh, h)
    js = [1728 * mp.kleinj(mp.mpc(mp.mpf(-b) / (2 * a), mp.sqrt(-d) / (2 * a))) for (a, b, c) in F]  # exact real part (a float here cost 1e-16, first run)
    if h > 1 and min(abs(js[i] - js[k]) for i in range(h) for k in range(i + 1, h)) < 1: bad_dist += 1
    co = [mp.mpc(1)]
    for j in js:
        co = [co[0]] + [co[i] - j * co[i - 1] for i in range(1, len(co))] + [-j * co[-1]]
    ints = [int(mp.nint(mp.re(c))) for c in co]
    dev = max(abs(c - n) for c, n in zip(co, ints))
    if dev > max(mp.mpf(10) ** -30, max(abs(c) for c in co) * mp.mpf(10) ** -140):  # relative to dps 150
        bad_int += 1
        if bad_int <= 3: print('  H1 debug d', d, 'h', h, 'dev', mp.nstr(dev, 5), 'maxcoef', mp.nstr(max(abs(c) for c in co), 5))
        continue
    P = Poly(sum(n * X ** (h - i) for i, n in enumerate(ints)), X, domain=ZZ)
    polys[d] = P
    fl = factor_list(P.as_expr())[1]
    if not (len(fl) == 1 and fl[0][1] == 1 and Poly(fl[0][0], X).degree() == h): bad_irrQ += 1
rep("H1 class polynomials H_d integral (coeffs within max(1e-30, 1e-140*max|coef|) of Z)", bad_int == 0, f"{len(ds)} fundamental d in [-399,-3], h <= {maxh}, failures {bad_int}")
rep("H2 j(a) = j(b) iff same class (distinct j over classes)", bad_dist == 0, f"{len(ds)} fields, failures {bad_dist}")
rep("H3 [Q(j):Q] = h (H_d irreducible over Q)", bad_irrQ == 0, f"{len(polys)} fields, failures {bad_irrQ}")
# irreducible over K = Q(sqrt d) for small h (finrank K K(j) = h)
bad = 0; n = 0
for d, P in polys.items():
    if 2 <= P.degree() <= 4 and -d < 300:
        fl = factor_list(P.as_expr(), extension=ssqrt(d))[1]; n += 1
        if len(fl) != 1: bad += 1
rep("H3b H_d irreducible over K (h in 2..4, |d|<300)", bad == 0, f"{n} fields, failures {bad}")
known = {-3: [1, 0], -4: [1, -1728], -7: [1, 3375], -8: [1, -8000], -163: [1, 640320 ** 3],
         -23: [1, 3491750, -5151296875, 12771880859375]}
ok = all(polys[d].all_coeffs() == c for d, c in known.items())
rep("H3c known values (j(rho)=0, j(i)=1728, -3375, 8000, -640320^3, H_-23)", ok, "")
# H4 CFT: p split in K, p not dividing disc(H_d): H_d mod p has h roots in F_p  <=>  p = principal form value
bad = 0; tested = 0
for d in [-23, -31, -39, -47, -55, -56, -59, -71, -79, -84, -87, -95, -104, -116]:
    if d not in polys: continue
    P = polys[d]; h = P.degree(); dc = int(discriminant(P))
    a0, b0, c0 = forms(d)[0]  # principal form (1, b, c)
    pr = set()
    for x in range(-120, 121):
        for y in range(0, 120):
            v = a0 * x * x + b0 * x * y + c0 * y * y
            if v < 3000: pr.add(v)
    for p in primerange(3, 3000):
        if d % p == 0 or dc % p == 0: continue
        if pow(d % p, (p - 1) // 2, p) != 1: continue  # p split in K
        co = [int(c) % p for c in P.all_coeffs()]
        nroots = sum(1 for t in range(p) if sum(cc * pow(t, h - i, p) for i, cc in enumerate(co)) % p == 0)
        tested += 1
        if (nroots == h) != (p in pr): bad += 1
rep("H4 split p: H_d splits completely mod p iff p = principal form", bad == 0, f"14 discriminants, {tested} (d,p) pairs p<3000, failures {bad}")
# R1 Radford: shuffle algebra = polynomial algebra on Lyndon words -> graded dimension identity
def lyndon(w): return all(w < w[i:] + w[:i] for i in range(1, len(w)))
ok = True
for a in (2, 3):
    Lk = [0] + [sum(1 for w in itertools.product(range(a), repeat=k) if lyndon(w)) for k in range(1, 11)]
    ser = [1] + [0] * 10
    for k in range(1, 11):
        for _ in range(Lk[k]):
            for t in range(k, 11): ser[t] += ser[t - k]
    ok &= ser == [a ** n for n in range(11)]
rep("R1 RadfordShufflePolynomial graded dims (|A|=2,3, length<=10)", ok, "prod over Lyndon words (1-t^|w|)^-1 = 1/(1-|A|t)")
