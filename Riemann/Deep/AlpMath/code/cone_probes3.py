# Probes M1 (AbMulNSt: #E[n] = n^2 via separable division polynomials), P1 (PartialZetaFunctionalEquation, h = 1 fields),
# U1 (ChevalleyUnits: congruence units are m-th powers, real quadratic instances)
import sys, math, random
import mpmath as mp
from sympy import Poly, symbols, GF, gcd as sgcd, primerange, factorint, isprime
random.seed(5)
def rep(tag, ok, msg): print(f"{tag}: {'PASS' if ok else 'FAIL'} {msg}"); sys.stdout.flush()
x = symbols('x')
# M1: division polynomials psi_n for y^2 = x^3 + a x + b, n odd: psi_n is a polynomial in x of degree (n^2-1)/2
def divpolys(a, b, nmax):
    # standard recursion with psi_2 = 2y; work with f = x^3+ax+b and y^2 = f; for odd n psi_n in Z[x]
    f = x**3 + a*x + b
    psi = {0: 0, 1: 1, 2: 2, 3: 3*x**4 + 6*a*x**2 + 12*b*x - a**2,
           4: 4*(x**6 + 5*a*x**4 + 20*b*x**3 - 5*a**2*x**2 - 4*a*b*x - 8*b**2 - a**3)}  # psi_2, psi_4 divided by y
    # represent psi_n = y * g_n for even n (store g_n), psi_n = g_n for odd n
    for n in range(5, nmax + 1):
        m = n // 2
        if n % 2 == 1:
            # psi_{2m+1} = psi_{m+2} psi_m^3 - psi_{m-1} psi_{m+1}^3 ; y^4 = f^2 where even factors appear
            if m % 2 == 0:
                val = psi[m + 2] * psi[m] ** 3 * f**2 - psi[m - 1] * psi[m + 1] ** 3
            else:
                val = psi[m + 2] * psi[m] ** 3 - psi[m - 1] * psi[m + 1] ** 3 * f**2
            psi[n] = val
        else:
            # psi_{2m} = psi_m (psi_{m+2} psi_{m-1}^2 - psi_{m-2} psi_{m+1}^2) / (2y)
            val = psi[m] * (psi[m + 2] * psi[m - 1] ** 2 - psi[m - 2] * psi[m + 1] ** 2)
            psi[n] = Poly(val, x).quo(Poly(2, x)).as_expr()
    return psi
bad = 0; cnt = 0
for p in primerange(5, 60):
    for _ in range(4):
        a, b = random.randrange(p), random.randrange(p)
        if (4 * a**3 + 27 * b**2) % p == 0: continue
        ps = divpolys(a, b, 7)
        for n in (3, 5, 7):
            if p % n == 0: continue
            P = Poly(ps[n], x, modulus=p)
            deg_ok = P.degree() == (n * n - 1) // 2
            sqf = Poly(sgcd(P, P.diff(x)), x, modulus=p).degree() == 0
            cnt += 1
            if not (deg_ok and sqf): bad += 1
rep("M1 AbMulNSt (g=1): #E[n] = n^2, n = 3,5,7", bad == 0, f"{cnt} (curve, p, n) cases, p < 60: psi_n mod p separable of degree (n^2-1)/2, failures {bad}")
# P1: completed partial zeta functional equation for F with h = 1 (C = C' = 1, different principal): d^(s/2) G_R(s)^n Z(s) = d^((1-s)/2) G_R(1-s)^n Z(1-s)
mp.mp.dps = 30
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
GR = lambda s: mp.pi ** (-s / 2) * mp.gamma(s / 2)
def zetaF(d, s):
    L = d ** (-s) * mp.fsum(kron(d, a) * mp.zeta(s, mp.mpf(a) / d) for a in range(1, d + 1))
    return mp.zeta(s) * L
err = 0
for d in (5, 8, 13, 17, 29):  # class number 1 real quadratic fields with a unit of norm -1
    for _ in range(4):
        s = mp.mpc(random.uniform(-1.5, 2.5), random.uniform(-8, 8))
        lhs = mp.mpf(d) ** (s / 2) * GR(s) ** 2 * zetaF(d, s)
        rhs = mp.mpf(d) ** ((1 - s) / 2) * GR(1 - s) ** 2 * zetaF(d, 1 - s)
        err = max(err, abs(lhs - rhs) / max(abs(lhs), mp.mpf(1e-30)))
for _ in range(4):  # F = Q
    s = mp.mpc(random.uniform(-1.5, 2.5), random.uniform(-8, 8))
    err = max(err, abs(GR(s) * mp.zeta(s) - GR(1 - s) * mp.zeta(1 - s)) / abs(GR(s) * mp.zeta(s)))
rep("P1 PartialZetaFunctionalEquation (h=1 instances)", err < 1e-20, f"F = Q, Q(sqrt d) d = 5,8,13,17,29 at 24 random s, max rel err {mp.nstr(err,3)}")
# U1 ChevalleyUnits: E = Q(sqrt D), units +-eps^k; find a coprime to b with: u = 1 mod a  =>  u = v^m
def unit_mod(D, a):
    # fundamental unit eps = (x + y sqrt D)/2 with x^2 - D y^2 = +-4 (D fundamental discriminant); work with u = (X + Y sqrt D)/2, integrality
    y = 1
    while True:
        for sg in (-4, 4):
            t = D * y * y + sg
            if t > 0 and math.isqrt(t) ** 2 == t: return (math.isqrt(t), y)
        y += 1
def mulq(u, v, D):  # (x1 + y1 r)/2 * (x2 + y2 r)/2 = ((x1x2 + D y1y2)/2 + (x1y2 + x2y1)/2 r)/2
    (x1, y1), (x2, y2) = u, v
    return ((x1 * x2 + D * y1 * y2) // 2, (x1 * y2 + x2 * y1) // 2)
def cong1(u, a, D):  # u - 1 in a O_E, O_E = Z[(D + sqrt D)/2]; u = (X + Y sqrt D)/2, u - 1 = ((X - 2) + Y sqrt D)/2
    X, Y = u
    # element (P + Q sqrt D)/2 is in a*O_E iff (P/a + (Q/a) sqrt D)/2 in O_E iff a | P, a | Q and (P/a - Q/a*D) even... use: (P + Q sqrt D)/2 in O_E iff P = Q D mod 2
    P, Q = X - 2, Y
    if P % a or Q % a: return False
    return (P // a - (Q // a) * D) % 2 == 0
ok = True; found = []
for D in (5, 8, 12, 13, 21):
    eps = unit_mod(D, 1)
    for m in (2, 3, 5):
        for b in (1, 2, 3, 5):
            hit = None
            for a in range(1, 400):
                if math.gcd(a, b) != 1: continue
                # units +-eps^k, k in [-K, K] (period of eps mod a divides |(O/a)^x| <= a^2); check k in 1..2a^2 and sign
                good = True; u = (2, 0)
                if cong1((-2, 0), a, D) and m % 2 == 0: good = False  # k = 0: -1 = 1 mod a but -1 not a square
                for k in range(1, 2 * a * a + 1):
                    u = mulq(u, eps, D)
                    for sgn in (1, -1):
                        uu = (sgn * u[0], sgn * u[1])
                        if cong1(uu, a, D):
                            # +-eps^k is an m-th power iff m | k and (sign +1 or m odd)
                            if not (k % m == 0 and (sgn == 1 or m % 2 == 1)): good = False; break
                    if not good: break
                if good: hit = a; break
            if hit is None: ok = False
            found.append((D, m, b, hit))
rep("U1 ChevalleyUnits (real quadratic instances)", ok, f"60 cases (D in 5,8,12,13,21; m in 2,3,5; b in 1,2,3,5): a found for all; e.g. {found[:4]}")
