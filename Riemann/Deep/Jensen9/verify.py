#!/usr/bin/env python3
"""verify.py -- re-runs the certificates of Riemann/Deep/Jensen9/REPORT.md (node N6, Jensen polynomials).

Part A (pure Python, exact): the degree bound
    D = floor(1/sin^2(2*arctan(a/T))) = floor((1+u^2)^2/(4u^2)),  u = a/T,
for T = 3,000,175,332,800 (Platt-Trudgian: RH true up to height T) and a = 3/8 (OpenAI 2026: no zero with
Re s > 7/8, hence |Re rho - 1/2| <= 3/8), plus the old bound a = 1/2 (trivial strip; GORTTW Thm 1.2 /
Kim-Lee 2021), the conditional half-planes, and the Xi-family bound 1 + T^2/a^2 (Kim-Lee Thm 4 with strip a).
Part B (exact): the sector geometry: every zero z = x+iy with |x| <= a, |y| > T gives w = z^2 with
    sin^2(angle(w, negative axis)) = Y^2/(X^2+Y^2) <= sin^2(2 arctan(a/T)).
Parts C-F need python-flint (pip install python-flint; Arb ball arithmetic):
C: FLINT/Arb recomputes N(T) at Platt-Trudgian's height (Turing's method): must be 12,363,153,437,138.
D: Obreschkoff/Chasse sector lemma spot checks with exact integer polynomials (certified root isolation):
   sharpness of d <= 1/sin^2(delta) for an edge pair; random polynomials with zeros in S(pi,delta) have
   hyperbolic Jensen polynomials up to degree floor(1/sin^2 delta); Gauss-Lucas closure of the sector.
E: consistency on the actual xi data: rigorous gamma(n) from Arb power series; J^{9,n} certified hyperbolic
   (9 sign changes at exact rational points) for 0 <= n <= 150, and J^{d,0} for d <= 40.
F: FLINT isolates zero #12,363,153,437,139 (the first above T) on the critical line; with it RH holds up to
   T' = 3000175332800.2865 and the bound rises to D' (secondary: depends on FLINT's Turing-method code).
Runtime (Apple M-series laptop, one thread): about 35 s, peak RSS about 76 MB (C: 11 s, F: 16 s).
"""
import sys, time, random
from fractions import Fraction as Fr

T_PT = 3000175332800
N_PT = 12363153437138

def D_bound(a, T):
    u = Fr(a) / T
    v = (1 + u * u) ** 2 / (4 * u * u)
    return v.numerator // v.denominator, v

def part_A():
    print('== A. degree bounds (exact rational arithmetic)')
    D_old, v_old = D_bound(Fr(1, 2), T_PT)
    D_new, v_new = D_bound(Fr(3, 8), T_PT)
    assert D_old == T_PT ** 2, 'strip-1/2 bound must equal T^2 (T integer)'
    assert (16 * T_PT ** 2 - 4) % 9 == 0 and D_new == (16 * T_PT ** 2 - 4) // 9
    print(f'  old (strip |Re z|<1/2, GORTTW Thm1.2 / Kim-Lee): D0 = T^2 = {D_old}  (~{D_old:.6e})')
    print(f'  new (7/8 half-plane, |Re z|<=3/8):               D  = (16T^2-4)/9 = {D_new}  (~{D_new:.6e})')
    print(f'  ratio D/D0 = {D_new / D_old:.9f}  (16/9 = {16/9:.9f})')
    for name, a in [('Liu 2610.12234 (0.874957..., conditional)', Fr(374957, 10 ** 6)),
                    ('Chang-Wang 2611.00924 (0.87446205769353, conditional)', Fr(37446205769353, 10 ** 14))]:
        print(f'  {name}: D = {D_bound(a, T_PT)[0]:.6e}')
    for name, a in [('Xi-family, strip 1/2 (Kim-Lee 2021)', Fr(1, 2)), ('Xi-family, a = 3/8 (new)', Fr(3, 8))]:
        v = 1 + Fr(T_PT ** 2) / (a * a)
        print(f'  {name}: d <= 1 + T^2/a^2 = {v.numerator // v.denominator}  (~{float(v):.6e})')
    return D_new

def part_B():
    print('== B. sector geometry (exact)')
    a, T = Fr(3, 8), Fr(T_PT)
    # extreme zero z = a + iT -> w = z^2 = X + iY
    X, Y = a * a - T * T, 2 * a * T
    s2 = Y * Y / (X * X + Y * Y)
    u = a / T
    assert s2 == (2 * u / (1 + u * u)) ** 2, 'sin^2 of the edge angle'
    assert X < 0 and Y > 0 and Y < -X, 'angle < pi/4 < pi/2'
    # monotonicity: for |x| <= a, |y| >= T the ratio |x|/|y| <= a/T, and t -> sin^2(2 arctan t) is increasing on [0,1]
    rnd = random.Random(1)
    for _ in range(2000):
        x = Fr(rnd.randint(-3000, 3000), 8000)            # |x| <= 3/8
        y = T + Fr(rnd.randint(0, 10 ** 9), rnd.randint(1, 1000))
        Xw, Yw = x * x - y * y, 2 * x * y
        assert Yw * Yw / (Xw * Xw + Yw * Yw) <= s2
    print('  edge zero 3/8 + iT: 1/sin^2 =', float(1 / s2), '; 2000 random zeros in the region: all inside the sector')

def part_C():
    from flint import arb, ctx
    print('== C. FLINT/Arb: N(T) at the Platt-Trudgian height (Turing method)')
    ctx.prec = 64
    t = time.time()
    N = arb(T_PT).zeta_nzeros()
    ok = (N == N_PT)
    print(f'  N({T_PT}) = {N}  expected {N_PT}  -> {"OK" if ok else "MISMATCH"}  ({time.time() - t:.1f} s)')
    assert ok

def jensen_int(q, d):
    """Jensen polynomial sum_j C(d,j) q^(j)(0) X^j of an integer polynomial q (list of coeffs, low first)."""
    from math import comb, factorial
    return [comb(d, j) * factorial(j) * (q[j] if j < len(q) else 0) for j in range(min(d, len(q) - 1) + 1)]

def n_real_roots(coeffs):
    from flint import fmpz_poly
    rr = fmpz_poly(coeffs).complex_roots()
    tot = sum(m for _, m in rr)
    real = sum(m for c, m in rr if c.imag == 0)     # Arb returns real roots with exactly zero imaginary part
    return real, tot

def polymul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i + j] += x * y
    return r

def part_D():
    from flint import fmpz_poly, ctx, arb
    ctx.prec = 128
    print('== D. Obreschkoff/Chasse sector lemma: spot checks (exact integer polynomials, certified roots)')
    # delta with tan(delta) = 1/5: sin^2 = 1/26, bound d <= 26.  Edge pair -5 +- i: q = w^2 + 10 w + 26.
    q = [26, 10, 1]
    for d, want in [(26, True), (27, False)]:
        J = jensen_int(q, d); r, t = n_real_roots(J)
        assert (r == t) == want
    print('  edge pair -5+-i (tan delta = 1/5): J^26 hyperbolic, J^27 not  -> bound 1/sin^2 = 26 is sharp')
    rnd = random.Random(2026)
    trials, worst = 0, None
    for _ in range(300):
        p = [1]
        for _k in range(rnd.randint(0, 6)):           # real negative zeros -b
            b = rnd.randint(1, 60)
            p = polymul(p, [b, 1])
        npairs = rnd.randint(1, 4)
        for _k in range(npairs):                      # pairs -u +- iv with v/u <= 1/5 (edge allowed)
            u = rnd.randint(1, 40)
            v = rnd.randint(0, u // 5) if rnd.random() < 0.6 else u // 5
            p = polymul(p, [u * u + v * v, 2 * u, 1])
        for d in [26, rnd.randint(2, 26)]:
            J = jensen_int(p, d); r, t = n_real_roots(J)
            assert r == t, ('counterexample?!', p, d)
            trials += 1
        # Gauss-Lucas: zeros of p' stay in the closed sector |arg(-w)| <= delta
        dp = [k * p[k] for k in range(1, len(p))]
        if len(dp) > 1:
            for c, m in fmpz_poly(dp).complex_roots():
                inside = -c.real - 5 * abs(c.imag)              # >= 0  <=>  |arg(-w)| <= delta
                assert inside + arb(10) ** -20 > 0, ('Gauss-Lucas violated?!', p)
    print(f'  {trials} Jensen polynomials (300 random polys, zeros in S(pi,delta), d <= 26): all hyperbolic;'
          ' derivative zeros stay in the sector')

def gammas(N, prec):
    from flint import arb, arb_series, ctx
    ctx.prec = prec; ctx.cap = 2 * N + 2
    s = arb_series([arb(1) / 2, 1], prec=2 * N + 2)               # s = 1/2 + z
    xi = (s * (s - 1) / 2) * (-(s / 2) * arb.pi().log()).exp() * arb_series.gamma(s / 2) * arb_series.zeta(s)
    c = xi.coeffs(); out = []; fact = arb(1)
    for n in range(N + 1):
        if n: fact *= n
        out.append(8 * fact * c[2 * n])                            # (-1+4z^2) Lambda(1/2+z) = 8 xi(1/2+z)
    return out

def certify_hyperbolic(coef):
    """coef: arb list, degree d, positive leading coeff. Certify d distinct real roots via strict sign
    alternation at d-1 exact dyadic points between approximate roots (plus the signs at -inf/+inf).
    The approximate roots only place the test points; the certificate is the ball-arithmetic sign check."""
    from flint import arb, acb, acb_poly, arb_poly
    d = len(coef) - 1
    roots = acb_poly([acb(c.mid()) for c in coef]).roots()       # approximate roots, full working precision
    rs = sorted((r.real.mid() for r in roots), key=lambda x: float(x))
    P = arb_poly(coef)
    signs = []
    for i in range(d - 1):
        p = arb(((arb(rs[i]) + arb(rs[i + 1])) / 2).mid())       # exact dyadic test point
        v = P(p)
        if v > 0: signs.append(1)
        elif v < 0: signs.append(-1)
        else: return False
    seq = [(-1) ** d] + signs + [1]                              # sign at -inf, interior points, +inf
    return all(seq[i] == -seq[i + 1] for i in range(len(seq) - 1))

def part_E():
    from flint import arb, ctx
    from math import comb
    print('== E. actual xi data: rigorous gamma(n), certified hyperbolicity of J^{9,n} and J^{d,0}')
    t = time.time()
    N = 160
    g = gammas(N, 4400)                                           # cancellation: needs ~4000 bits at n=160
    assert abs(float(g[0].mid()) - 3.9769662255065128793) < 1e-15
    ok9 = []
    for n in range(0, 151):
        coef = [comb(9, j) * g[n + j] / g[n] for j in range(10)]
        ok9.append(certify_hyperbolic(coef))
    assert all(ok9), [n for n, o in enumerate(ok9) if not o]
    okd = all(certify_hyperbolic([comb(d, j) * g[j] / g[0] for j in range(d + 1)]) for d in range(2, 41))
    assert okd
    print(f'  gamma(0) = {g[0].str(20)}; J^(9,n) hyperbolic for 0<=n<=150 and J^(d,0) for 2<=d<=40'
          f' (certified)  ({time.time() - t:.1f} s)')

def part_F():
    from flint import acb, arb, ctx
    print('== F. one zero above the Platt-Trudgian height (FLINT zeta_zeros: Turing-method isolation)')
    ctx.prec = 64
    t = time.time()
    z = acb.zeta_zeros(N_PT + 1, 1)[0]
    assert z.real == arb(1) / 2 and z.imag > T_PT
    lo = Fr(30001753328002865, 10 ** 4)
    assert z.imag > arb(30001753328002865) / 10 ** 4          # rigorous lower bound used below
    Dp = D_bound(Fr(3, 8), lo)[0]
    print(f'  zero #{N_PT + 1} = 1/2 + i{z.imag.str(18)} on the critical line ({time.time() - t:.1f} s)')
    print(f"  RH up to T' = 3000175332800.2865 => D' = {Dp} (D' - D = {Dp - D_bound(Fr(3, 8), T_PT)[0]})")

if __name__ == '__main__':
    t0 = time.time()
    part_A(); part_B()
    try:
        import flint  # noqa: F401
    except ImportError:
        print('python-flint not installed: parts C-F skipped (pip install python-flint)'); sys.exit(0)
    part_C(); part_D(); part_E(); part_F()
    print(f'ALL CHECKS PASSED in {time.time() - t0:.1f} s')
