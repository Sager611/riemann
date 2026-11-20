"""C4 model test: do neighbouring zeros heal an off-line pair in the x-family Jensen polynomials J^{d,0}?
F(x) = prod over zeros zeta of (1 + x/zeta^2) (off-line zeta come with their conjugates), J^{d,0}(X) = sum_j C(d,j) g_j X^j,
g_j = j! [x^j]F.  d_max = largest d with J^{d,0} real-rooted (monotone in d: J*^{d,0} = D J*^{d+1,0}/(d+1), Rolle).
Units: gamma0 = 10^4, Delta = 2pi/log(gamma0/2pi) ~ natural spacing, beta = 1.6 Delta (xi's ratio at T = 3e12: 0.373/0.234),
so Im w/Delta_w = 1.6 and the threshold heat time 1/(2d) ~ 1.3 Delta_w^2, as for xi at T.  Exact rationals, flint isolation."""
from fractions import Fraction as Fr
import flint, math, sys, time
G0 = Fr(10**4); DEL = Fr(85, 100); BETA = Fr(136, 100)     # 2pi/log(1e4/2pi) = 0.8526; beta = 1.6*0.85
def poly_mul(a, b):
    c = [Fr(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b): c[i+j] += x*y
    return c
def build(real_heights, pair_heights, beta=BETA):
    P = [Fr(1)]
    for g in real_heights: P = poly_mul(P, [Fr(1), 1/(g*g)])
    for g in pair_heights:            # (1 + x/z^2)(1 + x/zbar^2) = 1 + 2Re(1/z^2) x + x^2/|z|^4, z = g + i beta
        re = g*g - beta*beta; im = 2*g*beta; ab2 = re*re + im*im           # z^2 = re + i im
        P = poly_mul(P, [Fr(1), 2*re/ab2, 1/ab2])
    return [P[j]*math.factorial(j) for j in range(len(P))]
def real_rooted(g, d):
    n = min(d, len(g)-1)
    co = [math.comb(d, j)*g[j] for j in range(n+1)]
    den = 1
    for v in co: den = den*v.denominator//math.gcd(den, v.denominator)
    p = flint.fmpz_poly([int(v*den) for v in co])
    return sum(m for r, m in p.complex_roots() if r.imag == 0) == p.degree()
def dmax(g, hi):
    lo = 1
    assert not real_rooted(g, hi)
    while hi - lo > max(1, lo//20000):
        mid = (lo+hi)//2
        if real_rooted(g, mid): lo = mid
        else: hi = mid
    return lo
if __name__ == '__main__':
    M = int(sys.argv[1]); K = int(sys.argv[2])
    th = 2*math.atan(float(BETA/G0)); dA = 1/math.sin(th)**2
    print('isolated-pair prediction floor(1/sin^2 theta) =', math.floor(dA))
    H = Fr(742, 100)   # S(T)-worst gap: count 2 (the pair) vs 2h expected needs 2h - 2 <= 2*6.42 (Platt-Trudgian |S| at T = 3e12)
    cfg = {
      'A isolated pair': ([], [G0]),
      'B pair, real neighbours both sides (spacing Delta)': ([G0 + m*DEL for m in range(-M, M+1) if m != 0], [G0]),
      'C lattice of off-line pairs (spacing 2 Delta: 2 zeros each)': ([], [G0 + 2*k*DEL for k in range(-K, K+1)]),
      'D pair in S(T)-worst gap (half-width 7.42 Delta), real beyond': ([G0 + s*(H + m)*DEL for s in (-1, 1) for m in range(M)], [G0]),
      'D2 same gap, off-line pairs beyond (spacing 2 Delta)': ([], [G0] + [G0 + s*(H + 2*m)*DEL for s in (-1, 1) for m in range(K)]),
      'E real below (spacing Delta), pairs from G0 up (2 Delta)': ([G0 - m*DEL for m in range(1, M+1)], [G0 + 2*k*DEL for k in range(K)]),
    }
    for name, (rh, ph) in cfg.items():
        t0 = time.time(); g = build(rh, ph)
        dm = dmax(g, int(40*dA))
        print('%-55s zeros %3d  d_max = %d  ratio to isolated %.3f  (%.0f s)' % (name, len(rh)+2*len(ph), dm, dm/dA, time.time()-t0), flush=True)
