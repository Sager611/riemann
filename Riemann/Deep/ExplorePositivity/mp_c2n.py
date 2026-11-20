# mp_c2n.py — Romik's MP coefficients c_2n of Xi (zeta) [or DH with arg 'dh'] by Gauss-Legendre on precomputed nodes,
# two rules (deg 6 @ dps 45, deg 7 @ dps 60) for a cross-check; validation of the expansion against Xi.
import sys, time, mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre
from mp_core import *
NMAX = int(sys.argv[1]); out = sys.argv[2]
EDGES = [0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3.0, 3.3]
def coeffs(deg, dps):
    mp.mp.dps = dps
    gl = GaussLegendre(mp.mp)
    base = gl.calc_nodes(deg, mp.mp.prec)
    acc = [mp.mpf(0)] * (NMAX + 1)
    for a, b in zip(EDGES[:-1], EDGES[1:]):
        a, b = mp.mpf(a), mp.mpf(b)
        for x, w in base:
            u = (b - a) / 2 * x + (a + b) / 2
            W = w * (b - a) / 2
            e = mp.exp(2*u)
            g = omega(e) * (e + 1)**(-mp.mpf(3)/2) * 2 * e * 4 * mp.sqrt(2) * W
            th2 = mp.tanh(u)**2
            p = g
            for n in range(NMAX + 1):
                acc[n] += p; p *= th2
    return acc
t0 = time.time()
HI = len(sys.argv) > 3 and sys.argv[3] == 'hi'
cA, cB = (coeffs(7, 100), coeffs(8, 130)) if HI else (coeffs(6, 45), coeffs(7, 60))
mp.mp.dps = 130 if HI else 60
rel = max(abs(cA[n] - cB[n]) / cB[n] for n in range(NMAX + 1))
c = cB
with open(out, 'w') as fh:
    for n in range(NMAX + 1):
        fh.write('%d %s\n' % (2*n, mp.nstr(c[n], 120 if HI else 55)))
print('NMAX', NMAX, 'max rel diff', '(deg7,dps100) vs (deg8,dps130):' if HI else '(deg6,dps45) vs (deg7,dps60):', mp.nstr(rel, 3), 'time %.1fs' % (time.time() - t0))
print('c0 c2 c4 c6 =', [mp.nstr(c[k], 12) for k in range(4)])
asym = lambda n: 16*mp.sqrt(2)*mp.pi**1.5*mp.sqrt(n)*mp.exp(-4*mp.sqrt(mp.pi*n))
print('c_2n / Romik asymptotic, n = 10, 50, 100, NMAX:', [mp.nstr(c[n]/asym(n), 6) for n in (10, 50, 100, NMAX) if n <= NMAX])
for t in [0, 5, 14.134725141734693790457, 21.022039638771554993, 30, 40]:
    t = mp.mpf(t)
    fv = f_eval(2*NMAX, t/2)
    S = sum((-1)**n * c[n] * fv[2*n] for n in range(NMAX + 1))
    print('t=%s  Xi=%s  Xi_N=%s' % (mp.nstr(t, 8), mp.nstr(Xi(t), 12), mp.nstr(S, 12)))
