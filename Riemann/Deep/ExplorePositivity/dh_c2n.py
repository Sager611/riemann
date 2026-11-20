# dh_c2n.py — Davenport–Heilbronn analogue: F(s) = (5/pi)^((s+1)/2) Gamma((s+1)/2) f(s), f = sum a(n) n^-s,
# a(n) periodic mod 5 (1, k, -k, -1, 0), k = (sqrt(10-2 sqrt5)-2)/(sqrt5-1).  F(1/2+it) = int Phi_DH(u) e^{itu} du,
# Phi_DH(u) = 2 theta_a(e^{2u}) e^{3u/2}, theta_a(x) = sum n a(n) e^{-pi n^2 x/5}.  MP coefficients (up to one constant)
# c_2n = int_0^inf Phi_DH(u) cosh(u)^(-3/2) tanh(u)^(2n) du ; checks: Phi_DH even, sign, Xi_N/F constant, F at Spira's zero.
import sys, time, mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre
from mp_core import f_eval
NMAX = int(sys.argv[1]); out = sys.argv[2]
def setup():
    global K, A
    K = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
    A = [0, 1, K, -K, -1]
def theta(x, M):
    return sum(n * A[n % 5] * mp.exp(-mp.pi * n**2 * x / 5) for n in range(1, M + 1))
def Phi(u, M=80):
    return 2 * theta(mp.exp(2*u), M) * mp.exp(mp.mpf(3)*u/2)
def F(t):
    s = mp.mpf(1)/2 + 1j*t
    f = 5**(-s) * sum(A[r] * mp.zeta(s, mp.mpf(r)/5) for r in range(1, 5))
    return (5/mp.pi)**((s+1)/2) * mp.gamma((s+1)/2) * f
EDGES = [0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3.0, 3.3, 3.7]
def coeffs(deg, dps):
    mp.mp.dps = dps; setup()
    base = GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec)
    acc = [mp.mpf(0)] * (NMAX + 1)
    for a, b in zip(EDGES[:-1], EDGES[1:]):
        a, b = mp.mpf(a), mp.mpf(b)
        for x, w in base:
            u = (b - a)/2 * x + (a + b)/2
            g = Phi(u, 30 if u > 0.5 else 80) * mp.cosh(u)**(-mp.mpf(3)/2) * w * (b - a)/2
            th2 = mp.tanh(u)**2; p = g
            for n in range(NMAX + 1):
                acc[n] += p; p *= th2
    return acc
t0 = time.time()
mp.mp.dps = 40; setup()
print('kappa =', mp.nstr(K, 15))
print('evenness Phi(u)/Phi(-u) at u = 0.2, 0.5, 0.9:', [mp.nstr(Phi(mp.mpf(u))/Phi(-mp.mpf(u), 400), 15) for u in (0.2, 0.5, 0.9)])
grid = [mp.mpf(k)/100 for k in range(0, 361)]
vals = [Phi(u) for u in grid]
print('Phi_DH on [0, 3.6] step 0.01: min = %s at u = %s; Phi(0) = %s' % (mp.nstr(min(vals), 6), mp.nstr(grid[vals.index(min(vals))], 4), mp.nstr(vals[0], 8)))
print('F(1/2+it) imaginary/real at t = 3, 40:', [mp.nstr(abs(F(t).imag / F(t).real), 3) for t in (3, 40)])
z0 = mp.findroot(lambda t: F(t), mp.mpc('85.699348', '-0.308517'))
print('Spira zero: t =', mp.nstr(z0, 12), ' -> rho = 1/2 + i t =', mp.nstr(mp.mpf(1)/2 + 1j*z0, 12))
cA = coeffs(6, 45); cB = coeffs(7, 60)
mp.mp.dps = 60; setup()
rel = max(abs(cA[n] - cB[n]) / abs(cB[n]) for n in range(NMAX + 1))
c = cB
with open(out, 'w') as fh:
    for n in range(NMAX + 1):
        fh.write('%d %s\n' % (2*n, mp.nstr(c[n], 55)))
print('NMAX', NMAX, 'max rel diff of the two rules:', mp.nstr(rel, 3), ' min c_2n =', mp.nstr(min(c), 6), ' time %.1fs' % (time.time() - t0))
print('c0 c2 c4 =', [mp.nstr(c[k], 12) for k in range(3)], ' C1 ratio c2^2/(c0 c4) =', mp.nstr(c[1]**2/(c[0]*c[2]), 6))
for t in [0, 3, 10, 20, 30, 40]:
    t = mp.mpf(t)
    fv = f_eval(2*NMAX, t/2)
    S = sum((-1)**n * c[n] * fv[2*n] for n in range(NMAX + 1))
    print('t=%s  F=%s  Xi_N=%s  ratio=%s' % (mp.nstr(t, 4), mp.nstr(F(t).real, 12), mp.nstr(S, 12), mp.nstr(F(t).real / S, 10)))
