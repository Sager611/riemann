# hahn_toys.py — continuous Hahn coefficients d_2n for LP toys, with Romik's kernel applied to omega_Phi(x) = Phi(log(x)/2) x^(-1/4)
# (so that Xi_Phi(t) = int Phi(u) e^{itu} du; for Phi = omega(e^{2u}) e^{u/2} this is Romik's d_n).  Toys: atom u0 (cos(t u0)), gauss s.
import sys, mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre
kind, par, NMAX, out = sys.argv[1], mp.mpf(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
mp.mp.dps = 45
def kern(u, n):   # integrand factor per unit Phi(u) du (x = e^{2u} > 1), incl. Romik's prefactor and the factor 2 (even n)
    e = mp.exp(2*u); al = mp.tanh(u); m = 2*n
    return 2 * mp.rf(mp.mpf(3)/2, m) / (2**(m - mp.mpf(3)/2) * mp.factorial(m)) * e**(-mp.mpf(1)/4) * (e + 1)**(-mp.mpf(3)/2) * al**m * mp.hyp2f1(mp.mpf(m)/2 + mp.mpf(3)/4, mp.mpf(m)/2 + mp.mpf(5)/4, m + 2, al**2) * 2*e
if kind == 'atom':
    d = [kern(par, n) for n in range(NMAX + 1)]
    exact = lambda t: mp.cos(t*par)
else:
    base = GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec)
    E = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5]
    d = [mp.mpf(0)] * (NMAX + 1)
    for a, b in zip(E[:-1], E[1:]):
        a, b = mp.mpf(a), mp.mpf(b)
        for x, w in base:
            u = (b - a)/2*x + (a + b)/2; W = w*(b - a)/2 * mp.exp(-u**2/(2*par**2))
            for n in range(NMAX + 1): d[n] += W * kern(u, n)
    exact = lambda t: mp.quad(lambda u: mp.exp(-u**2/(2*par**2)) * mp.cos(t*u), [0, 2*par, 6*par, 12*par])
with open(out, 'w') as fh:
    for n in range(NMAX + 1): fh.write('%d %s\n' % (2*n, mp.nstr(d[n], 40)))
def gvals(nmax, x):
    v = [mp.mpf(1), 8*x/3]
    for n in range(1, nmax): v.append((8*x*v[n] - (2*n + 1)*v[n-1]) / (2*n + 3))
    return v
for t in [0, 1, 3]:
    t = mp.mpf(t); gv = gvals(2*NMAX, t/2)
    S = sum((-1)**n * d[n] * gv[2*n] for n in range(NMAX + 1))
    print('%s %s: t=%s exact=%s  CH sum=%s  ratio=%s' % (kind, mp.nstr(par, 3), mp.nstr(t, 3), mp.nstr(exact(t), 10), mp.nstr(S, 10), mp.nstr(exact(t)/S, 8)))
