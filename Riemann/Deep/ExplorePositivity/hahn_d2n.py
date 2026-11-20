# hahn_d2n.py — Romik's continuous Hahn expansion Xi(t) = sum (-1)^n d_2n g_2n(t/2),
# g_n = i^n (n+1) 3F2(-n, n+2, 3/4+ix; 3/2, 3/2; 1),  (2n+3) g_{n+1} = 8x g_n - (2n+1) g_{n-1};
# d_n = (3/2)_n/(2^{n-3/2} n!) int_0^inf omega(x)(x+1)^{-3/2} a^n 2F1(n/2+3/4, n/2+5/4; n+2; a^2) dx, a = (x-1)/(x+1)
# (even n: twice the x > 1 part).  Two Gauss-Legendre rules; validation against Xi.
import sys, time, mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre
from mp_core import omega, Xi
NMAX = int(sys.argv[1]); out = sys.argv[2]
EDGES = [0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3.0, 3.3]
def coeffs(deg, dps, ns=None):
    mp.mp.dps = dps
    ns = list(range(NMAX + 1)) if ns is None else ns
    base = GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec)
    acc = [mp.mpf(0)] * (NMAX + 1)
    for a, b in zip(EDGES[:-1], EDGES[1:]):
        a, b = mp.mpf(a), mp.mpf(b)
        for x, w in base:
            u = (b - a)/2*x + (a + b)/2
            e = mp.exp(2*u); al = mp.tanh(u)
            g = omega(e) * (e + 1)**(-mp.mpf(3)/2) * 2*e * w*(b - a)/2
            for n in ns:
                m = 2*n
                acc[n] += g * al**m * mp.hyp2f1(mp.mpf(m)/2 + mp.mpf(3)/4, mp.mpf(m)/2 + mp.mpf(5)/4, m + 2, al**2)
    return [2 * mp.rf(mp.mpf(3)/2, 2*n) / (2**(2*n - mp.mpf(3)/2) * mp.factorial(2*n)) * acc[n] if n in ns else None for n in range(NMAX + 1)]
t0 = time.time()
chk = sorted(set([0, 1, NMAX // 2, NMAX]))
dA = coeffs(5, 40); dB = coeffs(6, 55, chk)
mp.mp.dps = 55
rel = max(abs(dA[n] - dB[n]) / dB[n] for n in chk)
with open(out, 'w') as fh:
    for n in range(NMAX + 1): fh.write('%d %s\n' % (2*n, mp.nstr(dA[n], 38)))
print('NMAX', NMAX, 'rel diff (deg5,dps40) vs (deg6,dps55) at n =', chk, ':', mp.nstr(rel, 3), 'time %.0fs' % (time.time() - t0))
d = dA
asym = lambda n: 128 * mp.cbrt(2) * mp.pi**(mp.mpf(2)/3) * mp.exp(-2*mp.pi/3) / mp.sqrt(3) * n**(mp.mpf(4)/3) * mp.exp(-3 * mp.cbrt(4*mp.pi) * n**(mp.mpf(2)/3))
print('d0 d2 d4 =', [mp.nstr(v, 12) for v in d[:3]], ' d_2n/asym at n = 10, NMAX:', [mp.nstr(d[n]/asym(n), 5) for n in (10, NMAX) if n <= NMAX])
def gvals(nmax, x):
    v = [mp.mpf(1), 8*x/3]
    for n in range(1, nmax): v.append((8*x*v[n] - (2*n + 1)*v[n-1]) / (2*n + 3))
    return v
for t in [0, 5, 14.134725141734693790457, 30, 50]:
    t = mp.mpf(t); gv = gvals(2*NMAX, t/2)
    S = sum((-1)**n * d[n] * gv[2*n] for n in range(NMAX + 1))
    print('t=%s Xi=%s sum=%s ratio=%s' % (mp.nstr(t, 8), mp.nstr(Xi(t), 12), mp.nstr(S, 12), mp.nstr(Xi(t)/S, 10)))
