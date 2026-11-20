# herm_route.py — Turan's Hermite route (Polya 1915 / Turan 1959 Lemma II): Xi(t) = sum (-1)^n b_2n H_2n(t),
# b_n = int Phi(x) x^n e^{-x^2/4} dx / (2^n n!), Phi = Fourier density of Xi.  Pre-image p_N(z) = sum_{n<=N} (-1)^n b_2n z^2n
# real-rooted => Hermite approximant real-rooted.  Roots of p_N and of the approximants isolated by Arb.
import sys, mpmath as mp
from flint import arb, acb, acb_poly, ctx
from mp_core import Xi
mp.mp.dps = 60; ctx.prec = 700
NMAX = int(sys.argv[1])
def PhiT(u, M=12):
    u = abs(u)
    return sum((2*mp.pi**2*n**4*mp.exp(9*u/2) - 3*mp.pi*n**2*mp.exp(5*u/2)) * mp.exp(-mp.pi*n**2*mp.exp(2*u)) for n in range(1, M+1))
E = [0, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.3, 1.7, 2.2]
print('Xi(0) =', mp.nstr(Xi(0), 15), ' 2 int_0^inf PhiT =', mp.nstr(2*mp.quad(PhiT, E), 15))
b = [2*mp.quad(lambda x: PhiT(x) * x**(2*n) * mp.exp(-x**2/4), E) / (4**n * mp.factorial(2*n)) for n in range(NMAX + 1)]
print('b0 b2 b4 =', [mp.nstr(v, 12) for v in b[:3]], ' Newton k=1 ratio b2^2/(b0 b4) =', mp.nstr(b[1]**2/(b[0]*b[2]), 6), '(needed >= 2N/(N-1) for real p_N)')
H = [[1], [0, 2]]
for n in range(1, 2*NMAX):
    a, c = H[n], H[n-1]
    h = [0]*(n+2)
    for k, v in enumerate(a): h[k+1] += 2*v
    for k, v in enumerate(c): h[k] -= 2*n*v
    H.append(h)
for t in [0, 5, 14.134725141734693790457]:
    t = mp.mpf(t)
    S = sum((-1)**n * b[n] * sum(cf * t**k for k, cf in enumerate(H[2*n])) for n in range(NMAX + 1))
    print('t=%s Xi=%s Hermite sum=%s' % (mp.nstr(t, 8), mp.nstr(Xi(t), 12), mp.nstr(S, 12)))
def count(coef):
    P = acb_poly([acb(arb(mp.nstr(c, 55))) for c in coef])
    rts = P.roots(maxprec=2800)
    re = [r for r in rts if r.imag.contains(0)]; nr = [r for r in rts if not r.imag.contains(0)]
    return len(re), len(nr), (min(float(abs(r.imag.mid())) for r in nr) if nr else None), sorted(float(abs(r.real.mid())) for r in re)[::2][:4]
for N in [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40]:
    if N > NMAX: break
    pre = [mp.mpf(0)]*(2*N+1)
    for n in range(N+1): pre[2*n] = (-1)**n * b[n]
    app = [mp.mpf(0)]*(2*N+1)
    for n in range(N+1):
        for k, cf in enumerate(H[2*n]): app[k] += (-1)**n * b[n] * cf
    r1 = count(pre); r2 = count(app)
    print('N=%2d  pre-image p_N: real %d nonreal %d (min|Im| %s)   Hermite approximant: real %d nonreal %d (min|Im| %s) first real t>0 %s' % (
        N, r1[0], r1[1], r1[2] and '%.3g' % r1[2], r2[0], r2[1], r2[2] and '%.3g' % r2[2], ['%.4f' % v for v in r2[3]]), flush=True)
