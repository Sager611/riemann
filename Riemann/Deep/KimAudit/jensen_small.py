#!/usr/bin/env python3
"""kimaudit T7: (a) normalisation of gamma(m) (Riemann's phi, (I1)) against mpmath's xi; (b) hyperbolicity of J^{d,n} for small n
at d = floor(D_n), floor(D_n)+1, 2 floor(D_n) (D_n from the theorem's formula; the theorem itself is claimed only for n >= n_*).
Hyperbolicity is certified by d+1 sign alternations of J^{d,n} on (-oo, 0) at points between the polyroots approximations.
Run: bash .scratch/memcap.sh 290 900 uv run --no-project --with mpmath python3 jensen_small.py"""
import time
import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre
mp.mp.dps = 110
mpf, pi = mp.mpf, mp.pi
a = mpf('0.37338196661320752'); t0 = time.time()
MMAX = 110
def phi(t):
    return sum(pi*k**2*mp.exp(5*t/2)*(2*pi*k**2*mp.exp(2*t) - 3)*mp.exp(-pi*k**2*mp.exp(2*t)) for k in range(1, 14))
X = GaussLegendre(mp.mp).calc_nodes(7, mp.mp.prec)
Mt = [mpf(0)]*(MMAX+1); p = mpf(0)
while p < mpf('3.2'):
    q = p + mpf('0.1'); hw, c = (q-p)/2, (q+p)/2
    for (x, w) in X:
        t = hw*x + c; val = hw*w*phi(t); t2 = t*t
        for m in range(MMAX+1):
            Mt[m] += val; val *= t2
    p = q
gam = [32*mp.factorial(m)*Mt[m]/mp.factorial(2*m) for m in range(MMAX+1)]
xi = lambda s: s*(s-1)*pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)/2
print('(a) gamma(m) = 32 m! Mt(m)/(2m)! vs 8 m! xi^(2m)(1/2)/(2m)!  [%.0f s]' % (time.time()-t0))
mp.mp.dps = 40
for m in range(0, 3):
    ref = 8*mp.factorial(m)*mp.diff(xi, mpf(1)/2, 2*m)/mp.factorial(2*m)
    print('    m = %d: moments %s, xi %s, ratio-1 = %s' % (m, mp.nstr(gam[m], 15), mp.nstr(ref, 15), mp.nstr(gam[m]/ref-1, 3)))
mp.mp.dps = 110
Th = lambda w: sum(gam[m]*w**m/mp.factorial(m) for m in range(MMAX+1))
for w in [mpf(-50) + 10j, -mp.zetazero(1).imag**2]:
    print('    Theta(w) = 8 xi(1/2 + sqrt w) at w = %s: Taylor %s, xi %s' % (mp.nstr(w, 8), mp.nstr(Th(w), 12), mp.nstr(8*xi(mpf(1)/2 + mp.sqrt(w)), 12)))
def Z(u): return u*(pi*mp.exp(2*u) - mpf(9)/4)
def uinv(x):
    lo, hi = mpf('1e-9'), mpf(50)
    for _ in range(200):
        mid = (lo+hi)/2
        if Z(mid) < x: lo = mid
        else: hi = mid
    return (lo+hi)/2
print('(b) hyperbolicity of J^{d,n} (all roots real iff d+1 sign alternations)')
for n in [3, 5, 8, 10]:
    u = uinv(mpf(n)); R = mpf(9)/10*(n-mpf(1)/2)**2/(u+1)**2; D = R**2/(4*a**2*(R+a**2)); d0 = int(mp.floor(D))
    for d in [d0, d0+1, 2*d0]:
        if n + d > MMAX: continue
        co = [mp.binomial(d, j)*gam[n+j] for j in range(d+1)]
        J = lambda x: mp.polyval(co[::-1], x)
        rts, err = mp.polyroots(co[::-1], maxsteps=2000, extraprec=600, error=True)
        re_ = sorted(mp.re(r) for r in rts); maxim = max(abs(mp.im(r))/abs(r) for r in rts)
        pts = [re_[0]*2] + [(re_[i]+re_[i+1])/2 for i in range(d-1)] + [re_[-1]/2]
        sg = [mp.sign(J(x)) for x in pts]
        alt = sum(1 for i in range(d) if sg[i]*sg[i+1] < 0)
        print('    n = %2d, D_n = %8s, d = %3d: max |Im r|/|r| = %s, sign alternations %d of %d -> %s  [%.0f s]' % (n, mp.nstr(D, 6), d, mp.nstr(maxim, 3), alt, d, 'HYPERBOLIC' if alt == d else 'NOT certified', time.time()-t0))
