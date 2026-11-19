#!/usr/bin/env python3
"""kimaudit T9: the mixture identity of D Sec. 4 at moderate n, Theta^{(n)}(w)/gamma(n) = int p_n(t) 0F1(; n+1/2; t^2 w/4) dt,
against the Taylor series sum_m gamma(n+m) w^m/(gamma(n) m!) (moments of phi by Gauss-Legendre), at points of W_{R_n}.
Run: bash .scratch/memcap.sh 290 900 uv run --no-project --with mpmath python3 mixture.py"""
import time
import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre
mpf, pi = mp.mpf, mp.pi
a = mpf('0.37338196661320752'); t0 = time.time()
def phi(t): return sum(pi*k**2*mp.exp(5*t/2)*(2*pi*k**2*mp.exp(2*t) - 3)*mp.exp(-pi*k**2*mp.exp(2*t)) for k in range(1, 12))
def Z(u): return u*(pi*mp.exp(2*u) - mpf(9)/4)
def uinv(x):
    lo, hi = mpf('1e-9'), mpf(50)
    for _ in range(200):
        mid = (lo+hi)/2
        if Z(mid) < x: lo = mid
        else: hi = mid
    return (lo+hi)/2
for n, dps, M, tlo, thi in [(60, 30, 220, '0.5', '2.3')]:   # (a first version with mp.quad at n = 60, 200 hit the 900 s cap)
    mp.mp.dps = dps; tlo, thi = mpf(tlo), mpf(thi)
    X = GaussLegendre(mp.mp).calc_nodes(7, mp.mp.prec)
    Mt = [mpf(0)]*(M+1); p = tlo
    while p < thi:
        q = p + mpf('0.1'); hw, c = (q-p)/2, (q+p)/2
        for (x, w) in X:
            t = hw*x + c; val = hw*w*t**(2*n)*phi(t); t2 = t*t
            for m in range(M+1):
                Mt[m] += val; val *= t2
        p = q
    cm = [mpf(1)]; r = mpf(1)
    for m in range(1, M+1):
        k = n + m - 1; r *= Mt[m]/(2*(2*k+1)*Mt[m-1]); cm.append(r/mp.factorial(m))
    Fw = lambda w: mp.polyval(cm[::-1], w)
    u = uinv(mpf(n)); R = mpf(9)/10*(n-mpf(1)/2)**2/(u+1)**2; vm = 2*a*mp.sqrt(R + a**2)
    b = n + mpf(1)/2
    X6 = GaussLegendre(mp.mp).calc_nodes(6, mp.mp.prec)
    def mix(w):   # same Gauss-Legendre panels (deg 6), integrand t^{2n} phi(t) 0F1(b; t^2 w/4)
        s = mpf(0); p = tlo
        while p < thi:
            q = p + mpf('0.1'); hw, c = (q-p)/2, (q+p)/2
            s += sum(hw*wt*(hw*x+c)**(2*n)*phi(hw*x+c)*mp.hyp0f1(b, (hw*x+c)**2*w/4) for (x, wt) in X6); p = q
        return s/Mt[0]
    for w in [-R + 1j*vm, -R/2 + 1j*mp.mpf('0.7')*2*a*mp.sqrt(R/2 + a**2)]:
        T, Mx = Fw(w), mix(w)
        print('n = %d, w = %s: Taylor %s | mixture %s | rel diff %s  [%.0f s]' % (n, mp.nstr(w, 8), mp.nstr(T, 12), mp.nstr(Mx, 12), mp.nstr(abs(T-Mx)/abs(T), 3), time.time()-t0))
