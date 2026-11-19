#!/usr/bin/env python3
"""kimaudit T6: zeros of Theta^{(n)} near the claimed zero-free region W_R = {(Im w)^2 <= 4a^2(a^2 - Re w), Re w >= -R}
at moderate n, from the Taylor coefficients gamma(n+m) (moments of Riemann's phi by Gauss-Legendre panels), NOT from
the mixture identity. R = R_n := (9/10)(n-1/2)^2/(u(n)+1)^2 (the theorem's formula, claimed only for n >= n_*).
(1) argument principle on the boundary of W_{R_n}; (2) first sign change of Theta^{(n)}(-s), s > R_n;
(3) argument principle on W_{0.999 s1}. Run: bash .scratch/memcap.sh 290 1500 uv run --no-project --with mpmath python3 theta_n.py n dps M"""
import sys, time
import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre
n, dps, M = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
mp.mp.dps = dps
mpf, pi = mp.mpf, mp.pi
a = mpf('0.37338196661320752')
t0 = time.time()
def Z(u): return u*(pi*mp.exp(2*u) - mpf(9)/4)
def uinv(x):
    lo, hi = mpf('1e-9'), mpf(50)
    for _ in range(200):
        mid = (lo+hi)/2
        if Z(mid) < x: lo = mid
        else: hi = mid
    return (lo+hi)/2
def logphi(t):
    s = mpf(0)
    for k in range(1, 9):
        s += pi*k**2*mp.exp(5*t/2)*(2*pi*k**2*mp.exp(2*t) - 3)*mp.exp(-pi*k**2*mp.exp(2*t))
    return mp.log(s)
def logint(m, t): return 2*m*mp.log(t) + logphi(t)
# integration range: drop where log integrand < peak - (dps*2.31 + 40) for m in {n, n+M}
cut = dps*mpf('2.31') + 40
def peak(m): return logint(m, uinv(mpf(m)))
tlo = mpf('0.05')
while logint(n, tlo + mpf('0.05')) < peak(n) - cut: tlo += mpf('0.05')
thi = mpf(5)
while logint(n+M, thi - mpf('0.05')) < peak(n+M) - cut: thi -= mpf('0.05')
def moments(deg, width):
    X = GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec)
    Mt = [mpf(0)]*(M+1)
    p = tlo
    while p < thi:
        q = min(p + width, thi); hw, c = (q-p)/2, (q+p)/2
        for (x, w) in X:
            t = hw*x + c
            val = hw*w*mp.exp(logint(n, t)); t2 = t*t
            for m in range(M+1):
                Mt[m] += val; val *= t2
        p = q
    return Mt
Mt = moments(7, mpf('0.1'))
def moment_check(deg, width, m):
    X = GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec); s = mpf(0); p = tlo
    while p < thi:
        q = min(p + width, thi); hw, cc = (q-p)/2, (q+p)/2
        s += sum(hw*w*mp.exp(logint(n+m, hw*x+cc)) for (x, w) in X); p = q
    return s
err = max(abs(Mt[m]/moment_check(6, mpf('0.1'), m) - 1) for m in (0, M//2, M))
print('n = %d, dps = %d, M = %d, t-range [%s, %s], GL deg 7 vs 6 max rel diff of moments: %s (%.1f s)' % (n, dps, M, mp.nstr(tlo, 4), mp.nstr(thi, 4), mp.nstr(err, 3), time.time()-t0))
# c_m = gamma(n+m)/(gamma(n) m!), gamma(k) = 32 k! Mt(k)/(2k)!  =>  gamma(k+1)/gamma(k) = Mt(k+1)/(2(2k+1) Mt(k))
c = [mpf(1)]; r = mpf(1)
for m in range(1, M+1):
    k = n + m - 1
    r *= Mt[m]/(2*(2*k+1)*Mt[m-1])
    c.append(r/mp.factorial(m))
def F(w):
    s = mpf(0)
    for cm in reversed(c): s = s*w + cm
    return s
u = uinv(mpf(n)); R = mpf(9)/10*(n-mpf(1)/2)**2/(u+1)**2
print('u(n) = %s, R_n = %s, D_n = %s' % (mp.nstr(u, 10), mp.nstr(R, 10), mp.nstr(R**2/(4*a**2*(R+a**2)), 8)))
def tail_ok(wabs):  # last terms negligible relative to largest term at |w| = wabs
    terms = [abs(cm)*wabs**m for m, cm in enumerate(c)]
    return mp.nstr(max(terms[-5:])/max(terms), 3), mp.nstr(mp.log10(max(terms)), 4)
def winding(Rr):
    """number of zeros of Theta^{(n)} in W_{Rr}: arg change along upper boundary (vertex -> arc -> segment -> -Rr) over pi"""
    vm = 2*a*mp.sqrt(Rr + a**2)
    path = lambda tau: (a**2 - (tau*vm)**2/(4*a**2) + 1j*tau*vm) if tau <= 1 else (-Rr + 1j*vm*(2-tau))
    taus = [mpf(j)/400 for j in range(0, 801)]
    tot = mpf(0); maxstep = mpf(0); fprev = F(path(taus[0])); i = 0; refin = 0
    stack = list(zip(taus[:-1], taus[1:]))[::-1]
    while stack:
        t1, t2 = stack.pop()
        f1, f2 = F(path(t1)), F(path(t2))
        st = mp.arg(f2/f1)
        if abs(st) > mpf('0.3') and t2 - t1 > mpf(10)**-12:
            mid = (t1+t2)/2; stack.append((mid, t2)); stack.append((t1, mid)); refin += 1; continue
        tot += st; maxstep = max(maxstep, abs(st))
    return tot/pi, maxstep, refin, vm
KS = [int(k) for k in sys.argv[4].split(',')] if len(sys.argv) > 4 else [1]
for Rr, label in [(k*R, 'W_{%d R_n}' % k) for k in KS]:
    N, ms, rf, vm = winding(Rr)
    print('%s: log10|Theta^(n)/gamma(n)| at w = -Rr: %s;' % (label, mp.nstr(mp.log10(abs(F(-Rr))), 5)))
    print('%s: vmax = %s, tail/max term at |w|max: %s (log10 max term %s); winding (zeros inside) = %s, max arg step %s, refinements %d (%.1f s)' % (label, mp.nstr(vm, 6), *tail_ok(mp.sqrt(Rr**2 + vm**2)), mp.nstr(N, 6), mp.nstr(ms, 3), rf, time.time()-t0))
mn = min(F(-R*mpf(j)/200) for j in range(0, 201)) ; print('min of Theta^{(n)}/gamma(n) on [-R_n, 0] (sampled 201 pts) =', mp.nstr(mn, 6))
prev = F(-R); s1 = None
SMAX = max(KS)
for j in range(1, 200*(SMAX-1)+1):
    s = R*(1 + mpf(j)/200)
    v = F(-s)
    if v*prev < 0:
        lo, hi = R*(1 + mpf(j-1)/200), s
        for _ in range(40):
            mid = (lo+hi)/2
            if F(-mid)*prev > 0: lo = mid
            else: hi = mid
        s1 = (lo+hi)/2; break
    prev = v
if s1 is None:
    print('no sign change of Theta^{(n)}(-s) for s in [R_n, %d R_n] (step R_n/200)' % SMAX)
else:
    print('first sign change at s1 = %s = %s R_n  (j_{nu,1}^2/u^2 = %s); tail at |w|=s1: %s' % (mp.nstr(s1, 10), mp.nstr(s1/R, 6), mp.nstr(mp.besseljzero(n-mpf(1)/2, 1)**2/u**2, 8), tail_ok(s1)[0]))
    N, ms, rf, vm = winding(s1*mpf('0.999'))
    print('W_{0.999 s1}: winding = %s, max arg step %s, refinements %d (%.1f s)' % (mp.nstr(N, 6), mp.nstr(ms, 3), rf, time.time()-t0))
    N, ms, rf, vm = winding(s1*mpf('1.001'))
    print('W_{1.001 s1} (contains the real zero s1): winding = %s, refinements %d (%.1f s)' % (mp.nstr(N, 6), rf, time.time()-t0))
