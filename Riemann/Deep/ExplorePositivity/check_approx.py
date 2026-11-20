# check_approx.py — independent check of mp_approx_roots.py: (i) mpmath polyroots at 300 digits (N = 10, 20),
# (ii) real-zero count of Xi_N by sign changes of the recurrence-evaluated sum on a fine grid (N = 20, 50),
# (iii) residual of Xi_N at the first reported non-real root.
import mpmath as mp
from mp_core import f_coeffs, f_eval
mp.mp.dps = 300
cs = [mp.mpf(l.split()[1]) for l in open('logs/c2n_zeta.txt') if l.strip()]
F = f_coeffs(120)
def coeffs(N):
    co = [mp.mpf(0)] * (2*N + 1)
    for n in range(N + 1):
        for j, q in enumerate(F[2*n]):
            if q: co[j] += (-1)**n * cs[n] * mp.mpf(q.numerator) / q.denominator
    return co
for N in (10, 20):
    rr = mp.polyroots(coeffs(N)[::-1], maxsteps=2000, extraprec=2000)
    nre = sum(1 for z in rr if abs(mp.im(z)) < mp.mpf(10)**-30)
    print('N=%d polyroots(300 digits): real %d nonreal %d; nonreal t = 2x:' % (N, nre, len(rr) - nre),
          sorted(set((mp.nstr(abs(2*mp.re(z)), 6), mp.nstr(abs(2*mp.im(z)), 6)) for z in rr if abs(mp.im(z)) >= mp.mpf(10)**-30)))
def XiN(N, t):
    v = f_eval(2*N, t/2)
    return sum((-1)**n * cs[n] * v[2*n] for n in range(N + 1))
mp.mp.dps = 60
for N in (20, 50):
    T = 4.2*N; h = mp.mpf('0.01'); k = 0; prev = XiN(N, mp.mpf('0.005')); ch = 0
    t = mp.mpf('0.005')
    while t < T:
        t += h; cur = XiN(N, t)
        if cur * prev < 0: ch += 1
        prev = cur
    print('N=%d sign changes on (0, %.0f], step 0.01: %d  -> real zeros %d of %d' % (N, T, ch, 2*ch, 2*N))
mp.mp.dps = 300
t0 = 2 * mp.findroot(lambda x: sum(c * x**j for j, c in enumerate(coeffs(20))), mp.mpc(-14.4513, 18.6656) / 2)
print('N=20 root refined t =', mp.nstr(t0, 15), ' |Xi_N(t0)| =', mp.nstr(abs(XiN(20, t0)), 5), ' |Xi_N(t0+0.1)| =', mp.nstr(abs(XiN(20, t0 + mp.mpf('0.1'))), 5))
