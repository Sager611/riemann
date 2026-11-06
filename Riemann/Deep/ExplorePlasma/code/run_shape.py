import sys, time; sys.path.insert(0, __file__.rsplit('/', 1)[0])
# T2 + adversarial (a): shapes of Phi and PhiDH; single-hump threshold t_b and Penrose threshold t_P for K e^{s u^2} (Polymath t = 4 s)
from kern import *
mp.mp.dps = 30
t0 = time.time()
def shape(T, name, U=3.0, h=0.005):
    n = int(U/h); minF = mp.inf; maxd1 = -mp.inf; maxlc = -mp.inf; worst = -mp.inf; a1 = a2 = None
    for i in range(1, n+1):
        u = i*h; f0 = ev(T, u); f1 = ev(T, u, 1); f2 = ev(T, u, 2)
        lc = (f0*f2 - f1**2)/f0**2
        minF = min(minF, f0)
        if f1/f0 > maxd1: maxd1, a1 = f1/f0, u
        if lc > maxlc: maxlc, a2 = lc, u
        worst = max(worst, (f0*f2 - f1**2 - f0*f1/u)/(4*u**2*f0**2))   # (log K(sqrt w))'' at w = u^2
    print(name, 'on (0,%g]: min K %s | max K\'/K %s at %s | max (log K)\'\' %s at %s | max (log K(sqrt w))\'\' %s' % (U, mp.nstr(minF, 4), mp.nstr(maxd1, 6), a1, mp.nstr(maxlc, 6), a2, mp.nstr(worst, 6)))
shape(TZ, 'zeta Phi ')
shape(TD, 'DH PhiDH ')
for T, nm in [(TZ, 'zeta'), (TD, 'DH')]:
    f0 = ev(T, 0); f2 = ev(T, 0, 2); sb = -f2/(2*f0); U = mp.mpf(5)
    def P(s):
        def g(u):
            if u < mp.mpf('1e-6'): return f2/2 + s*f0
            return (ev(T, u)*mp.e**(s*u**2) - f0)/u**2
        return 2*(mp.quad(g, mp.linspace(0, U, 21), method='gauss-legendre') - f0/U)
    vals = [(s, P(s)) for s in [0, sb/2, sb, 1.5*sb, 2*sb, 3*sb]]
    print(nm, 'Penrose integral at the centre: ' + ' '.join('P(%s)=%s' % (mp.nstr(s, 5), mp.nstr(v, 6)) for s, v in vals))
    lo = hi = None
    for (s1, v1), (s2, v2) in zip(vals, vals[1:]):
        if v1 < 0 <= v2: lo, hi = s1, s2
    for it in range(45):
        mid = (lo + hi)/2
        if P(mid) < 0: lo = mid
        else: hi = mid
    sP = (lo + hi)/2
    hh = mp.mpf('1e-8'); dP = (P(hh) - P(-hh))/(2*hh); I = 2*mp.quad(lambda u: ev(T, u), mp.linspace(0, U, 21), method='gauss-legendre')
    print(nm, "K(0)=%s (log K)''(0)=%s s_b=%s t_b=%s s_P=%s t_P=%s | P'(0)=%s vs int K=%s" % tuple(mp.nstr(v, 10) for v in (f0, f2/f0, sb, 4*sb, sP, 4*sP, dP, I)))
print('time %.1f' % (time.time() - t0))
