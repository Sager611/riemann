# Independent high-precision check (mpmath) of C_fam(f_lam) = Q_{F1}(f_lam) and of the certified values at lam = 5/4, 7/6, 12/11, 1
from mpmath import mp, mpf, cos, pi, quad
mp.dps = 30
a = [mpf(v) / 10**6 for v in (122290, -45012, 16241, 1887, 2711, 3731, -5307, 2982, -5153, 5358)]
def f(u):
    return 1 + sum(aj * cos(2 * pi * (j + 1) * u) for j, aj in enumerate(a))
def R(d):
    return quad(lambda s: f(s) * f(s + d), [-mpf(1) / 2, mpf(1) / 2 - d])
I0 = quad(lambda u: f(u) ** 2, [-mpf(1) / 2, 0, mpf(1) / 2])
x = mpf(8722) / 10**6; b = mpf(4389) / 10**6
def C(lam):
    k = 1 / lam
    pts = [0, k, 1] if k < 1 else [0, 1]
    I1 = 2 * quad(lambda d: min(lam * d, 1) * R(d), pts)
    return I0 / lam + I1
for lam, J in [(mpf(1), '0.672501 (MT)'), (mpf(5) / 4, '0.797213'), (mpf(7) / 6, '0.764149'), (mpf(12) / 11, '0.727484')]:
    c = C(lam)
    s = (2 - c - lam * b) / (1 - x)
    print('lam=%s  C=%s  s=%s  Johns %s' % (mp.nstr(lam, 8), mp.nstr(c, 15), mp.nstr(s, 12), J), flush=True)
