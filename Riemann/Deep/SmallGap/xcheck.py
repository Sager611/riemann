# xcheck.py: independent mpmath cross-check of cert.json conventions (Fourier transform by quadrature, functional by quadrature)
# usage: uv run --no-project --with mpmath python3 xcheck.py cert.json
import sys, json
from fractions import Fraction
import mpmath as mp
mp.mp.dps = 60
C = json.load(open(sys.argv[1]))
d = C['d']; s = mp.mpf(Fraction(C['s']).numerator)/Fraction(C['s']).denominator
lam = mp.mpf(C['lam']); th = mp.mpf(399993)/700000
B = [mp.mpf(b) for b in C['B']]
def P(u):   # sum B_k H_{2k}(u), physicists Hermite via mpmath
    return mp.fsum(B[k]*mp.hermite(2*k, u) for k in range(d+1))
def Qp(v):
    return mp.fsum((B[k] if k % 2 == 0 else -B[k])*mp.hermite(2*k, v) for k in range(d+1))
c = mp.sqrt(2*mp.pi)
f = lambda x: P(c*x/s)*mp.exp(-mp.pi*x**2/s**2)
fhat_formula = lambda xi: s*Qp(c*s*xi)*mp.exp(-mp.pi*s**2*xi**2)
fhat_quad = lambda xi: 2*mp.quad(lambda x: f(x)*mp.cos(2*mp.pi*x*xi), mp.linspace(0, 14, 29))
worst = 0
for xi in [mp.mpf(0), mp.mpf('0.3'), lam, mp.mpf('0.9'), mp.mpf('1.4'), mp.mpf('2.2'), mp.mpf('3.1')]:
    a, b = fhat_formula(xi), fhat_quad(xi)
    worst = max(worst, abs(a-b))
    print('xi=%s  formula %s  quadrature %s' % (mp.nstr(xi, 5), mp.nstr(a, 15), mp.nstr(b, 15)))
print('max |formula - quadrature| =', mp.nstr(worst, 3))
fh0 = fhat_quad(0)
J = f(0) + 2*mp.quad(lambda a: f(a)*a, [0, th])
print('J/fhat(0) by quadrature =', mp.nstr(J/fh0, 20))
# coarse float grid scans (sanity only; the proof is cert.py)
mp.mp.dps = 30
fmin = min(f(mp.mpf(i)/100) for i in range(0, 1001))
fhmax = max(fhat_formula(lam + mp.mpf(i)/100) for i in range(0, 801))
print('min f on [0,10] grid =', mp.nstr(fmin, 5), '; max fhat on [lam, lam+8] grid =', mp.nstr(fhmax, 5))
