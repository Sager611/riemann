# Float (mpmath, 50 digits) sanity check of the upper end 0.7851 (SmallGap cert.json), NOT a certificate:
# f(x) = P(sqrt(2 pi) x / s) exp(-pi x^2/s^2), P(u) = sum_k B[k] H_{2k}(u); f^(xi) = s sum_k (-1)^k B[k] H_{2k}(sqrt(2 pi) s xi) exp(-pi s^2 xi^2).
import json, resource
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 50
C = json.load(open('Riemann/Deep/SmallGap/cert.json'))
B = [mp.mpf(b) for b in C['B']]; s = mp.mpf(Fr(C['s']).numerator) / Fr(C['s']).denominator
lam = mp.mpf(C['lam']); th = mp.mpf(399993) / 700000
def Hs(u, n):
    h0, h1 = mp.mpf(1), 2 * u; out = [h0, h1]
    for m in range(1, n):
        h0, h1 = h1, 2 * u * h1 - 2 * m * h0; out.append(h1)
    return out
def f(x):
    u = mp.sqrt(2 * mp.pi) * x / s; H = Hs(u, 2 * len(B))
    return sum(B[k] * H[2 * k] for k in range(len(B))) * mp.e ** (-mp.pi * x * x / s / s)
def fh(xi):
    u = mp.sqrt(2 * mp.pi) * s * xi; H = Hs(u, 2 * len(B))
    return s * sum((-1) ** k * B[k] * H[2 * k] for k in range(len(B))) * mp.e ** (-mp.pi * s * s * xi * xi)
F0 = fh(0); J = f(0) + 2 * mp.quad(lambda a: f(a) * a, [0, th / 2, th])
print('fhat(0) =', mp.nstr(F0, 15), ' J_theta(f) =', mp.nstr(J, 15), ' J/fhat(0) =', mp.nstr(J / F0, 15), '(cert.json sdp_value', C['sdp_value'], ')')
xs = [mp.mpf(i) / 200 for i in range(0, 2001)]
fmin = min((f(x), x) for x in xs); print('min f on [0,10] (step 0.005):', mp.nstr(fmin[0], 6), 'at', mp.nstr(fmin[1], 5))
ps = [lam + mp.mpf(i) / 400 for i in range(0, 3601)]
fhmax = max((fh(x), x) for x in ps); print('max f^ on [0.7851, 9.79] (step 0.0025):', mp.nstr(fhmax[0], 6), 'at', mp.nstr(fhmax[1], 6))
print('f^ just below lam:', mp.nstr(fh(lam - mp.mpf('0.001')), 6), ' at 0.775:', mp.nstr(fh(mp.mpf('0.775')), 6))
print(f'maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB')
