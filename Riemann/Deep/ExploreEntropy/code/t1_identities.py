# t1_identities.py (explore32): Beckner's constant vs the Weil archimedean symbol; Levy-Khintchine form; Burnol normalisation.
import mpmath as mp, numpy as np
mp.mp.dps = 30
D = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi)
th1 = mp.siegeltheta(0, derivative=1)
phi = lambda t: -mp.log(mp.pi) + mp.re(mp.digamma(mp.mpf(1)/4 + 1j*mp.mpf(t)/2))
print('D = psi(1/4) - ln pi =', D)
print("2 theta'(0)          =", 2*th1)
print('phi(0)               =', phi(0))
print('closed form -gamma - pi/2 - 3 ln2 - ln pi =', -mp.euler - mp.pi/2 - 3*mp.log(2) - mp.log(mp.pi))
# odd-sector (Gamma((s+1)/2)) constant
print('odd constant psi(3/4) - ln pi =', mp.digamma(mp.mpf(3)/4) - mp.log(mp.pi))
# Levy-Khintchine: phi(t) - phi(0) = 2 int_0^inf k(x) (1 - cos t x) dx, k(x) = e^{x/2}/(e^x - e^{-x})
k = lambda x: mp.exp(x/2)/(mp.exp(x) - mp.exp(-x))
for t in [0.5, 1, 5, 14.134725]:
    lk = 2*mp.quadosc(lambda x: k(x)*(1 - mp.cos(t*x)), [0, mp.inf], omega=t)
    print('t=%-9s phi(t)-phi(0)=%s  LK=%s' % (t, mp.nstr(phi(t) - phi(0), 15), mp.nstr(lk, 15)))
# killing constant check: D = -ln pi - gamma + int_0^inf (2e^{-x} - 2e^{x/2})/(e^x - e^{-x}) dx
I0 = mp.quad(lambda x: (2*mp.exp(-x) - 2*mp.exp(x/2))/(mp.exp(x) - mp.exp(-x)), [0, 1, mp.inf])
print('-ln pi - gamma + I0 =', -mp.log(mp.pi) - mp.euler + I0)
# Beckner for the self-dual Gaussian f = e^{-pi x^2} (f^ = f with e^{2 pi i x xi}): LHS/||f||^2 = 2 E ln|x|, x ~ N(0, 1/(4pi))
g = 2*mp.quad(lambda x: 2*mp.log(x)*mp.exp(-2*mp.pi*x*x), [0, 1, mp.inf])/mp.quad(lambda x: mp.exp(-2*mp.pi*x*x), [-mp.inf, mp.inf])
print('Gaussian: LHS/||f||^2 =', g, ' (closed form ln(1/2pi) - gamma - 2ln2 =', -mp.log(2*mp.pi) - mp.euler - 2*mp.log(2), ') vs D', D)
