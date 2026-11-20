#!/usr/bin/env python3
# dbnhigh admis_check.py: numeric checks of the derivative-design formulas used in lp.py and REPORT section 6.
import math
from mpmath import mp, mpf, diff, quad, cosh, exp, re, mpc, factorial, inf
mp.dps = 40
def K(h, d, b): return (h - b)/(d*d + (h - b)**2) + (h + b)/(d*d + (h + b)**2)
worst = 0
for (h0, k, d, b) in [(4, 3, 0.7, 0.5), (3, 6, 1.9, 0.3), (4, 10, 0.2, 0.74), (2.5, 8, 3.0, 0.0)]:
    num = (-1)**k*diff(lambda h: K(h, mpf(d), mpf(b)), mpf(h0), k)/factorial(k)
    cf = re((mpc(h0 - b, d))**(-k - 1) + (mpc(h0 + b, d))**(-k - 1))
    worst = max(worst, abs(num - cf)/abs(cf))
print('R_k closed form vs (1/k!)(-d/dh)^k K_h: max rel err %.1e' % worst)
worst = 0
for (h0, k, y) in [(4, 5, 0.5), (3, 12, 0.74), (4, 20, 0.12)]:
    lap = quad(lambda u: (u/2)**k*exp(-mpf(h0)*u/2)/factorial(k)*cosh(mpf(y)*u/2), [0, inf])
    cf = (mpf(h0) - y)**(-k - 1) + (mpf(h0) + y)**(-k - 1)
    worst = max(worst, abs(lap - cf)/cf)
print('O_k = int v_k cosh(yu/2) du vs closed form: max rel err %.1e' % worst)
# far field d^2 R_k(d,0) at large d
for k in (1, 2, 3):
    d = mpf(10)**6
    print('k=%d far field d^2 R_k(1e6,0) = %s' % (k, mp.nstr(d*d*re(2*(mpc(4, d))**(-k - 1)), 6)))
# Raymond's jet physical charge (TeX eq. src:jet-physical) vs the own-pair effect of one order-k term, p=4, r=2, delta=e^-15
delta = math.exp(-15); eta = delta/(1 - delta)
for y in (0.74, 0.5, 0.12):
    row = []
    for k in (10, 15, 20, 25, 30, 44):
        charge = eta/(4*2)*(k + 1)*(4/2)**k        # |b_k| = 1, Raymond's normalisation (1/p)(b p^k/k!) d^k
        effect = (1/4)*4**k*((4 - y)**(-k - 1) + (4 + y)**(-k - 1))   # own pair of the same term
        row.append('k=%d %.2g' % (k, charge/effect))
    print('y=%.2f charge/own-pair-effect:' % y, ', '.join(row))
