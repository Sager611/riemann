# B2' (reviewer fix 3): zeta(1-k, x) = -B_k(x)/k checked through Hurwitz's formula
# zeta(1-k, x) = 2 Gamma(k)/(2 pi)^k sum_{n>=1} cos(pi k/2 - 2 pi n x)/n^k, summed exactly over residues for x = a/b:
# sum_n f(n)/n^k = sum_{r=1}^{b} cos(pi k/2 - 2 pi r a/b) b^-k zeta(k, r/b)  (Hurwitz zeta at positive k only; no Bernoulli input)
import mpmath as mp
from fractions import Fraction
from sympy import Symbol, Poly, bernoulli as bp
X = Symbol('x'); mp.mp.dps = 40
worst = {}
for k in range(2, 13):
    e = 0
    for xv in [Fraction(1, 10), Fraction(37, 100), Fraction(1, 2), Fraction(93, 100), Fraction(1)]:
        a, b = xv.numerator, xv.denominator
        s = mp.fsum(mp.cos(mp.pi * k / 2 - 2 * mp.pi * r * a / b) * mp.mpf(b) ** (-k) * mp.zeta(k, mp.mpf(r) / b) for r in range(1, b + 1))
        ser = 2 * mp.gamma(k) / (2 * mp.pi) ** k * s
        ex = -Fraction(Poly(bp(k, X), X).eval(xv)) / k
        e = max(e, abs(ser - mp.mpf(ex.numerator) / ex.denominator))
    worst[k] = e
m = max(worst.values())
print({k: mp.nstr(v, 3) for k, v in worst.items()})
print(f"B2' Hurwitz zeta at 1-k by Hurwitz's Fourier formula (positive-k Hurwitz zeta over residues) vs -B_k(x)/k exact: {'PASS' if m < 1e-30 else 'FAIL'} k = 2..12, x in (1/10, 37/100, 1/2, 93/100, 1), max err {mp.nstr(m,3)}")
