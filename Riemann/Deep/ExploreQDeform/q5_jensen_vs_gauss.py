#!/usr/bin/env python3
"""q5_jensen_vs_gauss.py -- explore16: can the Jensen truncation be factored through a Gaussian (PSW sector-reducing) step?
The Jensen multiplier lambda_j = d!/((d-j)! d^j) = e^{-j(j-1)/(2d)} * mu_j.  If mu were a multiplier sequence, PSW's theorem
(q^{k^2}: sector S(t) -> S(g), cos g = cos t / q) plus mu would give J^d hyperbolicity from zero angles, and the cylinder-flow
zero ODE v_k' = d + sum_j coth((v_k - v_j)/2) (exact for Gaussian damping) would make explore15's screening rigorous.
Test (Polya-Schur algebraic characterisation): mu is a multiplier sequence iff sum_j C(m,j) mu_j X^j is real-rooted for all m;
we test m = d and also the Gaussian e^{-j(j-1)/(2D)} with D = c*d for c = 2, 4, 8 (residual lambda_j e^{j(j-1)/(2D)}).
Run: uv run --no-project --with python-flint python3 -u q5_jensen_vs_gauss.py
"""
from math import comb, factorial
from flint import arb, ctx
from q1_margins import real_count

def main():
    ctx.prec = 600
    print(' d   c=1(mu): real/d   c=2       c=4       c=8')
    for d in (3, 4, 6, 8, 12, 16, 24, 32, 40):
        row = []
        for c in (1, 2, 4, 8):
            D = c * d
            coef = [arb(comb(d, j)) * arb(factorial(d)) / (arb(factorial(d - j)) * arb(d) ** j) * (arb(j * (j - 1)) / (2 * D)).exp() for j in range(d + 1)]
            row.append(f'{real_count(coef)}/{d}')
        print(f'{d:3d}   ' + '   '.join(f'{x:>7}' for x in row))

if __name__ == '__main__':
    main()
