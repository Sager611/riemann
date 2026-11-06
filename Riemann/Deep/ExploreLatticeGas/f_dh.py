#!/usr/bin/env python3
"""f_dh.py -- explore15, adversarial test (iii-a): the Davenport-Heilbronn function through the same Jensen tests.

f(s) = ((1 - i k)/2) L(s, chi) + ((1 + i k)/2) L(s, conj chi), chi mod 5 with chi(2) = i,
k = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1); f(s) = 5^-s [z(s,1/5) - z(s,4/5) + k (z(s,2/5) - z(s,3/5))] (Hurwitz).
xi_DH(s) = (5/pi)^((s+1)/2) Gamma((s+1)/2) f(s) = xi_DH(1 - s), real on the line; it has off-line zeros
(first 0.808517 + 85.699348 i).  Theta_DH(w) = xi_DH(1/2 + sqrt w) = sum gamma_DH(n) w^n/n!.
Tests: odd Taylor coefficients vanish (functional equation), signs of gamma_DH, certified hyperbolicity of
J^{d,0} (d <= 60) and J^{9,n} (n <= 40) -- the data any Jensen/Pascal-recursion argument would use.
Run: uv run --no-project --with python-flint python3 -u f_dh.py   (about 1 min)
"""
import sys, time, math
from math import comb
from flint import arb, acb, arb_series, arb_poly, acb_poly, ctx
sys.path.insert(0, '.')
from c_reach import certified_real_count

def main():
    t0 = time.time()
    N = 110
    ctx.prec = 4000; ctx.cap = 2 * N + 2
    s = arb_series([arb(1) / 2, 1], prec=2 * N + 2)
    r5 = arb(5).sqrt()
    kap = ((10 - 2 * r5).sqrt() - 2) / (r5 - 1)
    Z = lambda a: arb_series.zeta(s, a)
    f = (-s * arb(5).log()).exp() * (Z(arb(1) / 5) - Z(arb(4) / 5) + kap * (Z(arb(2) / 5) - Z(arb(3) / 5)))
    xi = (((s + 1) / 2) * (arb(5) / arb.pi()).log()).exp() * arb_series.gamma((s + 1) / 2) * f
    c = xi.coeffs()
    odd = max(abs(c[2 * k + 1] / c[0]) for k in range(N))
    g = []; fact = arb(1)
    for n in range(N + 1):
        if n: fact *= n
        g.append(fact * c[2 * n])
    print(f'kappa = {kap.str(15)}; xi_DH(1/2) = {c[0].str(15)}; max |odd coeff|/c0 = {odd.str(3)} ({time.time() - t0:.0f} s)')
    signs = ''.join('+' if x > 0 else ('-' if x < 0 else '?') for x in g)
    print('signs of gamma_DH(0..%d): %s' % (N, signs))
    # normalise to positive leading behaviour if needed
    okd = []
    for d in range(2, 61):
        coef = [comb(d, j) * g[j] for j in range(d + 1)]
        okd.append((d, certified_real_count(coef)))
    bad = [(d, k) for d, k in okd if k != d]
    print(f'J^(d,0)_DH, d = 2..60: certified hyperbolic for {sum(1 for d, k in okd if k == d)}/59 degrees; failures: {bad[:10]}')
    ok9 = []
    for n in range(0, 41):
        coef = [comb(9, j) * g[n + j] for j in range(10)]
        ok9.append(certified_real_count(coef) == 9)
    print(f'J^(9,n)_DH, n = 0..40: hyperbolic for {sum(ok9)}/41 shifts')
    # where Chasse's sector bound puts the first failure (first off-line zero 0.808517 + 85.699348 i)
    x, y = 0.808517 - 0.5, 85.699348
    D = int(1 / math.sin(2 * math.atan(x / y)) ** 2)
    print(f'first off-line DH zero 1/2 + {x:.6f} + {y} i: 1/sin^2(2 atan(x/y)) = {D} (toy: d* = D_sec + 2..23)')
    print(f'done ({time.time() - t0:.0f} s)')

if __name__ == '__main__':
    main()
