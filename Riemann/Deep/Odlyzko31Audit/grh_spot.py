#!/usr/bin/env python3
"""odlaudit spot check of the data file odlyzko1976_tables12.csv (GRH columns) with our own formula at w = 0, kernel G(x/b):
max over b of (r1 log A + 2 r2 log B - E)/n. Usage: uv run --no-project --with mpmath python3 grh_spot.py > logs/grh_spot.log"""
import mpmath as mp
mp.mp.prec = 80
pi, EG = mp.pi, mp.euler
def G(t): return (1 - t / 2) * mp.cos(pi * t / 2) + mp.sin(pi * t / 2) / pi
def val(n, r1, b):
    X = 2 * b; pts = mp.linspace(0, X, int(X) + 2)
    Ic = mp.quad(lambda x: (1 - G(x / b)) / (2 * mp.sinh(x / 2)), pts) - mp.log(mp.tanh(X / 4))
    Ir = mp.quad(lambda x: (1 - G(x / b)) / (2 * mp.cosh(x / 2)), pts) + pi / 2 - 2 * mp.atan(mp.tanh(X / 4))
    E = 4 * mp.quad(lambda x: G(x / b) * mp.cosh(x / 2), pts)
    lB = EG + mp.log(8 * pi) - Ic; lA = lB + pi / 2 - Ir
    return mp.exp((r1 * lA + (n - r1) * lB - E) / n)
for n, r1, ref in ((100, 0, 19.268), (100, 100, 59.141), (1000, 0, 29.094), (8, 0, 5.734)):
    best = max((val(n, r1, mp.mpf(b) / 20), b / 20) for b in range(10, 400, 6))
    b0 = best[1]
    best = max((val(n, r1, mp.mpf(b0) + d), b0 + d) for d in [k / 50 for k in range(-15, 16)])
    print(f'GRH w=0, kernel G: n={n} r1={r1}: rd > {mp.nstr(best[0], 8)} at b={best[1]:.3f}; Odlyzko 1976 Table 1 (data file): {ref}')
