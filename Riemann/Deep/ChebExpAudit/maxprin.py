#!/usr/bin/env python3
"""chebaudit: numerical check of the maximum principle (DERIVATION 3.4) on a dense 2-D grid of the strip
delta <= beta <= 1 - delta, 0 <= gamma <= 60, for 5 kernels at their certificate x: interior max of R must not exceed
the boundary max M_x = max(A1, x^(2 delta - 1) A2)."""
import json, math, cmath
D = '/Users/adriansager/playground/riemann/Riemann/Deep/ChebExp/'
ch = {c['name']: c for c in json.load(open(D + 'chains.json'))}
B2 = (777 - math.sqrt(1081)) / 852
DEL = {'1/32': 1 / 32, '1-(777-1081**0.5)/852': 1 - B2}
def R(beta, g, h, ss, s0, lx, delta):
    z = complex(beta - ss, g)
    ph = h ** 3 if abs(z) < 1e-12 else abs(((1 - cmath.exp(-h * z)) / z) ** 3)
    return math.exp((beta - (1 - delta)) * lx) * ph * ((s0 - beta) ** 2 + g * g) / (s0 - beta)
for name, k in [('CTu', 0), ('CTu', 1), ('Lu', 0), ('Lu', 19), ('CT10000', 0), ('L10000', 0)]:
    c = ch[name]; iv = c['intervals'][k]; delta = DEL[c['delta_expr']]
    h, ss, s0 = iv['h'], iv['ss'], iv['s0']; lx = math.log(c['c'] * (iv['Lam_a'] + math.log(3))) / delta
    step = (2 * math.pi / h) / 60; ng = int(60 / step) + 1; nb = 41
    bmax = 0.0; imax = 0.0; iarg = None
    for ib in range(nb):
        beta = delta + (1 - 2 * delta) * ib / (nb - 1)
        m = max(R(beta, ig * step, h, ss, s0, lx, delta) for ig in range(ng))
        if ib in (0, nb - 1): bmax = max(bmax, m)
        elif m > imax: imax, iarg = m, beta
    print('%-8s %3d  h=%.3f s*=%.3f s0=%.4f log x=%.2f  boundary max %.6g  interior max %.6g (beta=%.4f)  ratio %.6f %s' %
          (name, k, h, ss, s0, lx, bmax, imax, iarg, imax / bmax, 'OK' if imax <= bmax * (1 + 1e-9) else 'VIOLATION'), flush=True)
