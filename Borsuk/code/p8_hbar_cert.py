#!/usr/bin/env python3
"""Barrier certificate for Tolmachev's outer body H(u) (arXiv:2611.16775), u = (0.8, 0, 0, 0): a finite set Y inside H(u)
whose graph of pairs at distance >= t is not 5-colourable, so every partition of H(u) into 5 parts has a part of
diameter >= t.  H(u) = B(0, sqrt rho) ∩ ⋂_{i<=3} B(P_i, 1) ∩ {x4 <= 1/2} with (Borsuk/code/p5_T.py H_params, u1 = u2 = u3 = 0):
t0 = 0.8 sqrt(1/15) - 0.2 sqrt(1/12) > 0, rho = 1/3 + t0^2, tau = 5 rho / 6, h = sqrt(rho / 6),
P_1 = (sqrt tau, 0, h, 0), P_{2,3} = (-sqrt(tau)/2, ±sqrt(3 tau)/2, h, 0)  (equilateral triangle, circumradius^2 tau).
Checks: exact rational points; membership with 80-digit decimals and an explicit margin; exact squared distances;
DIMACS CNF; kissat UNSAT; drat-trim VERIFIED.   Usage: python p8_hbar_cert.py POINTS.txt T_NUM T_DEN"""
import sys, subprocess
from fractions import Fraction as Fr
from decimal import Decimal, getcontext
getcontext().prec = 80
pts = [tuple(Fr(c) for c in line.split()) for line in open(sys.argv[1]) if line.strip() and not line.startswith('#')]
t2 = Fr(int(sys.argv[2]), int(sys.argv[3])) ** 2
D = lambda q: Decimal(q.numerator) / Decimal(q.denominator)
t0 = Decimal('0.8') * (Decimal(1) / Decimal(15)).sqrt() - Decimal('0.2') * (Decimal(1) / Decimal(12)).sqrt(); assert t0 > 0
rho = Decimal(1) / Decimal(3) + t0 * t0; tau = Decimal(5) * rho / Decimal(6); h = (rho / Decimal(6)).sqrt(); st = tau.sqrt()
P = [(st, Decimal(0), h, Decimal(0)), (-st / 2, (Decimal(3) * tau).sqrt() / 2, h, Decimal(0)), (-st / 2, -(Decimal(3) * tau).sqrt() / 2, h, Decimal(0))]
worst = Decimal(10)
for y in pts:
    yd = [D(c) for c in y]
    worst = min(worst, rho - sum(c * c for c in yd)); assert y[3] <= Fr(1, 2)
    for p in P: worst = min(worst, 1 - sum((a - b) ** 2 for a, b in zip(yd, p)))
assert worst > Decimal('1e-30'), worst
n = len(pts); k = 5; E = []
for i in range(n):
    for j in range(i + 1, n):
        if sum((a - b) ** 2 for a, b in zip(pts[i], pts[j])) >= t2: E.append((i, j))
with open('/tmp/hbar.cnf', 'w') as f:
    f.write('p cnf %d %d\n' % (n * k, n + k * len(E)))
    for v in range(n): f.write(' '.join(str(v * k + c + 1) for c in range(k)) + ' 0\n')
    for i, j in E:
        for c in range(k): f.write('-%d -%d 0\n' % (i * k + c + 1, j * k + c + 1))
T = 'Borsuk/.scratch/tools'
r = subprocess.run([T + '/kissat/build/kissat', '-q', '/tmp/hbar.cnf', '/tmp/hbar.drat', '--no-binary'], capture_output=True, text=True)
print('points', n, 'edges', len(E), 'rho', str(rho)[:12], 'min membership slack %.3e' % worst, 'kissat', {10: 'SAT', 20: 'UNSAT'}.get(r.returncode, r.returncode), flush=True)
if r.returncode == 20:
    v = subprocess.run([T + '/drat-trim/drat-trim', '/tmp/hbar.cnf', '/tmp/hbar.drat'], capture_output=True, text=True).stdout
    print([l for l in v.splitlines() if l.startswith('s ')])
