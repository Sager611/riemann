#!/usr/bin/env python3
"""P5 barrier certificate: a finite set Y inside Tolmachev's outer body H(u*) (arXiv:2611.16775, Def. H) whose
far graph (pairs at distance >= t) is not 5-colourable, so every partition of H(u*) into 5 parts has a part of
diameter >= t.  u* = (1,0,0,0): rho = 2/5, contact triangle = unit equilateral triangle at height h = sqrt(1/15)
(three vertices of a regular unit simplex inscribed in the Jung sphere).  Checks: exact rational points;
membership |y|^2 <= 2/5 (exact), x4 <= 1/2 (exact), |y - P_i|^2 <= 1 (60-digit decimals, margin reported);
exact squared distances; DIMACS CNF; kissat UNSAT; drat-trim VERIFIED.
Usage: python p5_cert.py POINTS.txt T_NUM T_DEN   (writes /tmp/p5cert.cnf, /tmp/p5cert.drat)
"""
import sys, subprocess
from fractions import Fraction as Fr
from decimal import Decimal, getcontext
getcontext().prec = 60
pts = [tuple(Fr(c) for c in line.split()) for line in open(sys.argv[1]) if line.strip() and not line.startswith('#')]
t2 = Fr(int(sys.argv[2]), int(sys.argv[3])) ** 2
D = lambda q: Decimal(q.numerator) / Decimal(q.denominator)
s3, s15 = Decimal(3).sqrt(), Decimal(15).sqrt()
P = [(1 / s3, Decimal(0), 1 / s15, Decimal(0)), (-1 / (2 * s3), Decimal('0.5'), 1 / s15, Decimal(0)), (-1 / (2 * s3), Decimal('-0.5'), 1 / s15, Decimal(0))]
rho = Fr(2, 5); worst = Decimal(10)
for y in pts:
    n2 = sum(c * c for c in y); assert n2 <= rho and y[3] <= Fr(1, 2)
    for p in P:
        d = sum((D(c) - q) ** 2 for c, q in zip(y, p)); worst = min(worst, 1 - d); assert 1 - d > Decimal('1e-40')
for p in P: assert abs(sum(q * q for q in p) - D(rho)) < Decimal('1e-50')
for i in range(3):
    for j in range(i + 1, 3): assert abs(sum((a - b) ** 2 for a, b in zip(P[i], P[j])) - 1) < Decimal('1e-50')
n = len(pts); k = 5; E = []
for i in range(n):
    for j in range(i + 1, n):
        if sum((a - b) ** 2 for a, b in zip(pts[i], pts[j])) >= t2: E.append((i, j))
with open('/tmp/p5cert.cnf', 'w') as f:
    f.write('p cnf %d %d\n' % (n * k, n + k * len(E)))
    for v in range(n): f.write(' '.join(str(v * k + c + 1) for c in range(k)) + ' 0\n')
    for i, j in E:
        for c in range(k): f.write('-%d -%d 0\n' % (i * k + c + 1, j * k + c + 1))
T = 'Borsuk/.scratch/tools'
r = subprocess.run([T + '/kissat/build/kissat', '-q', '/tmp/p5cert.cnf', '/tmp/p5cert.drat', '--no-binary'], capture_output=True, text=True)
print('points', n, 'edges', len(E), 'min membership slack', '%.3e' % worst, 'kissat', {10: 'SAT', 20: 'UNSAT'}.get(r.returncode, r.returncode))
if r.returncode == 20:
    v = subprocess.run([T + '/drat-trim/drat-trim', '/tmp/p5cert.cnf', '/tmp/p5cert.drat'], capture_output=True, text=True).stdout
    print([l for l in v.splitlines() if l.startswith('s ')])
