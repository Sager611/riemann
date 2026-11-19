#!/usr/bin/env python3
"""q5audit: true Arb ball radii of the q5cert.py quantities (reuses q5cert.py functions read-only from ../Degree2Mol).
Run: uv run --no-project --with python-flint --with mpmath python3 q5cert_radii.py"""
import os
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Degree2Mol')
src = open(os.path.join(D, 'q5cert.py')).read().split("nodes = json.load")[0]
exec(src)
nodes = json.load(open(os.path.join(D, 'q5nodes.json')))
s35 = arb(fmpq(3, 5)).sqrt(); s15 = arb(15).sqrt(); a0 = 2 / s15; b0 = arb(fmpq(1, 2)) - s35 / 3 - 2 / s15; rhs = (-(arb(2) / 3) * s35).exp()
for nd in nodes:
    R = F(nd['R']); h = [F(s) for s in nd['h']]; h[-1] = fmpq(1, 2) - sum(h[:-1], fmpq(0)); Q = Q_cheb(h)
    I = I_exact(Q, R); lhs = a0 * (I + arb(-R).exp())
    if nd['Rb'] is not None: lhs = lhs + b0 * arb(-F(nd['Rb'])).exp()
    print(nd['R'], "rad(I) =", I.rad().str(3), " rad(lhs) =", lhs.rad().str(3), " rhs-lhs =", (rhs - lhs).str(8))
