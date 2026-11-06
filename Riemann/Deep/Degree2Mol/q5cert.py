#!/usr/bin/env python3
"""Question 5 of arXiv 2508.11108 (l.958-960): kappa(theta) > 2 theta/3 for every theta in (0, 1/2], where
kappa(theta) = 1 - log c/R, R = sqrt(3/5)/theta, P(x) = x and Q = Q_R the minimiser (paper's S_R).
For P = x, beta = 1 (Q(1)=0) and theta R = sqrt(3/5), every admissible S (S(0)=1/2, S(R)=0, S in H^1) gives
    c(S) = 1/2 - sqrt(3/5)/3 + (2/sqrt15) e^R I_R(S),  I_R(S) = int_0^R e^t(S^2+S'^2) + e^{-t}((1-S)^2+S'^2) dt
(DERIVATION.md section 4), and c(Q_R) <= c(S) because Q_R minimises this convex functional.
kappa > 2theta/3  <=>  c < e^{R - (2/3)sqrt(3/5)}.
Zero extension: S_a optimal-ish at R_a, extended by 0 on [R_a, R]: I_R = I_a + e^{-R_a} - e^{-R}. Hence for R in [R_a, R_b]
    (2/sqrt15)(I_a + e^{-R_a}) + (1/2 - sqrt(3/5)/3 - 2/sqrt15) e^{-R_b} < e^{-(2/3)sqrt(3/5)}
suffices (the coefficient of e^{-R} is negative, so R = R_b is the worst case; R_b = infinity drops the term).
Run: uv run --no-project --with python-flint --with mpmath python3 q5cert.py
"""
import json, os
from fractions import Fraction
from flint import fmpq, fmpq_poly, arb, ctx
ctx.prec = 4000   # exact sums reach ~1e300 before cancelling to O(1) at small R
HERE = os.path.dirname(os.path.abspath(__file__))
def F(s):
    fr = Fraction(s); return fmpq(fr.numerator, fr.denominator)
X = fmpq_poly([0, 1]); ONE = fmpq_poly([1])
def cheb(m):
    Ts = [ONE, X]
    for k in range(2, m + 1): Ts.append(2 * X * Ts[-1] - Ts[-2])
    return Ts[m]
def Q_cheb(h):
    z = 1 - 2 * X; Q = ONE
    for n, hn in enumerate(h, start=1): Q += hn * (cheb(2 * n - 1)(z) - 1)
    return Q
def int_exp_half(G, a):   # int_{1/2}^{1} e^{a y} G(y) dy = e^{a} Xc - e^{a/2} Yc (exact Xc, Yc)
    Xc = fmpq(0); Yc = fmpq(0); D = G; k = 0; ap = a
    while not D.is_zero():
        s = 1 if k % 2 == 0 else -1
        Xc += s * D(fmpq(1)) / ap; Yc += s * D(fmpq(1, 2)) / ap
        D = D.derivative(); ap *= a; k += 1
    return arb(Xc) * arb(a).exp() - arb(Yc) * arb(a / 2).exp()
def I_exact(Q, R):
    """I_R(S), S(t) = Q(1/2 + t/(2R)); t = 2R(y - 1/2): I = 2R int_{1/2}^1 [e^{-R} e^{2Ry}(Q^2 + Q'^2/(4R^2)) + e^{R} e^{-2Ry}((1-Q)^2 + Q'^2/(4R^2))] dy"""
    dQ = Q.derivative(); d2 = dQ * dQ / (4 * R * R)
    A = int_exp_half(Q * Q + d2, 2 * R); B = int_exp_half((1 - Q) * (1 - Q) + d2, -2 * R)
    return 2 * arb(R) * (A * arb(-R).exp() + B * arb(R).exp())
nodes = json.load(open(os.path.join(HERE, 'q5nodes.json')))
s35 = arb(fmpq(3, 5)).sqrt(); s15 = arb(15).sqrt()
a0 = 2 / s15; b0 = arb(fmpq(1, 2)) - s35 / 3 - 2 / s15; rhs = (-(arb(2) / 3) * s35).exp()
R0 = F(nodes[0]['R'])
print("R_0 =", R0, "< 2 sqrt(3/5) =", (2 * s35).str(12), ":", bool(arb(R0) < 2 * s35), " (so theta <= 1/2 is covered)")
ok = True; prev_b = None
for nd in nodes:
    R = F(nd['R']); h = [F(s) for s in nd['h']]
    h[-1] = fmpq(1, 2) - sum(h[:-1], fmpq(0))        # Q(1) = 0 exactly (beta = 1)
    Q = Q_cheb(h)
    dQ = Q.derivative()
    assert Q(fmpq(0)) == 1 and Q(fmpq(1)) == 0 and (dQ - dQ(1 - X)).is_zero()
    assert prev_b is None or prev_b == R, "chain must be contiguous"
    I = I_exact(Q, R)
    lhs = a0 * (I + arb(-R).exp())
    if nd['Rb'] is not None:
        Rb = F(nd['Rb']); lhs = lhs + b0 * arb(-Rb).exp(); rng = f"[{float(R)}, {float(Rb)}]"; prev_b = Rb
    else:
        rng = f"[{float(R)}, infinity)"; prev_b = None
    good = bool(lhs < rhs); ok &= good
    # kappa - 2 theta/3 at R itself (from c(S_a)), for information
    c = arb(fmpq(1, 2)) - s35 / 3 + a0 * arb(R).exp() * I
    kap = 1 - c.log() / arb(R); th = s35 / arb(R)
    print(f"node R_a = {R} (theta = {float(th.mid()):.5f}): I_a in {I.str(15)}; covers R in {rng}: lhs {lhs.str(10)} < rhs {rhs.str(10)}: {good};"
          f" kappa(theta_a) >= {kap.lower().str(8, radius=False)} vs 2theta/3 = {(2*th/3).str(8)}")
print("Question 5 (kappa(theta) > 2 theta/3 for all theta in (0, 1/2]):", "PROVED" if ok and prev_b is None else "NOT PROVED")
