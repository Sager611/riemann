#!/usr/bin/env python3
# alprows P7: internal feasibility of the inert gain row LinkG_at (Formalism/SharpRestriction/GlobalTheoremInputs.lean:129-134)
# per (lift j, t), with vonMangoldtClassSum (Formalism/LocalRingCardinality/DataParams2.lean:240-243) exactly as defined
# (finite sum over i < Nat.log p a + 2, odd i, (a/p^i - wPart)^+), and the partner weight it needs if the gain per lift is
# the theta-data slope a/(p^2-1) (AlpToy: v_p(D^a theta(0)) = a/(p^2-1) + O(1); Lean rho, FnumEqRhoMulSumOfDemandRatioLe).
from fractions import Fraction as Fr
from sympy import primerange
def natlog(p, a):
    e = 0
    while p ** (e + 1) <= a: e += 1
    return e
worst = Fr(0)
for a in (10**3, 10**6, 10**9):
    for p in primerange(2, 31):
        dem = sum(Fr(a, p ** i) for i in range(natlog(p, a) + 2) if i % 2 == 1)   # wPart = 0
        sup = Fr(4 * a, p)                                                        # per-lift cap with |TT_w| = |TT|
        rho = Fr(a, p * p - 1)
        need_w = dem - rho                                                        # weights needed if gain = rho*a
        worst = max(worst, dem / sup)
        if a == 10**6 and p in (2, 3, 5, 7, 29):
            print('a=%d p=%2d demand(w=0)=%.6g cap=%.6g demand/cap=%.4f  theta-slope gain=%.6g  needed wPart=%.6g  a/(p+1)=%.6g'
                  % (a, p, float(dem), float(sup), float(dem / sup), float(rho), float(need_w), a / (p + 1)))
print('max demand/cap over a in {1e3,1e6,1e9}, p<=29:', float(worst), '(<= 1/3: the row is internally consistent)')

