# Exact exponent bookkeeping of Baluyot's thesis, ch. 3 (sympy).
import sympy as sp

A, N, U, V, B, T, y = sp.symbols('A N U V B T y', positive=True)
th, de, et, e0, sg, ph = sp.symbols('theta delta eta epsilon0 sigma phi', positive=True)

# 1. Lemma 3.1.6: X = (UNV)^{1/2}{(A^{-1}VU)^{1/2} + (N+A)^{1/4}[A^{-1}VU(N+A)(V+A^2) + N U^2]^{1/4}}
#    X^4 << (UNV)^2 {(A^{-1}VU)^2 + (N+A)[A^{-1}VU(N+A)(V+A^2) + N U^2]}, with (N+A)^k -> N^k + A^k
inner = (A**-1 * V * U)**2 + (N + A) * (A**-1 * V * U * (N + A) * (V + A**2) + N * U**2)
X4 = sp.expand((U * N * V)**2 * inner)
mons = set()
for term in sp.Add.make_args(X4):
    pw = sp.Poly(term * A**10, A, N, U, V).monoms()[0]
    mons.add((pw[0] - 10, pw[1], pw[2], pw[3]))
E_thesis = {(-2, 2, 4, 4), (1, 2, 3, 3), (3, 2, 3, 3), (-1, 4, 3, 4), (1, 4, 3, 3), (0, 4, 4, 2), (1, 3, 4, 2)}
print("monomials of X^4 (a,n,u,v):", sorted(mons))
print("in thesis E but not in X^4:", sorted(E_thesis - mons))
print("in X^4 but not in thesis E:", sorted(mons - E_thesis))
# the cross term N*A^2... comes with coefficient 2 from (N+A)^2; dominated by N^2,A^2 terms: check every monomial
# of X^4 is <= max of the E-set monomials with (1,2,3,3) replaced by (1,2,3,4)
E_fix = (E_thesis - {(1, 2, 3, 3)}) | {(1, 2, 3, 4)}
extra = sorted(mons - E_fix)
print("monomials of X^4 outside corrected E (need AM-GM domination):", extra)

# 2. (3.4.32)-(3.4.34): A^a N^n (AB)^u V^v <= A^a max{TN,ABV}^4 T^{-n} y^{2n+u+v-8},
#    y^{1/4} <= A << y^{1/2}: exponent of y is a/2+2n+u+v-8 (a>=0) or a/4+2n+u+v-8 (a<0). y = T^phi.
def yexp(a, n, u, v):
    return (sp.Rational(a, 2) if a >= 0 else sp.Rational(a, 4)) + 2 * n + u + v - 8

for name, E in (("thesis E", E_thesis), ("corrected E", E_fix)):
    rows = []
    for (a, n, u, v) in sorted(E):
        assert u + n - 4 >= 0 and v + n - 4 >= 0 and 0 <= n <= 4
        rows.append(((a, n, u, v), -n, yexp(a, n, u, v)))
    print(name)
    for r in rows:
        print("   ", r[0], " T^%s y^%s" % (r[1], r[2]))
    # max over E of T^{-n} y^{e} for y = T^phi, phi in [0, 8/13]: compare with T^{-2} y^{7/2}
    worst = []
    for r in rows:
        d = sp.simplify((r[1] + r[2] * ph) - (-2 + sp.Rational(7, 2) * ph))
        sol = sp.solve(sp.Eq(d, 0), ph)
        worst.append((r[0], sp.expand(d), sol))
    for w in worst:
        print("    exponent minus (-2 + 7phi/2):", w[0], w[1], " zero at phi =", w[2])
    ok = all(sp.expand(w[1]).subs(ph, sp.Rational(8, 13)) <= 0 and sp.expand(w[1]).subs(ph, 0) <= 0 for w in worst)
    print("    max is T^-2 y^(7/2) on 0<=phi<=8/13:", ok)

# 3. Lemma 3.1.4 with lambda = 1/T, Delta = T^{1-delta}: lambda^{-eta-5/2} (Delta)^{-(eta+7/2)} (logs dropped)
first = (et + sp.Rational(5, 2)) - (1 - de) * (et + sp.Rational(7, 2))
second = -(2 * sg - sp.Rational(7, 2)) - (1 - de) * (et + sp.Rational(5, 2))
print("Lemma 3.1.4 first assertion T-exponent:", sp.expand(first), " at eta=0:", sp.expand(first.subs(et, 0)))
print("Lemma 3.1.4 second assertion T-exponent:", sp.expand(second), " at eta=0:", sp.expand(second.subs(et, 0)))

# 4. error terms as T-exponents with y = T^theta (eta, eps -> 0)
E11 = sp.Rational(1, 2) + sp.Rational(7, 8) * th + first.subs(et, 0)
E12 = -sp.Rational(1, 2) + sp.Rational(7, 8) * th + second.subs(et, 0)
E2 = sp.Rational(5, 4) * th + first.subs(et, 0)
E32 = sp.Rational(7, 4) * th + first.subs(et, 0)
for nm, ex in (("E11", E11), ("E12", E12), ("E2", E2), ("E32", E32)):
    crit = sp.solve(sp.Eq(ex.subs({de: 0, sg: sp.Rational(1, 2)}), 0), th)
    print(nm, "T-exponent:", sp.expand(ex), "; at delta=0, sigma=1/2 negative iff theta <", crit)

# 5. section 3.6: given theta<4/7 choose delta, eps0 with (3.6.3)-(3.6.5); check (3.6.6),(3.6.7),(3.6.10)
import random
from fractions import Fraction as F
random.seed(0)
bad = 0
for _ in range(20000):
    t = F(random.randint(1, 570), 1000)              # theta in (0, 0.57]
    assert t < F(4, 7)
    d = (1 - F(7, 4) * t) / 14                       # 7 delta + 7 theta/4 = 1/2 + 7theta/8 < 1
    d = min(d, F(1, 4))
    # eps0 from (3.6.4),(3.6.5): 4 theta eps0 < delta and 7 delta + (7/4 + 4 eps0) theta < 1
    e = min(d / (8 * t), (1 - 7 * d - F(7, 4) * t) / (8 * t), F(1, 5))
    c1 = 7 * d + F(7, 4) * t < 1
    c2 = 4 * t * e < d
    c3 = 7 * d + (F(7, 4) + 4 * e) * t < 1
    # (3.6.6): -1/2 + 7d/2 + 7t/8 < -2 e t  ;  (3.6.7): -1 + 7d/2 + 7t/4 < -2 e t ; (3.6.10): delta/2 >= 2 e0 theta
    c4 = -F(1, 2) + F(7, 2) * d + F(7, 8) * t < -2 * e * t
    c5 = -1 + F(7, 2) * d + F(7, 4) * t < -2 * e * t
    c6 = d / 2 >= 2 * e * t
    if not all((c1, c2, c3, c4, c5, c6)):
        bad += 1
print("section 3.6 parameter choice: violations in 20000 random theta < 4/7:", bad)
# implications as linear algebra: (3.6.5) => (3.6.6) and (3.6.7)
lhs5 = 7 * de + (sp.Rational(7, 4) + 4 * e0) * th - 1          # < 0
lhs6 = (-sp.Rational(1, 2) + sp.Rational(7, 2) * de + sp.Rational(7, 8) * th) - (-2 * e0 * th)
lhs7 = (-1 + sp.Rational(7, 2) * de + sp.Rational(7, 4) * th) - (-2 * e0 * th)
print("(3.6.6) lhs = (3.6.5 lhs)/2 :", sp.simplify(lhs6 - lhs5 / 2) == 0)
print("(3.6.7) lhs - (3.6.5 lhs) =", sp.simplify(lhs7 - lhs5), "(<= 0 since delta, eps0, theta > 0)")

# 6. Ingham range: 3(1-s)/(2-s) < 1 + theta(1-2s) on (1/2, 1] for theta < 2/3 (concavity)
s = sp.symbols('s')
g = 3 * (1 - s) / (2 - s) - (1 + th * (1 - 2 * s))
print("Ingham gap g(1/2) =", sp.simplify(g.subs(s, sp.Rational(1, 2))), "; g'(1/2) =", sp.simplify(sp.diff(g, s).subs(s, sp.Rational(1, 2))),
      "; g'' =", sp.simplify(sp.diff(g, s, 2)), "; g(1) =", sp.simplify(g.subs(s, 1)))

# 7. Lemma 3.4.3 (type I) bookkeeping: J << (yABVq)^eps N^-eps (q A V^{3/2} + A B V^{1/2}) with
#    A <= y^{1/4}, AB <= y, q <= y, V <= y/q (exponents of y: a, b, qq, v >= 0). Maximise both exponents (LP by vertices).
from itertools import product
best1 = best2 = None
grid = [sp.Rational(k, 8) for k in range(0, 9)]
for a_, b_, q_, v_ in product([sp.Rational(k, 16) for k in range(0, 5)], grid, grid, grid):
    if a_ + b_ <= 1 and v_ <= 1 - q_:
        e1 = q_ + a_ + sp.Rational(3, 2) * v_
        e2 = a_ + b_ + v_ / 2
        best1 = e1 if best1 is None or e1 > best1 else best1
        best2 = e2 if best2 is None or e2 > best2 else best2
print("Lemma 3.4.3: max exponent of q A V^{3/2} =", best1, "; of A B V^{1/2} =", best2, "(thesis: 7/4 and 3/2)")
# S_2 / E_2: (y^{1/4} V)^{1+eta} summed over V <= y/q with weight 1/q -> y^{5/4}
print("E_2: y-exponent 1/4 + 1 =", sp.Rational(1, 4) + 1)
