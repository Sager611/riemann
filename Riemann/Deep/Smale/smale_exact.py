#!/usr/bin/env python3
"""Exact check of the univalence certificate of arXiv 2610.21472, Thm 6.1 (pair x=9/16, y=3/8).
Recomputes H = Q/B_a, 1/B(0)^2, the Schur-Cohn data of the Dieudonne polynomials P_c and
Sturm root counts on [-1,1]; special values c = 0, +-1/sqrt2.  Needs sympy.  ~15 s, < 100 MB."""
import sympy as sp
I = sp.I
x, y = sp.Rational(9, 16), sp.Rational(3, 8)
z, c = sp.symbols('z c')
u = -x - I*y
Q = sp.expand((z - 1)**4 - u**4)
a, b = 1 + u, 1 - I*u                      # the two zeros of Q in the disk
q, r = sp.div(sp.Poly(Q, z), sp.Poly(z - a, z)); assert r.is_zero
H = sp.Poly(sp.expand((1 - sp.conjugate(a)*z)*q.as_expr()), z)
h = [sp.nsimplify(sp.expand(H.coeff_monomial(z**k))) for k in range(5)]
paper = [-sp.Rational(6877, 4096) - sp.Rational(2253, 2048)*I, sp.Rational(302959, 65536) + sp.Rational(37401, 16384)*I,
         -sp.Rational(20499, 4096) - sp.Rational(5121, 2048)*I, sp.Rational(619, 256) + sp.Rational(3, 2)*I,
         -sp.Rational(7, 16) - sp.Rational(3, 8)*I]
print('H coefficients equal the paper:', all(sp.simplify(h[k] - paper[k]) == 0 for k in range(5)))
J = 1 + sp.expand((1 - a*sp.conjugate(a))*(1 - b*sp.conjugate(b)))/sp.expand((a - b)*sp.conjugate(a - b))
print('1/B(0)^2 =', sp.nsimplify(sp.simplify(J)), ', 1/B(0) =', sp.N(sp.sqrt(sp.nsimplify(sp.simplify(J))), 12))
U = [sp.Integer(1), 2*c, 4*c**2 - 1, 8*c**3 - 4*c]          # Chebyshev U_0..U_3
al = [sp.expand(h[k + 1]*U[k]) for k in range(4)]             # Dieudonne polynomial P_c, coefficients of w^0..w^3
cj = lambda e: sp.expand(sp.conjugate(e).subs(sp.conjugate(c), c))   # c real
def T(co):                                                    # Schur-Cohn transform
    n = len(co) - 1
    out = [sp.expand(cj(co[0])*co[k] - co[n]*cj(co[n - k])) for k in range(n + 1)]
    assert sp.simplify(out[n]) == 0
    return out[:n]
T1 = T(al); T2 = T(T1); T3 = T(T2)
def sturm(name, e):
    assert sp.expand(e - cj(e)) == 0
    P = sp.Poly(sp.expand(e), c)
    print('%-12s degree %2d, roots in [-1,1]: %d, value at c=1: %.6g' % (name, P.degree(), P.count_roots(-1, 1), float(P.eval(1))))
for name, e in [('delta1', T1[0]), ('delta2', T2[0]), ('delta3', T3[0]),
                ('|lead T1|^2', T1[-1]*cj(T1[-1])), ('|lead T2|^2', T2[-1]*cj(T2[-1]))]:
    sturm(name, e)
print('delta1 * 2^32 == paper D1:', sp.expand(T1[0]*2**32 - (114165512497 - 22817013760*c**2*(1 - 2*c**2)**2)) == 0)
print('only root of |lead T2|^2 in [-1,1]:', sorted(set(sp.Poly(sp.expand(T2[-1]*cj(T2[-1])), c).real_roots()) & {0}))
# degenerate values of c (leading coefficient U_3(c) = 0)
a1, a2, a3 = h[1], h[2], h[3]
print('c=0: P_0 = a1 - a3 w^2, |a1|^2-|a3|^2 =', sp.nsimplify(sp.expand(a1*sp.conjugate(a1) - a3*sp.conjugate(a3))))
for s in (1, -1):
    co = [a1, s*sp.sqrt(2)*a2, a3]
    t1 = [sp.expand(sp.conjugate(co[0])*co[k] - co[2]*sp.conjugate(co[2 - k])) for k in range(2)]
    d2 = sp.simplify(t1[0]**2 - t1[1]*sp.conjugate(t1[1]))
    print('c=%s/sqrt2: delta1 = %s, lead != 0: %s, delta2 = %s' % (s, sp.nsimplify(sp.simplify(t1[0])), sp.simplify(t1[1]) != 0, sp.nsimplify(d2)))
