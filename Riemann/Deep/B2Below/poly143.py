# Where PR 143's degree-9 polynomial comes from: closed form of the continuum gain + KKT contact of the scalar system.
import sympy as sp
c, d, l, t, b, rho = sp.symbols('c delta ell t b rho', real=True)
# (1) closed form of the ODE solution r z' = z - T, T = (14-18r-12z)/(30(3r-2)), z(7/9)=0
z_inf = (33*rho - 14)/48 - sp.Rational(35, 144)*(9*rho/7)**sp.Rational(4, 5)*(9*rho - 6)**sp.Rational(1, 5)
T = lambda r, z: (14 - 18*r - 12*z)/(30*(3*r - 2))
res = sp.simplify(rho*sp.diff(z_inf, rho) - (z_inf - T(rho, z_inf)))
print("ODE residual of closed form:", sp.nsimplify(sp.N(res.subs(rho, sp.Rational(3, 4)), 30)), " z(7/9) =", sp.simplify(z_inf.subs(rho, sp.Rational(7, 9))))
# (2) contact at x = 1/2: balance of the count <=> delta = thr(c, z(c)); b-coefficient of g vanishes <=> R = 5/9 + delta/6
x = sp.Rational(1, 2)
zc = sp.Symbol('zc')
V = sp.Rational(5, 6) - d
R = 1 - d + V*(sp.Rational(3, 2)*c - 1)
bal = sp.Eq((1 - d*(c + 2*x*((1 - c)/2 + zc))), R)
dsol_bal = sp.solve(bal, d)[0]
print("balance: delta =", sp.factor(dsol_bal), "  (= thr(c, z) = (5c/4-5/6)/(c-1/2-z))")
h = (1 + b + 3*l)/2; Tt = (1 + b - l) - h
g = sp.Rational(1, 4) - b/6 - 5*l/4 - d/2 - d*x*l + 3*h*(1 - R)/4 + d*x*Tt/4
gb = sp.expand(sp.diff(g, b))
dsol_kkt = sp.solve(sp.Eq(gb, 0), d)[0]
print("b-coefficient zero: delta =", sp.factor(dsol_kkt))
zsol = sp.solve(sp.Eq(dsol_bal, dsol_kkt), zc)[0]
print("=> z(c) =", sp.factor(zsol))
lsol = sp.solve(sp.Eq(g.subs(d, dsol_kkt), 0), l)[0]
tc = sp.factor(sp.simplify((sp.Rational(11, 3) - lsol)/4))
print("theta(c) =", tc)
csol = sp.solve(sp.Eq(t, tc), c)
print("c(theta) =", [sp.factor(s) for s in csol])
cs = csol[0]
# (3) fifth power: (144/35)((33c-14)/48 - z(c)) = (9c/7)^(4/5) (9c-6)^(1/5)
lhs = (sp.Rational(144, 35)*((33*c - 14)/48 - zsol))**5
rhs = (9*c/7)**4*(9*c - 6)
P = sp.factor(sp.numer(sp.together((lhs - rhs).subs(c, cs))))
print("eliminated polynomial in theta:", P)
Q = 3888*sp.Integer(19)**5*(15*t - 13)**4*(1 - t)**5 - 7*(7 - 6*t)**4*(9*t - 7)**5
print("ratio to PR143 polynomial:", sp.simplify(sp.expand(P)/sp.expand(Q)))
root = sp.nsolve(Q, t, 0.8733458)
print("root:", sp.N(root, 25), " c* =", sp.N(cs.subs(t, root), 15), " delta* =", sp.N(dsol_kkt.subs(c, cs.subs(t, root)), 15))
