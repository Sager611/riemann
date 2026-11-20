# explore55 follow-up: general p = P(1) recipe algebra (referee's g(x) = e^{-x l}(p - x lam)^2, u = p/2 + theta int P).
import sympy as sp
a, b, l, lam, p, t, I, Jp = sp.symbols('alpha beta ell lam p theta I Jp')
g = lambda s: sp.exp(-s*l)*(p - s*lam)**2
D = (b**2*g(a) - a**2*g(b))/(a*b*(b - a))
dd = sp.expand(sp.limit(sp.limit(sp.diff(sp.series(sp.series(D, a, 0, 3).removeO(), b, 0, 3).removeO(), a, b), a, 0), b, 0))
print('d_a d_b D at 0 =', dd)
g1 = lambda s: sp.exp(-s*l)*(p - s*lam)/s
print('first moment =', sp.expand(sp.series(g1(a) - p/a, a, 0, 2).removeO().coeff(a, 1)))
u = p/2 + t*I; B = dd.subs({l: 1, lam: t*I}); A = u**2 + p**2/12 + Jp/(12*t)
print('B - u^2 =', sp.simplify(B - u**2), ';  A expanded =', sp.expand(A))
print('Var Re Y = (A+B)/2 - u^2 =', sp.simplify((A + B)/2 - u**2), '; Var Im Y = (A-B)/2 =', sp.simplify((A - B)/2))
x = sp.symbols('x'); P = (4*x - x**2)/3; th = sp.Rational(1, 2)
uu = P.subs(x, 1)/2 + th*sp.integrate(P, (x, 0, 1)); JJ = sp.integrate(sp.diff(P, x)**2, (x, 0, 1))
print('theta=1/2, P=(4x-x^2)/3: u =', uu, ' int P\'^2 =', JJ, ' kappa2 =', sp.nsimplify(uu**2/(uu**2 + JJ/(24*th))))
