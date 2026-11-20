# explore55: exact (sympy, rational) check of the recipe algebra and of the 2D ceiling 2tQ/(1+2tQ), 7/8 at theta = 1/2.
import sympy as sp
a, b, x, l, lam, t, z1, z2 = sp.symbols('alpha beta x ell lam theta z1 z2')
g = lambda s: sp.exp(-s*l)*(1 - s*lam)**2
D = (b**2*g(a) - a**2*g(b))/(a*b*(b - a))
dd = sp.limit(sp.limit(sp.diff(sp.series(sp.series(D, a, 0, 3).removeO(), b, 0, 3).removeO(), a, b), a, 0), b, 0)
print('d_a d_b D at 0 =', sp.expand(dd), '  (expect l^3/6 + l^2 lam + l lam^2)')
g1 = lambda s: sp.exp(-s*l)*(1 - s*lam)/s                       # first moment: swap term, should give l^2/2 + l lam
print('first moment d_a at 0 =', sp.expand(sp.series(g1(a) - 1/a + l + lam, a, 0, 2).removeO().coeff(a, 1)), '  (expect l^2/2 + l lam)')
L, I, Jp = sp.symbols('L I Jp', positive=True)
u = sp.Rational(1, 2) + t*I; A = u**2 + sp.Rational(1, 12) + Jp/(12*t); B = sp.expand((l**3/6 + l**2*lam + l*lam**2).subs({l: 1, lam: t*I}))
print('B - u^2 =', sp.simplify(B - u**2))
k2 = sp.simplify(2*u**2/(A + B)); print('kappa_2D =', k2, ' ; kappa_CGG =', sp.simplify(u**2/A))
R = (sp.Rational(1, 2) + t)*x - t*x**2/2                      # Riesz representer of l(P) = P(1)/2 + theta int P on {P(0)=0}, <P,Q> = int P'Q'
P = sp.Function('P')
print('Riesz check: R(0) =', R.subs(x, 0), ", -R'' =", -sp.diff(R, x, 2), ", R'(1) =", sp.diff(R, x).subs(x, 1))
F = sp.integrate(sp.diff(R, x)**2, (x, 0, 1)); Q = sp.factor(12*F); print('sup u^2/J\' = int R\'^2 =', sp.factor(F), ' Q =', sp.expand(Q))
kap = sp.simplify(F/(F + 1/(24*t))); print('kappa_2D ceiling =', sp.factor(kap))
for th in [sp.Rational(1, 4), sp.Rational(1, 2), 1]:
    PR = sp.expand(R/R.subs(x, 1)); Iv = sp.integrate(PR, (x, 0, 1)).subs(t, th); Jv = sp.integrate(sp.diff(PR, x)**2, (x, 0, 1)).subs(t, th)
    print(f'theta={th}: P_2D = {PR.subs(t, th)}, kappa_2D = {sp.nsimplify(kap.subs(t, th))}, check = {sp.simplify(k2.subs({t: th, I: Iv, Jp: Jv}))}')
# CGG/BHB optimum for comparison (P = -theta x^2 + (1+theta) x at theta = 1/2 gives 19/27)
PB = -t*x**2 + (1 + t)*x; IB = sp.integrate(PB, (x, 0, 1)); JB = sp.integrate(sp.diff(PB, x)**2, (x, 0, 1))
print('BHB at theta=1/2:', sp.simplify((u**2/A).subs({I: IB, Jp: JB}).subs(t, sp.Rational(1, 2))))
# CGG optimum over quadratic P = a x^2 + (1-a) x (P(1) = 1): claim max = 1 - (1+theta)^-3 at a = -theta
aa = sp.symbols('a'); Pq = aa*x**2 + (1 - aa)*x; Iq = sp.integrate(Pq, (x, 0, 1)); Jq = sp.integrate(sp.diff(Pq, x)**2, (x, 0, 1))
kc = sp.simplify((u**2/A).subs({I: Iq, Jp: Jq})); crit = sp.solve(sp.diff(kc, aa), aa)
print('CGG critical a:', crit, ' value at a=-theta:', sp.factor(sp.simplify(kc.subs(aa, -t))), ' = 1-(1+theta)^-3 ?', sp.simplify(kc.subs(aa, -t) - (1 - (1 + t)**-3)) == 0)
print('2D ceiling = 1-(1+2theta)^-3 ?', sp.simplify(kap - (1 - (1 + 2*t)**-3)) == 0)
