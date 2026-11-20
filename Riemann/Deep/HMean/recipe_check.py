# hmean: sympy check of the recipe algebra (independent of explore55's exact.py).
import sympy as sp
a, b, z, l, lam, kap, mu = sp.symbols('alpha beta z ell lambda kappa mu')
g = lambda x: sp.exp(-x*l)*(1 - x*lam)**2
D = (b**2*g(a) - a**2*g(b))/(a*b*(b - a))
# (i) H main term: d_a d_b D at 0 (limit along a generic ray)
s = sp.series(sp.series(D, a, 0, 3).removeO(), b, 0, 3).removeO()
cab = sp.expand(s).coeff(a, 1).coeff(b, 1)
print('(i) d_a d_b D|0 =', sp.factor(sp.simplify(cab)), ' expected l^3/6+l^2 lam+l lam^2:', sp.simplify(cab - (l**3/6 + l**2*lam + l*lam**2)) == 0)
# (ii) Phi(beta) = coefficient of alpha^1 in D + (zeta'/zeta)(1+alpha) ~ -1/alpha + ... (only the alpha-linear part of D matters)
Da = sp.series(D, a, 0, 2).removeO()
Phi = sp.simplify(sp.expand(Da).coeff(a, 1))
g1, g2 = -(l + 2*lam), l**2/2 + 2*lam*l + lam**2
print('(ii) Phi - [1/b^2 + g1/b + g2 - g(b)/b^2] =', sp.simplify(Phi - (1/b**2 + g1/b + g2 - g(b)/b**2)))
print('     Phi regular at 0:', sp.limit(Phi, b, 0), '(finite)')
Psi = l*Phi + sp.diff(Phi, b)
# Perron with K = e^{kappa l}: polynomial part -> residue at 0 of K^z/z * poly(z); exponential part (K/tau)^z -> 0 for kappa<1
poly = sp.expand(l*(1/b**2 + g1/b + g2) + sp.diff(1/b**2 + g1/b + g2, b))
K = sp.exp(kap*l*z)
res = sp.residue(sp.expand(K/z*poly.subs(b, z)), z, 0)
f = sp.simplify((res/l**3).subs(lam, mu*l))
print('(iii) f(kappa) =', sp.expand(f))
print('     f - [1/6+mu+mu^2+mu(1-k)^2+(1-k)^3/3] =', sp.simplify(f - (sp.Rational(1,6) + mu + mu**2 + mu*(1-kap)**2 + (1-kap)**3/3)))
print('     f(1) =', sp.factor(f.subs(kap, 1)), '; u^2-1/12 with u=1/2+mu:', sp.expand((sp.Rational(1,2)+mu)**2 - sp.Rational(1,12)))
print('     f(0) =', sp.expand(f.subs(kap, 0)))
# (iv) BHB S1 from the same recipe: l + lam + [e^{-a l}(1-a lam) - 1]/a, d_a at 0
S1 = sp.series((sp.exp(-a*l)*(1 - a*lam) - 1)/a, a, 0, 2).removeO().coeff(a, 1)
print('(iv) S1/(T/2pi) =', sp.expand(S1), ' = L^2(1/2 + mu) with lam = mu L')
