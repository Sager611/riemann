# (a) exact: z_inf(r) = (33r-14)/48 - (35/144)(9r/7)^(4/5)(9r-6)^(1/5) solves r z' = z - T, T = (14-18r-12z)/(30(3r-2)), z(7/9) = 0,
#     and z_inf(r) < 7/6 - 3r/2 (T > 0, thr < 1/2) <=> r < 7/9 on (2/3, 7/9].
# (b) hypothetical twisted-family saving Z^(gamma tau): T = max(T_0, gamma/(2k)) (linear-fractional in rho' = tau/u', optimum at rho' in {0,1}); start 7/9.
import sympy as sp
r = sp.symbols('r', positive=True)
z = sp.Rational(33,48)*r - sp.Rational(14,48) - sp.Rational(35,144)*(9*r/7)**sp.Rational(4,5)*(9*r-6)**sp.Rational(1,5)
T = (14 - 18*r - 12*z)/(30*(3*r-2))
res = sp.simplify(r*sp.diff(z, r) - z + T)
print("(a) ODE residual simplified:", res, "| z(7/9) =", sp.nsimplify(z.subs(r, sp.Rational(7,9))))
for x in (0.70, 0.7335, 0.7523, 0.77, 0.7777):
    print("    r=%.4f residual=%.2e  z=%.6e  7/6-3r/2=%.6e" % (x, float(res.subs(r, x)) if res != 0 else 0.0, float(z.subs(r, x)), 7/6-1.5*x))
# margin: (7/6 - 3r/2) - z = (35/48)(3r-2)... reduce: z < 7/6-3r/2  <=>  3(3r-2) < 9r/7 (shown by hand); check the equivalence numerically
import math
ok = all(((7/6-1.5*x) - float(z.subs(r, x)) > 0) == (3*(3*x-2) < 9*x/7) for x in [0.667 + i*0.0001 for i in range(1, 1110)])
print("    sign of margin == sign of [3(3r-2) < 9r/7] on (0.6671, 0.7780):", ok)
