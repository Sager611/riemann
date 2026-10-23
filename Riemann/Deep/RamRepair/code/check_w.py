
# check_w.py -- independent checks of the BTY test-function ingredients (sympy exact + Arb balls)
# (1) w = sec^4(theta) * (g * g)(u), g(v) = cos(k v - theta) - cos(theta) on [0, D], k = tan(theta), D = 2 theta cot(theta):
#     autocorrelation of a nonnegative concave g => w >= 0, w even-extension positive definite (Re W(iy) >= 0), w' <= 0 on [0,D].
# (2) P(x) = |sum c_k e^{ikx}|^2/14912370: a_0 = 1, a_1 = 865534/497079, a = 2919857/828465, all a_k >= 0.
# (3) kappa = sum kappa_m = 433/859, p(x) = sum kappa_m x^m > 0 on [0, 1].
import sympy as sp
u, v, th = sp.symbols('u v theta', positive=True)
k = sp.tan(th); D = 2*th/sp.tan(th)
g = lambda x: sp.cos(k*x - th) - sp.cos(th)
h = sp.integrate(sp.expand_trig(g(v)*g(v + u)), (v, 0, D - u))
s2 = 1/sp.cos(th)**2
w = s2*(s2*(th/sp.tan(th) - u/2)*sp.cos(u*k) + 2*th/sp.tan(th) - u + sp.sin(2*th - u*k)/sp.sin(2*th) - 2*(1 + sp.sin(th - u*k)/sp.sin(th)))
diff = sp.simplify(sp.expand_trig(s2**2*h - w))
print('w - sec^4 (g*g) simplified:', diff)
import random
for _ in range(5):
    tv = random.uniform(0.3, 1.4); uv = random.uniform(0, 1)*float(D.subs(th, tv))
    print('  numeric check theta=%.3f u=%.3f: %.3e' % (tv, uv, float((s2**2*h - w).subs({th: tv, u: uv}).evalf(30))))
cc = [4,-8,2,20,-9,-34,27,91,-27,-201,32,895,1949,2389,1896,949,239]
from fractions import Fraction as Fr
ak = [Fr(sum(cc[j]*cc[j+kk] for j in range(17-kk))*(1 if kk == 0 else 2), 14912370) for kk in range(17)]
print('a0 =', ak[0], ' a1 =', ak[1], '==865534/497079:', ak[1] == Fr(865534, 497079), ' a =', sum(ak[1:]), '==2919857/828465:', sum(ak[1:]) == Fr(2919857, 828465))
print('all a_k >= 0:', all(x >= 0 for x in ak), [float(x) for x in ak])
km = [Fr(1), Fr(-851, 859), Fr(780, 859), Fr(-525, 859), Fr(171, 859), Fr(28, 859), Fr(-29, 859)]
print('kappa =', sum(km), '== 433/859:', sum(km) == Fr(433, 859), ' sum|kappa_m| =', float(sum(abs(x) for x in km)))
x = sp.symbols('x'); p = sum(sp.Rational(q.numerator, q.denominator)*x**m for m, q in enumerate(km))
print('real roots of p in [0,1]:', sp.Poly(p, x).count_roots(0, 1), ' p(0)=', p.subs(x, 0), ' p(1)=', p.subs(x, 1))
bm = [km[0]] + [km[m] + km[m-1] for m in range(1, 7)] + [km[6]]
print('b_m (1+x)p(x):', [str(b) for b in bm])
