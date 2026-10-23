# zfr_lemma5.py -- Lemma 5 (positivity Re F(z) + Re F(2 sigma - 1 - conj z) >= 0 on 0 <= Re z <= 2 sigma - 1), own check:
# on Re z = 0, with c = 2 sigma - 1, y' = y/c:  (c/(eta w0)) sum_{m>=1} b_m Re W((cm + iy)/eta) >= B(y') =
#   sum_{m=1}^{7} b_m (m - sgn(b_m) eps0/m)/(m^2 + y'^2),  eps0 = Cw(c/eta0) eta0^2/(w0 c^2)  (|W0(z)| <= Cw|z|^-3, |z| >= cm/eta),
# the m = 0 term Re W(iy/eta) >= 0 (w = sec^4 * autocorrelation, check_w3.py).  B decreasing in eps0; certify B > 0 at eps0 = 1/500.
import sympy as sp
from fractions import Fraction as Fr
KM = [Fr(1), Fr(-851, 859), Fr(780, 859), Fr(-525, 859), Fr(171, 859), Fr(28, 859), Fr(-29, 859)]
b = [KM[0]] + [KM[m] + KM[m-1] for m in range(1, 7)] + [KM[6]]
y = sp.symbols('y', real=True)
for eps in (sp.Rational(1, 500), sp.Rational(1, 400)):
    B = sum(sp.Rational(b[m].numerator, b[m].denominator)*(m - (1 if b[m] > 0 else -1)*eps/m)/(m*m + y**2) for m in range(1, 8))
    num, den = sp.fraction(sp.together(B))
    P = sp.Poly(sp.expand(num), y)
    print('eps0 =', eps, ' real roots of numerator:', P.count_roots(), ' B(0) =', float(B.subs(y, 0)), ' sign of leading coeff', sp.sign(P.LC()),
          ' min over grid', min(float(B.subs(y, t/10)) for t in range(0, 400)))
