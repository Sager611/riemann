"""gapcertaudit: independent derivation of Lemmas 3 and 4 (J = 5, 6 per-gap reductions) with sympy, from the coefficient maps
  f ~ b (sine), f' ~ w b (cos), f'' ~ -w^2 b (sine), f''' ~ d - w^3 b (cos, const dA/g), f'''' ~ w^4 b - w d (sine),
  f^(5) ~ e + w^5 b - w^2 d (cos, const dB/g), f^(6) ~ -(w e + w^6 b - w^3 d) (sine),
d = (2/g)((-1)^k A1 - A0), e = (2/g)((-1)^k B1 - B0), A = f'', B = f'''' at 0, g; bracket [F_p F_q]_0^g = const_p dQ + (g/2) sum coef_p q_k.
Checks: (i) per-mode coefficients vs exact6.py's qd, qe, gdd, gde, gee (as polynomials in x = w^2); J = 5 by a6 = phi52 = phi54 = 0
vs exact5.py; (ii) numer6 / numer5 output vs the numerators recomputed here, exactly, for the certified J = 6 form; (iii) degrees.
usage (from a copy of ExploreShape): lemma46.py logs/cert5_J6m80_3.367.json"""
import sys, json, os
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import sympy as sp
from exact6 import numer6
from exact5 import numer5

w, b, d, e, x = sp.symbols('w b d e x')
a = sp.symbols('a0:7'); f12, f14, f32, f34, f52, f54 = sp.symbols('f12 f14 f32 f34 f52 f54')
c = [b, w * b, -w ** 2 * b, d - w ** 3 * b, w ** 4 * b - w * d, e + w ** 5 * b - w ** 2 * d, -(w * e + w ** 6 * b - w ** 3 * d)]
Q = sp.expand(sum(a[j] * c[j] ** 2 for j in range(7)) + f12 * c[1] * d + f14 * c[1] * e + f32 * c[3] * d + f34 * c[3] * e
              + f52 * c[5] * d + f54 * c[5] * e)
P = Q.coeff(b, 2); B1 = Q.coeff(b, 1); Bd = B1.coeff(d, 1); Be = B1.coeff(e, 1)
Q0 = Q.coeff(b, 0); Gdd = Q0.coeff(d, 2); Gee = Q0.coeff(e, 2); Gde = Q0.coeff(d, 1).coeff(e, 1)
assert sp.expand(B1 - Bd * d - Be * e) == 0 and sp.expand(Q0 - Gdd * d ** 2 - Gde * d * e - Gee * e ** 2) == 0
X = lambda expr: sp.expand(expr.subs(w, sp.sqrt(x)))
qd = sp.expand(X(Bd / w)); qe = sp.expand(X(Be / w))
ref = {'p': sum(a[j] * x ** j for j in range(7)),
       'qd': f12 - (2 * a[3] + f32) * x + (f52 - 2 * a[4]) * x ** 2 - 2 * a[5] * x ** 3 - 2 * a[6] * x ** 4,
       'qe': f14 - f34 * x + (2 * a[5] + f54) * x ** 2 + 2 * a[6] * x ** 3,
       'gdd': a[3] + f32 + (a[4] - f52) * x + a[5] * x ** 2 + a[6] * x ** 3,
       'gde': f34 + f52 - (2 * a[5] + f54) * x - 2 * a[6] * x ** 2, 'gee': a[5] + f54 + a[6] * x}
mine = {'p': X(P), 'qd': qd, 'qe': qe, 'gdd': X(Gdd), 'gde': X(Gde), 'gee': X(Gee)}
for k in ref:
    print('(i) %s: own derivation == exact6 formula: %s' % (k, sp.expand(mine[k] - ref[k]) == 0))
# constant (zero-mode) part: a3|f'''|^2, a5|f^(5)|^2 constants and the bracket constants const_p * dQ
dA, dB = sp.symbols('dA dB')
const = a[3] * dA ** 2 + a[5] * dB ** 2 + f32 * dA * dA + f34 * dA * dB + f52 * dB * dA + f54 * dB * dB
Cref = sp.Matrix([[a[3] + f32, (f34 + f52) / 2], [(f34 + f52) / 2, a[5] + f54]])
print('(i) C: own == exact6:', sp.expand(const - (sp.Matrix([[dA, dB]]) * Cref * sp.Matrix([dA, dB]))[0]) == 0)
N11 = sp.expand(4 * mine['p'] * mine['gdd'] - x * qd ** 2); N12 = sp.expand(2 * mine['p'] * mine['gde'] - x * qd * qe)
N22 = sp.expand(4 * mine['p'] * mine['gee'] - x * qe ** 2); DET = sp.expand(N11 * N22 - N12 ** 2)
js = json.load(open(sys.argv[1]))
av = [Fr(s) for s in js['a']]; ph = {k: Fr(v) for k, v in js['phi'].items()}
sub = {a[j]: sp.Rational(av[j].numerator, av[j].denominator) for j in range(7)}
for k, s in (('12', f12), ('14', f14), ('32', f32), ('34', f34), ('52', f52), ('54', f54)):
    v = ph.get(k, Fr(0)); sub[s] = sp.Rational(v.numerator, v.denominator)
pE, N11E, N12E, N22E, DETE, CE = numer6(av, ph)
for name, mineP, theirs in (('N11', N11, N11E), ('N12', N12, N12E), ('N22', N22, N22E), ('DET', DET, DETE)):
    pm = sp.Poly(mineP.subs(sub), x); cm = [Fr(int(sp.fraction(cf)[0]), int(sp.fraction(cf)[1])) for cf in reversed(pm.all_coeffs())]
    print('(ii) %s: degree own %d, exact6 %d, identical coefficients: %s' % (name, len(cm) - 1, len(theirs) - 1, cm == [Fr(t) for t in theirs]))
# J = 5: a6 = phi52 = phi54 = 0 against exact5.numer5 (same form with a6, phi52, phi54 dropped)
sub5 = dict(sub); sub5[a[6]] = 0; sub5[f52] = 0; sub5[f54] = 0
ph5 = {k: v for k, v in ph.items() if k not in ('52', '54')}
out5 = numer5(av[:6], ph5)
for name, mineP, theirs in (('N11', N11, out5[1]), ('N12', N12, out5[2]), ('N22', N22, out5[3])):
    pm = sp.Poly(mineP.subs(sub5), x); cm = [Fr(int(sp.fraction(cf)[0]), int(sp.fraction(cf)[1])) for cf in reversed(pm.all_coeffs())]
    print('(ii) J=5 %s: identical to exact5: %s' % (name, cm == [Fr(t) for t in theirs]))
print('a6 =', float(av[6]), '; C =', [[float(CE[i][j]) for j in range(2)] for i in range(2)])
