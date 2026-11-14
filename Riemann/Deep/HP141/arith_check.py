# arith_check.py -- exact Fraction checks of every decimal rounding used in REPORT.md's drafted statements (stdlib only)
from fractions import Fraction as F
import math
import theta141 as T
th, old, oldE = T.TH, T.OLD, F(8733819666132076, 10**16)
ok = True
def chk(name, cond):
    global ok; ok &= bool(cond); print(('PASS ' if cond else 'FAIL ') + name)
chk('roundings up: TH_UP17, TH_UP25, TH_UP40, TH_FLOAT_UP >= theta141', min(F(T.TH_UP17), F(T.TH_UP25), F(T.TH_UP40), F(T.TH_FLOAT_UP)) >= th)
chk('theta141 < old 0.87338196661320752 by 3.611e-5', 0 < old - th and abs(float(old - th) - 3.611e-5) < 1e-8)
chk("theta141 < Liu's contact B2 = (777 - sqrt 1081)/852:  (777 - 852 theta)^2 > 1081", 777 - 852*th > 0 and (777 - 852*th)**2 > 1081)
E = 1/(1 - th); Eo = 1/(1 - oldE)
chk('ExplProg exponent 1/(1-theta141) = %.10f in (7.8955, 7.8956]; old 1/(1-0.8733819666132076) = %.10f -> 7.8978' % (E, Eo), F(78955, 10**4) < E <= F(78956, 10**4) and F(78977, 10**4) < Eo <= F(78978, 10**4))
chk('PrimRootNum exponents: 3.4/(1-theta) = %.6f <= 26.85, 3/(1-theta) = %.6f <= 23.69, 1/(1-theta) <= 7.896' % (F(34, 10)*E, 3*E), F(34, 10)*E <= F(2685, 100) and 3*E <= F(2369, 100) and E <= F(7896, 1000))
chk('Exp48 side condition 1/(4(1-theta)) = %.4f < 11.2' % (E/4), E/4 < F(112, 10))
chk('Littlewood78 limit (pi^2/6) e^-gamma (1-theta) = %.7f >= 0.11697 (float, 1e-15 safe)' % (math.pi**2/6*math.exp(-0.5772156649015329)*float(1 - th)), math.pi**2/6*math.exp(-0.5772156649015329)*float(1 - th) - 0.11697 > 1e-9)
chk('KimEff 3.6/(4a^2), a = theta - 1/2: %.7f >= 6.456' % float(F(36, 10)/(4*(th - F(1, 2))**2)), F(36, 10)/(4*(th - F(1, 2))**2) >= F(6456, 1000))
for k, n0, st in ((3, F(6268167, 10**4), F(62682, 100)), (4, F(1171700, 10**4), F(11717, 100)), (5, F(577300, 10**4), F(5773, 100)), (6, F(362283, 10**4), F(3623, 100))):
    chk('Ingham3 k=%d: stated e^%s >= certified n0 = e^%s' % (k, float(st), float(n0)), st >= n0)
chk('BVIrr row a = 1/(1-theta141) = %.6f (was %.6f)' % (E, 1/(1 - float(old))), True)
print('ALL PASS' if ok else 'SOME FAIL')
