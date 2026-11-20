# mpcheck.py -- mpmath cross-check (own code): |1/zeta| at the refined peaks of falsify.py, at zero ordinates, and at height T' = 3e12.
import math, time, json
import mpmath
from mpmath import mp, mpf, mpc, zeta
mp.dps = 20
BP = 0.87334585610671858; D = 0.875 - BP
lC = {'ML092g8_e05_Bp': (21.49078307, 21.46999417, 0.10644111, 0.92 - D), 'ML093g5_e10_Bp': (15.00367661, 14.9866546, 0.13267488, 0.93 - D)}
lCL, aL = 10.048373, 0.29142857
t0 = time.time()
for nm, (l1000, lTp, a, s) in lC.items():
    for t in [14.1217 if nm.startswith('ML092') else 14.1208, 1083.2225 if nm.startswith('ML092') else 1083.2247]:
        v = 1/abs(zeta(mpc(s, t)))
        print(nm, 's', round(s, 8), 't', t, '|1/zeta| =', mpmath.nstr(v, 10), 'ratio to C(1000) t^a =', '%.4e' % (float(v)/(math.exp(l1000)*t**a)) if t >= 1000 else '', flush=True)
v = 1/abs(zeta(mpc(0.85, 1083.2099))); print('sigma 0.85 t 1083.2099 |1/zeta| =', mpmath.nstr(v, 10), 'ratio to C_L t^a_L = %.4e' % (float(v)/(math.exp(lCL)*1083.2099**aL)), flush=True)
# at the ordinates of the first 30 zeros (where |1/zeta| peaks), the line of the final certificate
s = 0.92 - D; mx = 0
for n in range(1, 31):
    g = mpmath.zetazero(n).imag; v = float(1/abs(zeta(mpc(s, g)))); mx = max(mx, v)
print('max |1/zeta(%.6f + i gamma_n)|, n <= 30:' % s, round(mx, 6), round(time.time() - t0, 1), 's', flush=True)
# at height T' = 3e12 (Riemann-Siegel in mpmath): the three Gamma pieces meet here
TP = 3*10**12
for sig, lab, bnd in [(0.85, 'RH strip, C_L t^a_L', math.exp(lCL)*TP**aL), (0.918345856106719, 'final cert line, C_inf(T\') t^a', math.exp(lC['ML092g8_e05_Bp'][1])*TP**lC['ML092g8_e05_Bp'][2])]:
    for dt in [0.0, 0.37]:
        t1 = time.time(); v = float(1/abs(zeta(mpc(sig, TP + dt))))
        print('t = 3e12 +', dt, 'sigma', sig, '|1/zeta| =', round(v, 6), 'bound (%s) = %.4e' % (lab, bnd), 'ratio %.3e' % (v/bnd), round(time.time() - t1, 1), 's', flush=True)
