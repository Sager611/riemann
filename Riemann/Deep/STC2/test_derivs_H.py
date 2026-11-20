# test_derivs_H.py -- check certK2.py value, gradient and Hessian (vg, hess) of H_s against mpmath finite differences (L20)
import sys, json
sys.argv = ['x', 'L20', 'edge']
exec(open('Riemann/Deep/STC2/certK2.py').read().split("t0 = time.time()")[0])
from flint import arb, acb
import mpmath as mp
mp.mp.dps = 30
kk = KJ['L20']; dd = mp.mpf(6)/5
def Hf(s, b, t):
    tot = 0
    for x in (dd+b, dd-b):
        z = mp.mpc(x, t)
        for (X, j), c in zip(kk['B'], kk['c']):
            cs = c if (s == 1 or X == 'R') else -c
            w = z**(-j); tot += cs*(w.real if X == 'R' else -w.imag)
    f = mp.pi - sum(a*mp.atan(t/y) for a, y in ((2, dd+b), (2, dd-b), (-1, 2*dd+b), (-1, 2*dd-b)))
    return tot - s*f
for s in (1, -1):
    for (b, t) in ((0.3, 0.7), (0.5, 1.006), (0.1, 5.0)):
        H, Hb, Ht = vg(s, arb(b), arb(t))
        e = 1e-6
        nb = (Hf(s, b+e, t) - Hf(s, b-e, t))/(2*e); nt = (Hf(s, b, t+e) - Hf(s, b, t-e))/(2*e)
        hbb, hbt = hess(s, arb(b), arb(t))
        nbb = (Hf(s, b+e, t) - 2*Hf(s, b, t) + Hf(s, b-e, t))/e**2; ntt = (Hf(s, b, t+e) - 2*Hf(s, b, t) + Hf(s, b, t-e))/e**2
        nbt = (Hf(s, b+e, t+e) - Hf(s, b+e, t-e) - Hf(s, b-e, t+e) + Hf(s, b-e, t-e))/(4*e*e)
        print(s, b, t, 'H %.3e vs %.3e' % (float(H.mid()), Hf(s, b, t)), 'Hb %.6g/%.6g Ht %.6g/%.6g Hbb %.5g/%.5g Htt %.5g Hbt %.5g/%.5g' % (float(Hb.mid()), nb, float(Ht.mid()), nt, float(hbb.mid()), nbb, ntt, float(hbt.mid()), nbt))

