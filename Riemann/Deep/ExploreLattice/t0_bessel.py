# t0_bessel.py (explore48): timing/accuracy of arb K_{it}(x) (imaginary order, real argument) via python-flint.
import time
from flint import acb, arb, ctx
for t in (50.0, 150.0, 400.0):
    for fx in (0.5, 0.95, 1.0, 1.05, 2.0):
        x = fx*t
        for prec in (int(2.4*t)+80,):
            ctx.prec = prec
            t0 = time.time()
            v = acb(x).bessel_k(acb(0, t))
            dt = time.time()-t0
            # scale by e^{pi t/2}
            sc = v * (arb.pi()*t/2).exp()
            ctx.prec = 53
            print('t=%g x=%g prec=%d  e^{pi t/2}K = %s  time %.4f s' % (t, x, prec, sc.str(5, radius=True)[:60], dt))
