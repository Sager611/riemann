
import numpy as np
from flint import acb, arb, ctx, dirichlet_char
ctx.prec = 96
PI = acb(arb.pi())
ch = {n: dirichlet_char(13, n) for n in (12, 3, 9)}
def gam(s): return (acb(13)/PI)**(s/2) * (s/2).gamma()
def mid(z): return acb(z.real.mid(), z.imag.mid())
def kappa(n, t0=10.0):
    s0 = acb(0.5, t0); e0 = gam(s0); e0 = e0/abs(e0)
    return mid(acb(ch[n].hardy_z(acb(t0)).real) / (e0 * ch[n].l(s0)))
K = {n: kappa(n) for n in ch}
def f(c, s):
    cc, sc = acb(float(np.cos(c))), acb(float(np.sin(c)))
    return cc*K[12]*ch[12].l(s) + sc*(K[3]*ch[3].l(s) + K[9]*ch[9].l(s))
h = acb(1e-12)
def newton(c, s, it=60):
    for _ in range(it):
        fs = mid(f(c, s)); d = mid((mid(f(c, s+h)) - mid(f(c, s-h)))/(2*h))
        step = fs/d
        if float(abs(step).mid()) > 0.05: step = step*(0.05/float(abs(step).mid()))
        s = mid(s - step)
        if float(abs(step).mid()) < 1e-14: break
    return s, float(abs(f(c, s)).mid())
for c, starts in ((-0.003, (597.597, 710.136)), (-0.0022, (597.597, 710.136)), (0.0131, (604.0,)), (0.0128, (604.0,))):
    for t in starts:
        best = None
        for d0 in (0.01, 0.03, 0.06):
            s, res = newton(c, acb(0.5 + d0, t))
            if res < 1e-20 and abs(float(s.imag.mid()) - t) < 0.5:
                if best is None or float(s.real.mid()) > best[0]: best = (float(s.real.mid()), float(s.imag.mid()), res)
        print("c = %+.4f start t=%.3f -> zero %s" % (c, t, ("Re s = %.7f, Im s = %.6f, |f| = %.1e" % best) if best else "none converged"))
