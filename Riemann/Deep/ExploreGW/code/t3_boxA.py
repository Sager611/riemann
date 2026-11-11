
# box test for the mod-13 DH-type block A, which has a zero with beta > 1 at 1.183 + 77.19i
import numpy as np
from flint import acb, arb, ctx, dirichlet_char
ctx.prec = 96
PI = acb(arb.pi())
ch = {n: dirichlet_char(13, n) for n in (3, 9)}
def gam(s): return (acb(13)/PI)**(s/2) * (s/2).gamma()
def mid(z): return acb(z.real.mid(), z.imag.mid())
def kappa(n, t0=10.0):
    s0 = acb(0.5, t0); e0 = gam(s0); e0 = e0/abs(e0)
    return mid(acb(ch[n].hardy_z(acb(t0)).real) / (e0 * ch[n].l(s0)))
K = {n: kappa(n) for n in ch}
def xiA(w): return gam(w) * (K[3]*ch[3].l(w) + K[9]*ch[9].l(w))
h = acb(1e-12)
us = np.arange(0.05, 160.0, 0.05)
g = []
for u in us:
    w = acb(1.0, float(u)); g.append(float((mid(xiA(w+h)) - mid(xiA(w-h))).real.mid()/(2e-12)) if False else 0.0)
# log-derivative real part
def ld(w):
    a = mid(xiA(w)); d = mid((mid(xiA(w+h)) - mid(xiA(w-h)))/(2*h)); return float((d/a).real.mid())
g = np.array([ld(acb(1.0, float(u))) for u in us])
neg = us[g < 0]
print("block A: min Re(xi'/xi)(1+iu) on (0,160) = %.4f at u = %.2f; negative on %d grid points, u-range %s" % (g.min(), us[np.argmin(g)], len(neg), (round(float(neg.min()),2), round(float(neg.max()),2)) if len(neg) else None))
ts = np.arange(0.02, 80.0, 0.02)
X = []
for t in ts:
    v = xiA(acb(1.0, 2*float(t))); X.append(complex(float(v.real.mid()), float(v.imag.mid())))
X = np.array(X); ph = np.unwrap(np.angle(X))
b = np.abs(X)*np.cos(ph)
sc = int(np.sum(np.abs(np.diff(np.signbit(b).astype(int)))))
print("box block A a=1 on (0,80]: sign changes %d ; phase increase/pi %.3f ; min dphase/dt %.3f at t=%.2f" % (sc, (ph[-1]-ph[0])/np.pi, np.min(np.diff(ph))/0.02, ts[np.argmin(np.diff(ph))]))
