
import json, numpy as np
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
res = {}
for c in (0.5, 1.571):
    ws = json.load(open("t2_ws_%.3f.json" % c))
    zs = []
    for t in ws:
        for d0 in (0.05, 0.15, 0.35):
            s, r = newton(c, acb(0.5 + d0, t))
            b, g = float(s.real.mid()), float(s.imag.mid())
            if r < 1e-20 and b > 0.5 + 1e-8 and abs(g - t) < 1.0:
                if all(abs(g - z[1]) > 1e-6 for z in zs): zs.append((b, g))
                break
    zs.sort(key=lambda z: z[1])
    res[str(c)] = zs
    bb = np.array([z[0] for z in zs]); gg = np.array([z[1] for z in zs])
    L = np.log(13*gg/(2*np.pi))
    print("c=%.3f: %d distinct off-line zeros (Re s > 1/2) found from %d wrong-sign extrema" % (c, len(zs), len(ws)))
    for lo, hi in ((0, 200), (200, 400), (400, 600), (600, 800)):
        m = (gg >= lo) & (gg < hi)
        if m.any():
            print("   t in [%d,%d): n=%d, mean(beta-1/2)=%.3f, max=%.3f, mean (beta-1/2)*log(13t/2pi)/(2pi) [in mean spacings] = %.3f, #beta>1: %d" % (lo, hi, m.sum(), (bb[m]-0.5).mean(), (bb[m]-0.5).max(), ((bb[m]-0.5)*L[m]/(2*np.pi)).mean(), int((bb[m] > 1).sum())))
json.dump(res, open("t2_offline.json", "w"))
