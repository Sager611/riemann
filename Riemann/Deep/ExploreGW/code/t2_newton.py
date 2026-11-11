
# T2 Newton: locate off-line zeros of f_c(s) = cos c * k12 L(s,chi12) + sin c * (k3 L(s,chi3) + k9 L(s,chi9)) (mod 13, even)
import sys, json, numpy as np
from flint import acb, arb, ctx, dirichlet_char
ctx.prec = 96
PI = acb(arb.pi())
ch = {n: dirichlet_char(13, n) for n in (12, 3, 9)}
def gam(s): return (acb(13)/PI)**(s/2) * (s/2).gamma()
def kappa(n, t0=10.0):
    s0 = acb(0.5, t0); e0 = gam(s0); e0 = e0/abs(e0)
    return acb(ch[n].hardy_z(acb(t0)).real) / (e0 * ch[n].l(s0))
def mid(z): return acb(z.real.mid(), z.imag.mid())
K = {n: mid(kappa(n)) for n in ch}
for n in ch:  # constancy check
    print("kappa_%d = %s ; drift at t0=60: %s" % (n, K[n].str(8), abs(K[n]-kappa(n, 60.0)).mid().str(3)))
def f(c, s):
    cc, sc = acb(float(np.cos(c))), acb(float(np.sin(c)))
    return cc*K[12]*ch[12].l(s) + sc*(K[3]*ch[3].l(s) + K[9]*ch[9].l(s))
_t = f(0.1, acb(0.58, 277.4)); print("probe f:", _t.str(6))
h = acb(1e-12)
def newton(c, s, it=40):
    for _ in range(it):
        fs = mid(f(c, s)); d = mid((mid(f(c, s+h)) - mid(f(c, s-h)))/(2*h))
        step = fs/d
        if abs(step).mid() > 0.3: step = step*(0.3/float(abs(step).mid()))
        s = mid(s - step)
        if float(abs(step).mid()) < 1e-14: break
    return s, float(abs(f(c, s)).mid())
crit = json.load(open("t2_crit.json"))
negmin = crit["negmin"]
out = {}
for c in [float(x) for x in sys.argv[1:]]:
    tc = np.tan(c)
    try:
        ws = json.load(open("t2_ws_%.3f.json" % c))
    except Exception:
        ws = []
    cands = [(t, -1.0) for t in ws[:5]]
    found = []
    for (t, a) in cands:
        s, res = newton(c, acb(0.5 + 0.08, t))
        b = float(s.real.mid()) - 0.5
        found.append((round(t, 3), round(a, 5), round(float(s.real.mid()), 6), round(float(s.imag.mid()), 6), "%.1e" % res))
    out[c] = found
    print("c = %.3f (tan %.4f): EP births below threshold (first 6) -> Newton zeros (t*, |r*|, Re s, Im s, |f|):" % (c, tc))
    for x in found: print("   ", x)
json.dump({str(k): v for k, v in out.items()}, open("t2_newton.json", "w"))
