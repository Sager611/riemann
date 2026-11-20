# c1check.py -- code audit of the (C1)/(C2) certificates used by the record (ML092g8, ML093g5; heights = B'-heights exactly).
# Part A (their procedure, read-only import of Mertens78/cert.py): do the float-midpoint balls arb((lo+hi)/2,(hi-lo)/2) cover [lo,hi]?
#        Same question for J0_direct's t-balls at sigma_L = 0.85 (the J_L0 of the record).
# Part B (independent, mpmath.iv only): exact-endpoint cover of [-Y, Y] by bisection + the tail bracket, outward rounding.
import json, sys, time, math
t0 = time.time()
M78 = '/Users/adriansager/playground/riemann/Riemann/Deep/Mertens78'
sys.path.insert(0, M78)
from flint import arb, acb, ctx
import cert as C78
ctx.prec = 120
def gaps_check_C(c):
    # replicate check_C's piece generation and acceptance, counting balls that miss an endpoint of their piece
    x1 = c.x1; cmax = max(abs(float(y.mid())) for _, y, _ in c.m)
    Y = None
    def bracket(Yv):
        Ya = arb(Yv); br = -x1
        for sig, y, lam in c.m:
            x = sig - C78.B; cc = abs(y)
            br += lam*x*Ya**2/((Ya + cc)**2 + x**2) if lam > 0 else lam*x*Ya**2/((Ya - cc)**2)
        return br
    for Yc in [2*cmax + 50, 1e3, 1e4, 1e5, 1e6, 1e7]:
        if Yc > cmax and bracket(Yc) > 0: Y = Yc; break
    Y0 = min(Y, 2*cmax + 50)
    stack = [(-Y0 + 2*Y0*k/4000, -Y0 + 2*Y0*(k+1)/4000, 0) for k in range(4000)]
    w = Y0
    while w < Y:
        w2 = min(w*1.002, Y); stack.append((w, w2, 0)); stack.append((-w2, -w, 0)); w = w2
    last_hi = -Y0 + 2*Y0*4000/4000
    n = miss = 0; maxgap = 0.0; minlow = None
    while stack:
        lo, hi, dep = stack.pop(); n += 1
        iv = arb((lo+hi)/2, (hi-lo)/2)
        ok = iv.contains(arb(lo)) and iv.contains(arb(hi))
        if not ok:
            miss += 1; g = max(float((arb(lo) - (iv.mid() - iv.rad())).upper()), float(((arb(hi) - (iv.mid() + iv.rad()))).upper()), 0.0); maxgap = max(maxgap, g)
        v = c.F(iv)
        if v > 0:
            lw = float(v.lower()); minlow = lw if minlow is None else min(minlow, lw); continue
        mid = (lo+hi)/2; stack.append((lo, mid, dep+1)); stack.append((mid, hi, dep+1))
    return dict(Y=Y, Y0=Y0, uniform_last_hi_equals_Y0=(last_hi == Y0), balls=n, balls_missing_an_endpoint=miss, max_gap=maxgap, min_lower_F=minlow)
for nm, fn, eta in [('ML092g8', 'c_ml092_g8.json', '0.05'), ('ML093g5', 'c_ml093_g5.json', '0.1')]:
    d = json.load(open(M78 + '/' + fn)); d['eta'] = eta; d['T1'] = '1000'
    c = C78.Cert(d)
    print('A', nm, gaps_check_C(c), round(time.time() - t0, 1), 's', flush=True)
# J0_direct ball cover at sigma_L = 0.85 (count only; no zeta calls)
t = 0.0; miss = 0; n = 0; maxgap = 0.0
while t < 1000.0:
    w = 0.02 if t < 100 else 0.01; t2 = min(t + w, 1000.0)
    iv = arb((t+t2)/2, (t2-t)/2); n += 1
    if not (iv.contains(arb(t)) and iv.contains(arb(t2))):
        miss += 1; maxgap = max(maxgap, float((arb(t) - (iv.mid() - iv.rad())).upper()), float((arb(t2) - (iv.mid() + iv.rad())).upper()))
    t = t2
print('A J0_direct balls', n, 'missing an endpoint', miss, 'max gap', maxgap, flush=True)
# ---------------- Part B: independent mpmath.iv check ----------------
from mpmath import iv
iv.prec = 100
assert (iv.mpf([-1, 2])**2).a == 0 and (iv.mpf([-1, 2])**2).b == 4
def run_iv(fn):
    d = json.load(open(M78 + '/' + fn))
    E = iv.mpf('0.875')
    ms = [(iv.mpf(a) - E, iv.mpf(b), iv.mpf(lam), float(b), float(lam)) for a, b, lam in d['masses']]
    x1 = iv.mpf(d['sigma1']) - E
    c2 = sum((m[2]*m[0] for m in ms), iv.mpf(0)) - x1
    assert c2.a > 0
    cmax = max(abs(m[3]) for m in ms)
    def F(I):
        tot = -x1/(I**2 + x1**2)
        for x, y, lam, _, _ in ms: tot += lam*x/((I - y)**2 + x**2)
        return tot
    def bracket(Yv):
        Y = iv.mpf(Yv); br = -x1
        for x, y, lam, yf, lf in ms:
            cc = abs(yf)
            br += lam*x*Y**2/((Y + cc)**2 + x**2) if lf > 0 else lam*x*Y**2/((Y - cc)**2)
        return br
    Y = next(Yc for Yc in [2*cmax + 10, 2*cmax + 50, 100.0, 300.0, 1e3, 3e3, 1e4, 1e5, 1e6, 1e7] if Yc > cmax and bracket(Yc).a > 0)
    # exact-endpoint pieces: 3000 uniform on [-min(Y,60), min(Y,60)], geometric ratio 1.01 beyond, all endpoints exact binary floats
    Y0 = min(Y, 60.0); pts = [-Y0 + 2*Y0*k/3000 for k in range(3001)]; pts[0] = -Y0; pts[-1] = Y0
    right = [Y0]
    while right[-1] < Y: right.append(min(right[-1]*1.01, Y))
    pieces = [(pts[k], pts[k+1]) for k in range(3000)] + [(right[k], right[k+1]) for k in range(len(right)-1)] + [(-right[k+1], -right[k]) for k in range(len(right)-1)]
    # cover check: sorted endpoints chain from -Y to Y with no gap
    sp = sorted(pieces); assert sp[0][0] == -Y and sp[-1][1] == Y and all(sp[k][1] == sp[k+1][0] for k in range(len(sp)-1))
    n = 0; minlow = None; stack = [(a, b, 0) for a, b in pieces]
    while stack:
        a, b, dep = stack.pop(); n += 1
        v = F(iv.mpf([a, b]))
        if v.a > 0:
            minlow = float(v.a) if minlow is None else min(minlow, float(v.a)); continue
        assert dep < 45, ("C1 fails near", a, b, float(v.a))
        m = (a + b)/2; stack.append((a, m, dep+1)); stack.append((m, b, dep+1))
    # float minimum of F on a fine grid (non-rigorous, shows the tightness)
    import numpy as np
    w = np.linspace(-Y, Y, 400001); Ff = -float(x1.mid)/(w**2 + float(x1.mid)**2)
    for x, y, lam, yf, lf in ms: Ff += lf*float(x.mid)/((w - yf)**2 + float(x.mid)**2)
    return dict(Y=Y, C2_lower=float(c2.a), intervals=n, min_lower_F=minlow, float_min_F=float(Ff.min()), at_w=float(w[Ff.argmin()]))
for nm, fn in [('ML092g8', 'c_ml092_g8.json'), ('ML093g5', 'c_ml093_g5.json')]:
    print('B', nm, run_iv(fn), round(time.time() - t0, 1), 's', flush=True)
