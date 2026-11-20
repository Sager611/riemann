# falsify.py -- numerical falsification of the key lemmas behind the record (own code; flint acb zeta at 64 bits, numpy).
# (R1) |1/zeta(s1+it)| <= Phi(t) (M78 Thm A3 with B'), all t in [0, 1e4]; (R2) <= C_inf(1000) t^a (Cor. A4), t in [1000, 1e4];
# (R3) |1/zeta(0.85+it)| <= C_L t^{a_L} (MCH Lemma H1), t in [1000, 1e4]; (R4) J_L0 by trapezoid; (R5) Lemma A2 itself with true zeta values.
import json, math, time, sys
import numpy as np
from flint import acb, ctx
import mpmath
ctx.prec = 64
M78 = '/Users/adriansager/playground/riemann/Riemann/Deep/Mertens78/'
BP = 0.87334585610671858; D = 0.875 - BP
zr = lambda x: float(mpmath.zeta(x))
def cert(fn, eta):
    d = json.load(open(M78 + fn)); s1 = float(d['sigma1']) - D
    ms = np.array([[float(a) - D, float(b), float(c)] for a, b, c in d['masses']])
    kap = (1 + eta - BP)/2; A = math.log(3*zr(1 + eta)) - kap*math.log(1 + BP)
    Lam = ms[:, 2].sum(); a = (Lam - 1)*kap
    z1 = np.array([zr(x) for x in ms[:, 0]]); z2 = np.array([zr(2*x) for x in ms[:, 0]])
    R = s1; K = (Lam - 1)*A
    for (sg, y, lam), q1, q2 in zip(ms, z1, z2):
        c = abs(y)
        if lam > 0: R += lam*((1 + kap)*(c + 1 + sg) + 2*c + 2); K += lam*math.log(q1/q2)
        else: R += -lam*((3 + 2*kap)*c + sg); K += -lam*math.log(q1)
    return dict(s1=s1, ms=ms, kap=kap, A=A, a=a, lC=K + R/1000, z1=z1, z2=z2)
def logPhi(c, t):
    s1, ms, kap, A = c['s1'], c['ms'], c['kap'], c['A']
    h = lambda sg, tt: A + (1 + kap)*0.5*np.log((1 + sg)**2 + tt**2) - 0.5*np.log(sg**2 + tt**2)
    tot = -h(s1, t) + 0.5*np.log((s1 - 1)**2 + t**2) - 0.5*np.log(s1**2 + t**2)
    for (sg, y, lam), q1, q2 in zip(ms, c['z1'], c['z2']):
        tt = t + y; lm = 0.5*np.log(sg**2 + tt**2) - 0.5*np.log((sg - 1)**2 + tt**2)
        if lam > 0: tot += lam*(h(sg, tt) + math.log(q1/q2) + lm)
        else: tot += -lam*(-h(sg, tt) + math.log(q1) - lm)
    return tot
def absz(s, t): return abs(complex(acb(s, t).zeta()))
CS = {'ML092g8_e05_Bp': cert('c_ml092_g8.json', 0.05), 'ML093g5_e10_Bp': cert('c_ml093_g5.json', 0.1)}
kp = 0.3; Ap = math.log(3*zr(1.1)) + kp*math.log(2/3); LZ = math.log(zr(1.02)/zr(2.04)); LamL = (2.04 - 0.85 - 0.5)/0.35
aL = (LamL - 1)*kp; lCL = (LamL - 1)*Ap + LamL*LZ + (0.85 + LamL*(1.3*2.02 + 2))/1000
print('C_inf(1000) logs', {k: round(v['lC'], 6) for k, v in CS.items()}, 'a', {k: round(v['a'], 6) for k, v in CS.items()}, 'logC_L', round(lCL, 6), 'a_L', round(aL, 6), flush=True)
t0 = time.time()
def refine(f, t, w=0.05):
    g = (math.sqrt(5) - 1)/2; lo, hi = t - w, t + w
    for _ in range(22):
        c1 = hi - g*(hi - lo); d1 = lo + g*(hi - lo)
        if f(c1) > f(d1): hi = d1
        else: lo = c1
    tm = (lo + hi)/2; return tm, f(tm)
def scan(sig, step, tmax, ratio_fns):
    t = np.arange(0.0, tmax + step/2, step); iz = np.array([1/absz(sig, x) for x in t])
    out = {}
    for nm, (fn, tmin, pointwise) in ratio_fns.items():
        m = t >= tmin; r = iz[m]/fn(t[m]); tt = t[m]
        idx = [i for i in range(1, len(r) - 1) if r[i] >= r[i-1] and r[i] >= r[i+1]]
        top = sorted(idx, key=lambda i: -r[i])[:12]
        best = (0, None)
        for i in top:
            tm, v = refine(lambda x: (1/absz(sig, x))/pointwise(x), tt[i])
            if v > best[0]: best = (v, tm)
        out[nm] = dict(grid_max=float(r.max()), grid_tmax=float(tt[r.argmax()]), refined_max=best[0], refined_t=best[1], n=int(m.sum()))
        print(nm, 'sigma', round(sig, 6), out[nm], round(time.time() - t0), 's', flush=True)
    return t, iz
for nm, c in CS.items():
    step = 0.05 if nm.startswith('ML092') else 0.1
    fns = {nm + ' R1 |1/z|/Phi (t>=0)': (lambda tt, c=c: np.exp(logPhi(c, tt)), 0.0, lambda x, c=c: float(np.exp(logPhi(c, np.array([x])))[0])),
           nm + ' R2 |1/z|/(C t^a) (t>=1000)': (lambda tt, c=c: np.exp(c['lC'])*tt**c['a'], 1000.0, lambda x, c=c: math.exp(c['lC'])*x**c['a'])}
    scan(c['s1'], step, 10000.0, fns)
fnsL = {'R3 |1/z(0.85+it)|/(C_L t^a_L) (t>=1000)': (lambda tt: np.exp(lCL)*tt**aL, 1000.0, lambda x: math.exp(lCL)*x**aL)}
t, iz = scan(0.85, 0.05, 10000.0, fnsL)
m = t <= 1000.0; f = iz[m]/np.sqrt(0.85**2 + t[m]**2); JL = 2*np.sum((f[1:] + f[:-1])/2*np.diff(t[m]))
print('R4 J_L0 trapezoid (step 0.05)', JL, 'vs certified upper 19.5963', 'max|1/z| on [0,1000]', iz[m].max(), flush=True)
# (R5) Lemma A2 directly: v(s1) <= sum lam_j v(s_j), v = h - log|Z|, true zeta values, 160 t incl. the first 40 zero ordinates
c = CS['ML092g8_e05_Bp']; s1, ms, kap, A = c['s1'], c['ms'], c['kap'], c['A']
def v(sg, tt):
    s = complex(sg, tt); z = complex(acb(sg, tt).zeta())
    return A + (1 + kap)*math.log(abs(1 + s)) - math.log(abs(s)) - math.log(abs((s - 1)*z/s))
ts = [float(mpmath.zetazero(n).imag) for n in range(1, 41)] + list(np.linspace(0.3, 9990, 120))
worst = -1e9
for tt in ts:
    lhs = v(s1, tt); rhs = sum(lam*v(sg, tt + y) for sg, y, lam in ms)
    worst = max(worst, lhs - rhs)
print('R5 max over', len(ts), 't of v(s1) - sum lam_j v(s_j) (must be <= 0):', worst, round(time.time() - t0), 's', flush=True)
