# Gram-point repulsion test: 2406.06548 conj:rep / 2310.14415 Conj (G-B-G):
# g_n bad, g_{n-1}, g_{n+1} good  ==>  |Z'(g_n)| > 4 |Z(g_n)|.
# At a Gram point zeta(1/2+i g_n) is real, Z(g_n) = (-1)^n zeta, |Z'(g_n)| = |Im zeta'(1/2+i g_n)|; bad <=> zeta(1/2+i g_n) < 0.
import sys, time, math, json
import numpy as np
from flint import acb, arb, ctx, acb_series
ctx.prec = 64
n0, n1 = int(sys.argv[1]), int(sys.argv[2]); out = sys.argv[3]; THR = float(sys.argv[4]) if len(sys.argv) > 4 else 6.0
def theta_np(t):
    return t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3) + 31/(80640*t**5) + 127/(430080*t**7)
def dtheta_np(t):
    return 0.5*np.log(t/(2*np.pi)) - 1/(48*t**2) - 7/(1920*t**4) - 31/(16128*t**6)
def gram(ns):
    ns = np.asarray(ns, dtype=np.float64)
    x = (ns + 0.125)/np.e
    w = np.log(x + 1.0)  # Lambert W start
    for _ in range(40): w = w - (w*np.exp(w) - x)/(np.exp(w)*(w + 1))
    t = 2*np.pi*np.e*np.exp(w)
    t = np.maximum(t, 9.0)
    for _ in range(8): t = t - (theta_np(t) - np.pi*ns)/dtheta_np(t)
    return t
def zeta_val(t):
    z = acb(0.5, t).zeta()
    return float(z.real.mid())
def zeta_der(t):
    s = acb_series([acb(0.5, t), acb(1)], prec=2)
    f = s.zeta()
    c = f.coeffs()
    return c[0], c[1]
CH = 20000
lo = max(n0 - 1, -1)
prev2 = None
nb = 0; niso = 0; minm = (1e9, None); cands = []; t_start = time.time()
f = open(out, 'w')
vals = {}
n = lo
while n < n1 + 1:
    ns = np.arange(n, min(n + CH, n1 + 2))
    ts = gram(ns)
    for k, t in zip(ns.tolist(), ts.tolist()):
        vals[k] = (t, zeta_val(t))
    # process isolated bad in this chunk where neighbours known
    for k in ns.tolist():
        if k - 1 not in vals or k + 1 not in vals or k < n0 or k >= n1: continue
        t, z = vals[k]
        if z < 0:
            nb += 1
            if vals[k-1][1] > 0 and vals[k+1][1] > 0:
                niso += 1
                c0, c1 = zeta_der(t)
                m = abs(float(c1.imag.mid())) / abs(float(c0.real.mid()))
                if m < minm[0]: minm = (m, k, t)
                if m < THR: cands.append((k, t, m)); f.write(f'CAND n={k} t={t:.10f} zeta={float(c0.real.mid()):.6e} Imzeta1={float(c1.imag.mid()):.6e} m={m:.6f} m_norm={m/(0.5*math.log(t/(2*math.pi))):.5f}\n'); f.flush()
    # drop old
    for k in list(vals):
        if k < ns[-1] - 2: del vals[k]
    n = int(ns[-1])  # overlap one for neighbour
    f.write(f'progress n={n} t={ts[-1]:.1f} bad={nb} isolated_bad={niso} min_m={minm[0]:.5f} at n={minm[1]} elapsed={time.time()-t_start:.0f}s\n'); f.flush()
    if n >= n1: break
f.write(f'DONE n in [{n0},{n1}) bad={nb} isolated_bad={niso} min_m={minm} cands(m<THR)={len(cands)}\n')
