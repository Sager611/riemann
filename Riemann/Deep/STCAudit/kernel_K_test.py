
# Falsification test of Casana-Giner's kernel inequality (K): U(b,t) - f(b,t) >= 0 on 0<=b<=1/2, t != 0
# for the six line kernels (exact binary64 coefficients from Appendix D). float64 + scipy local minimisation.
import json, math, sys
import numpy as np
from scipy.optimize import minimize
K = json.load(open('Riemann/Deep/STCAudit/kernels.json'))
def Hfun(k, b, t, scaled=False):
    d = k['d']; B = k['B']; c = k['c']
    b = np.asarray(b, float); t = np.asarray(t, float)
    sg = np.sign(t)
    # f in the form (4.17) for |t|<=1, series-stable form for |t|>1
    small = np.abs(t) <= 1.0
    f = np.empty(np.broadcast(b, t).shape)
    bb, tt = np.broadcast_arrays(b, t)
    ys = [(2.0, d+bb), (2.0, d-bb), (-1.0, 2*d+bb), (-1.0, 2*d-bb)]
    fs = sg*math.pi + sum(-a*np.arctan(tt/y) for a, y in ys)  # (4.17)
    u = 1.0/tt
    def atm(x):  # arctan(x) - x, accurate for small x
        out = np.arctan(x) - x
        s = np.abs(x) < 0.05
        xs = x[s]; x2 = xs*xs
        out[s] = xs*x2*(-1/3 + x2*(1/5 + x2*(-1/7 + x2*(1/9 - x2/11))))
        return out
    fl = sum(a*atm(y*u) for a, y in ys)   # exact since sum a*y = 0
    f = np.where(small, fs, fl)
    U = np.zeros_like(f)
    for (X, j), cc in zip(B, c):
        for x in (d+bb, d-bb):
            z = (x + 1j*tt)**(-j)
            U += cc*(z.real if X == 'R' else -z.imag)
    H = U - f
    if scaled:
        return H*np.maximum(1.0, tt*tt)
    return H
res = {}
b = np.linspace(0, 0.5, 201)
tpos = np.unique(np.concatenate([np.logspace(-9, 7, 4000), np.linspace(1e-9, 30, 6000)]))
tg = np.concatenate([-tpos[::-1], tpos])
for name, k in K.items():
    cand = []
    gmin = (np.inf, 0, 0)
    for ch in np.array_split(tg, 40):
        BB, TT = np.meshgrid(b, ch, indexing='ij')
        Hs = Hfun(k, BB, TT, scaled=True)
        ii = np.argsort(Hs, axis=None)[:40]
        for q in ii: cand.append((Hs.flat[q], BB.flat[q], TT.flat[q]))
        del BB, TT, Hs
    cand.sort(); cand = cand[:40]
    gmin = cand[0]
    best = []
    for hv, b0, t0 in cand:
        x0 = [b0, t0]
        sgn = np.sign(x0[1])
        # minimise scaled H over b in [0,0.5], t of fixed sign
        lo, hi = (1e-12, 1e8) if sgn > 0 else (-1e8, -1e-12)
        fun = lambda x: float(Hfun(k, np.array([x[0]]), np.array([x[1]]), scaled=True)[0])
        r = minimize(fun, x0, method='L-BFGS-B', bounds=[(0, 0.5), (lo, hi)])
        best.append((r.fun, r.x[0], r.x[1]))
    best.sort()
    d = k['d']; cR1 = k['c'][k['B'].index(['R',1])]; cR2 = k['c'][k['B'].index(['R',2])]
    # t -> 0 limits: U(b,0) - pi and U(b,0)+pi (t->0- gives f->-pi)
    U0 = Hfun(k, b, np.full_like(b, 1e-13)) # = U - pi approx
    print(name, 'grid min scaled H = %.3e at b=%.4f t=%.4g' % gmin,
          '| local min %.3e at b=%.6f t=%.6g' % best[0], '| t^2H(u=0)=2dcR1-2cR2=%.6f' % (2*d*cR1-2*cR2),
          '| min_b H(b,0+)=%.4f' % U0.min())
    # distinct local minima
    seen = []
    for v in best:
        if all(abs(v[1]-s[1]) + abs(v[2]-s[2]) > 1e-3 for s in seen): seen.append(v)
    print('   distinct local minima:', ['%.2e@(%.4f,%.4g)' % s for s in seen[:8]])
    sys.stdout.flush()
