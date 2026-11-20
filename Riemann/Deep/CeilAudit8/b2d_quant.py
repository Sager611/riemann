# ceilaudit8 B2d: does the quantizer itself (2-point Lloyd-Max uniqueness, Fleischer's conclusion) on p = |f^(m)|^alpha
# see a zeta pair moved to delta >= 1.25 gbar (invisible to log-concavity of |f| = L_1(f))?  f^(m) via b2c's stable Taylor form.
import numpy as np, math, sys, importlib.util
spec = importlib.util.spec_from_file_location('b2c', 'Riemann/Deep/CeilAudit8/b2c_turan_stable.py')
src = open('Riemann/Deep/CeilAudit8/b2c_turan_stable.py').read().split("run(6708")[0]
src = src.replace("range(0, far_real.size, 1000)", "range(0, far_real.size, 250)").replace("far_real[None, j:j+1000]", "far_real[None, j:j+250]")
b2c = {}; exec(src, b2c)   # sized: 6001 x 250 doubles = 12 MB per temporary
Z = b2c['Z']; M = b2c['M']
def fm(x, real, pair, gbar, x0, m):
    near = real[np.abs(real-x0) <= b2c['RN']]; far = real[np.abs(real-x0) > b2c['RN']]
    T = np.zeros((M+1, x.size)); T[0] = 1.0
    for r in near:
        a = (x-r)/gbar; T[1:] = a*T[1:] + T[:-1]; T[0] = a*T[0]
    g0, d = pair; w = (x-g0)/gbar; dd = d/gbar
    P = np.zeros((3, x.size)); P[0] = w*w+dd*dd; P[1] = 2*w; P[2] = 1.0; T = b2c['conv'](T, P)
    G = b2c['far_D'](x, far, [-g0+1j*d, -g0-1j*d], gbar); F = b2c['conv'](T, G)
    D1 = G[1]/gbar; lg = np.concatenate([[0], np.cumsum((D1[1:]+D1[:-1])/2*np.diff(x))])  # log g_far(x)/g_far(x[0])
    val = math.factorial(m)*F[m]*np.exp(lg - lg.max())
    a = math.factorial(m+1)**2*F[m+1]**2; b = math.factorial(m)*math.factorial(m+2)*F[m]*F[m+2]
    return val, (a-b)/(a+np.abs(b))
def nfix(x, p):  # number of 2-point Lloyd-Max stationary quantizers of density p on [x0, x-1]
    w = p*np.gradient(x); W = np.cumsum(w); X = np.cumsum(w*x)
    s = slice(5, x.size-5); c1 = X[s]/W[s]; c2 = (X[-1]-X[s])/(W[-1]-W[s])
    phi = (c1+c2)/2 - x[s]; return int(np.sum(np.sign(phi[1:]) != np.sign(phi[:-1])))
def test(ia, ib, cases, tag):
    ga, gb = Z[ia], Z[ib]; g0 = (ga+gb)/2; gbar = 2*math.pi/math.log(g0/(2*math.pi))
    keep = np.ones(Z.size, bool); keep[[ia, ib]] = False; zz = Z[keep]; real = np.concatenate([zz, -zz])
    x = np.linspace(g0-6, g0+6, 6001)+1.3e-6
    for (dg, m) in cases:
        v, q = fm(x, real, (g0, dg*gbar), gbar, g0, m)
        sc = np.where(np.sign(v[1:]) != np.sign(v[:-1]))[0]   # zeros of f^(m)
        j = np.argmin(q); lo = sc[sc < j].max() if np.any(sc < j) else 0; hi = sc[sc >= j].min()+1 if np.any(sc >= j) else x.size-1
        res = []
        for al in (1, 2):
            p = np.abs(v[lo:hi+1])**al; xx = x[lo:hi+1]; whole = nfix(xx, p)
            idx = np.linspace(0, xx.size-1, 13).astype(int); nu = 0; tot = 0
            for i1 in range(len(idx)):
                for i2 in range(i1+2, len(idx)):
                    if idx[i2]-idx[i1] < 40: continue
                    tot += 1; nu += nfix(xx[idx[i1]:idx[i2]+1], p[idx[i1]:idx[i2]+1]) > 1
            res.append(f"alpha={al}: whole {whole} fixed pts, non-unique on {nu}/{tot} sub-intervals")
        print(f"{tag} delta/gbar={dg} m={m}: dip interval [{x[lo]-g0:+.3f},{x[hi]-g0:+.3f}] (x-g0), min q={q.min():+.3f}; " + "; ".join(res)); sys.stdout.flush()
test(6708, 6709, [(1.0, 0), (1.25, 0), (1.25, 2), (1.5, 2), (1.5, 3), (2.0, 4), (2.5, 6), (3.0, 7)], 'Lehmer')
test(6606, 6607, [(1.25, 1), (1.25, 2), (2.0, 4), (3.0, 7)], 'generic')
