# gl2audit: global float64 search, jointly complex (P,Q) vs real, same degree; F_c = inf over complex P (B l.796 / AT l.996)
import sys, numpy as np
from numpy.polynomial import chebyshev as Ch
from scipy.optimize import minimize
xi = float(sys.argv[1])/float(sys.argv[2]); K = int(sys.argv[3]); KE = int(sys.argv[4]); seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
xs, ws = np.polynomial.legendre.leggauss(80); Y = (xs+1)/2; Wt = ws/2; T = 1-2*Y
def bas(n):
    c = np.zeros(n+1); c[n] = 1; return Ch.chebval(T, c)-1, -2*Ch.chebval(T, Ch.chebder(c))
PH = [bas(2*k-1) for k in range(1, K+1)]; CHI = [bas(2*k) for k in range(1, KE+1)]
def Fc(R, b, d):
    qr = 1+sum(bk*p[0] for bk, p in zip(b, PH)); dqr = sum(bk*p[1] for bk, p in zip(b, PH))
    qi = sum(dk*p[0] for dk, p in zip(d, CHI)); dqi = sum(dk*p[1] for dk, p in zip(d, CHI))
    e = Wt*np.exp(2*R*Y)
    A = e@(qr*qr+qi*qi); C = e@((R*qr+dqr)**2+(R*qi+dqi)**2); I = e@(qi*dqr-qr*dqi)
    q1 = 1-2*sum(b); X = np.sqrt(max(A*C-I*I, 1e-300))
    return (1+np.exp(2*R)*q1*q1)/2+X/np.tanh(xi*X/A)
def kap(v, R, cplx):
    b = v[:K]; d = v[K:] if cplx else np.zeros(KE)
    return -(1-np.log(Fc(R, b, d))/R)
rng = np.random.default_rng(seed)
for R in [float(r) for r in sys.argv[6].split(',')]:
    best_r = None
    for s in range(4):
        v0 = np.zeros(K) if s == 0 else best_r[1]+1e-3*rng.standard_normal(K)
        o = minimize(kap, v0, args=(R, False), method='BFGS', options={'gtol': 1e-12, 'maxiter': 5000})
        if best_r is None or o.fun < best_r[0]: best_r = (o.fun, o.x)
    best_c = (np.inf, None)
    for s in range(12):
        sc = [1e-3, 1e-2, 1e-1, 1.0][s % 4]
        v0 = np.concatenate([best_r[1], sc*rng.standard_normal(KE)])
        o = minimize(kap, v0, args=(R, True), method='BFGS', options={'gtol': 1e-12, 'maxiter': 5000})
        if o.fun < best_c[0]: best_c = (o.fun, o.x)
    print(f'R={R:.4f} K={K} KE={KE}: real {-best_r[0]:.12f}  complex {-best_c[0]:.12f}  diff {best_r[0]-best_c[0]:+.2e}  |d|max {np.abs(best_c[1][K:]).max():.1e}', flush=True)
