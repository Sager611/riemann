# ceilaudit B1: freedom outside Bernard's input = jointly complex (P, Q) (equivalently imaginary mollifier shifts with complex Q).
# Independent Ritz test of the GL(2) one-piece functional at xi = nu/2 = 25/256 (nu -> 25/128):
# c = 1 + (1/xi) int int e^{2Rv} |P'(u)Q(v) + xi P(u)(Q'(v) + R Q(v))|^2, P(0)=0, P(1)=1, Q(0)=1, Q(y)+Q(1-y) = const (complex allowed).
import numpy as np, math
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize
xg, wg = leggauss(80); x = (xg + 1) / 2; w = wg / 2
XI = 25 / 256
def PQ(par, m, K, cplx):
    if cplx:
        z = par[:m + K] + 1j * par[m + K:]
    else:
        z = par[:m + K].astype(complex)
    p, cq = z[:m], z[m:]
    P = x + x * (1 - x) * sum(p[k] * x ** k for k in range(m))
    dP = 1 + sum(p[k] * ((k + 1) * x ** k - (k + 2) * x ** (k + 1)) for k in range(m))
    beta2 = 1 - cq.sum()   # Q(0) = beta/2 + sum c_k = 1
    Q = beta2 + sum(cq[k - 1] * (1 - 2 * x) ** (2 * k - 1) for k in range(1, K + 1))
    dQ = sum(cq[k - 1] * (2 * k - 1) * (1 - 2 * x) ** (2 * k - 2) * (-2) for k in range(1, K + 1))
    return P, dP, Q, dQ
def cfun(par, m, K, R, cplx):
    P, dP, Q, dQ = PQ(par, m, K, cplx)
    e = np.exp(2 * R * x)
    A = np.sum(w * abs(dP) ** 2); B = np.sum(w * e * abs(Q) ** 2)
    C = np.sum(w * abs(P) ** 2); D = np.sum(w * e * abs(dQ + R * Q) ** 2)
    X = np.sum(w * dP * np.conj(P)); Y = np.sum(w * e * Q * np.conj(dQ + R * Q))
    return 1 + (A * B + XI ** 2 * C * D + 2 * XI * (X * Y).real) / XI
rng = np.random.default_rng(1)
for R in [9.5, 9.867773, 10.3]:
    for (m, K) in [(4, 3), (6, 5)]:
        n = m + K
        best_r = None
        for s in range(6):
            x0 = rng.normal(0, 0.3, n) if s else np.zeros(n)
            r = minimize(lambda v: math.log(cfun(v, m, K, R, False)), x0, method='BFGS', options={'gtol': 1e-11, 'maxiter': 4000})
            if best_r is None or r.fun < best_r.fun: best_r = r
        kr = 1 - best_r.fun / R
        best_c = None
        for s in range(12):
            x0 = np.concatenate([best_r.x, np.zeros(n)]) + (rng.normal(0, 0.3, 2 * n) if s else 0)
            r = minimize(lambda v: math.log(cfun(v, m, K, R, True)), x0, method='BFGS', options={'gtol': 1e-11, 'maxiter': 6000})
            if best_c is None or r.fun < best_c.fun: best_c = r
        kc = 1 - best_c.fun / R
        imn = np.abs(best_c.x[n:]).max()
        print(f"R={R} deg(P)~{m+1} K={K}: real kappa={kr:.10f}  complex kappa={kc:.10f}  diff={kc-kr:.2e}  max|Im coeff| at complex optimum={imn:.1e}", flush=True)
print("certified ceiling 0.07005633 (float 0.0700563167 at R = 9.867773)")
