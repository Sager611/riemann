# e1rate: optimal Levinson-Conrey detector kappa(nu) in SL's class (E1 machinery), and h(theta;kappa).
# c(P,Q,R,nu) = 1/2 + (CP/nu) W2 + nu BP W2',  W2 = int e^{2Rv}Q^2, W2' = int e^{2Rv}(RQ+Q')^2  (Q(0)=1, Q(y)+Q(1-y)=1 => Q(1)=0)
# min over P (P(0)=0,P(1)=1) of a*CP + b*BP = sqrt(ab) coth(sqrt(b/a)) (P = sinh(lx)/sinh l; polynomials approach it).
import numpy as np, sys, json
from numpy.polynomial import chebyshev as C
from scipy.optimize import minimize
xg, wg = np.polynomial.legendre.leggauss(160); xg = (xg+1)/2; wg = wg/2
def Qvals(xs, y):
    u = 1-2*y; q = (1+u)/2 + 0*y; dq = -1+0*y   # 1-y
    for j, x in enumerate(xs, start=1):
        c = np.zeros(2*j+2); c[2*j+1] = 1
        Tj = C.chebval(u, c); dTj = C.chebval(u, C.chebder(c))
        q = q + x*(Tj-u); dq = dq + x*(dTj-1)*(-2)
    return q, dq
def cval(par, nu):
    R = par[0]; xs = par[1:]
    q, dq = Qvals(xs, xg); e = np.exp(2*R*xg)
    W2 = np.sum(wg*e*q*q); W2p = np.sum(wg*e*(R*q+dq)**2)
    a = W2/nu; b = nu*W2p; l = np.sqrt(b/a)
    return 0.5 + np.sqrt(a*b)/np.tanh(l)
def kappa(par, nu): return 1-np.log(cval(par, nu))/par[0]
def opt(nu, deg=6, R0=None):
    best = None
    for R0_ in ([R0] if R0 else [1.0, 1.5, 2.5, 4.0, 8.0]):
        x0 = np.concatenate([[R0_], np.zeros(deg)])
        f = lambda p: -kappa(p, nu) if p[0] > 0.05 else 10.0
        r = minimize(f, x0, method='Nelder-Mead', options={'maxiter': 40000, 'maxfev': 40000, 'xatol': 1e-10, 'fatol': 1e-13})
        r = minimize(f, r.x, method='BFGS', options={'gtol': 1e-11})
        if best is None or r.fun < best.fun: best = r
    return -best.fun, best.x
s2 = np.sqrt(2)
def cth(th): x = th/s2; return 2-th/2-(np.cos(x)/np.sin(x))/s2
def h(th, k): return (3+k-np.sqrt((1-k)*(9-k-8*cth(th))))/4
def lb(th, k): c = cth(th); return max(0.0, c, (c+2*k)/3, h(th, k))
if __name__ == "__main__":
    out = {}
    for nu in [float(v) for v in sys.argv[1:]]:
        k1, p1 = opt(nu, deg=4); k2, p2 = opt(nu, deg=8, R0=p1[0])
        out[nu] = (k2, list(p2))
        print(f"nu={nu:.6f} kappa(deg4)={k1:.9f} kappa(deg8)={k2:.9f} R={p2[0]:.5f} x={np.round(p2[1:],6).tolist()}", flush=True)
    json.dump(out, open('kappa_opt.json', 'w'))
