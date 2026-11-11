# usage: OMP_NUM_THREADS=1 uv run --no-project --with scipy --with numpy python3 Riemann/Deep/Harvest6/lp_gap.py
# quick LP estimate for Lamzouri 2611.14486 Thm Main2 (not a certificate):
# maximize f(0) + 2 int_0^theta f(a) a da  s.t. f = sum c_k hat_k (even, c>=0), int f = 1, fhat(xi) <= 0 for |xi| >= lam (sampled)
import numpy as np, sys
from scipy.optimize import linprog
theta = 4/7 - 1e-5
def run(lam, A=2.4, d=0.02, Xi=60.0, nxi=3000):
    centers = np.arange(0.0, A + 1e-12, d)          # hat at +-u_k (u_0 = 0 single hat)
    K = len(centers)
    # hat(a) = max(0, 1-|a-u|/d); even pair H_k(a) = hat(a-u)+hat(a+u) (k>0), H_0 = hat(a)
    # FT of hat centred u: d*sinc^2(d xi) e^{-2pi i u xi}; numpy sinc(x)=sin(pi x)/(pi x)
    xi = np.concatenate([np.linspace(lam, 8.0, nxi), np.linspace(8.0, Xi, nxi)])
    s2 = d*np.sinc(d*xi)**2
    M = np.empty((len(xi), K))
    M[:, 0] = s2
    for k in range(1, K):
        M[:, k] = 2*s2*np.cos(2*np.pi*centers[k]*xi)
    integ = np.array([d] + [2*d]*(K-1))               # int H_k
    # objective: H_k(0) + 2 int_0^theta H_k(a) a da (numerical)
    aa = np.linspace(0, theta, 4001); w = np.full_like(aa, aa[1]-aa[0]); w[0]*=0.5; w[-1]*=0.5
    obj = np.empty(K)
    for k in range(K):
        u = centers[k]
        Hk = np.maximum(0, 1-np.abs(aa-u)/d) + (np.maximum(0, 1-np.abs(aa+u)/d) if k>0 else 0)
        obj[k] = (1.0 if k==0 else 0.0) + 2*np.sum(w*Hk*aa)
    res = linprog(-obj, A_ub=M, b_ub=np.zeros(len(xi)), A_eq=integ[None,:], b_eq=[1.0], bounds=[(0,None)]*K, method='highs')
    return (-res.fun if res.status==0 else float('nan')), res.status
for lam in (0.75, 0.7892, 0.85):
    v, st = run(lam)
    print(f'lam={lam}: LP status {st} (2 = infeasible), value {v}')
print('Every nonnegative uniform-knot piecewise-linear f has fhat = d sinc^2(d xi) P(xi) with P periodic (period 1/d) and P(1/d) = P(0) > 0, so fhat > 0 near xi = 1/d: infeasible for any lam.')
