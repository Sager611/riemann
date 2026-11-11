# Global PF_3 / PF_4 search for the de Bruijn-Newman kernel K(u) = Phi(|u|) (Michalowski 2602.20313, Open Problem 1:
# "prove or disprove K in PF_4 for all admissible configurations, not just Toeplitz ones").
# PF_r: det[K(x_i - y_j)]_{r x r} >= 0 for all x_1<...<x_r, y_1<...<y_r.  Double-precision random search + local refinement;
# negative hits are re-evaluated in mpmath (see pf_verify.py).
import sys, numpy as np
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
r = int(sys.argv[1]) if len(sys.argv) > 1 else 4
PI = np.pi
def K(u):
    u = np.abs(u); s = np.zeros_like(u)
    for n in range(1, 8):
        s += (2*PI**2*n**4*np.exp(9*u) - 3*PI*n**2*np.exp(5*u))*np.exp(-PI*n**2*np.exp(4*u))
    return s
def ratio(X, Y):          # X, Y: (B, r) sorted; returns det / prod diag (sign test)
    M = K(X[:, :, None] - Y[:, None, :])
    d = np.linalg.det(M)
    sc = np.prod(np.abs(np.diagonal(M, axis1=1, axis2=2)), axis=1)
    return d/sc, M
best = []
W = float(sys.argv[3]) if len(sys.argv) > 3 else 0.8
for it in range(40):
    B = 200000
    X = np.sort(rng.uniform(-W, W, (B, r)), axis=1); Y = np.sort(rng.uniform(-W, W, (B, r)), axis=1)
    q, M = ratio(X, Y)
    # condition guard: skip configurations whose matrices are numerically singular
    cond = np.linalg.cond(M)
    ok = cond < 1e10
    q = np.where(ok, q, np.inf)
    i = np.argmin(q)
    best.append((q[i], X[i].copy(), Y[i].copy(), cond[i]))
best.sort(key=lambda b: b[0])
print('r =', r, 'W =', W, 'min ratio over', 40*200000, 'random configs:', best[0][0], 'cond', best[0][3])
print('X', best[0][1].tolist()); print('Y', best[0][2].tolist())
# local refinement by random-walk descent on the normalized ratio
q0, X0, Y0, c0 = best[0]
step = 0.05
for k in range(20000):
    X1 = np.sort(X0 + rng.normal(0, step, r)); Y1 = np.sort(Y0 + rng.normal(0, step, r))
    q1, M1 = ratio(X1[None], Y1[None])
    if np.linalg.cond(M1[0]) < 1e10 and q1[0] < q0:
        q0, X0, Y0 = q1[0], X1, Y1
    if k % 5000 == 4999: step /= 2
print('refined min ratio', q0); print('X', X0.tolist()); print('Y', Y0.tolist())
