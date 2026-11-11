# PF_4 search, version 2: objective q = det[K(x_i - y_j)] / (V(x) V(y) prod_i K(x_i - y_i)), V = Vandermonde product,
# so that the coalescent limit stays finite (q -> c H_4(u)/K^4) and descent cannot collapse to 0+.  Random starts + coordinate descent.
import sys, numpy as np
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
r = 4; PI = np.pi
def K(u):
    u = np.abs(u); s = np.zeros_like(u)
    for n in range(1, 8):
        s += (2*PI**2*n**4*np.exp(9*u) - 3*PI*n**2*np.exp(5*u))*np.exp(-PI*n**2*np.exp(4*u))
    return s
def V(z):
    v = np.ones(z.shape[:-1])
    for i in range(r):
        for j in range(i+1, r): v = v*(z[..., j] - z[..., i])
    return v
def q(X, Y):
    M = K(X[..., :, None] - Y[..., None, :])
    d = np.linalg.det(M)
    return d/(np.abs(V(X)*V(Y))*np.prod(np.abs(np.diagonal(M, axis1=-2, axis2=-1)), axis=-1)), np.linalg.cond(M)
best = (np.inf, None, None)
for it in range(30):
    B = 100000
    c = rng.uniform(-0.6, 0.6, (B, 1)); sx = 10**rng.uniform(-2.5, -0.3, (B, 1)); sy = 10**rng.uniform(-2.5, -0.3, (B, 1))
    X = np.sort(c + sx*rng.normal(size=(B, r)), axis=1); Y = np.sort(sy*rng.normal(size=(B, r)), axis=1)
    Q, C = q(X, Y); Q = np.where(C < 1e11, Q, np.inf)
    i = np.argmin(Q)
    if Q[i] < best[0]: best = (Q[i], X[i].copy(), Y[i].copy())
print('random: min q =', best[0]); print(' X', best[1].tolist()); print(' Y', best[2].tolist())
q0, X0, Y0 = best
for step in [0.02, 0.01, 0.005, 0.002, 0.001]:
    for k in range(3000):
        X1 = np.sort(X0 + rng.normal(0, step, r)); Y1 = np.sort(Y0 + rng.normal(0, step, r))
        Q1, C1 = q(X1, Y1)
        if C1 < 1e11 and Q1 < q0: q0, X0, Y0 = Q1, X1, Y1
print('refined: min q =', q0); print(' X', X0.tolist()); print(' Y', Y0.tolist())
