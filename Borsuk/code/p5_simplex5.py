# Sets of diameter 1 in R^4 that contain a regular unit simplex W: they lie in Omega5 = B(0, sqrt(2/5)) ∩ ∩_i B(w_i, 1)
# (Jung equality: circumcentre = centre of W).  Max piece diameter of Omega5 under the 5 cones aligned with W.
import numpy as np, sys
from scipy.optimize import minimize
rng = np.random.default_rng(0)
E = np.eye(5) - 1 / 5; U, S, Vt = np.linalg.svd(E); B = U[:, :4] * S[:4]; B /= np.linalg.norm(B[0]); R = np.sqrt(0.4); W = R * B
print('simplex edges', np.round([np.linalg.norm(W[i] - W[j]) for i in range(5) for j in range(i + 1, 5)], 6).tolist())
def inside(x): return (x @ x <= 0.4 + 1e-12) and all(((x - w) ** 2).sum() <= 1 + 1e-12 for w in W)
def boundary(N):
    d = rng.normal(size=(N, 4)); d /= np.linalg.norm(d, axis=1, keepdims=True); s = np.full(N, R)
    for w in W:
        dw = d @ w; s = np.minimum(s, dw + np.sqrt(dw * dw + 1 - 0.4))
    return d * s[:, None]
X = boundary(60000); lab = np.argmax(X @ W.T, axis=1)
best = 0; pair = None
for j in range(5):
    Y = X[lab == j]
    n2 = (Y * Y).sum(1)
    for i0 in range(0, len(Y), 100):
        D = n2[i0:i0 + 100, None] + n2[None] - 2 * sum(Y[i0:i0 + 100, k, None] * Y[None, :, k] for k in range(4)); m = D.max()
        if m > best: best = m; a, b = np.unravel_index(np.argmax(D), D.shape); pair = (Y[i0 + a], Y[b], j)
print('sampled max piece diameter', round(np.sqrt(best), 5))
# local maximisation of |x-y|^2 over x, y in Omega5 ∩ cone_j from the best pair
x0, y0, j = pair
cons = []
for z in (0, 1):
    sl = slice(4 * z, 4 * z + 4)
    cons.append({'type': 'ineq', 'fun': lambda v, sl=sl: 0.4 - v[sl] @ v[sl]})
    for w in W: cons.append({'type': 'ineq', 'fun': lambda v, sl=sl, w=w: 1 - ((v[sl] - w) ** 2).sum()})
    for k in range(5):
        if k != j: cons.append({'type': 'ineq', 'fun': lambda v, sl=sl, k=k: v[sl] @ (W[j] - W[k])})
res = minimize(lambda v: -((v[:4] - v[4:]) ** 2).sum(), np.concatenate([x0, y0]), constraints=cons, method='SLSQP', options={'ftol': 1e-14, 'maxiter': 500})
print('local max piece diameter', round(np.sqrt(-res.fun), 6), 'success', res.success)
