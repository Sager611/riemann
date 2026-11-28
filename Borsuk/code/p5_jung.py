# Jung corner u = (1,0,0,0): H(u) contains the regular unit simplex P1,P2,P3,p4,p5.  Evaluate the 5 cones aligned
# with the simplex vertices, then locally optimise a power diagram (5 affine functions) by an evolution strategy.
import numpy as np, sys, json
sys.path.insert(0, 'Borsuk/code')
from p5_T import H_params, boundary_sample
def dot(A, B): return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)
def phi(X, G, c):
    lab = np.argmax(dot(X, G) + c[None], axis=1); w = 0
    for j in range(len(G)):
        Y = X[lab == j]
        for i0 in range(0, len(Y), 200):
            if len(Y) > 1: w = max(w, d2(Y[i0:i0 + 200], Y).max())
    return np.sqrt(w)
u = json.loads(sys.argv[1]); rng = np.random.default_rng(0); rho, P = H_params(u)
R = np.sqrt(rho); h = P[0, 2]
# the two further simplex vertices: (0,0,h - s, +-w)
s = (2 / 3 + h * h - rho) / (2 * h) if False else None
# solve: points at distance 1 from P1..P3 on |x|^2 = rho: x = (0,0,z,w) with z^2+w^2 = rho, (z-h)^2 + w^2 = 1 - tau
tau = rho - h * h; z = (rho - (1 - tau) + h * h) / (2 * h); w = np.sqrt(max(rho - z * z, 0))
W = np.vstack([P, [0, 0, z, w], [0, 0, z, -w]])
print('simplex check: pairwise dists', np.round(np.sqrt(d2(W, W)[np.triu_indices(5, 1)]), 4).tolist())
X = boundary_sample(rho, P, 12000, rng); Xd = boundary_sample(rho, P, 40000, rng)
G0 = W.copy(); c0 = np.zeros(5)
print('aligned cones: phi(12k) =', round(phi(X, G0, c0), 4), ' phi(40k) =', round(phi(Xd, G0, c0), 4))
best = (phi(X, G0, c0), np.concatenate([G0.ravel(), c0])); mean = best[1].copy(); sig = 0.05
Xo = boundary_sample(rho, P, 4000, rng)
for it in range(120):
    Z = mean[None] + sig * rng.normal(size=(16, 25))
    f = np.array([phi(Xo, z[:20].reshape(5, 4), z[20:]) for z in Z]); o = np.argsort(f)
    if f[o[0]] < best[0]: best = (f[o[0]], Z[o[0]].copy())
    mean = 0.5 * mean + 0.5 * Z[o[:4]].mean(0); sig *= 0.98
z = best[1]; print('ES best: phi(4k) =', round(best[0], 4), ' phi(40k) =', round(phi(Xd, z[:20].reshape(5, 4), z[20:]), 4))
