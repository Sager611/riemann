# Jung corner: type-5 contact sets near the regular simplex with rho in [rho_lo, 0.4]; aligned partition = Voronoi cones of
# the contact directions; sampled max piece diameter of the all-contact outer body Omega(P) = B(0,R) ∩ ⋂ B(p_i, 1).
import numpy as np, sys, json
from scipy.optimize import minimize
sys.path.insert(0, 'Borsuk/code')
from p6_survey import boundary, d2
from p6_cones import SIMP
rho_lo, n, seed, M = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
def make(rho):
    R = np.sqrt(rho)
    for _ in range(200):
        P = SIMP + rng.normal(scale=rng.uniform(0.02, 0.25), size=(5, 4)); P = R * P / np.linalg.norm(P, axis=1, keepdims=True)
        # push pairwise distances <= 1 by a few projection steps
        for it in range(200):
            D = np.sqrt(d2(P, P)); bad = np.argwhere(np.triu(D > 1, 1))
            if len(bad) == 0: break
            for i, j in bad:
                m = (P[i] + P[j]) / 2; v = (P[i] - P[j]) / D[i, j]; P[i] = m + v / 2 * 0.999; P[j] = m - v / 2 * 0.999
            P = R * P / np.linalg.norm(P, axis=1, keepdims=True)
        if np.sqrt(d2(P, P)).max() > 1 + 1e-9: continue
        A = np.vstack([P.T, np.ones(5)]); lam = np.linalg.lstsq(A, np.r_[0, 0, 0, 0, 1.0], rcond=None)[0]
        if np.allclose(A @ lam, np.r_[0, 0, 0, 0, 1.0], atol=1e-9) and lam.min() > 0: return P, lam
    return None, None
out = []
for t in range(n):
    rho = rng.uniform(rho_lo, 0.4); P, lam = make(rho)
    if P is None: continue
    X = boundary(rho, P, [], M, rng); U = P / np.linalg.norm(P, axis=1, keepdims=True); lab = np.argmax(X @ U.T, 1); m = 0.0
    for j in range(5):
        Y = X[lab == j]
        if len(Y) > 1: m = max(m, float(np.sqrt(d2(Y, Y).max())))
    E = np.sqrt(d2(P, P))[np.triu_indices(5, 1)]
    out.append(m); print(json.dumps({'rho': round(rho, 4), 'min_edge': round(float(E.min()), 3), 'min_lam': round(float(lam.min()), 3), 'aligned_max_piece': round(m, 4)}), flush=True)
print(json.dumps({'rho_lo': rho_lo, 'configs': len(out), 'worst': round(max(out), 4) if out else None}))
