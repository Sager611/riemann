# Canonical contact-aligned partitions for type-4/5 bodies: (A) anti-facet cones K_i = cone{-p_k : k != i} (contain p_i),
# (B) Voronoi cones of the contact directions.  Sampled piece diameters over random admissible configurations.
import numpy as np, sys, json
sys.path.insert(0, 'Borsuk/code')
from p6_survey import random_config, boundary, d2
typ, n, seed, M = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
def pieces_max(X, lab, k):
    out = []
    for j in range(k):
        Y = X[lab == j]
        out.append(0.0 if len(Y) < 2 else float(np.sqrt(max(d2(Y[i:i + 500], Y).max() for i in range(0, len(Y), 500)))))
    return out
for t in range(n):
    rho, P, cuts, lam = random_config(typ, rng)
    X = boundary(rho, P, cuts, M, rng)
    if typ == 5:
        # anti-facet cones: x = sum_k mu_k (-p_k) with mu >= 0 for k != i  <=>  coefficients of x in basis {-p_k, k != i} nonnegative
        labA = np.full(M, -1)
        for i in range(5):
            B = -np.delete(P, i, 0).T; mu = np.linalg.solve(B, X.T).T
            labA[np.all(mu >= -1e-12, 1) & (labA < 0)] = i
        labB = np.argmax(X @ (P / np.linalg.norm(P, axis=1, keepdims=True)).T, 1)
        print(json.dumps({'rho': round(rho, 4), 'minlam': round(float(lam.min()), 3), 'antifacet': [round(v, 3) for v in pieces_max(X, labA, 5)], 'voronoi': [round(v, 3) for v in pieces_max(X, labB, 5)]}), flush=True)
