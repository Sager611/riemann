# Split-vertex type-5 contact sets: unit triangle + two points on its 'unit circle' at angle phi (rho from 3/8 to 2/5),
# aligned Voronoi cones of the contact directions on the all-contact outer body; also the best rotated fan.
import numpy as np, sys, json
sys.path.insert(0, 'Borsuk/code')
from p6_survey import boundary, d2
from p6_cones import SIMP, rot
rng = np.random.default_rng(0)
for phi_deg in [5, 15, 25, 38.3, 50, 60, 75.5]:
    # triangle p2,p3,p4 in plane spanned by e1,e2 centred at c; points p0,p1 on the circle of radius sqrt(2/3) around c in the (e3,e4)-plane
    phi = np.radians(phi_deg); r = np.sqrt(2 / 3); t = 2 * np.pi * np.arange(3) / 3
    T = np.stack([np.cos(t), np.sin(t), 0 * t, 0 * t], 1) / np.sqrt(3)
    Q = np.stack([[0, 0, r * np.cos(a), r * np.sin(a)] for a in (-phi / 2, phi / 2)])
    Pts = np.vstack([Q, T])
    # circumcentre of the 5 points (solve |x - p|^2 equal): least squares on the 4-dim affine condition
    A = 2 * (Pts[1:] - Pts[0]); b = (Pts[1:] ** 2).sum(1) - (Pts[0] ** 2).sum(); cc = np.linalg.lstsq(A, b, rcond=None)[0]
    P = Pts - cc; R2 = float((P[0] ** 2).sum())
    Aa = np.vstack([P.T, np.ones(5)]); lam = np.linalg.lstsq(Aa, np.r_[0, 0, 0, 0, 1.0], rcond=None)[0]
    X = boundary(R2, P, [], 5000, rng); U = P / np.linalg.norm(P, axis=1, keepdims=True)
    def maxpiece(lab):
        m = 0.0
        for j in range(5):
            Y = X[lab == j]
            if len(Y) > 1: m = max(m, float(np.sqrt(d2(Y, Y).max())))
        return m
    al = maxpiece(np.argmax(X @ U.T, 1))
    best = 9
    for rs in range(3):
        z = np.r_[rng.uniform(-np.pi, np.pi, 6), np.zeros(4)]; f = maxpiece(np.argmax((X - z[6:]) @ (SIMP @ rot(z[:6]).T).T, 1)); sig = np.r_[np.full(6, 0.3), np.full(4, 0.03)]
        for it in range(150):
            Z = z[None] + sig[None] * rng.normal(size=(8, 10)); fs = [maxpiece(np.argmax((X - q[6:]) @ (SIMP @ rot(q[:6]).T).T, 1)) for q in Z]; o = int(np.argmin(fs))
            if fs[o] < f: z, f = Z[o], fs[o]
            else: sig *= 0.92
            if sig[0] < 2e-3: break
        best = min(best, f)
    print(json.dumps({'phi': phi_deg, 'rho': round(R2, 4), 'edge_p0p1': round(float(np.linalg.norm(Q[0] - Q[1])), 3), 'lam_min': round(float(lam.min()), 3), 'aligned': round(al, 4), 'best_fan': round(best, 4)}), flush=True)
