# Illumination test: approximate normal cones of conv(F) for a nearly complete diameter-1 set F (sampled directions u,
# x(u) = argmax_F <x, u>), then search a rotated regular simplex of directions w_j such that every point x is
# illuminated by some w_j (<n, w_j> < 0 for all sampled normals n at x).  Illumination by 5 directions implies a
# covering by 5 smaller homothets, hence a Borsuk partition.  Reports the fraction of points illuminated.
import numpy as np, sys, json
sys.path.insert(0, 'Borsuk/code')
from p6_cones import rot, SIMP
F = np.load(sys.argv[1]); nd = int(sys.argv[2]); rng = np.random.default_rng(0)
U = rng.normal(size=(nd, 4)); U /= np.linalg.norm(U, axis=1, keepdims=True)
arg = np.empty(nd, int)
for i0 in range(0, nd, 2000): arg[i0:i0 + 2000] = np.argmax(U[i0:i0 + 2000] @ F.T, 1)
owners = np.unique(arg); groups = [U[arg == o] for o in owners]
rad = []
for g in groups:
    c = g.mean(0); c /= np.linalg.norm(c); rad.append(np.degrees(np.arccos(np.clip((g @ c).min(), -1, 1))))
rad = np.array(rad)
print(json.dumps({'set': sys.argv[1], 'boundary_points_hit': len(owners), 'normal_cone_radius_deg_quantiles': np.round(np.percentile(rad, [50, 90, 99, 100]), 1).tolist()}), flush=True)
fat = [g for g, r in zip(groups, rad) if r > 5]
def illum(z):
    W = SIMP @ rot(z).T; ok = 0
    for g in fat:
        if ((g @ W.T).max(0) < 0).any(): ok += 1
    return ok
best = (-1, None)
for rs in range(6):
    z = rng.uniform(-np.pi, np.pi, 6); f = illum(z); sig = 0.3
    for it in range(200):
        Z = z[None] + sig * rng.normal(size=(8, 6)); fs = [illum(q) for q in Z]; o = int(np.argmax(fs))
        if fs[o] > f: z, f = Z[o], fs[o]
        else: sig *= 0.93
        if sig < 1e-3 or f == len(fat): break
    if f > best[0]: best = (f, z)
print(json.dumps({'fat_points(radius>5deg)': len(fat), 'best_illuminated': best[0], 'all': best[0] == len(fat)}), flush=True)
