# Structure of a 5-colouring of the extreme type-3 body (rho = 1/3, unit equilateral contact triangle, slab cuts).
import numpy as np, sys, json
sys.path.insert(0, 'Borsuk/code')
from p6_survey import boundary, cegar, d2
rho, M, seed, thr = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
R = np.sqrt(rho); rng = np.random.default_rng(seed)
th = 2 * np.pi * np.arange(3) / 3; P = R * np.stack([np.cos(th), np.sin(th), 0 * th, 0 * th], 1)
cuts = [(np.array([0, 0, 1.0, 0]), 0.5), (np.array([0, 0, -1.0, 0]), 0.5), (np.array([0, 0, 0, 1.0]), 0.5)]
C = boundary(rho, P, cuts, M, rng)
r = cegar(C, thr, TL=900, rng=rng, want=True)
print(r[0], r[1])
if r[0] == 'COL':
    Z, zc = r[4], r[5]
    np.save('Borsuk/.scratch/ext3_col.npy', np.c_[Z, zc])
    a = Z[:, :2]; y = Z[:, 2:]; ra = np.linalg.norm(a, axis=1); ry = np.linalg.norm(y, axis=1)
    ang_a = np.degrees(np.arctan2(a[:, 1], a[:, 0])) % 360; ang_y = np.degrees(np.arctan2(y[:, 1], y[:, 0])) % 360
    for j in range(5):
        m = zc == j; Y = Z[m]
        dm = np.sqrt(max(d2(Y[i:i + 500], Y).max() for i in range(0, len(Y), 500)))
        print(json.dumps({'col': j, 'n': int(m.sum()), 'diam': round(float(dm), 4), 'centroid': np.round(Y.mean(0), 3).tolist(),
                          'ra_mean': round(float(ra[m].mean()), 3), 'ry_mean': round(float(ry[m].mean()), 3),
                          'ang_a_hist': np.histogram(ang_a[m & (ra > 0.15)], bins=12, range=(0, 360))[0].tolist(),
                          'ang_y_hist': np.histogram(ang_y[m & (ry > 0.15)], bins=12, range=(0, 360))[0].tolist()}))
