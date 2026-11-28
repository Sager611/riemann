# Fraction of uniformly random fan rotations (apex = circumcentre 0) that split a saved diameter-1 set F into pieces with
# max diameter < t (t = 0.999, 0.99, 0.98).  Also the distribution of max piece diameters.
import numpy as np, sys, json
sys.path.insert(0, 'Borsuk/code')
from p6_survey import d2
from p6_cones import SIMP
F = np.load(sys.argv[1]); N = int(sys.argv[2]); rng = np.random.default_rng(0)
n = len(F); D2 = np.empty((n, n), np.float32)
for i0 in range(0, n, 500): D2[i0:i0 + 500] = d2(F[i0:i0 + 500], F)
vals = []
for k in range(N):
    A = rng.normal(size=(4, 4)); Q, Rr = np.linalg.qr(A); Q = Q * np.sign(np.diag(Rr))[None, :]
    lab = np.argmax(F @ (SIMP @ Q.T).T, 1); m = 0.0
    for j in range(5):
        idx = np.flatnonzero(lab == j)
        if len(idx) > 1: m = max(m, float(D2[np.ix_(idx, idx)].max()))
    vals.append(np.sqrt(m))
vals = np.array(vals)
print(json.dumps({'set': sys.argv[1], 'n': n, 'rotations': N, 'frac<0.999': round(float((vals < 0.999).mean()), 4), 'frac<0.99': round(float((vals < 0.99).mean()), 4), 'frac<0.98': round(float((vals < 0.98).mean()), 4), 'min': round(float(vals.min()), 4), 'median': round(float(np.median(vals)), 4)}))
