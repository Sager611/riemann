# P6: rotated simplex-cone partitions (apex z, rotation Q) of all-contact outer bodies; minimise max piece diameter.
import numpy as np, sys, json, time
from scipy.linalg import expm
sys.path.insert(0, 'Borsuk/code')
from p6_survey import random_config, boundary
def dot(A, B): return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)
E = np.eye(5) - 1 / 5; U, S_, Vt = np.linalg.svd(E); SIMP = U[:, :4] * S_[:4]; SIMP /= np.linalg.norm(SIMP[0], keepdims=True)
IU = np.triu_indices(4, 1)
def rot(w):
    A = np.zeros((4, 4)); A[IU] = w; A -= A.T; return expm(A)
def pieces_diam(X, Q, z):
    D = SIMP @ Q.T                     # 5 cone axes
    lab = np.argmax(dot(X - z, D), axis=1); m = 0.0
    for j in range(5):
        Y = X[lab == j]
        for i0 in range(0, len(Y), 300):
            if len(Y) > 1: m = max(m, float(d2(Y[i0:i0 + 300], Y).max()))
    return np.sqrt(max(m, 0))
def optimise(X, rng, iters=150, lam=14):
    best = (9, None)
    for start in range(4):
        mean = np.concatenate([rng.normal(size=6) * 2, np.zeros(4)]); sig = np.array([0.5] * 6 + [0.03] * 4)
        for it in range(iters):
            Z = mean[None] + sig[None] * rng.normal(size=(lam, 10))
            f = np.array([pieces_diam(X, rot(z[:6]), z[6:]) for z in Z]); o = np.argsort(f)
            if f[o[0]] < best[0]: best = (f[o[0]], Z[o[0]].copy())
            mean = Z[o[:4]].mean(0); sig *= 0.985
    return best
if __name__ == '__main__':
    typ, n, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); rng = np.random.default_rng(seed)
    for k in range(n):
        rho, P, cuts, lam = random_config(typ, rng); t0 = time.time()
        X = boundary(rho, P, cuts, 3000, rng); f, z = optimise(X, rng)
        Xd = boundary(rho, P, cuts, 30000, rng); fd = pieces_diam(Xd, rot(z[:6]), z[6:])
        print(json.dumps({'type': typ, 'rho': round(rho, 4), 'phi_opt': round(f, 4), 'phi_dense': round(fd, 4), 'sec': round(time.time() - t0)}), flush=True)
