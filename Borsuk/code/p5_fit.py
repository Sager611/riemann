# P5 Stage 1: 5-cell power diagram (5 affine functions, cell = argmax) for H(u) with all cell diameters < 1.
# Smooth tabu colouring of a sample (far pairs must differ, near pairs should agree) -> LP fit of affine functions
# -> evaluate on a dense boundary sample -> add violating points, repeat.
import numpy as np, sys, json, time
from scipy.optimize import linprog
sys.path.insert(0, 'Borsuk/code')
from p5_T import H_params, boundary_sample
from tabu import colour_smooth

def dot(A, B): return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)

def interior_sample(rho, P, N, rng):
    out = []
    while sum(len(o) for o in out) < N:
        x = rng.uniform(-1, 1, (8 * N, 4)) * np.sqrt(rho)
        ok = ((x * x).sum(1) <= rho) & (x[:, 3] <= 0.5)
        for p in P: ok &= ((x - p) ** 2).sum(1) <= 1
        out.append(x[ok])
    return np.concatenate(out)[:N]

from scipy.sparse import coo_matrix
def fit(S, col, k=5, margin=1e-3):
    n = len(S); nv = 5 * k            # variables: g_j (4), c_j (1) for each j; slack per point
    R = []; Cc = []; Vv = []; r = 0
    for i in range(n):
        j = int(col[i])
        for kk in range(k):
            if kk == j: continue
            # f_kk(x) - f_j(x) + margin <= s_i
            for q in range(4):
                R += [r, r]; Cc += [5 * kk + q, 5 * j + q]; Vv += [S[i, q], -S[i, q]]
            R += [r, r, r]; Cc += [5 * kk + 4, 5 * j + 4, nv + i]; Vv += [1.0, -1.0, -1.0]; r += 1
    A = coo_matrix((Vv, (R, Cc)), shape=(r, nv + n)).tocsr(); b = np.full(r, -margin)
    c = np.concatenate([np.zeros(nv), np.ones(n)])
    bounds = [(-50, 50)] * nv + [(0, None)] * n
    res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method='highs')
    z = res.x[:nv].reshape(k, 5); return z[:, :4], z[:, 4], res.x[nv:]

def phi(X, G, c):
    lab = np.argmax(X @ G.T + c[None], axis=1); worst = 0; arg = None
    for j in range(len(G)):
        Y = X[lab == j]
        if len(Y) < 2: continue
        for i0 in range(0, len(Y), 100):
            D = d2(Y[i0:i0 + 100], Y); m = D.max()
            if m > worst: worst = m; a, bb = np.unravel_index(np.argmax(D), D.shape); arg = (Y[i0 + a], Y[bb])
    return np.sqrt(worst), arg

def knn_edges(S, kn=8):
    U = []; V = []
    for i0 in range(0, len(S), 500):
        D = d2(S[i0:i0 + 500], S); np.fill_diagonal(D[:, i0:i0 + 500], np.inf) if False else None
        for r_ in range(D.shape[0]): D[r_, i0 + r_] = np.inf
        nb = np.argsort(D, axis=1)[:, :kn]
        U.append(np.repeat(np.arange(i0, i0 + D.shape[0]), kn)); V.append(nb.ravel())
    U = np.concatenate(U); V = np.concatenate(V); m = U < V
    return U[m], V[m]

def far_edges(S, thr):
    U = []; V = []
    for i0 in range(0, len(S), 200):
        D = d2(S[i0:i0 + 200], S); ii, jj = np.nonzero(D >= thr * thr); m = jj > ii + i0; U.append(ii[m] + i0); V.append(jj[m])
    return np.concatenate(U), np.concatenate(V)

if __name__ == '__main__':
    u = json.loads(sys.argv[1]); thr = float(sys.argv[2]); seed = int(sys.argv[3])
    rng = np.random.default_rng(seed); rho, P = H_params(u)
    S = np.concatenate([boundary_sample(rho, P, 2500, rng), interior_sample(rho, P, 700, rng)])
    D = boundary_sample(rho, P, 40000, rng); t0 = time.time(); col = None
    for rnd in range(12):
        FU, FV = far_edges(S, thr); NU, NV = knn_edges(S)
        init = None if col is None else np.concatenate([col, np.random.default_rng(rnd).integers(0, 5, len(S) - len(col))])
        col, fc, nd = colour_smooth(len(S), FU, FV, NU, NV, 5, iters=400000, seed=rnd + 1, W=1000, init=init)
        G, c, sl = fit(S, col)
        f_s, _ = phi(S, G, c); f_d, arg = phi(D, G, c)
        print(json.dumps({'u': u, 'round': rnd, 'S': len(S), 'far_conf': fc, 'near_dis': nd, 'lp_misfit': int((sl > 1e-6).sum()), 'phi_S': round(f_s, 4), 'phi_dense': round(f_d, 4), 'sec': round(time.time() - t0, 1)}), flush=True)
        if f_d < 1: np.save('Borsuk/.scratch/p5_pd_%s.npy' % '_'.join(map(str, u)), np.concatenate([G, c[:, None]], 1)); break
        # add dense points violating: farthest pairs in cells
        lab = np.argmax(D @ G.T + c[None], axis=1); add = []
        for j in range(5):
            Y = np.flatnonzero(lab == j)
            if len(Y) < 2: continue
            II = []; JJ = []
            for i0 in range(0, len(Y), 100):
                Dj = d2(D[Y[i0:i0 + 100]], D[Y]); ii, jj = np.nonzero(Dj >= thr * thr); II.append(ii + i0); JJ.append(jj)
            ii = np.concatenate(II); jj = np.concatenate(JJ)
            if len(ii): sel = rng.choice(len(ii), min(60, len(ii)), replace=False); add += Y[ii[sel]].tolist() + Y[jj[sel]].tolist()
        S = np.concatenate([S, D[np.unique(add)]])
