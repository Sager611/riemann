# P6/P8: convex 5-cell partitions (power diagrams) of an outer body, found by minimising a softmax surrogate of
# "far pairs in a common cell"; evaluated exactly on samples (hard argmax).  Usage:
#   python p6_pd.py BODY_JSON NSAMPLE THR SEED   with BODY_JSON = {"u": [...]} (Tolmachev H(u)) or {"rho":..,"P":[[..]],"cuts":[[n,c],..]}
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from p5_T import H_params
from p6_survey import boundary

def dot(A, B): return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)

def body(spec):
    if 'u' in spec:
        rho, P = H_params(spec['u']); return rho, P, [(np.array([0, 0, 0, 1.0]), 0.5)]
    return spec['rho'], np.array(spec['P']), [(np.array(n), c) for n, c in spec.get('cuts', [])]

def interior(rho, P, cuts, N, rng):
    out = []
    while sum(len(o) for o in out) < N:
        x = rng.uniform(-1, 1, (6 * N, 4)) * np.sqrt(rho); ok = (x * x).sum(1) <= rho
        for p in P: ok &= ((x - p) ** 2).sum(1) <= 1
        for n, c in cuts: ok &= x @ n <= c
        out.append(x[ok])
    return np.concatenate(out)[:N]

def phi(X, G, c):
    lab = np.argmax(dot(X, G) + c[None], axis=1); w = 0.0
    for j in range(len(G)):
        Y = X[lab == j]
        for i0 in range(0, len(Y), 200):
            if len(Y) > 1: w = max(w, float(d2(Y[i0:i0 + 200], Y).max()))
    return np.sqrt(max(w, 0))

def far_pairs(X, t, cap, rng):
    I = []; J = []
    for i0 in range(0, len(X), 200):
        D = d2(X[i0:i0 + 200], X); ii, jj = np.nonzero(D >= t * t); m = jj > ii + i0; I.append(ii[m] + i0); J.append(jj[m])
    I = np.concatenate(I); J = np.concatenate(J)
    if len(I) > cap: s = rng.choice(len(I), cap, replace=False); I, J = I[s], J[s]
    return I, J

def optimise(X, I, J, G, c, iters=600, lr=0.02, rng=None):
    w = np.exp((np.sqrt(((X[I] - X[J]) ** 2).sum(1)) - 1.0) / 0.01); w /= w.mean()
    z = np.concatenate([G.ravel(), c]); m1 = np.zeros_like(z); m2 = np.zeros_like(z)
    for it in range(iters):
        tau = 0.05 * (0.002 / 0.05) ** (it / (iters - 1))
        G = z[:20].reshape(5, 4); c = z[20:]
        F = (dot(X, G) + c[None]) / tau; F -= F.max(1, keepdims=True); E = np.exp(F); p = E / E.sum(1, keepdims=True)
        same = (p[I] * p[J]).sum(1)                                  # prob. same cell
        # d same / d F(x)_j = p_j(x) (p_j(y) - same) / 1  (softmax derivative), times 1/tau for d/d f
        gx = p[I] * (p[J] - same[:, None]) * w[:, None]; gy = p[J] * (p[I] - same[:, None]) * w[:, None]
        gF = np.zeros_like(p); np.add.at(gF, I, gx); np.add.at(gF, J, gy); gF /= tau * len(I)
        gG = gF.T @ X; gc = gF.sum(0); g = np.concatenate([gG.ravel(), gc])
        m1 = 0.9 * m1 + 0.1 * g; m2 = 0.999 * m2 + 0.001 * g * g
        z -= lr * (m1 / (1 - 0.9 ** (it + 1))) / (np.sqrt(m2 / (1 - 0.999 ** (it + 1))) + 1e-12)
    return z[:20].reshape(5, 4), z[20:]

if __name__ == '__main__':
    spec = json.loads(sys.argv[1]); N = int(sys.argv[2]); t = float(sys.argv[3]); seed = int(sys.argv[4])
    rng = np.random.default_rng(seed); rho, P, cuts = body(spec)
    X = np.concatenate([boundary(rho, P, cuts, N, rng), interior(rho, P, cuts, N // 4, rng)])
    V = np.concatenate([boundary(rho, P, cuts, 4 * N, rng), interior(rho, P, cuts, N, rng)])
    best = (9, None); t0 = time.time()
    for rs in range(int(sys.argv[5]) if len(sys.argv) > 5 else 6):
        G = rng.normal(size=(5, 4)); G /= np.linalg.norm(G, axis=1, keepdims=True); c = np.zeros(5)
        I, J = far_pairs(X, t, 60000, rng)
        for rnd in range(3):
            G, c = optimise(X, I, J, G, c, rng=rng)
            lab = np.argmax(dot(X, G) + c[None], axis=1); bad = lab[I] == lab[J]
            I2, J2 = far_pairs(X, t, 200000, rng); viol = lab[I2] == lab[J2]
            I = np.concatenate([I, I2[viol][:20000]]); J = np.concatenate([J, J2[viol][:20000]])
        f = phi(V, G, c)
        if f < best[0]: best = (f, np.concatenate([G, c[:, None]], 1))
        print(json.dumps({'restart': rs, 'phi_valid': round(f, 4), 'sec': round(time.time() - t0, 1)}), flush=True)
    print(json.dumps({'spec': spec, 'best_phi': round(best[0], 4)}), flush=True)
    np.save('Borsuk/.scratch/pd_best.npy', best[1])
