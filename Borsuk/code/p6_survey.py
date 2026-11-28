# P6 margin survey: all-contact outer bodies Omega = B(0, sqrt rho) ∩ ∩ B(p_i, 1) ∩ frame cuts, for random admissible
# minimal contact sets of type 3 (triangle through 0), 4 (tetrahedron through 0 in x4 = 0), 5 (4-simplex).
# Estimate beta5 of a boundary sample by CEGAR (tabu colourings, greedy joint extension) at several thresholds.
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from tabu import colour_tabu

def dot(A, B): return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)

def boundary(rho, P, cuts, N, rng):
    d = rng.normal(size=(N, 4)); d /= np.linalg.norm(d, axis=1, keepdims=True); s = np.full(N, np.sqrt(rho))
    for p in P:
        dp = d @ p; s = np.minimum(s, dp + np.sqrt(dp * dp + 1 - rho))
    for n, c in cuts:          # <x, n> <= c, c > 0
        dn = d @ n; pos = dn > 1e-12; s[pos] = np.minimum(s[pos], c / dn[pos])
    return d * s[:, None] * (1 - 1e-12)

def admissible(P, rho):
    D = d2(P, P); return np.all(D <= 1 + 1e-12)

def random_config(typ, rng):
    for _ in range(100000):
        rho = rng.uniform(0.30, {3: 1 / 3, 4: 3 / 8, 5: 0.4}[typ]); R = np.sqrt(rho)
        if typ == 3:
            th = np.sort(rng.uniform(0, 2 * np.pi, 3)); P = R * np.stack([np.cos(th), np.sin(th), 0 * th, 0 * th], 1)
            phi = rng.uniform(0, 2 * np.pi); v1 = np.array([0, 0, np.cos(phi), np.sin(phi)]); v2 = np.array([0, 0, -np.sin(phi), np.cos(phi)])
            cuts = [(v1, 0.5), (-v1, 0.5), (v2, 0.5)]
        elif typ == 4:
            v = rng.normal(size=(4, 3)); v /= np.linalg.norm(v, axis=1, keepdims=True); P = R * np.hstack([v, np.zeros((4, 1))])
            cuts = [(np.array([0, 0, 0, 1.0]), 0.5)]
        else:
            v = rng.normal(size=(5, 4)); v /= np.linalg.norm(v, axis=1, keepdims=True); P = R * v; cuts = []
        if not admissible(P, rho): continue
        # 0 in the relative interior of conv P: solve sum lambda_i p_i = 0, sum lambda = 1, lambda > 0
        A = np.vstack([P.T, np.ones(len(P))]); b = np.zeros(5); b[-1] = 1
        lam, *_ = np.linalg.lstsq(A, b, rcond=None)
        if np.allclose(A @ lam, b, atol=1e-9) and np.all(lam > 0.02): return rho, P, cuts, lam
    raise RuntimeError('no config')

def cegar(C, thr, B=100, TL=300, rng=None, want=False):
    M = len(C); Y = list(rng.choice(M, 300, replace=False)); t0 = time.time(); prev = None
    while time.time() - t0 < TL:
        X = C[Y]; N = len(X); UU = []; VV = []
        for i0 in range(0, N, 64):
            D = d2(X[i0:i0 + 64], X); ii, jj = np.nonzero(D >= thr * thr); m = jj > ii + i0; UU.append(ii[m] + i0); VV.append(jj[m])
        U = np.concatenate(UU); V = np.concatenate(VV); col = None
        init = np.random.default_rng(N).integers(0, 5, N)
        if prev is not None: init[:len(prev)] = prev
        for sd in range(1, 7):
            c_, conf = colour_tabu(N, U, V, 5, iters=400000 * sd, seed=N + sd, init=init if sd % 2 else None)
            if conf == 0: col = c_; break
        if col is None: return ('NOT', N)
        prev = col
        worst = np.full(M, np.inf)
        for j in range(5):
            Yj = X[col == j]
            if len(Yj) == 0: worst[:] = -np.inf; break
            f = np.zeros(M)
            for i0 in range(0, M, 1000): f[i0:i0 + 1000] = d2(C[i0:i0 + 1000], Yj).max(1)
            worst = np.minimum(worst, f)
        bad = np.flatnonzero(worst >= thr * thr)
        if len(bad) == 0:
            # greedy joint extension
            Zb = np.zeros((M + N, 4)); Zb[:N] = X; zc = np.full(M + N, -1); zc[:N] = col; nz = N; Zn = np.zeros(M + N); Zn[:N] = (X * X).sum(1)
            Df = np.empty((20, M + N)); fail = []
            for b0 in range(0, M, 20):
                idx = rng.permutation(M)[b0:b0 + 20] if False else np.arange(b0, min(M, b0 + 20)); Cb = C[idx]
                Df[:len(idx)] = Zn[None, :]
                for k in range(4): Df[:len(idx)] -= 2 * Cb[:, k, None] * Zb.T[k][None, :]
                Df[:len(idx)] += (Cb * Cb).sum(1)[:, None]; far = Df[:len(idx)] >= thr * thr
                allowed = np.ones((len(idx), 5), bool)
                for j in range(5): allowed[:, j] = ~np.any(far & (zc == j)[None, :], axis=1)
                Db = d2(Cb, Cb) >= thr * thr; cb = np.full(len(idx), -1)
                for t in range(len(idx)):
                    al = allowed[t].copy()
                    for s_ in np.flatnonzero(Db[t, :t]):
                        if cb[s_] >= 0: al[cb[s_]] = False
                    if al.any(): cb[t] = int(np.flatnonzero(al)[0])
                    else: fail.append(int(idx[t]))
                ok = cb >= 0; q = int(ok.sum()); Zb[nz:nz + q] = Cb[ok]; Zn[nz:nz + q] = (Cb[ok] ** 2).sum(1); zc[nz:nz + q] = cb[ok]; nz += q
            if not fail:
                if want:
                    full = np.full(M, -1); full[np.arange(M)] = -1
                    return ('COL', N, X, col, Zb[:nz], zc[:nz])
                return ('COL', N)
            Y.extend(fail[:B]); continue
        Y.extend(bad[np.argsort(-worst[bad])][:B].tolist())
    return ('TIME', len(Y))

if __name__ == '__main__':
    typ, n, seed, M = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    thrs = [float(x) for x in sys.argv[5].split(',')]
    rng = np.random.default_rng(seed)
    for k in range(n):
        rho, P, cuts, lam = random_config(typ, rng)
        C = boundary(rho, P, cuts, M, rng); res = {}
        for t in thrs:
            r = cegar(C, t, rng=rng); res[t] = r[0]
            if r[0] == 'COL': break      # colourable at t: higher thresholds too
        print(json.dumps({'type': typ, 'k': k, 'rho': round(rho, 4), 'minlam': round(float(lam.min()), 3), 'maxpair': round(float(np.sqrt(d2(P, P).max())), 3), 'res': res}), flush=True)
