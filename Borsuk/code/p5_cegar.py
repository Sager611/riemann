# P5 barrier search (CEGAR): is there a finite Y ⊂ H(u) whose far graph (pairs at distance >= 1) is not 5-colourable?
# Y grows by candidate boundary points that have a far neighbour in every colour class of the current colouring.
import numpy as np, sys, time, json
from pysat.solvers import Solver
from tabu import colour_tabu
sys.path.insert(0, 'Borsuk/code')
from p5_T import H_params, boundary_sample

THR = 1.0
def dot(A, B):   # A @ B.T without BLAS (Accelerate caches a buffer per shape)
    return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))

def d2(A, B):
    return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)

PREV = {}
NOSAT = False; SAVE = 'Borsuk/.scratch/p5_hardY.npy'
def colour(X, k=5):
    N = len(X); UU = []; VV = []
    for i0 in range(0, N, 64):
        D2 = d2(X[i0:i0 + 64], X)
        ii, jj = np.nonzero(D2 >= THR * THR); m = jj > ii + i0
        UU.append(ii[m] + i0); VV.append(jj[m])
    U = np.concatenate(UU); V = np.concatenate(VV)
    init = np.random.default_rng(N).integers(0, k, N)
    if 'col' in PREV: init[:len(PREV['col'])] = PREV['col'][:N]
    col, c = colour_tabu(N, U, V, k, iters=300000, seed=N, init=init)
    if c == 0: PREV['col'] = col; return col
    for sd in range(2, 10):
        col, c = colour_tabu(N, U, V, k, iters=1000000 * sd, seed=N + sd, init=init if sd % 2 else None)
        if c == 0: PREV['col'] = col; return col
    print('tabu failed (conflicts', c, ') at N', N, file=sys.stderr, flush=True)
    if NOSAT: np.save(SAVE, X); return None
    s = Solver(name='g4'); first = None
    for v in range(N): s.add_clause([v * k + c + 1 for c in range(k)])
    for i0 in range(0, N, 64):
        D2 = d2(X[i0:i0 + 64], X)
        ii, jj = np.nonzero(D2 >= THR * THR); m = jj > ii + i0
        for a, b in zip((ii[m] + i0).tolist(), jj[m].tolist()):
            if first is None: first = a
            for c in range(k): s.add_clause([-(a * k + c + 1), -(b * k + c + 1)])
    if first is not None: s.add_clause([first * k + 1])
    if not s.solve(): s.delete(); return None
    m = s.get_model(); s.delete()
    return np.array([next(c for c in range(k) if m[v * k + c] > 0) for v in range(N)])

def run(u, M, seed, B):
    rng = np.random.default_rng(seed); rho, P = H_params(u)
    C = boundary_sample(rho, P, M, rng)
    Y = list(rng.choice(M, 300, replace=False)); t0 = time.time()
    for it in range(400):
        if time.time() - t0 > TL: print(json.dumps({'u': u, 'result': 'time limit', 'Y': len(Y)}), flush=True); return
        X = C[Y]; col = colour(X)
        if col is None:
            print(json.dumps({'u': u, 'it': it, 'Y': len(Y), 'result': 'NOT 5-colourable'}), flush=True); np.save('Borsuk/.scratch/p5_barrier_%s.npy' % '_'.join(map(str, u)), X); return
        worst = np.full(M, np.inf)
        for j in range(5):
            Yj = X[col == j]
            if len(Yj) == 0: worst[:] = -np.inf; break
            f = np.zeros(M)
            for i0 in range(0, M, 1000):
                f[i0:i0 + 1000] = d2(C[i0:i0 + 1000], Yj).max(1)
            worst = np.minimum(worst, f)
        bad = np.flatnonzero(worst >= THR * THR)
        if it % 10 == 0: print(json.dumps({'u': u, 'it': it, 'Y': len(Y), 'bad': len(bad), 'max_minfar': round(float(worst.max()), 4), 'sec': round(time.time() - t0, 1)}), flush=True)
        if len(bad) == 0:
            # greedy joint extension over all candidates (random order, blocks)
            Zbuf = np.zeros((M + len(X), 4)); Zbuf[:len(X)] = X; nz = len(X)
            zc = np.full(M + len(X), -1); zc[:nz] = col; fail = []
            Zn = np.zeros(M + len(X)); Zn[:nz] = (X * X).sum(1)
            order = rng.permutation(M); Dfix = np.empty((20, M + len(X))); ZT = Zbuf.T  # fixed shapes: macOS malloc keeps growing-size blocks
            for b0 in range(0, M, 20):
                idx = order[b0:b0 + 20]; Cb = C[idx]
                Dfix[:len(idx)] = Zn[None, :]
                for kk in range(4): Dfix[:len(idx)] -= 2 * Cb[:, kk, None] * ZT[kk][None, :]
                Dfix[:len(idx)] += (Cb * Cb).sum(1)[:, None]
                far = Dfix[:len(idx)] >= THR * THR
                allowed = np.ones((len(idx), 5), bool)
                for j in range(5): allowed[:, j] = ~np.any(far & (zc == j)[None, :], axis=1)
                Db = d2(Cb, Cb) >= THR * THR
                cb = np.full(len(idx), -1)
                for t in range(len(idx)):
                    al = allowed[t].copy()
                    for s_ in np.flatnonzero(Db[t, :t]):
                        if cb[s_] >= 0: al[cb[s_]] = False
                    if al.any(): cb[t] = int(np.flatnonzero(al)[0])
                    else: fail.append(int(idx[t]))
                ok = cb >= 0; q = int(ok.sum()); Zbuf[nz:nz + q] = Cb[ok]; Zn[nz:nz + q] = (Cb[ok] ** 2).sum(1); zc[nz:nz + q] = cb[ok]; nz += q; ZT = Zbuf.T
            print(json.dumps({'it': it, 'greedy_fail': len(fail)}), flush=True)
            if not fail:
                print(json.dumps({'u': u, 'thr': THR, 'result': 'joint colouring of all candidates found', 'Y': len(Y), 'M': M}), flush=True); return
            Y.extend(fail[:B]); continue
        add = bad[np.argsort(-worst[bad])][:B]; Y.extend(add.tolist())

if __name__ == '__main__':
    us = json.loads(sys.argv[1]); M = int(sys.argv[2]); seed = int(sys.argv[3]); B = int(sys.argv[4])
    THR = float(sys.argv[5]); TL = float(sys.argv[6]); NOSAT = len(sys.argv) > 7 and sys.argv[7] == 'nosat'
    for u in us: run(u, M, seed, B)
