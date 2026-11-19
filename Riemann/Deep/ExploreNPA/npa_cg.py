#!/usr/bin/env python3
"""explore52 — column generation for the off-line-pair CUE(N) analogue (see offline_lp.py for the model).
Exact enumeration of all configurations on Z_M (on-line points + off-line pairs with depths from a grid), priced in
numpy blocks against the LP duals, so memory stays O(block) instead of O(#configs).
usage: npa_cg.py N M s1,s2,... SETS [world]   world: on (RH world: no pairs) | off (default)
"""
import sys, itertools, math, time, resource
import numpy as np
from scipy.optimize import linprog

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k): yield (k,) + rest

def zlam(lam):
    z = 1
    for j in set(lam):
        a = lam.count(j); z *= j ** a * math.factorial(a)
    return z

def rowspec(N, sets):
    spec = []
    for n in range(1, N + 1):
        parts = list(partitions(n))
        for i, lam in enumerate(parts):
            for mu in parts[i:]:
                L = len(lam) + len(mu)
                if sets == 'P' and not (len(lam) == 1 and len(mu) == 1): continue
                if sets == 'P3' and L > 3: continue
                if sets == 'L4' and L > 4: continue
                if sets == 'Q4' and not (len(lam) == len(mu) and len(lam) <= 2): continue
                # 'open top' variants: words of length >= 3 only strictly inside the RS range (|lam| <= N - 1); all pair words kept
                if sets == 'P3o' and (L > 3 or (L >= 3 and n == N)): continue
                if sets == 'L4o' and (L > 4 or (L >= 3 and n == N)): continue
                # 'closed' level 2 (an NPA principal block): every off-diagonal word comes with both diagonal partners -> drop (3,1) words
                if sets == 'L4c' and (L > 4 or len(lam) == 3 or len(mu) == 3): continue
                t = zlam(lam) if lam == mu else 0.0
                s = math.sqrt(zlam(lam) * zlam(mu))
                spec.append((lam, mu, t / s, s))
    return spec

def feats(P, spec):
    F = np.empty((len(spec), P.shape[0]))
    for r, (lam, mu, t, s) in enumerate(spec):
        a = np.prod([P[:, k] for k in lam], axis=0); b = np.prod([P[:, k] for k in mu], axis=0)
        F[r] = (a * np.conj(b)).real / s
    return F

def ks(N): return np.arange(N + 1)

def online_block(N, M, non):
    if non == 0:
        return np.zeros((1, N + 1), complex), np.zeros(1), [()]
    X = np.array(list(itertools.combinations_with_replacement(range(M), non)), dtype=np.int64)
    P = np.exp(2j * np.pi * X[:, :, None] * ks(N)[None, None, :] / M).sum(axis=1)
    srt = np.sort(X, axis=1)
    eq_prev = np.concatenate([np.zeros((len(X), 1), bool), srt[:, 1:] == srt[:, :-1]], axis=1)
    eq_next = np.concatenate([srt[:, 1:] == srt[:, :-1], np.zeros((len(X), 1), bool)], axis=1)
    ns = (~eq_prev & ~eq_next).sum(axis=1).astype(float)
    return P, ns, [tuple(x) for x in X]

def blocks(N, M, depths, world, chunk=int(__import__('os').environ.get('NPA_CHUNK', 60000))):
    """yield (P, ns, lab) with lab(i) -> label of row i: ALL configurations (rotation fixed: first pair at 0, pure on-line
    configurations as dihedral orbits; duplicates harmless). Vectorised over the remaining pairs."""
    S = len(depths); k = ks(N)
    from offline_lp import online_configs
    X = np.array(online_configs(N, M), dtype=np.int64)
    P = np.exp(2j * np.pi * X[:, :, None] * k[None, None, :] / M).sum(axis=1)
    ns = np.array([(np.unique(x, return_counts=True)[1] == 1).sum() for x in X], float)
    yield P, ns, (lambda i, X=X: (tuple(int(v) for v in X[i]), ()))
    if world == 'on': return
    W = np.array([2 * np.cosh(k * s) for s in depths])               # S x (N+1)
    E = np.exp(2j * np.pi * np.arange(M)[:, None] * k[None, :] / M)   # M x (N+1)
    items = (W[:, None, :] * E[None, :, :]).reshape(S * M, N + 1)     # item j*M + x
    for npair in range(1, N // 2 + 1):
        non = N - 2 * npair
        Pon, nson, labon = online_block(N, M, non)
        R = (np.array(list(itertools.combinations_with_replacement(range(S * M), npair - 1)), dtype=np.int64)
             if npair > 1 else np.zeros((1, 0), dtype=np.int64))
        rest = items[R].sum(axis=1) if npair > 1 else np.zeros((1, N + 1), complex)
        nR, nO = len(R), len(Pon); step = max(1, chunk // nO)
        for j0 in range(S):
            for c0 in range(0, nR, step):
                pp = W[j0][None, :] + rest[c0:c0 + step]
                Pb = (pp[:, None, :] + Pon[None, :, :]).reshape(-1, N + 1)
                nsb = np.tile(nson, len(pp))
                def lab(i, c0=c0, j0=j0, R=R, labon=labon, nO=nO):
                    r, o = divmod(i, nO); rr = R[c0 + r]
                    return (tuple(int(v) for v in labon[o]), ((0, depths[j0]),) + tuple((int(q % M), depths[q // M]) for q in rr))
                yield Pb, nsb, lab

def ns_fix(lab, ns):
    # on-line configs with a point fixed at 0: recount simple points including the 0
    out = np.empty(len(lab))
    for i, l in enumerate(lab):
        x = (0,) + l; u, c = np.unique(x, return_counts=True); out[i] = (c == 1).sum()
    return out

def rays(N, spec, nang=720, Deltas=(0, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.5, 4.0, np.inf)):
    """infinitely deep off-line pairs at infinitesimal weight: eps * e^{2Ns} fixed. Only the top words (|lam| = |mu| = N) survive:
    p_k ~ e^{ks} c_k, c_k = 1 + e^{-k Delta} w^k (two pairs, depth gap Delta, relative angle w; Delta = inf: one pair).
    Returns rows x rays (objective 0, normalisation 0)."""
    k = np.arange(N + 1); cols = []
    for D in Deltas:
        for j in range(nang if np.isfinite(D) else 1):
            ck = 1 + (np.exp(-k * D) if np.isfinite(D) else 0) * np.exp(2j * np.pi * k * j / nang)
            col = []
            for (lam, mu, t, s) in spec:
                if sum(lam) < N: col.append(0.0); continue
                a = np.prod([ck[q] for q in lam]) * np.conj(np.prod([ck[q] for q in mu])); col.append(a.real / s)
            cols.append(col)
    if N >= 6:   # three deep pairs: c_k = 1 + e^{-k D1} w1^k + e^{-k D2} w2^k
        for D1 in (0, 0.3, 1.0):
            for D2 in (0, 0.3, 1.0):
                R3 = int(__import__('os').environ.get('NPA_RAY3', 48))
                for j1 in range(R3):
                    for j2 in range(R3):
                        ck = 1 + np.exp(-k * D1) * np.exp(2j * np.pi * k * j1 / R3) + np.exp(-k * D2) * np.exp(2j * np.pi * k * j2 / R3)
                        col = []
                        for (lam, mu, t, s) in spec:
                            if sum(lam) < N: col.append(0.0); continue
                            a = np.prod([ck[q] for q in lam]) * np.conj(np.prod([ck[q] for q in mu])); col.append(a.real / s)
                        cols.append(col)
    return np.array(cols).T

RAYS = None
def solve(A, b, c):
    n = A.shape[1]
    if RAYS is not None:
        A = np.hstack([A, RAYS]); c = np.concatenate([c, np.zeros(RAYS.shape[1])])
    Aeq = np.vstack([A, np.concatenate([np.ones(n), np.zeros(A.shape[1] - n)])[None, :]]); beq = np.concatenate([b, [1.0]])
    cs = np.maximum(1.0, np.abs(Aeq).max(axis=0))
    r = linprog(c / cs, A_eq=Aeq / cs, b_eq=beq, bounds=(0, None), method='highs-ds')
    if r.status != 0: return None
    x = r.x / cs; y = r.eqlin.marginals
    resid = np.abs(Aeq @ x - beq).max()
    return r.fun, x[:n], y, resid

if __name__ == '__main__':
    N = int(sys.argv[1]); M = int(sys.argv[2]); depths = [float(t) for t in sys.argv[3].split(',')]
    sets = sys.argv[4]; world = sys.argv[5] if len(sys.argv) > 5 else 'off'
    spec = rowspec(N, sets); b = np.array([t for (_, _, t, _) in spec])
    t0 = time.time()
    if '-rays' in sys.argv and world == 'off':
        RAYS = rays(N, spec); print(f"rays: {RAYS.shape[1]} asymptotic directions (infinitely deep pairs)", flush=True)
    # initial columns: all pure on-line configurations
    g = blocks(N, M, depths, 'on'); P0, ns0, lab0 = next(g)
    A = feats(P0, spec); c = ns0 / N; labs = [lab0(i) for i in range(len(ns0))]
    it = 0
    while True:
        it += 1
        sol = solve(A, b, c)
        if sol is None: print("INFEASIBLE"); break
        val, x, y, res = sol
        # price all configurations
        best = []; nconf = 0; minrc = 0.0
        for P, ns, lab in blocks(N, M, depths, world):
            F = feats(P, spec); nconf += len(ns)
            cs = np.maximum(1.0, np.abs(np.vstack([F, np.ones((1, F.shape[1]))])).max(axis=0))
            rc = (ns / N - y[:-1] @ F - y[-1]) / cs
            idx = np.where(rc < -1e-10)[0]
            if len(idx):
                minrc = min(minrc, rc[idx].min())
                o = idx[np.argsort(rc[idx])[:200]]
                best += [(rc[i], F[:, i], ns[i] / N, lab(i)) for i in o]
        print(f"it={it} val={val:.6f} cols={A.shape[1]} resid={res:.1e} priced={nconf} minrc={minrc:.2e} t={time.time()-t0:.0f}s", flush=True)
        if not best: break
        best.sort(key=lambda z: z[0]); best = best[:400]
        A = np.hstack([A, np.array([z[1] for z in best]).T]); c = np.concatenate([c, [z[2] for z in best]]); labs += [z[3] for z in best]
        if it > 60: print("iteration cap"); break
    supp = np.where(x > 1e-12)[0]
    offw = sum(x[i] * 2 * len(labs[i][1]) / N for i in supp)
    print(f"RESULT N={N} M={M} sets={sets} world={world} depths={depths}: minSimpleOn={val:.6f} supp={len(supp)} offline-zero share={offw:.2e} "
          f"dual_obj={float(np.concatenate([b,[1.0]]) @ y):.6f} maxrss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20:.0f}MB")
    if '-v' in sys.argv:
        for i in supp[np.argsort(-x[supp])]: print(f"   w={x[i]:.3e} {labs[i]}")
