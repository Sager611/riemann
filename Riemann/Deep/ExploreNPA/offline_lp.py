#!/usr/bin/env python3
"""explore52: off-line pairs in the exact finite-N CUE analogue of the simple-on-line LP (unconditional world).
Adapted from Riemann/Deep/ExploreSchurWeyl/cue_lp.py (explore50; on-line multisets only = RH world).
A configuration has N zeros: on-line points z = w^x (w = e^{2 pi i/M}, count 1) and off-line pairs
{e^{s} w^x, e^{-s} w^x} (count 2; reflection rho -> 1 - conj rho), so p_k = sum z^k gets 2cosh(ks) w^{kx} from a pair
and p_{-k} = conj p_k (holomorphic moments = zeta23's bilinear Gabor traces; 2cosh(ks) = the huge entries of off-line zeros).
Law: invariant under rotation/reflection; constraints (Diaconis-Shahshahani values, RS range a+b <= N):
  P : E|p_k|^2 = k (1 <= k <= N)                  [bandwidth one]
  P3: P + E Re(p_a p_b conj p_{a+b}) = 0 (1 <= a <= b, a+b <= N)   [NPA level-2 words, unconditional for zeta by AF s7(e)]
Objective: min E[#simple on-line zeros]/N (and max E[sum over on-line distinct m^2]/N for the RH-world check).
usage: offline_lp.py N M s1,s2,... [P,P3]
"""
import sys, itertools, math, time
import numpy as np
from scipy.optimize import linprog

def canon(g):
    n = len(g); best = None
    for s in range(n):
        r = tuple(g[s:] + g[:s]); best = r if best is None or r < best else best
    gr = g[::-1]
    for s in range(n):
        r = tuple(gr[s:] + gr[:s]); best = r if r < best else best
    return best

def online_configs(N, M):
    seen = set(); out = []
    for rest in itertools.combinations_with_replacement(range(M), N - 1):
        x = (0,) + rest
        g = [x[i+1] - x[i] for i in range(N-1)] + [M - x[-1]]
        c = canon(g)
        if c in seen: continue
        seen.add(c); out.append(tuple(np.cumsum((0,) + c[:-1]).tolist()))
    return out

def mixed_configs(N, M, S):
    """first pair at x = 0 (rotation fixed); returns list of (online_tuple, pairs_tuple of (x, j))"""
    out = []
    items = [(x, j) for x in range(M) for j in range(S)]
    for npair in range(1, N // 2 + 1):
        non = N - 2 * npair
        for j0 in range(S):
            for rest in itertools.combinations_with_replacement(range(len(items)), npair - 1):
                prs = ((0, j0),) + tuple(items[i] for i in rest)
                for on in itertools.combinations_with_replacement(range(M), non):
                    out.append((on, prs))
    return out

def features(confs, N, M, depths):
    nC = len(confs); P = np.zeros((nC, N + 1), dtype=complex)
    ns = np.zeros(nC); m2 = np.zeros(nC); npairs = np.zeros(nC)
    ks = np.arange(N + 1)
    for i, (on, prs) in enumerate(confs):
        for x in on: P[i] += np.exp(2j * np.pi * ks * x / M)
        for (x, j) in prs: P[i] += 2 * np.cosh(ks * depths[j]) * np.exp(2j * np.pi * ks * x / M)
        if on:
            _, cnt = np.unique(on, return_counts=True)
            ns[i] = (cnt == 1).sum(); m2[i] = (cnt ** 2).sum()
        npairs[i] = len(prs)
    return P, ns / N, m2 / N, npairs

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

def constraints(P, N, sets):
    """Haar (Diaconis-Shahshahani) moments E p_lam conj p_mu = delta z_lam, |lam| = |mu| <= N (as in explore50's cue_lp.py).
    P: len(lam) = len(mu) = 1; P3: len(lam)+len(mu) <= 3; L4: <= 4; DS: all."""
    rows, targ, lab = [], [], []
    for n in range(1, N + 1):
        parts = list(partitions(n))
        pl = {lam: np.prod([P[:, k] for k in lam], axis=0) for lam in parts}
        for i, lam in enumerate(parts):
            for mu in parts[i:]:
                L = len(lam) + len(mu)
                if sets == 'P' and not (len(lam) == 1 and len(mu) == 1): continue
                if sets == 'P3' and L > 3: continue
                if sets == 'L4' and L > 4: continue
                if sets == 'Q4' and not (len(lam) == len(mu) and len(lam) <= 2): continue   # pairs + quartics, no odd words
                if sets == 'Q4e' and not (lam == mu and len(lam) <= 2): continue           # pairs + diagonal quartics E|p_a p_b|^2
                f = (pl[lam] * np.conj(pl[mu])).real
                t = zlam(lam) if lam == mu else 0.0
                s = math.sqrt(zlam(lam) * zlam(mu))
                rows.append(f / s); targ.append(t / s); lab.append((lam, mu))
    return np.array(rows), np.array(targ), lab

def solve(A, b, obj, sense):
    """column-scaled LP: u_i = w_i * c_i, c_i = max(1, max_j |A_ji|) (deep off-line pairs have huge moments)."""
    n = A.shape[1]
    Aeq = np.vstack([A, np.ones((1, n))]); beq = np.concatenate([b, [1.0]])
    cs = np.maximum(1.0, np.abs(Aeq).max(axis=0))
    c = (-obj if sense == 'max' else obj) / cs
    r = linprog(c, A_eq=Aeq / cs, b_eq=beq, bounds=(0, None), method='highs-ds')
    if r.status != 0: return None, r
    r.x = r.x / cs
    res = np.abs(Aeq @ r.x - beq).max(); r.resid = res
    return (-r.fun if sense == 'max' else r.fun), r

if __name__ == '__main__':
    N = int(sys.argv[1]); M = int(sys.argv[2]); depths = [float(t) for t in sys.argv[3].split(',')]
    sets = sys.argv[4].split(',') if len(sys.argv) > 4 else ['P', 'P3']
    t0 = time.time()
    on = [(c, ()) for c in online_configs(N, M)]
    mx = mixed_configs(N, M, len(depths))
    confs = on + mx; non = len(on)
    P, ns, m2, npairs = features(confs, N, M, depths)
    print(f"N={N} M={M} depths={depths} online={non} mixed={len(mx)} t={time.time()-t0:.1f}s", flush=True)
    for s in sets:
        A, b, lab = constraints(P, N, s)
        for world, idx in [('on ', np.arange(non)), ('off', np.arange(len(confs)))]:
            v, r = solve(A[:, idx], b, ns[idx], 'min')
            if v is None: print(f"  {s:3s} {world}: INFEAS({r.status})"); continue
            w = r.x; supp = np.where(w > 1e-9)[0]
            pairmass = float((w * npairs[idx]).sum()) * 2 / N
            msg = f"  {s:3s} {world} ncons={len(b)}: minSimpleOn={v:.6f} supp={len(supp)} offline-zero share={pairmass:.2e} resid={r.resid:.1e} minw={w.min():.1e}"
            if world == 'on ':
                v2, _ = solve(A[:, idx], b, m2[idx], 'max'); msg += f"  maxM2={v2:.6f}"
            print(msg, flush=True)
            if '-v' in sys.argv and len(supp) <= 30:
                for i in supp[np.argsort(-w[supp])]:
                    c = confs[idx[i]]; print(f"      w={w[i]:.5f} on={list(c[0])} pairs={[(x, depths[j]) for x, j in c[1]]}")
    import resource
    print(f"done t={time.time()-t0:.1f}s maxrss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20:.0f}MB")
