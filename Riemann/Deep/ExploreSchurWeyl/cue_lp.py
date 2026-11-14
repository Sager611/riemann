#!/usr/bin/env python3
"""Finite-N CUE analogue of the simple-zero LP (explore50).
Random N-multisets on Z_M (points z = w^x, w = e^{2 pi i/M}), rotation/reflection-invariant law.
Constraints: Haar (Diaconis-Shahshahani) moments E p_lam conj(p_mu) = delta z_lam for chosen (lam, mu), |lam| = |mu| <= N.
Objectives: max E[sum_rho m_rho^2]/N ; min E[#simple]/N ; min E[#distinct]/N.
Constraint sets: P (pair: lam = mu = (k)), P3 (+ triples (a,b)|(a+b)), L4 (all with len(lam)+len(mu) <= 4), DS (all).
usage: script.py N M [sets] """
import sys, itertools, math, time
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
        a = lam.count(j); z *= j**a * math.factorial(a)
    return z

def canon(g):
    # gaps g (cyclic, sum M); dihedral-canonical tuple
    n = len(g); best = None
    for s in range(n):
        r = tuple(g[s:] + g[:s]); best = r if best is None or r < best else best
    gr = g[::-1]
    for s in range(n):
        r = tuple(gr[s:] + gr[:s]); best = r if r < best else best
    return best

def configs(N, M):
    seen = set(); out = []
    for rest in itertools.combinations_with_replacement(range(M), N - 1):
        x = (0,) + rest
        g = [x[i+1] - x[i] for i in range(N-1)] + [M - x[-1]]
        c = canon(g)
        if c in seen: continue
        seen.add(c)
        pts = np.cumsum((0,) + c[:-1])
        out.append(pts)
    return np.array(out, dtype=np.int64)

def features(C, N, M, sets):
    w = np.exp(2j*np.pi/M)
    P = np.stack([np.exp(2j*np.pi*k*C/M).sum(axis=1) for k in range(N+1)], axis=1)  # P[:,k] = p_k
    rows, targ, lab = [], [], []
    for n in range(1, N+1):
        parts = list(partitions(n))
        pl = {lam: np.prod([P[:, k] for k in lam], axis=0) for lam in parts}
        for i, lam in enumerate(parts):
            for mu in parts[i:]:
                L = len(lam) + len(mu)
                if sets == 'P' and not (len(lam) == 1 and len(mu) == 1): continue
                if sets == 'P3' and L > 3: continue
                if sets == 'L4' and L > 4: continue
                f = (pl[lam] * np.conj(pl[mu])).real
                t = zlam(lam) if lam == mu else 0.0
                s = math.sqrt(zlam(lam) * zlam(mu))
                rows.append(f / s); targ.append(t / s); lab.append((lam, mu))
    return np.array(rows), np.array(targ), lab

def mult_stats(C, N):
    m2 = np.zeros(len(C)); ns = np.zeros(len(C)); nd = np.zeros(len(C))
    for i, c in enumerate(C):
        _, cnt = np.unique(c, return_counts=True)
        m2[i] = (cnt**2).sum(); ns[i] = (cnt == 1).sum(); nd[i] = len(cnt)
    return m2 / N, ns / N, nd / N

def solve(A, b, obj, sense):
    n = A.shape[1]
    Aeq = np.vstack([A, np.ones((1, n))]); beq = np.concatenate([b, [1.0]])
    c = -obj if sense == 'max' else obj
    r = linprog(c, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
    if r.status != 0: return None, r
    return (-r.fun if sense == 'max' else r.fun), r

if __name__ == '__main__':
    N = int(sys.argv[1]); M = int(sys.argv[2])
    sets = sys.argv[3].split(',') if len(sys.argv) > 3 else ['P', 'P3', 'L4', 'DS']
    t0 = time.time()
    C = configs(N, M)
    m2, ns, nd = mult_stats(C, N)
    print(f"N={N} M={M} configs={len(C)} t={time.time()-t0:.1f}s", flush=True)
    for s in sets:
        A, b, lab = features(C, N, M, s)
        res = []
        objs = [('maxM2', m2, 'max')] if '-m' in sys.argv else [('maxM2', m2, 'max'), ('minSimple', ns, 'min'), ('minDistinct', nd, 'min')]
        for name, obj, sense in objs:
            v, r = solve(A, b, obj, sense)
            if v is None: res.append(f"{name}=INFEAS({r.status})"); continue
            wsol = r.x; supp = np.where(wsol > 1e-9)[0]
            res.append(f"{name}={v:.6f}[supp {len(supp)}]")
            if name == 'maxM2' and len(supp) <= 40 and '-v' in sys.argv:
                for i in supp[np.argsort(-wsol[supp])]:
                    print(f"      w={wsol[i]:.5f} pts={C[i].tolist()} m2={m2[i]*N:.0f}")
        print(f"  {s:3s} ncons={len(b)}: " + "  ".join(res), flush=True)
    print(f"done t={time.time()-t0:.1f}s")
