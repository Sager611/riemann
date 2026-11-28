"""Shared CEGAR for 'is there a compatible S (all pairs at squared distance <= s, containing vertex 0) whose graph
H_s[S] (pairs at squared distance exactly s) is not 5-colourable?'.  Inputs: m, edge arrays U, V of H_s on L, and a
boolean matrix BAD (pairs at squared distance > s).  Blocking clauses come from near-colourings of all of L found by
tabu search with the current S frozen: if S avoids every conflict edge of such a colouring, the colouring is proper
on S.  Returns ('UNSAT', iters) or ('CANDIDATE', S) or ('TIMEOUT', iters)."""
import time
import numpy as np
from pysat.solvers import Solver
from tabu import colour_tabu

def prune(m, U, V, BAD, root=0):
    """drop vertices that cannot lie in a 6-critical compatible S: fewer than 5 compatible H_s-neighbours alive"""
    alive = np.ones(m, bool); adj = [[] for _ in range(m)]
    for a, b in zip(U.tolist(), V.tolist()):
        if not BAD[a, b]: adj[a].append(b); adj[b].append(a)
    changed = True
    while changed:
        changed = False
        for v in range(m):
            if alive[v] and sum(1 for w in adj[v] if alive[w]) < 5: alive[v] = False; changed = True
    return alive

def cegar(m, U, V, BAD, tl=1800, root=0):
    alive = prune(m, U, V, BAD)
    if not alive[root]: return ('UNSAT', 0, 'root pruned')
    idx = np.flatnonzero(alive); loc = -np.ones(m, int); loc[idx] = np.arange(len(idx))
    keep = alive[U] & alive[V]; U, V = loc[U[keep]].astype(np.int32), loc[V[keep]].astype(np.int32)
    BAD = BAD[np.ix_(idx, idx)]; m = len(idx); r0 = int(loc[root])
    A = Solver(name='g4'); A.add_clause([r0 + 1])
    bi, bj = np.nonzero(np.triu(BAD, 1))
    for a, b in zip(bi.tolist(), bj.tolist()): A.add_clause([-(a + 1), -(b + 1)])
    nxt = m + 1; t0 = time.time(); it = 0; col = None
    def block(c):
        nonlocal nxt
        mono = np.flatnonzero(c[U] == c[V]); cl = []
        for e in mono:
            y = nxt; nxt += 1; A.add_clause([-y, int(U[e]) + 1]); A.add_clause([-y, int(V[e]) + 1]); cl.append(y)
        A.add_clause(cl); return len(mono)
    c0, _ = colour_tabu(m, U, V, 5, iters=300000, seed=1); block(c0)
    while True:
        it += 1
        if time.time() - t0 > tl: return ('TIMEOUT', it, m)
        if not A.solve(): return ('UNSAT', it, m)
        mod = A.get_model(); inS = np.array([mod[v] > 0 for v in range(m)])
        compat = ~BAD
        for v in np.random.default_rng(it).permutation(m):
            if not inS[v] and compat[v, inS].all(): inS[v] = True
        S = np.flatnonzero(inS); sl = -np.ones(m, int); sl[S] = np.arange(len(S))
        se = inS[U] & inS[V]; Us, Vs = sl[U[se]].astype(np.int32), sl[V[se]].astype(np.int32)
        cs = None
        for sd in range(1, 5):
            c_, conf = colour_tabu(len(S), Us, Vs, 5, iters=300000 * sd, seed=sd)
            if conf == 0: cs = c_; break
        if cs is None: return ('CANDIDATE', S, m)
        init = np.random.default_rng(it).integers(0, 5, m); init[S] = cs
        c_, conf = colour_tabu(m, U, V, 5, iters=200000, seed=it, init=init, frozen=inS.astype(np.int32))
        assert (c_[S] == cs).all()
        block(c_)


def cegar_lens(m, U, V, BAD, tl=1800, root=0, reps=None):
    """Split by the first neighbour: a 6-critical H inside S may be moved to contain root, which then has a
    neighbour y in H; so S lies in the lens {v compatible with root and y}.  reps: candidate y's (default: all
    neighbours of root).  Returns ('UNSAT', number of lenses) or the first non-UNSAT result."""
    nb = sorted(set(V[U == root].tolist()) | set(U[V == root].tolist()))
    ys = nb if reps is None else [y for y in reps if y in set(nb)]
    t0 = time.time(); done = 0
    for y in ys:
        keep = ~BAD[root] & ~BAD[y]; keep[root] = keep[y] = True
        idx = np.flatnonzero(keep); loc = -np.ones(m, int); loc[idx] = np.arange(len(idx))
        e = keep[U] & keep[V]; Ul, Vl = loc[U[e]].astype(np.int32), loc[V[e]].astype(np.int32)
        res = cegar_forced(len(idx), Ul, Vl, BAD[np.ix_(idx, idx)], int(loc[root]), int(loc[y]), tl=max(60, tl - (time.time() - t0)))
        if res[0] != 'UNSAT': return (res[0], res[1] if res[0] != 'CANDIDATE' else idx[res[1]], y)
        done += 1
    return ('UNSAT', done, m)

def cegar_forced(m, U, V, BAD, r0, y0, tl=600):
    alive = prune(m, U, V, BAD)
    if not (alive[r0] and alive[y0]): return ('UNSAT', 0, m)
    idx = np.flatnonzero(alive); loc = -np.ones(m, int); loc[idx] = np.arange(len(idx))
    keep = alive[U] & alive[V]; U, V = loc[U[keep]].astype(np.int32), loc[V[keep]].astype(np.int32)
    BAD = BAD[np.ix_(idx, idx)]; m = len(idx); r0, y0 = int(loc[r0]), int(loc[y0])
    A = Solver(name='g4'); A.add_clause([r0 + 1]); A.add_clause([y0 + 1])
    bi, bj = np.nonzero(np.triu(BAD, 1))
    for a, b in zip(bi.tolist(), bj.tolist()): A.add_clause([-(a + 1), -(b + 1)])
    nxt = m + 1; t0 = time.time(); it = 0
    def block(c):
        nonlocal nxt
        cl = []
        for e in np.flatnonzero(c[U] == c[V]):
            y = nxt; nxt += 1; A.add_clause([-y, int(U[e]) + 1]); A.add_clause([-y, int(V[e]) + 1]); cl.append(y)
        A.add_clause(cl)
    c0, _ = colour_tabu(m, U, V, 5, iters=300000, seed=1); block(c0)
    compat = ~BAD
    while True:
        it += 1
        if time.time() - t0 > tl: return ('TIMEOUT', it, m)
        if not A.solve(): return ('UNSAT', it, m)
        mod = A.get_model(); inS = np.array([mod[v] > 0 for v in range(m)])
        for v in np.random.default_rng(it).permutation(m):
            if not inS[v] and compat[v, inS].all(): inS[v] = True
        S = np.flatnonzero(inS); sl = -np.ones(m, int); sl[S] = np.arange(len(S))
        se = inS[U] & inS[V]; Us, Vs = sl[U[se]].astype(np.int32), sl[V[se]].astype(np.int32)
        cs = None
        for sd in range(1, 5):
            c_, conf = colour_tabu(len(S), Us, Vs, 5, iters=300000 * sd, seed=sd)
            if conf == 0: cs = c_; break
        if cs is None: return ('CANDIDATE', idx[S], m)
        init = np.random.default_rng(it).integers(0, 5, m); init[S] = cs
        c_, conf = colour_tabu(m, U, V, 5, iters=200000, seed=it, init=init, frozen=inS.astype(np.int32))
        assert (c_[S] == cs).all()
        block(c_)
