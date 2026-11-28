#!/usr/bin/env python3
"""D1: subsets of uniform 4-polytope vertex sets and root systems (Borsuk in R^4).

For a finite set P and a squared distance s of P let H_s(P) = graph on P joining pairs at squared
distance exactly s.  If S is a subset of P with diam(S)^2 = s, its diameter graph is the induced subgraph
H_s(P)[S].  So if every H_s(P) is 5-colourable, no subset of P is a counterexample.  This script checks
that for every set and every s: max degree <= 4 (greedy), empty 5-core (greedy), or an explicit SAT
5-colouring checked edge by edge.  Thresholds where H_s(P) is NOT 5-colourable are reported (they need
the refined search over subsets S that avoid pairs at squared distance > s).
All arithmetic is exact: coordinates in Z[w] (w = phi or sqrt 2) scaled to integers, int64 with an
overflow guard.  Usage: python d1.py [group ...]   (groups: A4 B4 F4 H4 ROOTS)
"""
import sys, itertools, json, time
from fractions import Fraction as Fr
from math import lcm
import numpy as np
from pysat.solvers import Solver

FIELDS = {'Q': (0, 0), 'Q2': (0, 2), 'Q5': (1, 1)}   # w^2 = M1 w + M0
M1, M0, FIELD = 0, 0, 'Q'
WF = {'Q': 0.0, 'Q2': 2 ** 0.5, 'Q5': (1 + 5 ** 0.5) / 2}

def set_field(f):
    global M1, M0, FIELD
    M1, M0 = FIELDS[f]; FIELD = f

class K:
    __slots__ = ('a', 'b')
    def __init__(s, a, b=0):
        s.a = a if isinstance(a, Fr) else Fr(a); s.b = b if isinstance(b, Fr) else Fr(b)
    @staticmethod
    def c(o): return o if isinstance(o, K) else K(o)
    def __add__(s, o): o = K.c(o); return K(s.a + o.a, s.b + o.b)
    __radd__ = __add__
    def __sub__(s, o): o = K.c(o); return K(s.a - o.a, s.b - o.b)
    def __rsub__(s, o): return K.c(o) - s
    def __neg__(s): return K(-s.a, -s.b)
    def __mul__(s, o):
        o = K.c(o); return K(s.a * o.a + M0 * s.b * o.b, s.a * o.b + s.b * o.a + M1 * s.b * o.b)
    __rmul__ = __mul__
    def conj(s): return K(s.a + M1 * s.b, -s.b)
    def inv(s):
        c = s.conj(); n = s * c; assert n.b == 0 and n.a != 0; return K(c.a / n.a, c.b / n.a)
    def __truediv__(s, o): return s * K.c(o).inv()
    def sign(s): return sgn(s.a, s.b)
    def __eq__(s, o): o = K.c(o); return s.a == o.a and s.b == o.b
    def __hash__(s): return hash((s.a, s.b))

def sgn(A, B):
    """exact sign of A + B w"""
    if FIELD == 'Q5': p, q, r = 2 * A + B, B, 5
    elif FIELD == 'Q2': p, q, r = A, B, 2
    else: p, q, r = A, 0, 1
    sp = (p > 0) - (p < 0); sq = (q > 0) - (q < 0)
    if sq == 0 or sp == sq: return sp
    if sp == 0: return sq
    return sp if p * p - r * q * q > 0 else sq

def ip(x, y):
    t = K(0)
    for u, v in zip(x, y): t = t + u * v
    return t

def refl(x, r, rr):
    c = (ip(x, r) * 2) / rr
    return tuple(u - c * v for u, v in zip(x, r))

def orbit(seed, roots):
    rrs = [ip(r, r) for r in roots]
    seen = {seed}; todo = [seed]
    while todo:
        x = todo.pop()
        for r, rr in zip(roots, rrs):
            y = refl(x, r, rr)
            if y not in seen: seen.add(y); todo.append(y)
    return list(seen)

def vec(*c): return tuple(K.c(x) for x in c)
def add(x, y): return tuple(u + v for u, v in zip(x, y))
def scal(c, x): return tuple(K.c(c) * u for u in x)

def solve_dual(alpha):
    """omega_i with <omega_i, alpha_j> = delta_ij (Gauss-Jordan over K)."""
    n = len(alpha); A = [[alpha[j][k] for k in range(n)] for j in range(n)]   # rows alpha_j
    # want W with W A^T = I  ->  solve A X = I, W = X^T
    M = [A[i][:] + [K(1 if i == j else 0) for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = next(r for r in range(col, n) if not M[r][col] == K(0))
        M[col], M[piv] = M[piv], M[col]
        iv = M[col][col].inv(); M[col] = [x * iv for x in M[col]]
        for r in range(n):
            if r != col and not M[r][col] == K(0):
                f = M[r][col]; M[r] = [x - f * y for x, y in zip(M[r], M[col])]
    X = [row[n:] for row in M]            # X = A^{-1}
    return [tuple(X[k][i] for k in range(n)) for i in range(n)]   # columns of A^{-1}

def to_int(pts):
    L = 1
    for p in pts:
        for u in p: L = lcm(L, u.a.denominator, u.b.denominator)
    X = np.array([[[int(u.a * L), int(u.b * L)] for u in p] for p in pts], dtype=np.int64)
    return X

def ringings(omega, c):
    out = []
    for m in range(1, 16):
        R = [i for i in range(4) if m >> i & 1]
        x = vec(*([0] * len(omega[0])))
        for i in R: x = add(x, scal(c[i], omega[i]))
        out.append((''.join(str(i) for i in R), x))
    return out

def icosians():
    h = Fr(1, 2); P = []
    for i in range(4):
        for s in (1, -1):
            v = [0] * 4; v[i] = s; P.append(vec(*v))
    for sg in itertools.product((h, -h), repeat=4): P.append(vec(*sg))
    base = [K(0, h), K(h), K(-h, h), K(0)]          # phi/2, 1/2, 1/(2 phi), 0
    ev = [p for p in itertools.permutations(range(4))
          if sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j]) % 2 == 0]
    for p in ev:
        for sg in itertools.product((1, -1), repeat=3):
            v = [None] * 4
            vals = [base[0] * sg[0], base[1] * sg[1], base[2] * sg[2], base[3]]
            for k in range(4): v[p[k]] = vals[k]
            P.append(tuple(v))
    assert len(set(P)) == 120
    return P

def sets_for(group):
    """yield (name, field, points as tuples of K, vertex_transitive)"""
    if group == 'A4':
        set_field('Q')
        for m in range(1, 16):
            R = [i for i in range(4) if m >> i & 1]
            x = [sum(1 for i in R if k <= i) for k in range(5)]   # sum of omega_i, omega_i = (1^(i+1), 0..)
            pts = sorted(set(itertools.permutations(x)))
            yield ('A4:' + ''.join(map(str, R)), 'Q', [vec(*p) for p in pts], True)
    elif group == 'B4':
        set_field('Q2'); r2 = K(0, 1); h = Fr(1, 2)
        om = [vec(1, 0, 0, 0), vec(1, 1, 0, 0), vec(1, 1, 1, 0), vec(h, h, h, h)]
        al = [vec(1, -1, 0, 0), vec(0, 1, -1, 0), vec(0, 0, 1, -1), vec(0, 0, 0, 1)]
        for nm, x in ringings(om, [K(1), K(1), K(1), r2]):
            yield ('B4:' + nm, 'Q2', orbit(x, al), True)
    elif group == 'F4':
        set_field('Q2'); r2 = K(0, 1); h = Fr(1, 2)
        al = [vec(0, 1, -1, 0), vec(0, 0, 1, -1), vec(0, 0, 0, 1), vec(h, -h, -h, -h)]
        om = [vec(1, 1, 0, 0), vec(2, 1, 1, 0), vec(3 * h, h, h, h), vec(1, 0, 0, 0)]
        for i in range(4):
            for j in range(4):
                assert ip(om[i], al[j]) * 2 / ip(al[j], al[j]) == K(1 if i == j else 0)
        for nm, x in ringings(om, [K(1), K(1), r2, r2]):
            yield ('F4:' + nm, 'Q2', orbit(x, al), True)
    elif group == 'H4':
        set_field('Q5'); roots = icosians()
        g = vec(10, 3, Fr(17, 10), Fr(1, 5))
        pos = [r for r in roots if ip(r, g).sign() > 0]; assert len(pos) == 60
        simple = []
        for a in pos:
            ps = set(pos)
            if all(refl(b, a, K(1)) in ps for b in pos if b != a): simple.append(a)
        assert len(simple) == 4, len(simple)
        om = solve_dual(simple)
        for nm, x in ringings(om, [K(1)] * 4):
            yield ('H4:' + nm, 'Q5', orbit(x, simple), True)
    elif group == 'ROOTS':
        set_field('Q'); h = Fr(1, 2)
        long_ = [vec(*v) for v in set(itertools.permutations((1, 1, 0, 0))) for _ in [0]]
        D4 = set()
        for v in set(itertools.permutations((1, 1, 0, 0))):
            for sg in itertools.product((1, -1), repeat=4):
                D4.add(tuple(a * s for a, s in zip(v, sg)))
        E = set()
        for i in range(4):
            for s in (1, -1):
                v = [0] * 4; v[i] = s; E.add(tuple(v))
        Hc = set(itertools.product((h, -h), repeat=4))
        E2 = {tuple(2 * a for a in v) for v in E}
        yield ('F4roots', 'Q', [vec(*v) for v in D4 | E | Hc], False)
        yield ('B4roots', 'Q', [vec(*v) for v in D4 | E], False)
        yield ('C4roots', 'Q', [vec(*v) for v in D4 | E2], False)
        yield ('24cell+dual', 'Q', [vec(*v) for v in D4 | {tuple(2 * a for a in v) for v in E | Hc}], False)

# ---------------- graph checks ----------------
OFF, SH = 1 << 30, 1 << 31

def keys_of(DA, DB): return (DA + OFF) * SH + (DB + OFF)

def unkey(k): return int(k // SH) - OFF, int(k % SH) - OFF

def sqd_block(X, i0, i1):
    a = X[:, :, 0]; b = X[:, :, 1]
    IA = a[i0:i1] @ a.T + M0 * (b[i0:i1] @ b.T)
    IB = a[i0:i1] @ b.T + b[i0:i1] @ a.T + M1 * (b[i0:i1] @ b.T)
    NA = (a * a).sum(1) + M0 * (b * b).sum(1); NB = 2 * (a * b).sum(1) + M1 * (b * b).sum(1)
    DA = NA[i0:i1, None] + NA[None, :] - 2 * IA; DB = NB[i0:i1, None] + NB[None, :] - 2 * IB
    return DA, DB

def csr(n, U, V):
    W = np.concatenate([U, V]); Z = np.concatenate([V, U])
    o = np.argsort(W, kind='stable'); W = W[o]; Z = Z[o]
    ptr = np.searchsorted(W, np.arange(n + 1))
    return ptr, Z

def smallest_last(n, ptr, Z):
    """degeneracy order (bucket queue); returns (order to colour, degeneracy)"""
    deg = (ptr[1:] - ptr[:-1]).tolist(); Zl = Z.tolist(); P = ptr.tolist()
    maxd = max(deg) if n else 0
    buckets = [set() for _ in range(maxd + 1)]
    for v in range(n): buckets[deg[v]].add(v)
    alive = [True] * n; order = []; dg = 0; cur = 0
    for _ in range(n):
        cur = max(cur - 1, 0)
        while not buckets[cur]: cur += 1
        v = buckets[cur].pop(); dg = max(dg, cur); alive[v] = False; order.append(v)
        for w in Zl[P[v]:P[v + 1]]:
            if alive[w]:
                buckets[deg[w]].discard(w); deg[w] -= 1; buckets[deg[w]].add(w)
    return order[::-1], dg

def colour_check(n, U, V, k=5, budget=300_000):
    """('deg'|'greedy'|'sat', colouring) or ('unsat'|'unknown', None)"""
    ptr, Z = csr(n, U, V); Zl = Z.tolist(); P = ptr.tolist()
    order, dg = smallest_last(n, ptr, Z)
    col = [-1] * n
    for v in order:
        used = {col[w] for w in Zl[P[v]:P[v + 1]]}
        col[v] = next(c for c in range(n + 1) if c not in used)
    if max(col) < k: return ('deg' if dg < k else 'greedy'), col
    s = Solver(name='g4')
    for v in range(n): s.add_clause([v * k + c + 1 for c in range(k)])
    for u, v in zip(U.tolist(), V.tolist()):
        for c in range(k): s.add_clause([-(u * k + c + 1), -(v * k + c + 1)])
    if len(U): s.add_clause([int(U[0]) * k + 1])
    s.conf_budget(budget); r = s.solve_limited()
    if r is None: s.delete(); return 'unknown', None
    if not r: s.delete(); return 'unsat', None
    m = s.get_model(); s.delete()
    col = [next(c for c in range(k) if m[v * k + c] > 0) for v in range(n)]
    return 'sat', col

def verify(U, V, col, k=5):
    c = np.array(col); return bool(((c >= 0) & (c < k)).all() and (c[U] != c[V]).all())

def check_set(name, field, pts, vt, out):
    set_field(field); t0 = time.time()
    X = to_int(pts); n = len(X)
    mx = int(np.abs(X).max()); assert 16 * X.shape[1] * mx * mx < (1 << 29), ('overflow guard', name, mx)
    if vt:
        DA, DB = sqd_block(X, 0, 1); kk = keys_of(DA[0], DB[0]); kk = kk[np.arange(n) != 0]
        vals, cnt = np.unique(kk, return_counts=True); degmax = dict(zip(vals.tolist(), cnt.tolist()))
    else:
        degmax = {}
        for i0 in range(0, n, 256):
            DA, DB = sqd_block(X, i0, min(n, i0 + 256)); kk = keys_of(DA, DB)
            for r in range(kk.shape[0]):
                v, c = np.unique(np.delete(kk[r], i0 + r), return_counts=True)
                for a, b in zip(v.tolist(), c.tolist()): degmax[a] = max(degmax.get(a, 0), b)
    wf = WF[field]
    fv = lambda key: (lambda ab: ab[0] + ab[1] * wf)(unkey(key))
    allv = sorted(degmax, key=fv); smax = allv[-1]
    need = [v for v in allv if degmax[v] >= 5]
    res = {'set': name, 'n': n, 'ndist': len(allv), 'maxdeg<5': len(allv) - len(need)}
    stats = {}; bad = []
    B = max(8, int(3e6 / n))
    def run(batch):
        T = np.array(sorted(batch), dtype=np.int64); E = {v: ([], []) for v in batch}
        for i0 in range(0, n, 128):
            i1 = min(n, i0 + 128); DA, DB = sqd_block(X, i0, i1); kk = keys_of(DA, DB)
            ii, jj = np.nonzero(np.isin(kk, T)); kv = kk[ii, jj]; ii = ii + i0; m = jj > ii
            ii, jj, kv = ii[m], jj[m], kv[m]
            o = np.argsort(kv, kind='stable'); ii, jj, kv = ii[o], jj[o], kv[o]
            cuts = np.flatnonzero(np.diff(kv)) + 1
            for a, b, c in zip(np.split(ii, cuts), np.split(jj, cuts), np.split(kv, cuts)):
                if len(c): E[int(c[0])][0].append(a.astype(np.int32)); E[int(c[0])][1].append(b.astype(np.int32))
        for v in batch:
            U = np.concatenate(E[v][0]) if E[v][0] else np.zeros(0, np.int32)
            V = np.concatenate(E[v][1]) if E[v][1] else np.zeros(0, np.int32)
            del E[v]
            tt = time.time(); st, col = colour_check(n, U, V)
            print(f'  {name} s/smax={fv(v) / fv(smax):.5f} deg={degmax[v]} m={len(U)} {st} {time.time() - tt:.1f}s', file=sys.stderr, flush=True)
            if col is not None: assert verify(U, V, col), ('bad colouring', name, v)
            else: bad.append({'s/smax': round(fv(v) / fv(smax), 6), 'key': int(v), 'st': st, 'deg': degmax[v]})
            stats[st] = stats.get(st, 0) + 1
    batch = []; tot = 0
    for v in need:
        batch.append(v); tot += degmax[v]
        if tot >= B: run(batch); batch = []; tot = 0
    if batch: run(batch)
    res.update(stats); res['bad'] = bad; res['sec'] = round(time.time() - t0, 1)
    print(json.dumps(res), flush=True); out.write(json.dumps(res) + '\n'); out.flush()

if __name__ == '__main__':
    groups = sys.argv[1:] or ['ROOTS', 'A4', 'B4', 'F4', 'H4']
    with open('Borsuk/.scratch/d1_results.jsonl', 'a') as out:
        for g in groups:
            for name, field, pts, vt in sets_for(g):
                print(f'{name}: n={len(pts)}', file=sys.stderr, flush=True)
                check_set(name, field, pts, vt, out)
