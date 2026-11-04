"""Tests of OpenAI's "The crossing number of complete bipartite graphs" (Sep 23 2026), exact arithmetic.
 B1 Lemma 2.1 (subspace bound): sum_{i,j} dim(U_i cap V_j) <= floor(m^2/4) N when U_i cap V_i = 0
    (exhaustive over lines/planes from a pool in Q^2 and Q^3 for small m, plus hill-climbing in Q^3, Q^4)
 B2 Lemma 4.2 (rank detector): Phi kills additive and diagonal matrices, Phi(E_pq) rank 1, Phi(S_order) invertible
 B3 drawings of K_{m,n} (n odd), straight and bent: Lemma 3.2 additive relation, rank Phi(J^ii) = s^2,
    sum_{i!=k} rank Phi(J^ik) >= 2 d_m s^2, and a hill-climb on |T|/2 (pairs with nonzero algebraic crossing)
usage: python3 kmn_check.py SEED"""
import random, sys, itertools
from fractions import Fraction as Fr
from kn_check import orient, seg_cross, rank

def dim_span(vecs):
    return rank([list(map(Fr, v)) for v in vecs]) if vecs else 0

def dim_cap(U, V):
    return dim_span(U) + dim_span(V) - dim_span(U + V)

def B1(rng):
    # exhaustive: N=2 with subspaces from {0, 4 lines, Q^2}; N=3 with {0, 7 lines, 7 planes, Q^3} over {0,1} vectors
    for N, m_max in ((2, 5), (3, 3)):
        vecs = [v for v in itertools.product((0, 1, -1), repeat=N) if any(v)]
        subs = {}
        for k in range(0, N+1):
            for comb in itertools.combinations(vecs, k):
                if dim_span(list(comb)) == k:
                    # canonical key: row-reduced basis
                    key = tuple(sorted(tuple(x) for x in comb))
                    subs.setdefault(frozenset(vv for vv in vecs if dim_span(list(comb)+[vv]) == k) if k else frozenset(), list(comb))
        S = list(subs.values())
        capd = [[dim_cap(A, B) for B in S] for A in S]
        pairs = [(a, b) for a in range(len(S)) for b in range(len(S)) if capd[a][b] == 0]
        for m in range(2, m_max+1):
            best = 0
            # maximise sum over choices of m admissible pairs (multiset, order irrelevant)
            for choice in itertools.combinations_with_replacement(range(len(pairs)), m):
                tot = sum(capd[pairs[i][0]][pairs[j][1]] for i in choice for j in choice)
                if tot > best: best = tot
            bound = (m*m//4)*N
            assert best <= bound, (N, m, best, bound)
            print(f"B1: N={N} m={m}: {len(S)} subspaces, {len(pairs)} admissible pairs, max sum={best}, bound={bound}", flush=True)
    # hill climbing in Q^3, Q^4 with spans of small integer vectors
    for N in (3, 4):
        for m in (3, 4, 5):
            best_all = 0
            for rep in range(6):
                def rsub():
                    k = rng.randint(0, N)
                    return [tuple(rng.randint(-1, 1) for _ in range(N)) for _ in range(k)]
                fam = []
                while len(fam) < m:
                    U, V = rsub(), rsub()
                    if dim_cap(U, V) == 0: fam.append((U, V))
                def score(f): return sum(dim_cap(f[i][0], f[j][1]) for i in range(m) for j in range(m))
                cur = score(fam)
                for it in range(400):
                    f2 = list(fam); i = rng.randrange(m)
                    U, V = rsub(), rsub()
                    if rng.random() < 0.5: U = list(fam[rng.randrange(m)][1]) if rng.random() < 0.7 else U
                    else: V = list(fam[rng.randrange(m)][0]) if rng.random() < 0.7 else V
                    if dim_cap(U, V): continue
                    f2[i] = (U, V); sc = score(f2)
                    if sc >= cur: fam, cur = f2, sc
                best_all = max(best_all, cur)
            assert best_all <= (m*m//4)*N
            print(f"B1: climb N={N} m={m}: best sum {best_all} <= bound {(m*m//4)*N}", flush=True)

def phi_setup(xs):
    n = len(xs); s = (n-1)//2
    Fp = []
    for p in range(n):
        d = Fr(1)
        for q in range(n):
            if q != p: d *= (xs[p]-xs[q])
        Fp.append(d)
    def h(p, q): return [xs[p]**a * xs[q]**b for a in range(s) for b in range(s)]
    B = {}
    for p in range(n):
        for q in range(n):
            hv = h(p, q); c = (xs[p]-xs[q])/(Fp[p]*Fp[q])
            B[(p, q)] = [[c*x*y for y in hv] for x in hv]
    def Phi(H):
        d = s*s
        out = [[Fr(0)]*d for _ in range(d)]
        for (p, q), M in B.items():
            if H[p][q]:
                for a in range(d):
                    for b in range(d):
                        out[a][b] += H[p][q]*M[a][b]
        return out
    return Phi, s

def B2(rng):
    for s in range(1, 5):
        n = 2*s+1
        xs = [Fr(x) for x in rng.sample(range(-40, 40), n)]
        Phi, _ = phi_setup(xs)
        r = [rng.randint(-3, 3) for _ in range(n)]; c = [rng.randint(-3, 3) for _ in range(n)]
        Hadd = [[r[p]+c[q] for q in range(n)] for p in range(n)]
        assert all(x == 0 for row in Phi(Hadd) for x in row)
        Hd = [[(p+3 if p == q else 0) for q in range(n)] for p in range(n)]
        assert all(x == 0 for row in Phi(Hd) for x in row)
        for p in range(n):
            for q in range(n):
                Epq = [[1 if (a, b) == (p, q) else 0 for b in range(n)] for a in range(n)]
                assert rank(Phi(Epq)) == (1 if p != q else 0)
        for t in range(5):
            order = list(range(n)); rng.shuffle(order); pos = {v: k for k, v in enumerate(order)}
            S = [[(Fr(1, 2) if pos[p] < pos[q] else Fr(-1, 2)) if p != q else 0 for q in range(n)] for p in range(n)]
            assert rank(Phi(S)) == s*s
    print("B2: rank detector properties (i)-(iii) ok for s<=4", flush=True)

def kmn_drawing(m, n, rng, bentprob, scale):
    R = 10**6
    A = [(rng.randint(-R, R), rng.randint(-R, R)) for _ in range(m)]
    Bv = [(rng.randint(-R, R), rng.randint(-R, R)) for _ in range(n)]
    poly = {}
    for i in range(m):
        for p in range(n):
            path = [A[i]]
            if rng.random() < bentprob:
                (x1, y1), (x2, y2) = A[i], Bv[p]; L = max(abs(x2-x1), abs(y2-y1))
                path.append(((x1+x2)//2 + int(rng.uniform(-scale, scale)*L), (y1+y2)//2 + int(rng.uniform(-scale, scale)*L)))
            path.append(Bv[p]); poly[(i, p)] = path
    return A, Bv, poly

def kmn_I(m, n, poly):
    I = {}; tot = 0
    for e, f in itertools.combinations(sorted(poly), 2):
        s = 0; cnt = 0
        for a, b in zip(poly[e], poly[e][1:]):
            for c, d in zip(poly[f], poly[f][1:]):
                x = seg_cross(a, b, c, d)
                if x: s += x; cnt += 1
        if (e[0] == f[0] or e[1] == f[1]) and cnt: raise ValueError("adjacent")
        I[(e, f)] = s; I[(f, e)] = -s; tot += cnt
    return I, tot

def kmn_J(m, n, A, Bv, poly, I, rng):
    import math
    def order_at(center, edges, first):
        cut = rng.uniform(0, 2*math.pi)
        keys = []
        for e in edges:
            path = poly[e] if first else poly[e][::-1]
            q = path[1]
            ang = math.atan2(q[1]-center[1], q[0]-center[0]) - cut
            while ang < 0: ang += 2*math.pi
            keys.append((ang, e))
        keys.sort(); return {e: k for k, (_, e) in enumerate(keys)}
    ordA = [order_at(A[i], [(i, p) for p in range(n)], True) for i in range(m)]
    ordB = [order_at(Bv[p], [(i, p) for i in range(m)], False) for p in range(n)]
    def J(e, f):
        if e == f: return Fr(0)
        v = Fr(I[(e, f)])
        if e[0] == f[0]: v += Fr(1 if ordA[e[0]][e] < ordA[e[0]][f] else -1, 2)          # both out of a_i
        if e[1] == f[1]: v += Fr(1 if ordB[e[1]][e] < ordB[e[1]][f] else -1, 2)          # both into b_p: (-1)(-1)
        return v
    return J

def B3_check(m, n, poly, J, I):
    s = (n-1)//2; d = s*s
    blocks = {(i, k): [[J((i, p), (k, q)) for q in range(n)] for p in range(n)] for i in range(m) for k in range(m)}
    for i in range(m):
        for k in range(m):
            H = [[blocks[(i, k)][p][q] - blocks[(i, 0)][p][q] - blocks[(0, k)][p][q] + blocks[(0, 0)][p][q] for q in range(n)] for p in range(n)]
            for p in range(n):
                for q in range(n):
                    assert H[p][q] - H[p][0] - H[0][q] + H[0][0] == 0, "additive relation fails"
    xs = [Fr(x) for x in random.Random(5).sample(range(-40, 40), n)]
    Phi, _ = phi_setup(xs)
    tot = 0
    for i in range(m):
        assert rank(Phi(blocks[(i, i)])) == d
        for k in range(m):
            if i != k:
                rk = rank(Phi(blocks[(i, k)]))
                nz = sum(1 for p in range(n) for q in range(n) if p != q and I[((i, p), (k, q))] != 0)
                assert rk <= nz
                tot += rk
    dm = (m//2)*((m-1)//2)
    assert tot >= 2*dm*d, ("rank sum", tot, 2*dm*d)
    return tot

def tpairs(m, n, I):
    return sum(1 for (e, f), v in I.items() if e < f and v != 0 and e[0] != f[0] and e[1] != f[1])

def B3(rng, reps=4, iters=1500):
    R = 10**6
    for (m, n) in ((3, 5), (4, 5), (5, 5), (3, 7)):
        dm = (m//2)*((m-1)//2); dn = (n//2)*((n-1)//2)
        res = []
        for rep in range(reps):
            while True:
                try:
                    A, Bv, poly = kmn_drawing(m, n, rng, 0.0, 0.0); I, tot = kmn_I(m, n, poly); break
                except ValueError: pass
            best = tpairs(m, n, I)
            for it in range(iters):
                A2, B2v = list(A), list(Bv); poly2 = {e: list(v) for e, v in poly.items()}
                r = rng.random()
                if r < 0.3:
                    side = rng.randrange(2); v = rng.randrange(m if side == 0 else n)
                    if side == 0:
                        A2[v] = (A[v][0]+rng.randint(-R//5, R//5), A[v][1]+rng.randint(-R//5, R//5))
                        for p in range(n): poly2[(v, p)][0] = A2[v]
                    else:
                        B2v[v] = (Bv[v][0]+rng.randint(-R//5, R//5), Bv[v][1]+rng.randint(-R//5, R//5))
                        for i in range(m): poly2[(i, v)][-1] = B2v[v]
                else:
                    e = rng.choice(sorted(poly2)); (x1, y1), (x2, y2) = poly2[e][0], poly2[e][-1]
                    L = max(abs(x2-x1), abs(y2-y1))
                    poly2[e] = [poly2[e][0], poly2[e][-1]] if rng.random() < 0.2 else [poly2[e][0], ((x1+x2)//2+int(rng.uniform(-.8, .8)*L), (y1+y2)//2+int(rng.uniform(-.8, .8)*L)), poly2[e][-1]]
                try:
                    I2, tot2 = kmn_I(m, n, poly2)
                except ValueError:
                    continue
                c2 = tpairs(m, n, I2)
                if c2 <= best: A, Bv, poly, I, best = A2, B2v, poly2, I2, c2
            J = kmn_J(m, n, A, Bv, poly, I, rng)
            rs = B3_check(m, n, poly, J, I)
            res.append(best)
        print(f"B3: K_({m},{n}): additive relation, detector ranks, rank sum ok; hill-climb min |T|/2 = {min(res)} (d_m d_n = {dm*dn}); runs {res}", flush=True)
        assert min(res) >= dm*dn

if __name__ == "__main__":
    rng = random.Random(int(sys.argv[1]))
    B2(rng); B1(rng); B3(rng)
