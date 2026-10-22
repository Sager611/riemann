#!/usr/bin/env python3
"""verify.py -- computer checks for NOTE.md (evolution algebras A(G) and A_RW(G)).

Default run: about 25 s (python3 verify.py; hard alarm 60 s).  python3 verify.py --full adds the slow checks (6-8 min).
Needs: python3, sympy, networkx.  Memory < 150 MB.

[N] Theorem N over finite fields: every natural basis of A(G) and of A_RW(G) is enumerated by brute force
    (up to order and scaling); the number equals the product over twin classes predicted by Theorem N.
[P] Theorem P (char p): A(G) ~= A_RW(G) decided by brute force (via all natural bases of A_RW(G)) equals the
    criterion of Theorem P; K_{1,1,3} and the tree T6 are not isomorphic over F_3 but are over F_9.
[G] Full GL_n Groebner systems over Q (no theory used): P3 forces p00^3 = 4, K_{1,3} forces p00^3 = 9,
    the diamond, P4 and the paw have no solution over C.
[L] Lemma 2: k(gamma^2 - alpha) = 1^T(P^T P - alpha I)1 - (P1 - gamma 1)^T(P1 + gamma 1) for k <= 6, and
    for k = 2 the elimination ideal of {P^T P = alpha I, P1 = gamma 1} is exactly (alpha - gamma^2).
[A] Theorem A: for all 995 connected graphs with 2..7 vertices (networkx graph atlas) the twin-class system
    (Theorem I + Lemma 2) is solved over C, R and Q with sympy and compared with "regular, or biregular with
    d1/d2 a cube"; every predicted isomorphism is built and checked exactly. Extra: K_{1,8} (iso over Q),
    the subdivided K_4 and the CRR (2,3)-biregular graph (iso over R, not over Q), the Petersen graph.
[H] Remarks: a non-zero non-injective homomorphism A_RW(G) -> A(G) over C for the diamond and the double star;
    BCS22 Cor. 2.7 needs |k| >= 4 (unique natural basis without (2LI) over F_2 and F_3); BCS22 Thm 2.4(i) "<="
    and Cor. 2.17 fail in characteristic 2 (F_2 examples); BCS22 Thm 2.11 "r = dim A^2" fails for A(P_5).
"""
import sys, time, signal, itertools
from fractions import Fraction
import networkx as nx
from networkx.generators.atlas import graph_atlas_g
import sympy as sp

FULL = '--full' in sys.argv
signal.alarm(900 if FULL else 60)          # hard timeout
T0 = time.time()
FAILS = []

def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print('  FAIL:', msg)

# ----------------------------------------------------------------------------------------------
# small finite fields: F_p (p prime) and F_9 = F_3[i]/(i^2 + 1); elements are 0..q-1
class GF:
    def __init__(self, q):
        self.q = q
        if q in (2, 3, 5, 7, 11, 13):
            self.p = q
            self.add = [[(a + b) % q for b in range(q)] for a in range(q)]
            self.mul = [[(a * b) % q for b in range(q)] for a in range(q)]
        elif q == 9:                                   # x = a + 3b  <->  a + b i
            self.p = 3
            enc = lambda a, b: (a % 3) + 3 * (b % 3)
            self.add = [[enc(x % 3 + y % 3, x // 3 + y // 3) for y in range(9)] for x in range(9)]
            self.mul = [[enc((x % 3) * (y % 3) - (x // 3) * (y // 3), (x % 3) * (y // 3) + (x // 3) * (y % 3))
                         for y in range(9)] for x in range(9)]
        else:
            raise ValueError(q)
        self.neg = [next(y for y in range(q) if self.add[x][y] == 0) for x in range(q)]
        self.inv = [None] + [next(y for y in range(q) if self.mul[x][y] == 1) for x in range(1, q)]
        self.I = 3 if q == 9 else None                 # sqrt(-1) in F_9
    def n(self, k):                                    # image of the integer k
        r = 0
        for _ in range(k % self.p): r = self.add[r][1]
        return r
    def sub(self, a, b): return self.add[a][self.neg[b]]

def kernel(F, rows, n):
    R = [list(r) for r in rows]; piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, len(R)) if R[i][c]), None)
        if pr is None: continue
        R[r], R[pr] = R[pr], R[r]
        iv = F.inv[R[r][c]]; R[r] = [F.mul[iv][x] for x in R[r]]
        for i in range(len(R)):
            if i != r and R[i][c]:
                f = R[i][c]; R[i] = [F.sub(R[i][j], F.mul[f][R[r][j]]) for j in range(n)]
        piv.append(c); r += 1
        if r == len(R): break
    basis = []
    for fc in [c for c in range(n) if c not in piv]:
        v = [0] * n; v[fc] = 1
        for i, pc in enumerate(piv): v[pc] = F.neg[R[i][fc]]
        basis.append(tuple(v))
    return basis

def span(F, basis, n):
    out = []
    for co in itertools.product(range(F.q), repeat=len(basis)):
        v = [0] * n
        for c, b in zip(co, basis):
            if c:
                for j in range(n): v[j] = F.add[v[j]][F.mul[c][b[j]]]
        out.append(tuple(v))
    return out

def normalized(v):
    for x in v:
        if x: return x == 1
    return False

def prod(F, M, f, g):
    """f*g in an evolution algebra with structure matrix M (rows = squares of the basis)."""
    m = len(M[0]); out = [0] * m
    for w in range(len(f)):
        c = F.mul[f[w]][g[w]]
        if c:
            for j in range(m): out[j] = F.add[out[j]][F.mul[c][M[w][j]]]
    return tuple(out)

def rank(F, vecs, n):
    return len(vecs) - len(kernel(F, [list(c) for c in zip(*vecs)], len(vecs))) if vecs else 0

def natural_bases(F, M):
    """all natural bases up to order and scaling (no theory used)."""
    n = len(M); order = {}
    res = []
    def rec(chosen, rows):
        if len(chosen) == n:
            res.append(tuple(chosen)); return
        K = kernel(F, rows, n) if rows else [tuple(int(i == j) for j in range(n)) for i in range(n)]
        if len(K) < n - len(chosen): return
        for v in span(F, K, n):
            if not normalized(v) or (chosen and v <= chosen[-1]): continue
            if not any(prod(F, M, v, v)): continue        # f^2 = 0 is impossible in a natural basis (ann = 0)
            if rank(F, chosen + [v], n) < len(chosen) + 1: continue
            rec(chosen + [v], rows + [[F.mul[v[w]][M[w][j]] for w in range(n)] for j in range(len(M[0]))])
    rec([], [])
    return res

def coords(F, B, y):
    n = len(B)
    A = [[B[a][j] for a in range(n)] + [y[j]] for j in range(n)]
    for c in range(n):
        pr = next(i for i in range(c, n) if A[i][c])
        A[c], A[pr] = A[pr], A[c]
        iv = F.inv[A[c][c]]; A[c] = [F.mul[iv][x] for x in A[c]]
        for i in range(n):
            if i != c and A[i][c]:
                f = A[i][c]; A[i] = [F.sub(A[i][j], F.mul[f][A[c][j]]) for j in range(n + 1)]
    return [A[i][n] for i in range(n)]

def iso_via_basis(F, M, B, G):
    """is there a rescaling/relabelling of the natural basis B of A_RW(G) with g_v^2 = sum_{u~v} g_u ?"""
    n = len(B)
    S = [coords(F, B, prod(F, M, B[a], B[a])) for a in range(n)]
    if any(S[a][a] for a in range(n)): return False
    H = nx.Graph(); H.add_nodes_from(range(n))
    for a in range(n):
        for b in range(n):
            if S[a][b]:
                if not S[b][a]: return False
                H.add_edge(a, b)
    if not nx.is_isomorphic(H, G): return False
    for l0 in range(1, F.q):                       # lambda_b = lambda_a^2 S[a][b] along the support
        lam = {0: l0}; ok = True; st = [0]
        while st and ok:
            a = st.pop()
            for b in H[a]:
                val = F.mul[F.mul[lam[a]][lam[a]]][S[a][b]]
                if b in lam:
                    if lam[b] != val: ok = False; break
                else:
                    lam[b] = val; st.append(b)
        if ok: return True
    return False

def adj_M(F, G):
    n = G.number_of_nodes()
    return [[1 if G.has_edge(i, j) else 0 for j in range(n)] for i in range(n)]

def rw_M(F, G):
    n = G.number_of_nodes()
    return [[F.inv[F.n(G.degree(i))] if G.has_edge(i, j) else 0 for j in range(n)] for i in range(n)]

def twin_classes(G):
    cl = {}
    for v in sorted(G.nodes()): cl.setdefault(frozenset(G[v]), []).append(v)
    return [tuple(c) for c in cl.values()]

def quotient(G):
    classes = twin_classes(G); idx = {v: i for i, c in enumerate(classes) for v in c}
    H = nx.Graph(); H.add_nodes_from(range(len(classes)))
    for u, v in G.edges(): H.add_edge(idx[u], idx[v])
    return classes, H

_bc = {}
def thmN_count(F, M, G):
    """Theorem N: #natural bases (up to order, scaling) = prod over twin classes T of
       #(orthogonal bases up to order, scaling) of (k^T, sum lambda_i x_i y_i)."""
    tot = 1
    for T in twin_classes(G):
        r0 = M[T[0]]; j0 = next(j for j in range(len(r0)) if r0[j])
        lams = tuple(sorted(F.mul[M[w][j0]][F.inv[r0[j0]]] for w in T))
        key = (F.q, lams)
        if key not in _bc:
            _bc[key] = len(natural_bases(F, [[l] for l in lams])) if len(lams) > 1 else 1
        tot *= _bc[key]
    return tot

def thmP_iso(F, G):
    """Theorem P (char p odd) / Theorem I + Lemma 2 (char 2): exists c: V(H) -> k^x, constant on each N_H(S)
       (value gamma_S), with gamma_S^2 = c_S d_S if p does not divide |S|, and c_S d_S in Sim_|S|(F_q) otherwise
       (Sim_k(F_q) = F_q^x for k even, squares for k odd)."""
    classes, H = quotient(G); m = len(classes)
    d = [F.n(G.degree(c[0])) for c in classes]; k = [len(c) for c in classes]
    sq = {F.mul[x][x] for x in range(1, F.q)}
    for c in itertools.product(range(1, F.q), repeat=m):
        ok = True
        for S in range(m):
            nb = sorted(H[S])
            if any(c[T] != c[nb[0]] for T in nb): ok = False; break
            g = c[nb[0]]; cd = F.mul[c[S]][d[S]]
            if k[S] % F.p or F.p == 2:
                if F.mul[g][g] != cd: ok = False; break
            elif k[S] % 2 == 1 and cd not in sq: ok = False; break
        if ok: return True
    return False

def check_iso_ff(F, G, phi):
    """phi[v] = image of e_v in A_RW(G) (coordinates); checks bijectivity and multiplicativity."""
    n = G.number_of_nodes(); M = rw_M(F, G)
    if rank(F, [phi[v] for v in range(n)], n) < n: return False
    for v in range(n):
        for u in range(n):
            lhs = prod(F, M, phi[v], phi[u])
            if u != v:
                if any(lhs): return False
            else:
                rhs = [0] * n
                for w in G[v]:
                    rhs = [F.add[x][y] for x, y in zip(rhs, phi[w])]
                if list(lhs) != rhs: return False
    return True

def sec(title):
    print('\n' + title + '   [t = %.1f s]' % (time.time() - T0)); sys.stdout.flush()

# ----------------------------------------------------------------------------------------------
sec('[N]+[P] finite fields: natural bases and isomorphisms by brute force')
small = [g for g in graph_atlas_g() if 2 <= g.number_of_nodes() <= 5 and nx.is_connected(g)]
plan = [(3, 5), (5, 4), (7, 4)] + ([(5, 5), (7, 5)] if FULL else [])
for q, nmax in plan:
    F = GF(q); tested = isos = 0
    for G in small:
        n = G.number_of_nodes()
        if n > nmax or any(G.degree(v) % F.p == 0 for v in G): continue
        tested += 1
        Mrw, Ma = rw_M(F, G), adj_M(F, G)
        NBrw, NBa = natural_bases(F, Mrw), natural_bases(F, Ma)
        check(len(NBrw) == thmN_count(F, Mrw, G) and len(NBa) == thmN_count(F, Ma, G),
              'Theorem N count q=%d %s' % (q, list(G.edges())))
        brute = any(iso_via_basis(F, Mrw, B, G) for B in NBrw)
        isos += brute
        check(brute == thmP_iso(F, G), 'Theorem P q=%d %s' % (q, list(G.edges())))
    print('  F_%d, n <= %d: %d graphs (all degrees prime to p), %d with A(G) ~= A_RW(G); counts and criterion agree'
          % (q, nmax, tested, isos))

K113 = nx.Graph([(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4)])   # u=0, v=1, twins 2,3,4
T6 = nx.Graph([(0, 1), (0, 2), (0, 3), (0, 4), (4, 5)])                     # claw 0;1,2,3 plus path 0-4-5
for name, G, qs in [('K_{1,1,3}', K113, [3, 9] if FULL else [3]), ('T6', T6, [3, 9] if FULL else [3])]:
    for q in qs:
        F = GF(q); M = rw_M(F, G); NB = natural_bases(F, M)
        brute = any(iso_via_basis(F, M, B, G) for B in NB)
        check(len(NB) == thmN_count(F, M, G), 'count %s q=%d' % (name, q))
        check(brute == (q == 9) and thmP_iso(F, G) == (q == 9), '%s q=%d' % (name, q))
        print('  %s over F_%d: %d natural bases of A_RW (Theorem N predicts %d); isomorphic: %s (brute force), %s (Theorem P)'
              % (name, q, len(NB), thmN_count(F, M, G), brute, thmP_iso(F, G)))
# explicit isomorphisms over F_9 = F_3(i): P = [[0,-1,-1],[1,i,-i],[1,-i,i]] (P^T P = -I, P 1 = 1)
F = GF(9); i_, m1 = F.I, F.neg[1]; mi = F.neg[F.I]
Pcols = [(0, 1, 1), (m1, i_, mi), (m1, mi, i_)]          # columns of P
def emb(n, pos, vals):
    v = [0] * n
    for p_, x in zip(pos, vals): v[p_] = x
    return tuple(v)
phiK = {0: emb(5, [0], [1]), 1: emb(5, [1], [1])}
for a in range(3): phiK[2 + a] = emb(5, [2, 3, 4], Pcols[a])
check(check_iso_ff(F, K113, phiK), 'explicit iso K113 F9')
phiT = {0: emb(6, [0], [m1]), 4: emb(6, [4], [1]), 5: emb(6, [5], [m1])}
for a in range(3): phiT[1 + a] = emb(6, [1, 2, 3], Pcols[a])
check(check_iso_ff(F, T6, phiT), 'explicit iso T6 F9')
print('  explicit isomorphisms over F_9 (twins mixed by P = [[0,-1,-1],[1,i,-i],[1,-i,i]]) verified for K_{1,1,3} and T6')
# char 2: A_RW(G) = A(G) when all degrees are odd; the double star S(2,2) is neither regular nor biregular
S22 = nx.Graph([(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)])
F2 = GF(2); check(rw_M(F2, S22) == adj_M(F2, S22), 'char 2 equality')
print('  char 2: A_RW(S(2,2)) = A(S(2,2)) literally (degrees 3,3,1,1,1,1)')

# ----------------------------------------------------------------------------------------------
sec('[G] full GL_n systems (n^2 unknowns + 1/det), sympy Groebner over Q')
def full_system(n, edges):
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(edges)
    A = sp.Matrix(n, n, lambda i, j: 1 if G.has_edge(i, j) else 0)
    d = [G.degree(i) for i in range(n)]
    Mp = sp.Matrix(n, n, lambda i, j: sp.Rational(A[i, j], d[i]))    # A_RW: e'_i^2 = sum_j Mp[i,j] e'_j
    P = sp.Matrix(n, n, lambda w, v: sp.Symbol('p_%d_%d' % (w, v)))   # phi(e_v) = sum_w P[w,v] e'_w
    t = sp.Symbol('t'); eqs = []
    for v in range(n):
        for u in range(v + 1, n):
            for j in range(n): eqs.append(sp.expand(sum(P[w, v] * P[w, u] * Mp[w, j] for w in range(n))))
        for j in range(n):
            eqs.append(sp.expand(sum(P[w, v]**2 * Mp[w, j] for w in range(n)) - sum(P[j, u] for u in range(n) if A[v, u])))
    eqs.append(sp.expand(t * P.det() - 1))
    return [e for e in eqs if e != 0], P, t
eqs, P, t = full_system(3, [(0, 1), (0, 2)])
GB = sp.groebner(eqs, t, *[x for x in P if x != P[0, 0]], P[0, 0], order='lex')
uni = [g for g in GB.exprs if g.free_symbols == {P[0, 0]}]
check(uni == [P[0, 0]**3 - 4], 'P3 lex')
print('  P3: lex basis contains', uni, '-> no isomorphism over Q (4 is not a cube); exists over R, C')
tests = [('K_{1,3}', 4, [(0, 1), (0, 2), (0, 3)], 9), ('diamond', 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)], None)]
if FULL:
    tests += [('P4', 4, [(0, 1), (1, 2), (2, 3)], None), ('paw', 4, [(0, 1), (0, 2), (1, 2), (0, 3)], None),
              ('C4', 4, [(0, 1), (1, 2), (2, 3), (0, 3)], 0), ('K4', 4, list(itertools.combinations(range(4), 2)), 0)]
for name, n, edges, cube in tests:
    eqs, P, t = full_system(n, edges)
    GB = sp.groebner(eqs, *list(P), t, order='grevlex')
    if cube is None:
        check(list(GB.exprs) == [1], name)
        print('  %s: Groebner basis = [1] -> not isomorphic over C, hence over no field of characteristic 0' % name)
    elif cube == 0:
        check(list(GB.exprs) != [1], name)
        print('  %s (regular): Groebner basis != [1], consistent with phi = deg * id' % name)
    else:
        check(list(GB.exprs) != [1] and GB.reduce(P[0, 0]**3 - cube)[1] == 0, name)
        print('  %s: p00^3 - %d lies in the ideal -> no isomorphism over Q; exists over R, C' % (name, cube))

# ----------------------------------------------------------------------------------------------
sec('[L] Lemma 2 (norm identity) and the k = 2 elimination')
al, ga, w = sp.symbols('alpha gamma w')
for k in range(1, 7):
    P = sp.Matrix(k, k, lambda i, j: sp.Symbol('q%d_%d' % (i, j))); one = sp.ones(k, 1)
    lhs = k * (ga**2 - al)
    rhs = (one.T * (P.T * P - al * sp.eye(k)) * one)[0] - ((P * one - ga * one).T * (P * one + ga * one))[0]
    check(sp.expand(lhs - rhs) == 0, 'Lemma 2 identity k=%d' % k)
P = sp.Matrix(2, 2, lambda i, j: sp.Symbol('q%d_%d' % (i, j)))
M2 = P.T * P - al * sp.eye(2)
eqs = [M2[0, 0], M2[0, 1], M2[1, 1]] + list(P * sp.ones(2, 1) - ga * sp.ones(2, 1)) + [w * al * ga - 1]
GB = sp.groebner([sp.expand(e) for e in eqs], *list(P), w, al, ga, order='lex')
elim = [g for g in GB.exprs if g.free_symbols <= {al, ga}]
check(elim == [al - ga**2], 'k=2 elimination')
print('  identity holds for k = 1..6; k = 2 elimination ideal:', elim)

# ----------------------------------------------------------------------------------------------
sec('[A] characteristic 0: all connected graphs with 2..7 vertices, over C, R, Q')
def is_rat_cube(fr):
    def ic(n):
        r = round(abs(n) ** (1 / 3))
        return any(c**3 == abs(n) for c in (r - 1, r, r + 1))
    return ic(fr.numerator) and ic(fr.denominator)

def criterion(G):
    dg = dict(G.degree())
    if len(set(dg.values())) == 1: return 'regular', None
    if nx.is_bipartite(G):
        X, Y = nx.bipartite.sets(G)
        dX, dY = {dg[v] for v in X}, {dg[v] for v in Y}
        if len(dX) == 1 and len(dY) == 1: return 'biregular', (X, Y, dX.pop(), dY.pop())
    return 'no', None

def predicted(G):
    kind, data = criterion(G)
    if kind == 'regular': return dict(C=True, R=True, Q=True)
    if kind == 'biregular': return dict(C=True, R=True, Q=is_rat_cube(Fraction(data[2], data[3])))
    return dict(C=False, R=False, Q=False)

def c_system_decide(G):
    """Theorem I + Lemma 2: c_S (S twin class), c constant on N_H(S), (common value)^2 = c_S d_S."""
    classes, H = quotient(G); m = len(classes)
    c = sp.symbols('c0:%d' % m); z = sp.Symbol('z')
    d = [G.degree(cl[0]) for cl in classes]; eqs = []
    for S in range(m):
        nb = sorted(H[S])
        eqs += [c[T] - c[nb[0]] for T in nb[1:]] + [c[nb[0]]**2 - c[S] * d[S]]
    eqs.append(z * sp.Mul(*c) - 1)
    GB = sp.groebner(eqs, *c, z, order='lex')
    if list(GB.exprs) == [1]: return dict(C=False, R=False, Q=False)
    sols = sp.solve(list(GB.exprs), list(c) + [z], dict=True)
    real = any(all(sp.simplify(v).is_real for v in s.values()) for s in sols)
    rat = any(all(sp.simplify(v).is_rational for v in s.values()) for s in sols)
    return dict(C=len(sols) > 0, R=real, Q=rat)

def explicit_iso_ok(G, field):
    """build the diagonal isomorphism of Theorem A and check every product exactly."""
    kind, data = criterion(G); n = G.number_of_nodes(); nodes = sorted(G.nodes())
    if kind == 'regular':
        dd = G.degree(nodes[0]); scale = {v: sp.Integer(dd) for v in nodes}
    else:
        X, Y, d1, d2 = data
        tt = sp.Rational(d1, d2) ** sp.Rational(1, 3)          # real cube root
        if field == 'Q':
            tt = sp.nsimplify(tt)
            if not tt.is_rational: return False
        scale = {v: (tt**2 * d2 if v in X else tt * d2) for v in nodes}
    deg = dict(G.degree())
    for v in nodes:     # phi(e_v)^2 = scale_v^2/deg(v) * sum_{w~v} e'_w  must equal  sum_{w~v} scale_w e'_w
        for w_ in G[v]:
            if sp.simplify(scale[v]**2 / deg[v] - scale[w_]) != 0: return False
    return all(sp.simplify(s).is_real for s in scale.values()) and (field != 'Q' or all(s.is_rational for s in scale.values()))

graphs = [g for g in graph_atlas_g() if g.number_of_nodes() >= 2 and nx.is_connected(g)]
extra = [('K_{1,8}', nx.complete_bipartite_graph(1, 8)),
         ('subdivided K_4', nx.Graph([(a, 4 + i) for i, (a, b) in enumerate(itertools.combinations(range(4), 2))] +
                                     [(b, 4 + i) for i, (a, b) in enumerate(itertools.combinations(range(4), 2))])),
         ('subdivided K_5', nx.Graph([(a, 5 + i) for i, (a, b) in enumerate(itertools.combinations(range(5), 2))] +
                                     [(b, 5 + i) for i, (a, b) in enumerate(itertools.combinations(range(5), 2))])),
         ('Petersen', nx.petersen_graph()), ('K_{2,16}', nx.complete_bipartite_graph(2, 16))]
counts = {'regular': 0, 'biregular': 0, 'no': 0, 'C': 0, 'R': 0, 'Q': 0}
def run_A(g):
    p = predicted(g); dcs = c_system_decide(g)
    check(p == dcs, 'Theorem A vs c-system %s' % list(g.edges()))
    for fld in 'CRQ':
        if p[fld]: check(explicit_iso_ok(g, fld), 'explicit iso %s over %s' % (list(g.edges()), fld))
    return dcs
for g in graphs:
    dcs = run_A(g); counts[criterion(g)[0]] += 1
    for fld in 'CRQ': counts[fld] += dcs[fld]
print('  %d atlas graphs: %d regular, %d biregular non-regular, %d other; A(G) ~= A_RW(G) over C: %d, R: %d, Q: %d'
      % (len(graphs), counts['regular'], counts['biregular'], counts['no'], counts['C'], counts['R'], counts['Q']))
for name, g in extra:
    print('  %-16s %-9s A(G) ~= A_RW(G) over C, R, Q: %s' % (name, criterion(g)[0], run_A(g)))
print('  every predicted isomorphism (diagonal, Theorem A) checked exactly; c-system == criterion for all graphs')

# ----------------------------------------------------------------------------------------------
sec('[H] remarks: homomorphisms over C; BCS22 (Cor. 2.7, Thm 2.4, Cor. 2.17, Thm 2.11) in small cases')
I = sp.I
def hom_ok(G, fimg):
    """fimg[v] = image of e_v (dict vertex -> coefficient) of a linear map A_RW(G) -> A(G)."""
    nodes = sorted(G.nodes()); deg = dict(G.degree())
    def mulA(x, y):            # product in A(G)
        out = {}
        for v in nodes:
            c = sp.expand(x.get(v, 0) * y.get(v, 0))
            if c != 0:
                for w_ in G[v]: out[w_] = sp.expand(out.get(w_, 0) + c)
        return {k: v for k, v in out.items() if v != 0}
    def lin(coeffs):           # f(sum coeffs_v e_v)
        out = {}
        for v, c in coeffs.items():
            for w_, a in fimg.get(v, {}).items(): out[w_] = sp.expand(out.get(w_, 0) + c * a)
        return {k: v for k, v in out.items() if v != 0}
    for v in nodes:
        for u in nodes:
            prod_rw = {} if u != v else {w_: sp.Rational(1, deg[v]) for w_ in G[v]}
            if lin(prod_rw) != mulA(fimg.get(v, {}), fimg.get(u, {})): return False
    return True
diamond = nx.Graph([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)])        # twins 2, 3
check(hom_ok(diamond, {2: {2: 1, 3: I}, 3: {2: -1, 3: -I}}), 'hom diamond')
dstar = nx.Graph([(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)])          # C_{2,2}, twins 2, 3
check(hom_ok(dstar, {2: {2: 1, 3: I}, 3: {2: -1, 3: -I}}), 'hom double star')
print('  e_x -> e_x + i e_y, e_y -> -(e_x + i e_y), others -> 0 is a non-zero, non-injective homomorphism')
print('  A_RW(G) -> A(G) over C for the diamond and the double star C_{2,2} (neither is isomorphic)')
for q, Mx in [(3, [[0, 0, 1], [0, 0, 2], [1, 0, 0]]), (2, [[0, 0, 1], [0, 0, 1], [1, 0, 0]])]:
    F = GF(q); nb = natural_bases(F, Mx)
    check(len(nb) == 1, 'BCS F_%d' % q)
    print('  F_%d, e1^2 = e3, e2^2 = %s e3, e3^2 = e1: non-degenerate, not (2LI), natural bases up to order/scaling: %d'
          % (q, '-' if q == 3 else '', len(nb)))
# BCS22 Thm 2.4(i) "<=", Prop 2.16, Cor 2.17 in characteristic 2
F = GF(2)
M4 = [[0, 0, 0, 1], [0, 0, 0, 1], [0, 0, 0, 1], [1, 0, 0, 0]]          # e1^2 = e2^2 = e3^2 = e4, e4^2 = e1
NB4 = natural_bases(F, M4); u = (1, 1, 1, 0)
check(any(prod(F, M4, u, u)) and not any(u in B for B in NB4), 'BCS 2.4 char 2')
M5 = [[0, 0, 0, 0, 1]] * 4 + [[1, 0, 0, 0, 0]]                         # e1^2 = .. = e4^2 = e5, e5^2 = e1
NB5 = natural_bases(F, M5); u1, u2 = (1, 1, 1, 0, 0), (0, 0, 0, 1, 0)
check(any(u1 in B for B in NB5) and any(u2 in B for B in NB5) and not any(u1 in B and u2 in B for B in NB5)
      and not any(prod(F, M5, u1, u2)), 'BCS 2.17 char 2')
print('  F_2: u = e1+e2+e3 has u^2 = e4 != 0 and rank-1 support squares but lies in no natural basis (BCS22 Thm 2.4(i));')
print('       e1+e2+e3 and e4 (5-dim example) are orthogonal natural vectors in no common natural basis (BCS22 Cor 2.17)')
P5A = sp.Matrix(5, 5, lambda i, j: 1 if abs(i - j) == 1 else 0)
check(P5A.rank() == 4 and len(twin_classes(nx.path_graph(5))) == 5, 'P5')
print('  A(P_5) over Q: dim A^2 = rank = %d, but 5 twin classes (BCS22 Thm 2.11 says r = dim A^2)' % P5A.rank())

print('\n%s   [total %.1f s]' % ('ALL CHECKS PASSED' if not FAILS else 'FAILURES: %d' % len(FAILS), time.time() - T0))
sys.exit(1 if FAILS else 0)
