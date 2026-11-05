#!/usr/bin/env python3
"""sym.py - search 3-colourings of K_n with a prescribed automorphism group such that
every K_k uses all 3 colours (i.e. each union of two colour classes is K_k-free).
Edge orbits of the group are the SAT variables; hyperedges come from ./hyp (C)."""
import sys, os, subprocess, itertools, time, json
HERE = os.path.dirname(os.path.abspath(__file__))

class UF:
    def __init__(s, n): s.p = list(range(n))
    def f(s, x):
        while s.p[x] != x:
            s.p[x] = s.p[s.p[x]]; x = s.p[x]
        return x
    def u(s, a, b):
        a, b = s.f(a), s.f(b)
        if a != b: s.p[max(a, b)] = min(a, b)

def orbits_from_gens(n, gens):
    """gens: list of permutations (lists). Returns class matrix C (n x n, -1 diag), m, vertex roots."""
    idx = {}
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    for i, p in enumerate(pairs): idx[p] = i
    uf = UF(len(pairs))
    for g in gens:
        for i, (a, b) in enumerate(pairs):
            x, y = g[a], g[b]
            if x > y: x, y = y, x
            uf.u(i, idx[(x, y)])
    rep = {}
    C = [[-1] * n for _ in range(n)]
    for i, (a, b) in enumerate(pairs):
        r = uf.f(i)
        if r not in rep: rep[r] = len(rep)
        C[a][b] = C[b][a] = rep[r]
    vuf = UF(n)
    for g in gens:
        for a in range(n): vuf.u(a, g[a])
    roots = sorted({vuf.f(a) for a in range(n)})
    return C, len(rep), roots

def hyperedges(n, C, m, roots, k=6):
    inp = f"{n} {m} {k} {len(roots)}\n" + " ".join(map(str, roots)) + "\n" + "\n".join(" ".join(map(str, row)) for row in C) + "\n"
    out = subprocess.run([os.path.join(HERE, "hyp")], input=inp, capture_output=True, text=True, check=True)
    lines = out.stdout.split("\n")
    m2, ne = map(int, lines[0].split())
    E = [list(map(int, l.split())) for l in lines[1:1 + ne]]
    return E, out.stderr.strip()

def solve(m, E, solver_name="cadical153", time_budget=None, extra=None):
    from pysat.solvers import Solver
    v = lambda c, col: 3 * c + col + 1
    cl = []
    for c in range(m):
        cl.append([v(c, 0), v(c, 1), v(c, 2)])
        for a, b in itertools.combinations(range(3), 2): cl.append([-v(c, a), -v(c, b)])
    for T in E:
        if len(T) < 3: return None, "trivially UNSAT (hyperedge size %d)" % len(T)
        for col in range(3): cl.append([v(c, col) for c in T])
    # colour-permutation symmetry breaking
    cl.append([v(0, 0)])
    if m > 1: cl.append([-v(1, 2)])
    if extra: cl += extra
    with Solver(name=solver_name, bootstrap_with=cl) as s:
        ok = s.solve()
        if not ok: return None, "UNSAT"
        model = set(x for x in s.get_model() if x > 0)
        col = [next(cc for cc in range(3) if v(c, cc) in model) for c in range(m)]
        return col, "SAT"

def full_colouring(n, C, col):
    return [[-1 if a == b else col[C[a][b]] for b in range(n)] for a in range(n)]

def check_fast(M, k=6):
    """count k-cliques (up to 1) in each two-colour union via bitset backtracking"""
    n = len(M)
    for miss in range(3):
        adj = [0] * n
        for a in range(n):
            for b in range(n):
                if a != b and M[a][b] != miss: adj[a] |= 1 << b
        def ext(cand, size):
            if size == k: return True
            while cand:
                v = cand.bit_length() - 1
                cand &= ~(1 << v)
                if size + 1 + bin(cand & adj[v]).count("1") >= k - 0 and ext(cand & adj[v], size + 1): return True
                if size + bin(cand).count("1") < k: return False
            return False
        if ext((1 << n) - 1, 0): return False, miss
    return True, None

def write_colouring(fn, M, note=""):
    n = len(M)
    with open(fn, "w") as f:
        f.write(f"# 3-colouring of K_{n}; row i, column j = colour of edge ij (0,1,2); '.' on diagonal. {note}\n")
        f.write(f"{n}\n")
        for a in range(n): f.write("".join("." if a == b else str(M[a][b]) for b in range(n)) + "\n")

# ---------------- group constructions ----------------
def left_regular(mul):
    n = len(mul)
    return [list(mul[g]) for g in range(n)]  # x -> g*x  (mul[g][x])

def cyc_gen(n): return [[(x + 1) % n for x in range(n)]]

def block_circ(mod, r, fixed=0, mults=()):
    """Z_mod acting on r orbits (vertex (x,i) -> i*mod+x) plus 'fixed' fixed points; optional multipliers a acting x->a x on all orbits"""
    n = r * mod + fixed
    g = list(range(n))
    for i in range(r):
        for x in range(mod): g[i * mod + x] = i * mod + (x + 1) % mod
    gens = [g]
    for a in mults:
        h = list(range(n))
        for i in range(r):
            for x in range(mod): h[i * mod + x] = i * mod + (a * x) % mod
        gens.append(h)
    return n, gens

if __name__ == "__main__":
    pass

