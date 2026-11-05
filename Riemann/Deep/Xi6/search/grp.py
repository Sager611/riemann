#!/usr/bin/env python3
"""grp.py - small groups as multiplication tables (element 0 = identity)."""
import itertools
def table(els, mul):
    idx = {e: i for i, e in enumerate(els)}
    T = [[idx[mul(a, b)] for b in els] for a in els]
    # move identity to 0
    e = next(i for i in range(len(els)) if all(T[i][j] == j for j in range(len(els))))
    if e != 0:
        perm = list(range(len(els))); perm[0], perm[e] = perm[e], perm[0]
        inv = {perm[i]: i for i in range(len(els))}
        T = [[inv[T[perm[a]][perm[b]]] for b in range(len(els))] for a in range(len(els))]
    return T
def check_group(T):
    n = len(T)
    assert all(T[0][x] == x and T[x][0] == x for x in range(n))
    for a in range(n):
        assert sorted(T[a]) == list(range(n))
    for a, b, c in itertools.product(range(n), repeat=3):
        if T[T[a][b]][c] != T[a][T[b][c]]: return False
    return True
def abelian(ns):
    els = list(itertools.product(*[range(x) for x in ns]))
    return table(els, lambda a, b: tuple((x + y) % q for x, y, q in zip(a, b, ns)))
def product(T1, T2):
    n1, n2 = len(T1), len(T2)
    els = [(a, b) for a in range(n1) for b in range(n2)]
    return table(els, lambda u, v: (T1[u[0]][v[0]], T2[u[1]][v[1]]))
def semidirect(ns, Mat, d):
    """N = Z_ns (abelian), automorphism given by integer matrix Mat (acts on column vectors), cyclic top Z_d.
    (x,i)(y,j) = (x + Mat^i y, i+j). Mat^d must be identity on N."""
    r = len(ns)
    def app(Mt, y): return tuple(sum(Mt[a][b] * y[b] for b in range(r)) % ns[a] for a in range(r))
    pw = [[[int(a == b) for b in range(r)] for a in range(r)]]
    for i in range(1, d + 1):
        P = pw[-1]
        pw.append([[sum(Mat[a][c] * P[c][b] for c in range(r)) for b in range(r)] for a in range(r)])
    Nels = list(itertools.product(*[range(x) for x in ns]))
    for y in Nels: assert app(pw[d], y) == y, "Mat^d != id"
    els = [(x, i) for i in range(d) for x in Nels]
    return table(els, lambda u, v: (tuple((a + b) % q for a, b, q in zip(u[0], app(pw[u[1]], v[0]), ns)), (u[1] + v[1]) % d))
def cyclic_sd(nn, a, d):  # Z_nn x| Z_d via x -> a x
    return semidirect([nn], [[a]], d)
def dihedral(k): return cyclic_sd(k, k - 1, 2)
def gen_dihedral(ns):  # Dih(A): inversion on abelian A
    r = len(ns); M = [[-1 if a == b else 0 for b in range(r)] for a in range(r)]
    return semidirect(ns, M, 2)
def dicyclic(k):  # order 4k: <a,x | a^{2k}=1, x^2=a^k, x a x^-1 = a^-1>; elements (i,j)=a^i x^j
    els = [(i, j) for j in range(2) for i in range(2 * k)]
    def mul(u, v):
        i1, j1 = u; i2, j2 = v
        if j1 == 0: return ((i1 + i2) % (2 * k), j2)
        # a^i1 x a^i2 x^j2 = a^{i1 - i2} x^{1+j2}
        if j2 == 0: return ((i1 - i2) % (2 * k), 1)
        return ((i1 - i2 + k) % (2 * k), 0)  # x^2 = a^k
    return table(els, mul)
def heis(p):  # Heisenberg group mod p (order p^3)
    els = list(itertools.product(range(p), repeat=3))
    return table(els, lambda u, v: ((u[0] + v[0]) % p, (u[1] + v[1]) % p, (u[2] + v[2] + u[0] * v[1]) % p))
def A4():
    return semidirect([2, 2], [[0, 1], [1, 1]], 3)
def S3(): return dihedral(3)

