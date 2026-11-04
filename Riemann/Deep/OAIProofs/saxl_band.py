"""Lemma-level checks of OpenAI's "A Cyclic Polytabloid Proof of Saxl's Conjecture"
(Sections 3-5), recomputed from the definitions with exact arithmetic.

 C  band high set: the only 0-1 filling of B_m with the staircase high margins is B_h (m<=9),
    and P_s(w_m) = w~_h (x) u_{m,s} with all cross signs cancelling (full expansion, m<=5)
 D  path contraction <u_{m,2}, X_1 (x)...(x) X_K> = [X_1 adj(X_2) X_3 ... X_K]_{11} (random integer X, m=3..7),
    the sign u = (-1)^(m-1) u_{m,2}, the segment table A_{r,p}, Example 5.4
 D3 for every eta |- 2m-1 with <=4 rows (m=3..8): an explicit invertible Q with <u, t_N> != 0 (N_i = Q M_i Q^-1),
    pairing computed by brute-force expansion, not by the matrix formula
 H  Q(v_R (x) v_R) = (-1)^{sum C(r,2)} prod r! t_m   (full expansion, m<=4)
"""
import itertools, random, sys
from fractions import Fraction
import mn

def perm_sign(p):
    s = 1; p = list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s

def triangle(m):
    return [(i, j) for i in range(1, m+1) for j in range(1, m+2-i)]

def rows_of(cells):
    d = {}
    for (i, j) in cells: d.setdefault(i, []).append((i, j))
    return [sorted(v, key=lambda c: c[1]) for k, v in sorted(d.items())]

def cols_of(cells):
    d = {}
    for (i, j) in cells: d.setdefault(j, []).append((i, j))
    return [sorted(v, key=lambda c: c[0]) for k, v in sorted(d.items())]

def alt_blocks(blocks, letter_offset=0):
    """list of (assignment dict cell->letter, coeff) for the product of alternating blocks"""
    terms = [({}, 1)]
    for blk in blocks:
        r = len(blk); new = []
        for perm in itertools.permutations(range(r)):
            sg = perm_sign(perm)
            for (asg, c) in terms:
                a2 = dict(asg)
                for k, cell in enumerate(blk): a2[cell] = perm[k] + 1 + letter_offset
                new.append((a2, c*sg))
        terms = new
    return terms

# ---------------- C: high-set uniqueness and factorization ----------------
def count_fillings(m, s):
    h = m - s
    cells = triangle(m)
    rsum = {i: max(h+1-i, 0) for i in range(1, m+1)}
    csum = {j: max(h+1-j, 0) for j in range(1, m+1)}
    sols = []
    def rec(i, colleft, chosen):
        if i > m:
            if all(v == 0 for v in colleft.values()): sols.append(list(chosen))
            return
        rowcells = [c for c in cells if c[0] == i]
        for sub in itertools.combinations(rowcells, rsum[i]):
            ok = True; cl = dict(colleft)
            for (_, j) in sub:
                cl[j] -= 1
                if cl[j] < 0: ok = False; break
            if ok: rec(i+1, cl, chosen + list(sub))
    rec(1, dict(csum), [])
    return sols, h

def check_factorization(m, s):
    h = m - s
    B = triangle(m); D = set(triangle(h)); F = [c for c in B if c not in D]
    vR = alt_blocks(rows_of(B)); vC = alt_blocks(cols_of(B))
    def highset(asg): return frozenset(c for c, a in asg.items() if a > s)
    from collections import defaultdict
    gR = defaultdict(list); gC = defaultdict(list)
    for asg, c in vR: gR[highset(asg)].append((asg, c))
    for asg, c in vC: gC[highset(asg)].append((asg, c))
    common = set(gR) & set(gC)
    assert common == {frozenset(D)}, ("high sets", m, s, len(common))
    def as_dict(terms):
        return {tuple(sorted(a.items())): c for a, c in terms}
    def shifted(terms):
        return [({cell: a+s for cell, a in asg.items()}, c) for asg, c in terms]
    tot_sign = 1
    for name, blocks_full, blocksD, blocksF, grp in (("rows", rows_of(B), rows_of(sorted(D)), rows_of(F), gR), ("cols", cols_of(B), cols_of(sorted(D)), cols_of(F), gC)):
        lhs = as_dict(grp[frozenset(D)])
        hi = shifted(alt_blocks(blocksD)) if D else [({}, 1)]
        lo = alt_blocks(blocksF)
        prod = {}
        for a1, c1 in hi:
            for a2, c2 in lo:
                a = dict(a1); a.update(a2); prod[tuple(sorted(a.items()))] = c1*c2
        assert set(prod) == set(lhs), (name, "support")
        ratios = {lhs[k]*prod[k] for k in lhs}
        assert len(ratios) == 1, (name, "nonuniform sign")
        sg = ratios.pop()
        pred = 1
        for blk in blocks_full:
            r = len(blk); a = max(r-s, 0); b = min(r, s); pred *= (-1)**(a*b)
        assert sg == pred, (name, sg, pred)
        tot_sign *= sg
    assert tot_sign == 1
    return len(vR)

# ---------------- D: path contraction ----------------
def band_cells(m):
    return [c for c in triangle(m) if c[0]+c[1] >= m]

def path_positions(m):
    K = 2*m-1; p = {}
    for j in range(1, m+1): p[2*j-1] = (m+1-j, j)
    for j in range(1, m): p[2*j] = (m-j, j)
    return [p[q] for q in range(1, K+1)]

def u_terms(m, path_order=False):
    """u_{m,2} (paper orders) or the path-ordered u, as list of (dict cell->(a,b), coeff)"""
    F = band_cells(m); pos = path_positions(m); idx = {c: k for k, c in enumerate(pos)}
    rb = rows_of(F); cb = cols_of(F)
    if path_order:
        rb = [sorted(b, key=lambda c: idx[c]) for b in rb]
        cb = [sorted(b, key=lambda c: idx[c]) for b in cb]
    R = alt_blocks(rb); C = alt_blocks(cb)
    out = []
    for a1, c1 in R:
        for a2, c2 in C:
            out.append(({cell: (a1[cell], a2[cell]) for cell in F}, c1*c2))
    return out

def matmul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def adj(X):
    return [[X[1][1], -X[0][1]], [-X[1][0], X[0][0]]]

def pairing(terms, mats_by_cell):
    tot = 0
    for asg, c in terms:
        pr = c
        for cell, (a, b) in asg.items():
            pr *= mats_by_cell[cell][a-1][b-1]
            if pr == 0: break
        tot += pr
    return tot

def check_contraction(m, trials, rng):
    pos = path_positions(m); K = len(pos)
    U = u_terms(m); Up = u_terms(m, path_order=True)
    dU = {tuple(sorted(a.items())): c for a, c in U}; dUp = {tuple(sorted(a.items())): c for a, c in Up}
    assert set(dU) == set(dUp) and all(dUp[k] == (-1)**(m-1)*dU[k] for k in dU), "sign u = (-1)^(m-1) u_{m,2}"
    for _ in range(trials):
        X = [[[rng.randint(-3, 3) for _ in range(2)] for _ in range(2)] for _ in range(K)]
        lhs = pairing(U, {pos[q]: X[q] for q in range(K)})
        P = [[1, 0], [0, 1]]
        for q in range(K):
            P = matmul(P, X[q] if q % 2 == 0 else adj(X[q]))
        assert lhs == P[0][0], (m, lhs, P[0][0])
    return len(U)

I2 = [[1, 0], [0, 1]]; Zm = [[1, 0], [0, -1]]; Tm = [[0, 1], [1, 0]]; Jm = [[0, 1], [-1, 0]]
BASIS = [I2, Zm, Tm, Jm]

def segment(r, p, basis=BASIS):
    tot = [[0, 0], [0, 0]]
    for perm in itertools.permutations(range(r)):
        sg = perm_sign(perm); P = I2
        for k in range(r):
            M = basis[perm[k]]
            P = matmul(P, adj(M) if (p+k) % 2 == 0 else M)
        tot = [[tot[i][j] + sg*P[i][j] for j in range(2)] for i in range(2)]
    return tot

def check_segments():
    def sc(c, M): return [[c*x for x in row] for row in M]
    table = {(1, 1): I2, (1, 0): I2, (2, 1): sc(-2, Zm), (2, 0): sc(2, Zm), (3, 1): sc(-6, Jm), (3, 0): sc(-6, Jm), (4, 1): sc(-24, I2), (4, 0): sc(24, I2)}
    for r in range(1, 5):
        for p in (1, 2, 3, 4):
            assert segment(r, p) == table[(r, p % 2)], (r, p, segment(r, p))
    assert matmul(sc(-6, Jm), sc(2, Zm)) == sc(12, Tm)
    Q = [[1, 1], [1, 2]]; Qi = [[2, -1], [-1, 1]]
    assert matmul(Q, Qi) == I2 and matmul(matmul(Q, Tm), Qi) == [[1, 0], [3, -1]]

def conjugated_basis(Q):
    det = Q[0][0]*Q[1][1] - Q[0][1]*Q[1][0]
    Qi = [[Fraction(Q[1][1], det), Fraction(-Q[0][1], det)], [Fraction(-Q[1][0], det), Fraction(Q[0][0], det)]]
    return [matmul(matmul(Q, M), Qi) for M in BASIS]

def tN_value(word_by_pos, segments, N):
    """value of the dual polytabloid t_N on a pair-word, by direct expansion"""
    val = 1
    for seg in segments:
        r = len(seg); s = 0
        for perm in itertools.permutations(range(r)):
            pr = perm_sign(perm)
            for k, q in enumerate(seg):
                a, b = word_by_pos[q]
                pr *= N[perm[k]][a-1][b-1]
                if pr == 0: break
            s += pr
        val *= s
        if val == 0: return 0
    return val

def check_two_band(m, rng):
    pos = path_positions(m); K = len(pos); U = u_terms(m)
    words = [([asg[pos[q]] for q in range(K)], c) for asg, c in U]
    results = {}
    for eta in mn.partitions(K):
        if len(eta) > 4: continue
        cols = mn.conj(eta)
        segs = []; st = 0
        for r in cols: segs.append(list(range(st, st+r))); st += r
        found = None
        cands = [[[1, 0], [0, 1]], [[1, 1], [1, 2]], [[1, 2], [0, 1]], [[2, 1], [1, 1]]] + [[[rng.randint(-4, 4) for _ in range(2)] for _ in range(2)] for _ in range(40)]
        for Q in cands:
            if Q[0][0]*Q[1][1] - Q[0][1]*Q[1][0] == 0: continue
            N = conjugated_basis(Q)
            val = sum(c*tN_value(w, segs, N) for w, c in words)
            if val != 0: found = (Q, val); break
        results[eta] = found
        assert found is not None, ("no nonzero pairing found", m, eta)
    return results

# ---------------- H: dominated quotient map ----------------
def check_dominated(m):
    B = triangle(m); vR = alt_blocks(rows_of(B))
    from collections import Counter
    out = Counter()
    for a1, c1 in vR:
        for a2, c2 in vR:
            word = tuple(a1[c] + a2[c] - 1 for c in B)
            cnt = Counter(word)
            if all(cnt.get(r, 0) == r for r in range(1, m+1)) and sum(cnt.values()) == len(B):
                out[word] += c1*c2
    out = {w: c for w, c in out.items() if c != 0}
    tm = tuple(m+1-c[0] for c in B)
    from math import factorial, comb, prod
    pred = (-1)**sum(comb(r, 2) for r in range(1, m+1)) * prod(factorial(r) for r in range(1, m+1))
    assert out == {tm: pred}, (m, out)

if __name__ == "__main__":
    rng = random.Random(7)
    for m in range(2, 10):
        for s in (1, 2):
            if m - s < 1: continue
            sols, h = count_fillings(m, s)
            assert len(sols) == 1 and set(sols[0]) == set(triangle(h)), (m, s, len(sols))
    print("C: unique high set B_h for all 2<=m<=9, s=1,2")
    for m in range(2, 6):
        for s in (1, 2):
            if m - s < 1: continue
            n = check_factorization(m, s)
    print("C: P_s(w_m) = w~_h (x) u_{m,s}, cross signs (-1)^{ab} as stated, product 1, m<=5")
    for m in range(3, 8):
        nt = check_contraction(m, 60 if m < 7 else 15, rng)
    print("D: path contraction formula and u = (-1)^(m-1) u_{m,2} verified, m=3..7")
    check_segments(); print("D: segment table and Example 5.4 verified")
    for m in range(3, 9):
        res = check_two_band(m, rng)
        print(f"D3: m={m}: all {len(res)} shapes eta |- {2*m-1} with <=4 rows have a nonzero pairing <u, t_N>")
    for m in range(1, 5):
        check_dominated(m)
    print("H: Q(v_R (x) v_R) = (-1)^sum C(r,2) prod r! t_m for m<=4")
