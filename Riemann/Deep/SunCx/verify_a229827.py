"""Z.-W. Sun, OEIS A229827 (2013), Conjecture (i), general part: for every finite field F_q and every P(x) of degree < q-1
(equivalently sum_x P(x) = 0) not of the form c - x, some circular order (a_1, ..., a_q) of F_q has the q sums
P(a_1)+a_2, ..., P(a_q)+a_1 pairwise distinct.  Earlier monomial form (OEIS revisions #16-#36; NMBRTHRY list, 2013-10-01):
P = x^k, 0 <= k <= q-2, except 2 | q and k = 1.

Certified here (stdlib only; about 1 minute, < 200 MB):
 [1] Reduction.  For bijective P, an order works iff u o P is a q-cycle for some complete mapping u of (F_q,+)
     (u and u + id bijective; take u = pi o P^-1, pi: a_i -> a_{i+1}).  A q-cycle has sign (-1)^(q-1).
 [2] q = 4:  all 8 complete mappings are even (they are w x + c, w^2 x + c), a 4-cycle is odd, so every even P fails:
     besides the excluded x + c these are P = w x + c and w^2 x + c (w^2 = w + 1).  Every odd P works.
 [3] q = 8:  all 384 complete mappings are affine, hence even (AGL(3,2) <= A_8), an 8-cycle is odd, so every even P fails
     (20160 functions, 8 of them x + c), e.g. x^2, x^4, lambda x (lambda not in F_2), x^5+x^3+x, x^6+x^3+x^2.  Every odd P works.
 [4] q = 9:  9-cycles are even.  Of the 2241 complete mappings exactly 108 are odd, namely x -> Mx + c with M of order 8;
     AGL(2,3) contains no 9-cycle, so every odd affine P fails (x^3, c x^3 for c in {+-1,+-i}, lambda x with lambda
     primitive).  Exactly 648 odd P fail (5 orbits, incl. x^6+x^4+x^3+x^2 and x^6+x^4-x^3+x^2); all other odd P work.
     (Even P other than c - x and all non-bijective P work: classify_a229827.c.)
 [5] q <= 8, all P with sum 0 (also non-bijective): orbits under P -> gPg^-1 + c (g in AGL), each orbit decided by
     search: the failures are exactly [2], [3] (and c - x); q = 2, 3, 5, 7 have none.  q = 9: classify_a229827.c.
 [6] Witness orders (a229827_witnesses.txt) for every x^k, q in {4,8,9,16,25,27,32,49,64,81}, except (8,2), (8,4), (9,3),
     and for lambda x (lambda primitive, one per Galois orbit, q = 16..81).
 [7] (ii) fails on {0,1,2} with P = 2x (known: Sun, arXiv 1309.1679, Remark after Thm 1.3).
Run: python3 verify_a229827.py
"""
import itertools, os, sys, time
T0 = time.time()
def say(*a): print('[%5.1fs]' % (time.time() - T0), *a, flush=True)

# ---------- fields F_{p^k}: elements = base-p digit vectors, t root of lexicographically first monic irreducible ----------
def find_irred(p, k):
    def irreducible(c):          # c: coefficients low..high, monic, degree k
        for d in range(1, k // 2 + 1):
            for f in itertools.product(range(p), repeat=d):
                g = list(f) + [1]; r = c[:]
                for i in range(len(r) - 1, d - 1, -1):
                    co = r[i] % p
                    if co:
                        for j in range(d + 1): r[i - d + j] = (r[i - d + j] - co * g[j]) % p
                if all(x % p == 0 for x in r[:d]): return False
        return True
    for low in itertools.product(range(p), repeat=k):
        if low[0] and irreducible(list(low) + [1]): return low
def field(p, k):
    q = p ** k; low = find_irred(p, k)
    dig = lambda x: [(x // p ** i) % p for i in range(k)]
    enc = lambda d: sum((d[i] % p) * p ** i for i in range(k))
    def mul(a, b):
        A, B = dig(a), dig(b); r = [0] * (2 * k)
        for i in range(k):
            for j in range(k): r[i + j] += A[i] * B[j]
        for d in range(2 * k - 1, k - 1, -1):
            c = r[d] % p
            if c:
                r[d] = 0
                for i in range(k): r[d - k + i] -= c * low[i]
        return enc([x % p for x in r[:k]])
    ADD = [[enc([x + y for x, y in zip(dig(a), dig(b))]) for b in range(q)] for a in range(q)]
    MUL = [[mul(a, b) for b in range(q)] for a in range(q)]
    NEG = [enc([-x for x in dig(a)]) for a in range(q)]
    return q, ADD, MUL, NEG, low
def power(MUL, x, e):
    r = 1
    for _ in range(e): r = MUL[r][x]
    return r
def sign(P):
    seen = [0] * len(P); s = 1
    for i in range(len(P)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]: seen[j] = 1; j = P[j]; L += 1
            if L % 2 == 0: s = -s
    return s
def is_cycle(P):
    j, L = P[0], 1
    while j != 0: j = P[j]; L += 1
    return L == len(P)
def works(P, ADD, q, want_max=False):
    """search circular orders with a_1 = 0; returns (True, order) or (False, max number of distinct sums)"""
    A = [0] * q; best = [0]
    def go(d, x, usedV, usedC, cnt):
        if d == q:
            c = ADD[P[x]][0]; t = cnt + (not usedC >> c & 1)
            best[0] = max(best[0], t)
            return t == q
        if want_max and cnt + (q - d) + 1 <= best[0]: return False
        for y in range(q):
            if usedV >> y & 1: continue
            c = ADD[P[x]][y]; new = not usedC >> c & 1
            if not new and not want_max: continue
            A[d] = y
            if go(d + 1, y, usedV | 1 << y, usedC | 1 << c, cnt + new): return True
        return False
    ok = go(1, 0, 1, 0, 0)
    return (True, A[:]) if ok else (False, best[0])
def complete_mappings(q, ADD):
    R = range(q)
    return [u for u in itertools.permutations(R) if len({ADD[u[x]][x] for x in R}) == q]
def affine(u, q, ADD):
    return all(ADD[u[ADD[x][y]]][u[0]] == ADD[u[x]][u[y]] for x in range(q) for y in range(q))
def poly(P, q, ADD, MUL, NEG, names):
    """coefficients c_j (j = q-1..0) of the interpolating polynomial: c_0 = P(0), c_j = -sum_a P(a) a^(q-1-j)"""
    c = [P[0]] + [0] * (q - 1)
    for j in range(1, q):
        s = 0
        for a in range(q): s = ADD[s][MUL[P[a]][power(MUL, a, q - 1 - j)]]
        c[j] = NEG[s]
    terms = [('' if c[j] == 1 and j else names[c[j]] + ('*' if j else '')) + ('x^%d' % j if j > 1 else 'x' if j == 1 else '') for j in range(q - 1, -1, -1) if c[j]]
    return ' + '.join(terms) or '0'

# ---------- [1] + [2] q = 4 ----------
q, ADD, MUL, NEG, low = field(2, 2); w = 2                      # F_4 = F_2[t]/(t^2+t+1), w = t
assert low == (1, 1)
cm = complete_mappings(q, ADD)
assert len(cm) == 8 and all(sign(u) == 1 for u in cm)
bad = [P for P in itertools.permutations(range(q)) if not works(P, ADD, q)[0]]
assert sorted(bad) == sorted(P for P in itertools.permutations(range(q)) if sign(P) == 1)
lin = [tuple(ADD[MUL[a][x]][c] for x in range(q)) for a in (2, 3) for c in range(q)]
assert sorted(bad) == sorted(lin + [tuple(ADD[x][c] for x in range(q)) for c in range(q)])
say('q=4: complete mappings 8, all even; failing bijective P = the 12 even permutations = x+c (excluded) and w x + c, w^2 x + c')

# ---------- [3] q = 8 ----------
q, ADD, MUL, NEG, low = field(2, 3); assert low == (1, 0, 1)   # F_8 = F_2[t]/(t^3 + t^2 + 1)
cm = complete_mappings(q, ADD)
assert len(cm) == 384 and all(affine(u, q, ADD) and sign(u) == 1 for u in cm)
cycles = []
for rest in itertools.permutations(range(1, q)):
    s = [0] * q; seq = (0,) + rest
    for i in range(q): s[seq[i]] = seq[(i + 1) % q]
    cycles.append(s)
good = set()
for u in cm:
    ui = [0] * q
    for x in range(q): ui[u[x]] = x
    for s in cycles: good.add(tuple(ui[s[x]] for x in range(q)))       # P = u^-1 o sigma
odd = [P for P in itertools.permutations(range(q)) if sign(P) == -1]
assert len(good) == 20160 and all(P in good for P in odd)
for name, P in [('x^2', [power(MUL, x, 2) for x in range(q)]), ('x^4', [power(MUL, x, 4) for x in range(q)]),
                ('t x', [MUL[2][x] for x in range(q)]), ('x^5+x^3+x', [ADD[ADD[power(MUL, x, 5)][power(MUL, x, 3)]][x] for x in range(q)]),
                ('x^6+x^3+x^2', [ADD[ADD[power(MUL, x, 6)][power(MUL, x, 3)]][power(MUL, x, 2)] for x in range(q)])]:
    assert len(set(P)) == q and sign(P) == 1
    ok, mx = works(P, ADD, q, want_max=True)
    assert not ok
    say('q=8: P = %-12s even permutation, no circular order (max %d distinct sums of 8)' % (name, mx))
say('q=8: complete mappings 384, all affine and even; failing bijective P = exactly the 20160 even permutations')

# ---------- [4] q = 9 ----------
q, ADD, MUL, NEG, low = field(3, 2); assert low == (1, 0)      # F_9 = F_3[i]/(i^2+1); element a + b i = a + 3b
cm = complete_mappings(q, ADD)
oddcm = [u for u in cm if sign(u) == -1]
assert len(cm) == 2241 and len(oddcm) == 108
def linpart(u):                       # u affine: x -> L(x) + u(0)
    return tuple(ADD[u[x]][NEG[u[0]]] for x in range(q))
def order(L):
    M = list(L); n = 1
    while M != list(range(q)): M = [L[M[x]] for x in range(q)]; n += 1
    return n
assert all(affine(u, q, ADD) and order(linpart(u)) == 8 for u in oddcm)
AGL = []
for m in itertools.product(range(3), repeat=4):
    L = [((m[0] * (x % 3) + m[1] * (x // 3)) % 3) + 3 * ((m[2] * (x % 3) + m[3] * (x // 3)) % 3) for x in range(q)]
    if len(set(L)) == q: AGL += [[ADD[L[x]][c] for x in range(q)] for c in range(q)]
assert len(AGL) == 432 and not any(is_cycle(f) for f in AGL)
cycles = []
for rest in itertools.permutations(range(1, q)):
    s = [0] * q; seq = (0,) + rest
    for i in range(q): s[seq[i]] = seq[(i + 1) % q]
    cycles.append(s)
good = set()
for u in oddcm:
    ui = [0] * q
    for x in range(q): ui[u[x]] = x
    for s in cycles: good.add(tuple(ui[s[x]] for x in range(q)))
assert len(good) == 181440 - 648                           # odd P that work; 648 odd P fail
i = 3; names = {0: '0', 1: '1', 2: '-1', 3: 'i', 6: '-i', 4: '(1+i)', 5: '(-1+i)', 7: '(1-i)', 8: '(-1-i)'}
X = lambda e: [power(MUL, x, e) for x in range(q)]
reps = [('x^3', X(3)), ('-x^3', [NEG[v] for v in X(3)]), ('i x^3', [MUL[i][v] for v in X(3)]), ('(1+i)x', [MUL[4][x] for x in range(q)]),
        ('(-1+i)x', [MUL[5][x] for x in range(q)]),
        ('x^6+x^4+x^3+x^2', [ADD[ADD[ADD[a][b]][c]][d] for a, b, c, d in zip(X(6), X(4), X(3), X(2))]),
        ('x^6+x^4-x^3+x^2', [ADD[ADD[ADD[a][b]][NEG[c]]][d] for a, b, c, d in zip(X(6), X(4), X(3), X(2))])]
for name, P in reps:
    assert len(set(P)) == q and sign(P) == -1 and tuple(P) not in good
    ok, mx = works(P, ADD, q, want_max=True)
    assert not ok
    say('q=9: P = %-16s odd permutation, no circular order (max %d distinct sums of 9)' % (name, mx))
assert works(X(7), ADD, q)[0] and works(X(5), ADD, q)[0] and works([MUL[3][x] for x in range(q)], ADD, q)[0]
say('q=9: complete mappings 2241, odd ones = 108 affine (order-8 linear part); 648 odd P fail; x^5, x^7, i x work')

# ---------- [5] all P with sum 0, q <= 8, up to conjugation by AGL and adding constants ----------
def classify(p, k):
    q, ADD, MUL, NEG, low = field(p, k)
    G = []
    for m in itertools.product(range(p), repeat=k * k):
        dig = lambda x: [(x // p ** j) % p for j in range(k)]
        L = [sum(((sum(m[r * k + c] * dig(x)[c] for c in range(k))) % p) * p ** r for r in range(k)) for x in range(q)]
        if len(set(L)) == q:
            for v in range(q):
                g = [ADD[L[x]][v] for x in range(q)]; gi = [0] * q
                for x in range(q): gi[g[x]] = x
                G.append((g, gi))
    pw = [q ** x for x in range(q)]
    seen = bytearray(q ** q); orbits = 0; fails = []
    for head in itertools.product(range(q), repeat=q - 1):
        s = 0
        for v in head: s = ADD[s][v]
        P = head + (NEG[s],)
        code = sum(P[x] * pw[x] for x in range(q))
        if seen[code]: continue
        orbits += 1; size = 0
        for g, gi in G:
            Q = [g[P[gi[x]]] for x in range(q)]
            for c in range(q):
                cd = 0
                for x in range(q): cd += ADD[Q[x]][c] * pw[x]
                if not seen[cd]: seen[cd] = 1; size += 1
        if not works(P, ADD, q)[0]: fails.append((P, size))
    return q, orbits, fails
expect = {2: (1, 0, 0), 3: (3, 3, 0), 4: (6, 4, 8), 5: (16, 5, 0), 7: (477, 7, 0), 8: (312, 8, 20152)}
for p, k in [(2, 1), (3, 1), (2, 2), (5, 1), (7, 1), (2, 3)]:
    q, orbits, fails = classify(p, k)
    qq, ADD, MUL, NEG, low = field(p, k)
    excl = sum(sz for P, sz in fails if len({ADD[P[x]][x] for x in range(q)}) == 1)
    other = sum(sz for P, sz in fails if len({ADD[P[x]][x] for x in range(q)}) > 1)
    assert all(len(set(P)) == q and sign(P) == (1 if q in (4, 8) else -1 if q == 9 else sign(P)) for P, sz in fails if len({ADD[P[x]][x] for x in range(q)}) > 1)
    assert (orbits, excl, other) == expect[q], (q, orbits, excl, other)
    say('q=%d: %d orbits; failing functions: %d of the form c-x, %d others' % (q, orbits, excl, other))

# ---------- [6] witnesses ----------
here = os.path.dirname(os.path.abspath(__file__)); cache = {}; n = 0
for line in open(os.path.join(here, 'a229827_witnesses.txt')):
    if line.startswith('#') or not line.strip(): continue
    p, k, typ, par, cyc = line.split(); p, k, par = int(p), int(k), int(par)
    if (p, k) not in cache: cache[(p, k)] = field(p, k)
    q, ADD, MUL, NEG, low = cache[(p, k)]
    P = [power(MUL, x, par) for x in range(q)] if typ == 'mono' else [MUL[par][x] for x in range(q)]
    if typ == 'lin':
        o, y = 1, par
        while y != 1: y = MUL[y][par]; o += 1
        assert o == q - 1                                   # lambda primitive
    a = list(map(int, cyc.split(',')))
    assert sorted(a) == list(range(q)) and len({ADD[P[a[j]]][a[(j + 1) % q]] for j in range(q)}) == q
    n += 1
mono = {(p, k): set() for (p, k) in cache}
for line in open(os.path.join(here, 'a229827_witnesses.txt')):
    if line.startswith('#') or not line.strip(): continue
    s = line.split()
    if s[2] == 'mono': mono[(int(s[0]), int(s[1]))].add(int(s[3]))
for (p, k), es in mono.items():
    q = p ** k
    missing = sorted(set(range(q - 1)) - es - ({1} if p == 2 else set()))
    assert missing == {(2, 3): [2, 4], (3, 2): [3]}.get((p, k), []), (q, missing)
say('witnesses: %d circular orders verified; every x^k (k <= q-2, q in {4,8,9,16,25,27,32,49,64,81}) works except x^2, x^4 on F_8 and x^3 on F_9' % n)

# ---------- [7] (ii), known ----------
assert all(len({2 * c[j] + c[(j + 1) % 3] for j in range(3)}) < 3 for c in [(0, 1, 2), (0, 2, 1)])
say('(ii): {0,1,2}, P = 2x fails (known, Sun 1309.1679)')
