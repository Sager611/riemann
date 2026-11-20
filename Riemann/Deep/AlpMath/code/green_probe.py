# Probe: Green's parametrisation of cuspidal characters of GL_d(F_g), as stated in
# Development/AbelianVarietyTheory/PublishedTheorems/GaloisCohomology₂/GreenRegularCharacters.lean:13
# (on the cone via green_regular_characters_of_literature, EveryCharacterOfLiteratureOfVWin.lean:399,
#  GlobalTheory.GreenBridge.green_of_regular in SiegelEisensteinSecondTermHypotheses.lean:72).
# Independent check: Burnside's algorithm character table of GL_d(F_g); cuspidality via unipotent radicals U_a.
import sys, itertools, numpy as np, cmath
np.random.seed(1)

def run(p, n):
    N2 = n * n
    W = p ** np.arange(N2)
    # enumerate invertible matrices
    allm = np.array(list(itertools.product(range(p), repeat=N2)), dtype=np.int64).reshape(-1, n, n)
    det = np.round(np.linalg.det(allm.astype(float))).astype(np.int64) % p
    G = allm[det != 0]
    del allm
    order = len(G)
    keys = (G.reshape(order, N2) * W).sum(1)
    look = -np.ones(p ** N2, dtype=np.int64); look[keys] = np.arange(order)
    def idx(A):  # A: (..., n, n)
        return look[((A.reshape(-1, N2) % p) * W).sum(1)]
    # inverses: brute force via products with all? use solve: for each g, find h with g h = I
    eye_key = (np.eye(n, dtype=np.int64).reshape(N2) * W).sum()
    inv = np.empty(order, dtype=np.int64)
    # inverse by Gaussian elimination mod p (vectorised loop over elements in python, small)
    def minv(A):
        A = [list(map(int, r)) + [1 if i == j else 0 for j in range(n)] for i, r in enumerate(A)]
        for c in range(n):
            piv = next(r for r in range(c, n) if A[r][c] % p)
            A[c], A[piv] = A[piv], A[c]
            iv = pow(A[c][c], p - 2, p)
            A[c] = [(x * iv) % p for x in A[c]]
            for r in range(n):
                if r != c and A[r][c] % p:
                    f = A[r][c]; A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
        return np.array([row[n:] for row in A], dtype=np.int64)
    Ginv = np.empty_like(G)
    for t in range(order):
        Ginv[t] = minv(G[t])
    inv = idx(Ginv)
    # conjugacy classes
    cls = -np.ones(order, dtype=np.int64); reps = []; sizes = []
    for t in range(order):
        if cls[t] >= 0: continue
        orb = np.unique(idx(G @ G[t] @ Ginv))
        cls[orb] = len(reps); reps.append(t); sizes.append(len(orb))
    k = len(reps); sizes = np.array(sizes)
    # class multiplication coefficients c[i, j, kk] = #{x in C_i : x^{-1} z_kk in C_j}
    c = np.zeros((k, k, k), dtype=np.int64)
    for kk, r in enumerate(reps):
        y = idx(Ginv @ G[r])
        np.add.at(c, (cls, cls[y], kk), 1)
    # Burnside: w eigenvector of A_j with (A_j)[i, kk] = c[j, i, kk]
    A = np.tensordot(np.random.randn(k), c.astype(float), axes=(0, 0))
    ev, V = np.linalg.eig(A)
    id_cls = cls[look[eye_key]]
    chars = []
    for col in range(k):
        w = V[:, col] / V[id_cls, col]
        deg2 = order / np.sum(np.abs(w) ** 2 / sizes)
        deg = np.sqrt(deg2.real)
        chi = w * deg / sizes
        chars.append(chi)
    chars = np.array(chars)
    # sanity: orthonormality
    gram = (chars * sizes) @ chars.conj().T / order
    orth_err = np.abs(gram - np.eye(k)).max()
    degs = np.round(chars[:, id_cls].real).astype(int)
    # cuspidal: sum over U_a of chi = 0 for 0<a<n
    cusp = np.ones(k, dtype=bool)
    for a in range(1, n):
        Us = []
        for B in itertools.product(range(p), repeat=a * (n - a)):
            M = np.eye(n, dtype=np.int64)
            M[:a, a:] = np.array(B).reshape(a, n - a)
            Us.append(M)
        uc = cls[idx(np.array(Us))]
        s = chars[:, uc].sum(1)
        cusp &= np.abs(s) < 1e-6 * len(Us)
    # field K = F_{p^n}: find irreducible monic m of degree n, elements as coefficient tuples
    q = p ** n
    def polymulmod(a, b, m):
        r = [0] * (2 * n - 1)
        for i in range(n):
            for j in range(n):
                r[i + j] = (r[i + j] + a[i] * b[j]) % p
        for d in range(2 * n - 2, n - 1, -1):
            f = r[d]
            if f:
                for i in range(n + 1):
                    r[d - n + i] = (r[d - n + i] - f * m[i]) % p
        return tuple(r[:n])
    elems = list(itertools.product(range(p), repeat=n))
    def is_irred(m):
        # no root of m in any smaller extension: test x^(p^e) = x mod m only for e = n, via gcd-free brute force: degree <= 4
        # brute force: m irreducible iff it has no monic factor of degree 1..n//2
        for dfac in range(1, n // 2 + 1):
            for cf in itertools.product(range(p), repeat=dfac):
                f = list(cf) + [1]
                # polynomial division m / f
                r = list(m)
                for d in range(n, dfac - 1, -1):
                    co = r[d]
                    if co:
                        for i in range(dfac + 1):
                            r[d - dfac + i] = (r[d - dfac + i] - co * f[i]) % p
                if all(x == 0 for x in r[:dfac]):
                    return False
        return True
    m = next(list(cf) + [1] for cf in itertools.product(range(p), repeat=n) if is_irred(list(cf) + [1]))
    E = {e: i for i, e in enumerate(elems)}
    mul = [[E[polymulmod(a, b, m)] for b in elems] for a in elems]
    zero = E[tuple([0] * n)]; one = E[tuple([1] + [0] * (n - 1))]
    def kpow(x, e):
        r = one
        for _ in range(e): r = mul[r][x]
        return r
    # generator of K^x
    def korder(x):
        r = x; o = 1
        while r != one: r = mul[r][x]; o += 1
        return o
    gen = next(x for x in range(len(elems)) if x != zero and korder(x) == q - 1)
    dlog = {}; r = one
    for e in range(q - 1): dlog[r] = e; r = mul[r][gen]
    # embed F_p in K: scalar a -> element (a,0,...)
    sc = [E[tuple([a] + [0] * (n - 1))] for a in range(p)]
    add = [[E[tuple((x + y) % p for x, y in zip(a, b))] for b in elems] for a in elems]
    neg = [E[tuple((-x) % p for x in a)] for a in elems]
    def kdet(Mx):  # Leibniz in K
        tot = zero
        for perm in itertools.permutations(range(n)):
            sgn = 1
            for i in range(n):
                for j in range(i + 1, n):
                    if perm[i] > perm[j]: sgn = -sgn
            t = one
            for i in range(n): t = mul[t][Mx[i][perm[i]]]
            tot = add[tot][t if sgn == 1 else neg[t]]
        return tot
    subfield = set()
    for e in range(1, n):
        if n % e == 0:
            for x in range(len(elems)):
                if x != zero and kpow(x, p ** e) == x: subfield.add(x)
    ell = {}  # class index -> a root z (dlog) of charpoly, for elliptic regular classes
    for ci, r in enumerate(reps):
        y = G[r]
        roots = []
        for z in range(len(elems)):
            if z == zero: continue
            Mx = [[add[z if i == j else zero][neg[sc[int(y[i, j])]]] for j in range(n)] for i in range(n)]
            if kdet(Mx) == zero: roots.append(z)
        if len(roots) == n and not any(z in subfield for z in roots):
            ell[ci] = dlog[roots[0]]
    regs = [kk for kk in range(q - 1) if len({(kk * p ** i) % (q - 1) for i in range(n)}) == n]
    orbits = sorted({min((kk * p ** i) % (q - 1) for i in range(n)) for kk in regs})
    ecls = sorted(ell)
    pred = []
    for o in orbits:
        v = [(-1) ** (n - 1) * sum(cmath.exp(2j * cmath.pi * o * ell[ci] * p ** i / (q - 1)) for i in range(n)) for ci in ecls]
        pred.append(np.array(v))
    cus = [i for i in range(k) if cusp[i]]
    target_deg = int(np.prod([p ** i - 1 for i in range(1, n)]))
    deg_ok = all(degs[i] == target_deg for i in cus)
    # match predicted vectors with cuspidal characters' elliptic values
    used = set(); match_ok = True
    for v in pred:
        hit = [i for i in cus if i not in used and np.abs(chars[i, ecls] - v).max() < 1e-6]
        if not hit: match_ok = False; break
        used.add(hit[0])
    print(f"GL_{n}(F_{p}): |G|={order} classes={k} orth_err={orth_err:.1e} cuspidal={len(cus)} regular_orbits={len(orbits)} "
          f"degree_ok={deg_ok} (target {target_deg}) elliptic_classes={len(ecls)} values_match={match_ok} "
          f"PASS={len(cus)==len(orbits) and deg_ok and match_ok and orth_err<1e-6}")
    sys.stdout.flush()

for (p, n) in [(2, 2), (3, 2), (5, 2), (7, 2), (2, 3), (3, 3), (2, 4), (11, 2)]:
    run(p, n)
