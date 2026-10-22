
import itertools, sys, time, random
from fractions import Fraction as Fr
from math import lcm, comb
sys.path.insert(0, '.')
from qc_check import build_qc, is_mono, is_rainbow

PR = (1 << 61) - 1

def rank_mod_p(rows, p=PR):
    rows = [r[:] for r in rows]
    rk, ncol = 0, (len(rows[0]) if rows else 0)
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % p), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], p - 2, p)
        rows[rk] = [(v * inv) % p for v in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % p:
                f = rows[i][c]
                rows[i] = [(a - f * b) % p for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def nullspace_Q(rows, ncol):
    """rows: list of lists of ints (matrix A, len(rows) x ncol). Return basis of {h: A h = 0} over Q (as Fractions)."""
    A = [[Fr(v) for v in r] for r in rows]
    piv_cols, r = [], 0
    for c in range(ncol):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]; A[r] = [v / pv for v in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv_cols.append(c); r += 1
    free = [c for c in range(ncol) if c not in piv_cols]
    basis = []
    for fc in free:
        h = [Fr(0)] * ncol; h[fc] = Fr(1)
        for i, pc in enumerate(piv_cols): h[pc] = -A[i][fc]
        basis.append(h)
    return basis

def P_monomials(k, D, n):
    smax = int((1 - 1 / (k - 1)) * n + 1e-12)
    return [d for d in itertools.product(range(D), repeat=n) if sum(e >= k - 1 for e in d) <= smax]

def is_sunflower(vs):
    n = len(vs[0])
    return len(set(vs)) == len(vs) and all(is_mono([v[r] for v in vs]) or is_rainbow([v[r] for v in vs]) for r in range(n))

def greedy_sf(k, D, n, tries, seed=1):
    rnd = random.Random(seed)
    pts = list(itertools.product(range(1, D + 1), repeat=n))
    best = []
    for _ in range(tries):
        rnd.shuffle(pts)
        V = []
        for x in pts:
            if not any(is_sunflower((x,) + c) for c in itertools.combinations(V, k - 1)):
                V.append(x)
        if len(V) > len(best): best = V[:]
    return best

def exact_max_sf(k, D, n):
    pts = list(itertools.product(range(1, D + 1), repeat=n))
    idx = {p: i for i, p in enumerate(pts)}
    edges = [sum(1 << idx[p] for p in c) for c in itertools.combinations(pts, k) if is_sunflower(c)]
    best = 0
    for m in range(1 << len(pts)):
        if bin(m).count('1') > best and not any((m & e) == e for e in edges):
            best = bin(m).count('1')
    return best, len(edges)

def contraction_test(k, D, n, q, nh=3, seed=7):
    """V=[D]^n (not sunflower-free). For h in P^perp, S(a)=sum T(a,x1..x_{k-1}) prod h(x_i) must be 0 for every a."""
    L = 1
    for v in q.values(): L = lcm(L, v.denominator)
    qi = {a: int(v * L) for a, v in q.items()}
    pts = list(itertools.product(range(1, D + 1), repeat=n))
    mons = P_monomials(k, D, n)
    A = [[prod_pow(x, d) for x in pts] for d in mons]  # dimP x |V|
    basis = nullspace_Q(A, len(pts))
    rnd = random.Random(seed)
    worst = 0
    for t in range(nh):
        cf = [rnd.randint(-5, 5) for _ in basis]
        if not any(cf): cf[0] = 1
        h = [sum(c * b[i] for c, b in zip(cf, basis)) for i in range(len(pts))]
        assert all(sum(Fr(A[r][i]) * h[i] for i in range(len(pts))) == 0 for r in range(len(A)))
        Lh = 1
        for v in h: Lh = lcm(Lh, v.denominator)
        hi = [int(v * Lh) for v in h]
        nz = [(x, hv) for x, hv in zip(pts, hi) if hv]
        for a in pts:
            S = 0
            for combo in itertools.product(nz, repeat=k - 1):
                val = 1
                for r in range(n):
                    val *= qi[(a[r],) + tuple(c[0][r] for c in combo)]
                    if not val: break
                if val:
                    w = val
                    for c in combo: w *= c[1]
                    S += w
            if S != 0: worst += 1
    return len(basis), worst

def prod_pow(x, d):
    p = 1
    for xi, di in zip(x, d): p *= xi ** di
    return p

if __name__ == "__main__":
    t0 = time.time()
    k = 4
    print("== dim P (formula) vs rank of monomial evaluation matrix on [D]^n (mod 2^61-1; rank_p<=rank_Q<=count) ==")
    for D in (4, 5, 6):
        for n in (1, 2, 3):
            mons = P_monomials(k, D, n)
            pts = list(itertools.product(range(1, D + 1), repeat=n))
            smax = int((1 - 1 / (k - 1)) * n + 1e-12)
            formula = sum(comb(n, s) * (D - k + 1) ** s * (k - 1) ** (n - s) for s in range(smax + 1))
            rk = rank_mod_p([[prod_pow(x, d) % PR for x in pts] for d in mons])
            print(f"k=4 D={D} n={n}: formula={formula} count={len(mons)} rank={rk} D^n={D**n}")
    print("== exact V_4(4,2) by exhaustive search ==")
    best, ne = exact_max_sf(4, 4, 2)
    print(f"V_4(4,2)={best} (4-sunflowers in [4]^2: {ne}); dim P(4,2)=15")
    print("== greedy sunflower-free V, Lemma 2.1 conclusion rank(P|_V)=|V| <= dim P ==")
    for (D, n, tries) in ((4, 2, 50), (5, 2, 50), (4, 3, 6), (5, 3, 3)):
        V = greedy_sf(4, D, n, tries)
        assert not any(is_sunflower(c) for c in itertools.combinations(V, 4))
        mons = P_monomials(4, D, n)
        rk = rank_mod_p([[prod_pow(x, d) % PR for x in V] for d in mons])
        # diagonality of T_c on V^4 (all 4-tuples incl. repeats)
        q, _ = build_qc(4, D)
        diag = True
        for tup in itertools.product(V, repeat=4) if len(V) ** 4 <= 3_000_000 else []:
            val = 1
            for r in range(n):
                val *= q[tuple(t[r] for t in tup)]
                if val == 0: break
            if (val != 0) != (tup[0] == tup[1] == tup[2] == tup[3]): diag = False; break
        print(f"D={D} n={n}: |V|={len(V)} rank(P|_V)={rk} dimP={len(mons)} T_c diagonal on V^4={diag if len(V)**4<=3_000_000 else 'skipped'} [{time.time()-t0:.0f}s]")
    print("== contraction test on V=[D]^n: S(a)=sum_{x1..x3} T_c(a,x1,x2,x3) h(x1)h(x2)h(x3) for h in P^perp ==")
    for (D, n) in ((4, 2), (5, 2)):
        q, _ = build_qc(4, D)
        dperp, bad = contraction_test(4, D, n, q)
        naive = {a: Fr(1 if (is_mono(a) or is_rainbow(a)) else 0) for a in q}
        dperp2, bad2 = contraction_test(4, D, n, naive)
        print(f"D={D} n={n}: dim P^perp={dperp}; interpolated q_c: nonzero S(a) count={bad} (expect 0); naive q (c=1): nonzero S(a) count={bad2} (control) [{time.time()-t0:.0f}s]")

