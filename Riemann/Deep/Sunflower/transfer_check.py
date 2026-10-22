
import itertools, math, time
from fractions import Fraction as Fr
from math import comb, factorial

def set_sunflower(sets):
    sets = [frozenset(s) for s in sets]
    if len(set(sets)) < len(sets): return False
    inters = {a & b for a, b in itertools.combinations(sets, 2)}
    return len(inters) == 1

def vec_sunflower(vs):
    if len(set(vs)) < len(vs): return False
    return all(len({v[r] for v in vs}) in (1, len(vs)) for r in range(len(vs[0])))

def correspondence(k, D, t, rsize, extra=1):
    # ground set: blocks B_i = {(i,j)}, R = {('R',m)}, plus unused extra elements
    R = [('R', m) for m in range(rsize)]
    G = {}
    for v in itertools.product(range(D), repeat=t):
        G[v] = frozenset(R) | {(i, v[i]) for i in range(t)}
    vecs = list(G)
    bad_fwd = bad_bwd = cnt = 0
    for c in itertools.combinations(vecs, k):
        vs, ss = vec_sunflower(c), set_sunflower([G[v] for v in c])
        cnt += vs
        if vs and not ss: bad_fwd += 1
        if ss and not vs: bad_bwd += 1
    return len(vecs), cnt, bad_fwd, bad_bwd

def symmetry(N, D, s):
    t = min(s, (N - s) // (D - 1))
    counts = {S: 0 for S in itertools.combinations(range(N), s)}
    for perm in itertools.permutations(range(N)):
        blocks = [perm[i * D:(i + 1) * D] for i in range(t)]
        R = perm[D * t: D * t + (s - t)]
        for choice in itertools.product(*blocks):
            counts[tuple(sorted(R + choice))] += 1
    vals = set(counts.values())
    expect = Fr(factorial(N) * D ** t, comb(N, s))
    return t, vals, expect

if __name__ == "__main__":
    t0 = time.time()
    print("== (a) block transfer: vector k-sunflower in G <=> set k-sunflower (all k-subsets of G) ==")
    for (k, D, t, rs) in ((4, 4, 2, 1), (4, 5, 2, 0), (4, 4, 3, 2), (5, 5, 2, 1), (3, 4, 3, 1)):
        nG, cnt, bf, bb = correspondence(k, D, t, rs)
        print(f"k={k} D={D} t={t} |R|={rs}: |G|={nG} vector-sunflowers={cnt} vec-not-set={bf} set-not-vec={bb}")
    print("== (b) symmetry: #perms with S in G is the same for every s-set S, = N! D^t / C(N,s) ==")
    for (N, D, s) in ((7, 2, 3), (8, 3, 3), (7, 3, 2)):
        t, vals, expect = symmetry(N, D, s)
        print(f"N={N} D={D} s={s} t={t}: counts={sorted(vals)} expected={expect}")
    print("== (c) rho^t <= rho^s + rho^(-1+(N-s)/(D-1)), t=min(s,floor((N-s)/(D-1))); sum identity ==")
    worst = 0.0; idbad = 0
    for rho in (0.5, 0.63, 0.8, 0.9, 0.99):
        for D in range(2, 21):
            for N in range(1, 81):
                tot = 0.0
                for s in range(N + 1):
                    t = min(s, (N - s) // (D - 1))
                    lhs = rho ** t; rhs = rho ** s + rho ** (-1 + (N - s) / (D - 1))
                    worst = max(worst, lhs / rhs)
                    tot += comb(N, s) * rhs
                closed = (1 + rho) ** N + (1 / rho) * (1 + rho ** (1 / (D - 1))) ** N
                if abs(tot - closed) > 1e-9 * closed: idbad += 1
                if closed > (1 + 1 / rho) * (1 + rho ** (1 / (D - 1))) ** N * (1 + 1e-12): idbad += 1
    print(f"max lhs/rhs = {worst:.6f} (must be <= 1); identity/final-inequality failures = {idbad}")
    print("== (d) constants for D=16k^2, rho=4^(-1/(k-1)) (exact rational bounds) ==")
    LN4_LO = Fr(1386, 1000)  # ln 4 = 1.38629... > 1.386
    bad = 0
    for k in range(3, 401):
        D = 16 * k * k
        # rho in [1/2,1): 4^(-1/(k-1)) >= 4^(-1/2) = 1/2 for k>=3
        if not (1 + k * (D - 1) <= 16 * k ** 3): bad += 1
        if not (LN4_LO / (k - 1) >= Fr(1, k)): bad += 1
        x = LN4_LO / ((k - 1) * (D - 1))  # x <= true exponent; e^{-x_true} <= 1 - x + x^2/2
        if not (1 - x + x * x / 2 <= 1 - Fr(1, 16 * k ** 3)): bad += 1
    print(f"k=3..400: failures = {bad}")
    print("== (e) Chernoff step, exact integers: S^(k-1) <= (4k^2)^n D^((k-2)n), S = sum_(s<=floor((k-2)n/(k-1))) C(n,s)(D-k+1)^s(k-1)^(n-s) ==")
    for k in (3, 4, 5, 6, 8):
        for D in (4 * k * k + 1, 8 * k * k, 16 * k * k):
            fails = 0; tight = 0.0
            for n in range(1, 151):
                smax = (k - 2) * n // (k - 1)
                S = sum(comb(n, s) * (D - k + 1) ** s * (k - 1) ** (n - s) for s in range(smax + 1))
                lhs, rhs = S ** (k - 1), (4 * k * k) ** n * D ** ((k - 2) * n)
                if lhs > rhs: fails += 1
                tight = max(tight, math.exp((math.log(S) * (k - 1) - math.log(rhs)) / ((k - 1) * n)))
            print(f"k={k} D={D}: n=1..150 failures={fails}; max per-coordinate ratio (S/bound)^(1/n)={tight:.4f}")
    print("== (f) k=4 numbers ==")
    k = 4; D = 16 * k * k; C = (4 * k * k) ** (1 / (k - 1)) * D ** (1 - 1 / (k - 1)); rho = C / D
    print(f"D={D} C={C:.3f} rho={rho:.6f} 1+rho^(1/(D-1))={1 + rho ** (1 / (D - 1)):.8f} 2-1/(16k^3)={2 - 1 / (16 * k ** 3):.8f}")
    print(f"[{time.time() - t0:.0f}s]")

