r"""Random tests of the two reduction lemmas of OpenAI's Seymour paper (Section 1),
on arbitrary oriented graphs (the lemmas' internal steps do not need a counterexample).

Lemma 1.2 (subset deficit), per-step claims for any D, S, T subsetneq E and v in S
with Q = F(v) & (E \ T) nonempty, where D_T deletes all arcs S -> E \ T:
  (a) F_{D_T}^2(v) subset I | T | F(T)
  (b) newly counted exact second out-neighbours of v lie in Q & F(T)
  (c) C = F(Q) \ (I|T|Q|F(T)) subset N2_D(v) and C misses F_{D_T}^2(v)
  (d) |N2_{D_T}(v)| <= d2 - |C| + |Q & F(T)|
  (e) g(T|Q) - g(T) = |C| - |Q & F(T)|   (exact)
  (f) vertices that lose no out-neighbour have N2_{D_T} subset N2_D
  (g) F^2(S) \ F(S) = F(E) \ (I|E)
Lemma 1.3 (blow-up D[T_m]), for D with positive indegrees, any m >= 1, U:
  (h) |F_G U| = m|F(S)| + sum_{v in S\F(S)} |f(T_v)|
  (i) |U| - 2|F_G U| + |F_G^2 U| <= m(|F^2 S \ F S| - |F S \ S|) + |S \ F S|  (S proper)
  (j) S = V(D)  =>  F_G U = F_G^2 U = V(G)
"""
import random, sys

def pc(x): return bin(x).count('1')

def rand_oriented(n, p, rng):
    out = [0]*n
    for a in range(n):
        for b in range(a+1, n):
            if rng.random() < p:
                if rng.random() < 0.5: out[a] |= 1 << b
                else: out[b] |= 1 << a
    return out

def F(out, U):
    r = 0
    while U:
        lb = U & -U; v = lb.bit_length()-1; U ^= lb; r |= out[v]
    return r

def N2(out, v):
    f1 = out[v]; return F(out, f1) & ~f1 & ~(1 << v)

def test_deficit(trials, rng):
    cnt = 0
    for _ in range(trials):
        n = rng.randint(3, 9); out = rand_oriented(n, rng.uniform(0.2, 0.9), rng)
        full = (1 << n) - 1
        S = rng.randint(1, full-1)
        FS = F(out, S); I = FS & S; E = FS & ~S
        if E == 0: continue
        # (g)
        assert (F(out, FS) & ~FS) == (F(out, E) & ~(I | E)), "g"
        # random T subsetneq E
        while True:
            T = rng.randint(0, E) & E
            if T != E: break
        outT = [ (o & ~(E & ~T)) if (S >> v) & 1 else o for v, o in enumerate(out)]
        FT = F(out, T)
        g = lambda X: pc(F(out, X) & ~(I | X))
        for v in range(n):
            Q = out[v] & (E & ~T)
            if not ((S >> v) & 1) or Q == 0:
                # (f)
                assert (N2(outT, v) & ~N2(out, v)) == 0, "f"
                continue
            cnt += 1
            F2T = F(outT, outT[v])
            assert (F2T & ~(I | T | FT)) == 0, "a"
            new = N2(outT, v) & ~N2(out, v)
            assert (new & ~(Q & FT)) == 0, "b"
            C = F(out, Q) & ~(I | T | Q | FT)
            assert (C & ~N2(out, v)) == 0 and (C & F2T) == 0, "c"
            d2 = pc(N2(out, v))
            assert pc(N2(outT, v)) <= d2 - pc(C) + pc(Q & FT), "d"
            assert g(T | Q) - g(T) == pc(C) - pc(Q & FT), "e"
    return cnt

def test_blowup(trials, rng):
    cnt = 0
    for _ in range(trials):
        n = rng.randint(2, 5); out = rand_oriented(n, rng.uniform(0.4, 1.0), rng)
        indeg = [0]*n
        for a in range(n):
            for b in range(n):
                if (out[a] >> b) & 1: indeg[b] += 1
        if min(indeg) == 0: continue
        m = rng.randint(1, 4)
        N = n*m
        # vertex (v,k) -> index v*m+k ; internal arcs k<l ; external all arcs of D
        outG = [0]*N
        for v in range(n):
            for k in range(m):
                x = 0
                for l in range(k+1, m): x |= 1 << (v*m+l)
                for w in range(n):
                    if (out[v] >> w) & 1:
                        for l in range(m): x |= 1 << (w*m+l)
                outG[v*m+k] = x
        fullG = (1 << N) - 1
        for _ in range(30):
            U = rng.randint(1, fullG-1)
            Tv = [(U >> (v*m)) & ((1 << m)-1) for v in range(n)]
            S = 0
            for v in range(n):
                if Tv[v]: S |= 1 << v
            FS = F(out, S); F2S = F(out, FS)
            FU = F(outG, U); F2U = F(outG, FU)
            def f(T):
                if T == 0: return 0
                a = (T & -T).bit_length()-1
                return ((1 << m)-1) & ~((1 << (a+1))-1)
            pred = m*pc(FS) + sum(pc(f(Tv[v])) for v in range(n) if (S >> v) & 1 and not (FS >> v) & 1)
            assert pc(FU) == pred, "h"
            if S != (1 << n)-1:
                lhs = pc(U) - 2*pc(FU) + pc(F2U)
                rhs = m*(pc(F2S & ~FS) - pc(FS & ~S)) + pc(S & ~FS)
                assert lhs <= rhs, ("i", lhs, rhs)
            else:
                assert FU == fullG and F2U == fullG, "j"
            cnt += 1
    return cnt

if __name__ == "__main__":
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    rng = random.Random(seed)
    a = test_deficit(200000, rng)
    b = test_blowup(20000, rng)
    print(f"seed={seed} subset-deficit step instances={a} blowup instances={b}: all assertions (a)-(j) hold")
