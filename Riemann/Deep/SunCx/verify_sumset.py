"""Sun's restricted-sumset bounds: Wang-Sun, arXiv 2210.12044 = Electron. J. Combin. 33(2) (2026) #P2.64, Conjecture 1
(arXiv numbering: Conjecture 1.1) = OEIS A357130, Conjecture 2 (2022).  For finite A_1..A_n (n > 1, |A_i| > 1) in an additive
group G, with S = sum |A_i| and p(G) the least order of a nonzero element,
  L = {a_1+...+a_n : a_i in A_i, a_i != a_{i+1} (0<i<n)},          (L)  |L| >= min{p(G), S - 2n + 1 + (n mod 2)},
  C = {a_1+...+a_n : ... and also a_n != a_1},                      (C)  |C| >= min{p(G), S - 2n + (-1)^n (1 + (n mod 2))}.
Certified here (stdlib only, about 1 minute):
 [1] (L), n = 3, G = Z (or any group with p(G) >= 5): A = ({0,1,2},{0,2},{0,1,2}): L = {2,3,4}, 3 < 4.
     (Pattern (3,2,3) already stated as disproved: github.com/darkness-38/counter-hunter, paper/counter_hunter_verification.tex, 2026-10-17.)
 [2] (L), every odd n >= 3: n-3 copies of {0,1} followed by [1]: |L| = 3 < 4 (the prefix alternates, so L = (n-3)/2 + {2,3,4}).
 [3] (L), n = 3, G = Z/p, every prime p >= 5 (checked p <= 101): A = ({0,1}, {0,..,h}, {2} u {4,..,h}), h = (p+3)/2:
     S = p + 4, bound p, but 4 notin L (the integer sums lie in [2, p+4]; sum 4 and sum p+4 = 1+h+h are impossible), |L| = p-1.
 [4] (C), n = 4, G = Z/7: A = ({0,1},{0,1,2,4,5},{2,4,5},{0,1,4,5}): S = 14, bound 7, C = Z/7 minus {2}, |C| = 6 (integer sums miss 2, 9, 16).
     (Unique up to x -> ux+v, rotation and reversal among all 4-tuples in Z/7; none in Z/5 for n <= 6.)
 [5] exhaustive counts of violating tuples: Z/5 n=3: L 20, C 0; n=4: 0, 0; n=5: 470, 0; Z/7 n=3: 203, 0; n=4: 0, 336.
 Larger exhaustive runs (Z/5 n<=6, Z/7 n<=5, Z/11, Z/13 n=3, Z/11 n=4 with A_1={0,1}, groups with p(G)<=3, Z intervals,
 equal sizes): classify_sumset.c (cc -O2; table in NOVELTY.md).
Run: python3 verify_sumset.py
"""
import itertools, time
T0 = time.time()
def say(*a): print('[%5.1fs]' % (time.time() - T0), *a, flush=True)
def sums(As, mod=None, cyc=False):
    """L (cyc=False) or C (cyc=True) by dynamic programming over (first, last) elements"""
    n = len(As); out = set()
    firsts = As[0] if cyc else [None]
    for f in firsts:
        S = {a: {a} for a in As[0] if (f is None or a == f)}
        for i in range(1, n):
            T = {}
            for b in As[i]:
                acc = set()
                for a, s in S.items():
                    if a != b: acc |= {x + b for x in s}
                if acc: T[b] = acc
            S = T
        for b, s in S.items():
            if f is None or b != f: out |= s
    return {x % mod for x in out} if mod else out
def bound(As, cyc, p=None):
    n = len(As); s = sum(map(len, As))
    b = s - 2 * n + ((-1) ** n * (1 + n % 2) if cyc else 1 + n % 2)
    return min(p, b) if p else b
def brute(As, mod=None, cyc=False):
    n = len(As); out = set()
    for t in itertools.product(*As):
        if all(t[i] != t[i + 1] for i in range(n - 1)) and (not cyc or t[-1] != t[0]): out.add(sum(t) % mod if mod else sum(t))
    return out
# [1]
A = [{0, 1, 2}, {0, 2}, {0, 1, 2}]
assert sums(A) == brute(A) == {2, 3, 4} and bound(A, False) == 4
for p in [5, 7, 11, 13]: assert len(sums(A, p)) == 3 < bound(A, False, p) == 4
say('[1] Z and Z/p (p>=5): L({0,1,2},{0,2},{0,1,2}) = {2,3,4}, 3 < 4')
# [2]
for n in range(3, 22, 2):
    A = [{0, 1}] * (n - 3) + [{0, 1, 2}, {0, 2}, {0, 1, 2}]
    L = sums(A); assert L == {(n - 3) // 2 + v for v in (2, 3, 4)} and bound(A, False) == 4
    if n <= 9: assert brute(A) == L
say('[2] odd n = 3..21: n-3 copies of {0,1} then [1]: |L| = 3 < 4')
# [3]
primes = [p for p in range(5, 102) if all(p % d for d in range(2, p))]
for p in primes:
    h = (p + 3) // 2
    A = [{0, 1}, set(range(h + 1)), {2} | set(range(4, h + 1))]
    L = sums(A, p); assert sum(map(len, A)) == p + 4 and bound(A, False, p) == p and L == set(range(p)) - {4}
    if p <= 13: assert brute(A, p) == L
say('[3] Z/p, primes 5..101: ({0,1},{0..h},{2}u{4..h}), h=(p+3)/2: L = Z/p minus {4}, |L| = p-1 < p')
# [4]
A = [{0, 1}, {0, 1, 2, 4, 5}, {2, 4, 5}, {0, 1, 4, 5}]
C = sums(A, 7, cyc=True)
assert C == brute(A, 7, cyc=True) == set(range(7)) - {2} and bound(A, True, 7) == 7
assert sums(A, cyc=True) == brute(A, cyc=True) and len(sums(A, cyc=True)) >= bound(A, True)   # fine in Z: the failure is mod 7
say('[4] (C) Z/7, n=4: C({0,1},{0,1,2,4,5},{2,4,5},{0,1,4,5}) = Z/7 minus {2}, 6 < 7')
# [5]
def count(p, n):
    """number of n-tuples of subsets (|A_i| >= 2) of Z/p violating (L), resp. (C); A_1 runs over orbit representatives
    under x -> ux+v (u != 0), weighted by orbit size; sums kept as p-bit masks; last set handled by an OR over subsets"""
    full = (1 << p) - 1
    rot = lambda m, b: ((m << b) | (m >> (p - b))) & full
    subs = [m for m in range(1 << p) if bin(m).count('1') >= 2]
    elts = {m: [a for a in range(p) if m >> a & 1] for m in subs}
    seen, reps = set(), []
    for m in subs:
        if m in seen: continue
        orb = {sum(1 << ((u * a + v) % p) for a in elts[m]) for u in range(1, p) for v in range(p)}
        seen |= orb; reps.append((m, len(orb)))
    res = [0, 0]
    def step(prev, S, b):
        T = 0
        for a in elts[prev]:
            if a != b: T |= S[a]
        return rot(T, b)
    def rec(i, prev, S, SC, s, w):
        if i == n - 1:
            RL = [step(prev, S, b) for b in range(p)]
            RC = [0] * p
            for b in range(p):
                for f, Sf in SC.items():
                    if f != b: RC[b] |= step(prev, Sf, b)
            Lm = [0] * (1 << p); Cm = [0] * (1 << p)
            for m in range(1, 1 << p):
                lb = (m & -m).bit_length() - 1
                Lm[m] = Lm[m & (m - 1)] | RL[lb]; Cm[m] = Cm[m & (m - 1)] | RC[lb]
            for m in subs:
                t = s + len(elts[m])
                if bin(Lm[m]).count('1') < min(p, t - 2 * n + 1 + n % 2): res[0] += w
                if bin(Cm[m]).count('1') < min(p, t - 2 * n + (-1) ** n * (1 + n % 2)): res[1] += w
            return
        for m in subs:
            S2 = {b: step(prev, S, b) for b in elts[m]}
            SC2 = {f: {b: step(prev, Sf, b) for b in elts[m]} for f, Sf in SC.items()}
            rec(i + 1, m, S2, SC2, s + len(elts[m]), w)
    for m, w in reps:
        S = {a: 1 << a for a in elts[m]}; SC = {f: {a: (1 << a if a == f else 0) for a in elts[m]} for f in elts[m]}
        rec(1, m, S, SC, len(elts[m]), w)
    return tuple(res)
def count_slow(p, n):
    subs = [frozenset(c) for r in range(2, p + 1) for c in itertools.combinations(range(p), r)]
    bl = bc = 0
    for As in itertools.product(subs, repeat=n):
        if len(sums(As, p)) < bound(As, False, p): bl += 1
        if len(sums(As, p, cyc=True)) < bound(As, True, p): bc += 1
    return bl, bc
assert count_slow(5, 3) == count(5, 3) == (20, 0)
for p, n, exp in [(5, 3, (20, 0)), (5, 4, (0, 0)), (5, 5, (470, 0)), (7, 3, (203, 0)), (7, 4, (0, 336))]:
    got = count(p, n); assert got == exp, (p, n, got)
    say('[5] Z/%d n=%d: violations of (L) %d, of (C) %d' % (p, n, got[0], got[1]))
