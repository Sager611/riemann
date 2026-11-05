"""Z.-W. Sun, OEIS A229827 (2013), conjectures on circular permutations with sums P(a_i) + a_{i+1}.
(i) general part: for any finite field F_q and P of degree < q-1, not of the form c - x, some circular permutation of F_q has
    {P(a_1)+a_2, ..., P(a_q)+a_1} = F_q.  FAILS for P = Frobenius: q = 8, P = x^2 or x^4; q = 9, P = x^3 (exhaustive below).
(ii) any field, |A| = n > 2: fails for A = {0,1,2}, P(x) = 2x -- already noted by Sun (arXiv 1309.1679, Remark after Thm 1.3).
Run: python3 verify_a229827.py (stdlib only, a few seconds)."""
import itertools
def best_cycle(q, add, P):
    """max number of distinct sums over all circular permutations of range(q) with a_1 = 0, and number of complete ones"""
    best = 0; ok = 0
    for perm in itertools.permutations(range(1, q)):
        c = (0,) + perm
        s = len({add(P[c[i]], c[(i + 1) % q]) for i in range(q)})
        best = max(best, s); ok += (s == q)
    return best, ok
# F_8 = F_2[x]/(x^3+x+1): 3-bit ints, addition xor
def mul8(a, b):
    r = 0
    for i in range(3):
        if b >> i & 1: r ^= a << i
    for d in (4, 3):
        if r >> d & 1: r ^= 0b1011 << (d - 3)
    return r
def pw(mul, x, e):
    r = 1
    for _ in range(e): r = mul(r, x)
    return r
for e in (2, 4):
    P = [pw(mul8, x, e) for x in range(8)]
    best, ok = best_cycle(8, lambda a, b: a ^ b, P)
    assert ok == 0
    print("F_8, P(x) = x^%d: no complete circular permutation (max %d distinct sums of 8)" % (e, best))
# F_9 = F_3[t]/(t^2+1): element a + b t encoded 3a + b... use (a, b) -> 3*b + a
def dec9(n): return n % 3, n // 3
def enc9(a, b): return a % 3 + 3 * (b % 3)
def add9(x, y):
    (a, b), (c, d) = dec9(x), dec9(y); return enc9(a + c, b + d)
def mul9(x, y):
    (a, b), (c, d) = dec9(x), dec9(y); return enc9(a * c - b * d, a * d + b * c)
P = [pw(mul9, x, 3) for x in range(9)]
assert enc9(0, 0) == 0
best, ok = best_cycle(9, add9, P)
assert ok == 0
print("F_9, P(x) = x^3: no complete circular permutation (max %d distinct sums of 9)" % best)
# P is not of the form c - x on F_8 / F_9 (checked as functions):
assert not any(all(pw(mul8, x, 2) == c ^ x for x in range(8)) for c in range(8))
assert not any(all(P[x] == add9(c, mul9(enc9(2, 0), x)) for x in range(9)) for c in range(9))
# (ii), known: A = {0,1,2}, P(x) = 2x over Q
assert all(len({2 * c[i] + c[(i + 1) % 3] for i in range(3)}) < 3 for c in [(0, 1, 2), (0, 2, 1)])
print("(ii) {0,1,2}, P = 2x fails (known, Sun's Remark 1.3)")

