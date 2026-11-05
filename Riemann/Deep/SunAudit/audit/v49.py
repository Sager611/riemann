import itertools
def S(A, n, p=None):
    out = set()
    for seq in itertools.permutations(A, n):
        s = sum((k + 1) * a for k, a in enumerate(seq))
        out.add(s % p if p else s)
    return out
def bound(sz, n, p=None):
    b = (sz - n) * n * (n + 1) // 2 + n * (n * n - 1) // 6 + 1
    return min(p, b) if p else b
cases = [((0,1,2), 2, None), ((0,1,2), 2, 5), ((0,1,2,4), 2, 7), ((0,1,2,5), 3, 11), ((0,1,2,5), 4, 11), ((0,1,3,4), 4, 13), ((0,1,3,4), 3, 13), ((0,1,2,3,6), 3, 13)]
for A, n, p in cases:
    s = S(A, n, p)
    print(A, 'n=%d' % n, 'F_%s' % p if p else 'Q', '|S|=%d' % len(s), 'bound=%d' % bound(len(A), n, p), sorted(s))
# all failures in Q for n=2, A subset {0..12}, |A|<=6, normalized min 0 gcd 1
from math import gcd
fails = []
for sz in range(2, 7):
    for rest in itertools.combinations(range(1, 13), sz - 1):
        g = 0
        for x in rest: g = gcd(g, x)
        if g != 1: continue
        A = (0,) + rest
        if len(S(A, 2)) < bound(sz, 2): fails.append(A)
print('Q n=2 failures (A in {0..12}, |A|<=6):', fails)
