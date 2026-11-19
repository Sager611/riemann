# gen_adv.py (ImqAudit): falsification samples in the sieve range [2.8e10, 7.627e11] (seed fixed).
# A: d = 3 mod 8 (2 inert) with (-d|p) = -1 for every odd p <= 47 (the hardest classes: smallest M and h),
#    log-uniform over the range; B: random fundamental d of the six classes next to the range ends and chunk edges;
#    C: d with every p <= 13 split (large M). Output: one d per line, fundamental only.
import random, math, sys
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1119)
LO, HI = 28_000_000_000, 762_700_000_000
P = [3,5,7,11,13,17,19,23,29,31,37,41,43,47]
def jac(a, n):
    a %= n; t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: t = -t
        a %= n
    return t if n == 1 else 0
def sqf(u):
    if u % 4 == 0: return False
    p = 3
    while p * p <= u:
        if u % (p * p) == 0: return False
        p += 2
    return True
def fund(d):
    r = d % 16
    if r in (3, 7, 11, 15): return sqf(d)
    if r in (4, 8): return sqf(d // 4)
    return False
out = []
mod = 8; base = 3
for p in P[:9]:  # CRT up to 29
    nr = [x for x in range(1, p) if jac(-x, p) == -1]
    x = random.choice(nr)
    while (base - x) % p: base += mod
    mod *= p
A = []
while len(A) < 400:
    # fresh residues each time for variety
    m, b = 8, 3
    for p in P[:9]:
        x = random.choice([y for y in range(1, p) if jac(-y, p) == -1])
        while (b - x) % p: b += m
        m *= p
    t = math.exp(random.uniform(math.log(LO), math.log(HI)))
    d = b + m * int((t - b) // m)
    if d < LO: d += m
    if d > HI: continue
    if all(jac(-d, p) == -1 for p in P[9:]) and fund(d): A.append(d)
Bs = []
edges = [LO, HI, 30_999_999_999, 33_999_999_999, 47_244_640_255, 99_999_999_999, 100_000_000_000, 29_796_335_616, 762_600_393_220]
for e in edges:
    k = 0; d = e
    while k < 8 and LO <= d <= HI:
        if fund(d): Bs.append(d); k += 1
        d += 1 if e != HI else -1
    d = e - 1; k = 0
    while k < 4 and LO <= d <= HI:
        if fund(d): Bs.append(d); k += 1
        d -= 1
C = []
while len(C) < 40:
    d = random.randrange(LO, HI)
    if d % 8 == 7 and all(jac(-d, p) == 1 for p in P[:5]) and fund(d): C.append(d)
for d in sorted(set(A + Bs + C)): print(d)
print(len(A), len(set(Bs)), len(C), file=sys.stderr)
