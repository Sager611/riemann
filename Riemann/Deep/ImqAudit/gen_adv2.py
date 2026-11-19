# gen_adv2.py LO HI COUNT PMAX SEED: COUNT fundamental d in [LO, HI], d = 3 mod 8, (-d|p) = -1 for every odd p <= PMAX
# (CRT over p <= 29, rejection for the rest): the classes with the smallest M(d), for the top of the sieve range.
# Non-squarefree d are kept here; hcount reports them as non-fundamental.
import random, sys, math
LO, HI, CNT, PMAX, SEED = map(int, sys.argv[1:6]); random.seed(SEED)
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
def primes(n): return [p for p in range(3, n + 1) if all(p % q for q in range(2, int(p ** .5) + 1))]
P = primes(PMAX); P1 = [p for p in P if p <= 29]; P2 = [p for p in P if p > 29]
NR = {p: [y for y in range(1, p) if jac(-y, p) == -1] for p in P1}
def sqf(u):
    p = 3
    while p * p <= u:
        if u % (p * p) == 0: return False
        p += 2
    return True
out = set(); tries = 0
while len(out) < CNT:
    m, b = 8, 3
    for p in P1:
        x = random.choice(NR[p])
        while (b - x) % p: b += m
        m *= p
    for k in range((LO - b) // m + 1, (HI - b) // m + 1):
        d = b + m * k; tries += 1
        if LO <= d <= HI and all(jac(-d, p) == -1 for p in P2): out.add(d)   # squarefreeness: hcount checks
for d in sorted(out): print(d)
print('tries', tries, file=sys.stderr)
