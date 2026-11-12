# lemmaP_mine.py (rc5audit, own code): Lemma P of RCApps5 NOTE §2.
# LHS  Q'_D(b) = sum_{n,n'} b(n) conj b(n') sum_{c|(n,n')} (D/c) r_D(nn'/c^2),  r_D(M) = #{(T,y): T^2-Dy^2=4M, T=yD mod 2}/2
# RHS  sum_{e|f sqfree} (w_e/2) sum_C |S_C^{(D/e^2)}(b_e)|^2,  S_C = sum_{m,(m,f_e)=1} b_e(m) r_{F_C}(m)/w_e  (reduced forms F_C)
# b: random complex on odd squarefree n <= NMAX built from small odd primes (so that they divide conductors).
import math, random, itertools, sys
random.seed(7)
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
PR = [3, 5, 7, 11, 13]
def sqf_odd_from(pr, nmax):
    out = [1]
    for p in pr:
        out += [x * p for x in out if x * p <= nmax]
    return sorted(set(out))
NS = [n for n in sqf_odd_from(PR, NMAX) if n > 1]
def jac(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0
def kron(D, c):  # (D/c) for odd squarefree c: product of Legendre symbols
    r = 1
    for p in PR:
        if c % p == 0: r *= jac(D, p)
    return r
def isdisc(D): return D % 4 in (0, 1)
def conductor(D):
    f = 1; Dc = D
    for p in range(2, int(math.isqrt(-D)) + 2):
        while Dc % (p * p) == 0 and isdisc(Dc // (p * p)):
            Dc //= p * p; f *= p
    return f, Dc
def rcount(D, M):  # #{(T,y) : T^2 - D y^2 = 4M, T = yD mod 2}
    cnt = 0; ymax = math.isqrt(4 * M // (-D)) if D < 0 else 0
    for y in range(-ymax, ymax + 1):
        t2 = 4 * M + D * y * y
        if t2 < 0: continue
        T = math.isqrt(t2)
        if T * T != t2: continue
        for TT in ((T, -T) if T else (0,)):
            if (TT - y * D) % 2 == 0: cnt += 1
    return cnt
def forms(D):
    out = []; a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b - D) % 2: continue
            num = b * b - D
            if num % (4 * a): continue
            c = num // (4 * a)
            if c < a or (b < 0 and a == c): continue
            if math.gcd(math.gcd(a, b), c) != 1: continue
            out.append((a, b, c))
        a += 1
    return out
def reps(F, m):
    a, b, c = F; D = b * b - 4 * a * c; cnt = 0
    ymax = math.isqrt(4 * a * m // (-D)) + 1
    for y in range(-ymax, ymax + 1):
        disc = 4 * a * m + D * y * y
        if disc < 0: continue
        s = math.isqrt(disc)
        for x in range((-b * y - s) // (2 * a) - 1, (-b * y + s) // (2 * a) + 2):
            if a * x * x + b * x * y + c * y * y == m: cnt += 1
    return cnt
def lhs(D, b):
    tot = 0
    for n in NS:
        for n2 in NS:
            g = math.gcd(n, n2)
            for c in [d for d in NS + [1] if g % d == 0]:
                k = kron(D, c) if c > 1 else 1
                if k == 0: continue
                M = n * n2 // (c * c); r = rcount(D, M)
                if r: tot += b[n] * b[n2].conjugate() * k * r / 2
    return tot
def rhs(D, f, b):
    tot = 0; blocks = []
    fp = [p for p in PR if f % p == 0]
    for k in range(len(fp) + 1):
        for E in itertools.combinations(fp, k):
            e = math.prod(E); De = D // (e * e); fe = f // e
            we = 6 if De == -3 else (4 if De == -4 else 2)
            be = {}
            for n in NS:
                if n % e == 0 and math.gcd(n // e, f) == 1: be[n // e] = b[n]
            S = []
            for F in forms(De):
                s = 0
                for m, v in be.items():
                    if math.gcd(m, fe) == 1:
                        r = reps(F, m)
                        if r: s += v * r / we
                S.append(s)
            blk = we / 2 * sum(abs(s) ** 2 for s in S); blocks.append((e, blk)); tot += blk
    return tot, blocks
tests = 0; worst = 0; negblk = 0; nonmax = 0; nontriv = 0; small = 0
Ds = [D for D in range(-7, -2600, -1) if isdisc(D)]
cand = [D for D in Ds if any(conductor(D)[0] % p == 0 for p in PR)]
random.shuffle(cand)
for D in cand[:int(sys.argv[2]) if len(sys.argv) > 2 else 60]:
    f, DK = conductor(D)
    b = {n: complex(random.gauss(0, 1), random.gauss(0, 1)) for n in NS}
    L = lhs(D, b); R, blocks = rhs(D, f, b)
    dev = abs(L - R) / max(abs(L), 1e-12); worst = max(worst, dev); tests += 1; nonmax += f > 1
    nontriv += any(e > 1 and bl > 1e-9 for e, bl in blocks); small += any((D // (e * e)) in (-3, -4) for e, _ in blocks)
    negblk += any(bl < -1e-9 * max(1, abs(L)) for _, bl in blocks)
    if abs(L.imag) > 1e-8 * max(1, abs(L)): print('complex LHS?', D, L)
print(f"NS={len(NS)} (odd squarefree <= {NMAX} from {PR}); tests={tests} nonmaximal={nonmax} with nonzero e>1 block={nontriv} with Delta_e in (-3,-4)={small}; max rel dev={worst:.2e}; negative blocks={negblk}")
