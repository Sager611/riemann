import sys, math
from decimal import Decimal as D, getcontext
from fractions import Fraction as Fr
from math import comb, isqrt, gcd
getcontext().prec = 80
sys.path.insert(0, '.')
from lib import *
part = sys.argv[1]
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
def ln(x):
    if isinstance(x, Fr):
        return D(x.numerator).ln() - D(x.denominator).ln()
    return D(x).ln()
EPS = D(10) ** -60
def check_seq(name, a, root_dir=None, root_start=None, ratio_dir=None, ratio_start=None, plain_inc_start=None):
    """a: dict n -> positive value (int or Fraction) for n = 1..N.
       root_dir 'dec'/'inc': (a_n^{1/n}) strictly monotone for n >= root_start
       ratio_dir 'inc'/'dec': (a_{n+1}^{1/(n+1)}/a_n^{1/n}) strictly monotone for n >= ratio_start
       plain_inc_start: (a_{n+1}/a_n) strictly increasing for n >= start (n index as given)"""
    ns = sorted(a)
    l = {n: ln(a[n]) / n for n in ns if n > 0}
    L = {n: ln(a[n]) for n in ns}
    out = []
    close = []
    N = ns[-1]
    if root_dir:
        for n in range(root_start, N):
            d = l[n + 1] - l[n]
            if abs(d) < EPS: close.append(('root', n))
            if (root_dir == 'dec' and not d < 0) or (root_dir == 'inc' and not d > 0): out.append(('root', n))
    if ratio_dir:
        for n in range(ratio_start, N - 1):
            d2 = l[n + 2] - 2 * l[n + 1] + l[n]
            if abs(d2) < EPS: close.append(('ratio', n))
            if (ratio_dir == 'inc' and not d2 > 0) or (ratio_dir == 'dec' and not d2 < 0): out.append(('ratio', n))
    if plain_inc_start is not None:
        for n in range(plain_inc_start, N - 1):
            d2 = L[n + 2] - 2 * L[n + 1] + L[n]
            if not d2 > 0: out.append(('plain', n))
    # also report the last failure below the thresholds (sanity: threshold tightness)
    below = []
    if root_dir:
        for n in range(1 if 1 in l else ns[0], root_start):
            if n + 1 in l and n in l:
                d = l[n + 1] - l[n]
                if (root_dir == 'dec' and not d < 0) or (root_dir == 'inc' and not d > 0): below.append(('root', n))
    if ratio_dir:
        for n in range(1 if 1 in l else ns[0], ratio_start):
            if n in l and n + 2 in l:
                d2 = l[n + 2] - 2 * l[n + 1] + l[n]
                if (ratio_dir == 'inc' and not d2 > 0) or (ratio_dir == 'dec' and not d2 < 0): below.append(('ratio', n))
    print(name, 'n<=%d' % N, 'FAILS:' if out else 'ok', out[:10], '| close:', close[:5], '| below-threshold failures:', below[-4:], flush=True)

if part == 'comb':
    N = NMAX
    # 3.1 Fibonacci
    F = [0, 1]
    for i in range(N + 3): F.append(F[-1] + F[-2])
    check_seq('3.1 Fib', {n: F[n] for n in range(2, N + 1)}, 'inc', 2, 'dec', 4)
    # 3.2 Bell
    bell = [1]; row = [1]
    for i in range(N + 2):
        nr = [row[-1]]
        for x in row: nr.append(nr[-1] + x)
        row = nr; bell.append(row[0])
    check_seq('3.2 Bell', {n: bell[n] for n in range(1, N + 1)}, 'inc', 1, 'dec', 1)
    # 3.3 derangements
    d = [1, 0]
    for n in range(2, N + 3): d.append((n - 1) * (d[n-1] + d[n-2]))
    check_seq('3.3 Derang', {n: d[n] for n in range(2, N + 1)}, 'inc', 2, 'dec', 3)
    # 3.5 tangent numbers T(n): tan x = sum T(n) x^{2n-1}/(2n-1)!
    # via Bernoulli: T(n) = (-1)^{n-1} 2^{2n}(2^{2n}-1) B_{2n}/(2n)
    B = [Fr(1)]
    for m in range(1, 2 * N + 3):
        B.append(-sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1))
    T = {n: int((-1)**(n-1) * 2**(2*n) * (2**(2*n) - 1) * B[2*n] / (2*n)) for n in range(1, N + 1)}
    check_seq('3.5 Tangent', T, 'inc', 1, 'dec', 2, plain_inc_start=1)
    check_seq('2.15 Bernoulli', {n: (-1)**(n-1) * B[2*n] for n in range(1, N + 1)}, 'inc', 1, 'dec', 2)
    E = {0: 1}
    for n in range(1, N + 2):
        E[2*n] = -sum(comb(2*n, 2*k) * E[2*k] for k in range(n))
    check_seq('2.15 Euler', {n: (-1)**n * E[2*n] for n in range(1, N + 1)}, 'inc', 1, 'dec', 1)
    # 3.6 central trinomial
    Tc = {n: sum(comb(n, 2*k) * comb(2*k, k) for k in range(n//2 + 1)) for n in range(1, N + 1)}
    check_seq('3.6 CentralTrinomial', Tc, 'inc', 1, 'dec', 1)
    M = [1, 1]
    for n in range(2, N + 2): M.append(((2*n + 1) * M[n-1] + (3*n - 3) * M[n-2]) // (n + 2))
    check_seq('3.7 Motzkin', {n: M[n] for n in range(1, N + 1)}, 'inc', 1, 'dec', 1)
    # 3.9 g_n
    g = {n: sum(comb(n, k)**2 * comb(2*k, k) for k in range(n + 1)) for n in range(1, N + 1)}
    check_seq('3.9 g_n', g, 'inc', 1, 'dec', 1)
    # 3.11 large Schroder numbers 1,2,6,22,90
    S = [1, 2]
    for n in range(2, N + 2): S.append(((6*n - 3) * S[n-1] - (n - 2) * S[n-2]) // (n + 1))
    check_seq('3.11 Schroder', {n: S[n] for n in range(1, N + 1)}, 'inc', 1, 'dec', 1)
    Dm = {n: sum(comb(n, k)**2 * comb(2*k, k) * comb(2*(n-k), n-k) for k in range(n + 1)) for n in range(0, N + 1)}
    check_seq('3.12 Domb', {n: Dm[n] for n in range(1, N + 1)}, 'inc', 1, 'dec', 1)
    check_seq('3.12 Domb ratio D(n+1)/D(n) (n>=0)', Dm, plain_inc_start=0)
    Pn = {n: 2**n * sum(comb(n, 2*k) * comb(2*k, k)**2 * 4**(n - 2*k) for k in range(n//2 + 1)) for n in range(0, N + 1)}
    check_seq('3.13 P_n', {n: Pn[n] for n in range(1, N + 1)}, 'inc', 1, 'dec', 1)
    check_seq('3.13 P ratio (n>=0)', Pn, plain_inc_start=0)
    # Springer numbers: e.g.f. 1/(cos x - sin x); S_n: 1, 1, 3, 11, 57, 361
    # compute via boustrophedon of sequence 1,1,0,0,... ? use recurrence from e.g.f.: f*(cos-sin)=1
    from fractions import Fraction
    cs = []  # coefficients of cos x - sin x : c_k / k!
    for k in range(N + 3):
        if k % 2 == 0: cs.append((-1)**(k//2))
        else: cs.append(-((-1)**((k-1)//2)))
    Sp = [1]
    for n in range(1, N + 2):
        # sum_{k=0}^n C(n,k) S_{n-k} c_k = 0 for n>=1
        s = sum(comb(n, k) * Sp[n-k] * cs[k] for k in range(1, n + 1))
        Sp.append(-s // cs[0] if cs[0] == 1 else s)
    assert Sp[:6] == [1, 1, 3, 11, 57, 361], Sp[:6]
    check_seq('3.4 Springer', {n: Sp[n] for n in range(1, N + 1)}, None, None, 'dec', 1)
    check_seq('3.4 Springer ratio S_{n+1}/S_n (n>=0)', {n: Sp[n] for n in range(0, N + 1)}, plain_inc_start=0)
    # 3.8 Franel f^{(r)} with N(r); 3.10 Apery-like A^{(r)} with M(r)
    Nr = {2:1,3:1,4:1,5:1,6:1,7:3,8:3,9:3,10:5,11:5,12:7,13:7,14:9,15:9,16:9,17:11,18:11}
    for r, st in Nr.items():
        f = {n: sum(comb(n, k)**r for k in range(n + 1)) for n in range(1, min(N, 400) + 1)}
        check_seq('3.8 f^(%d) N=%d' % (r, st), f, None, None, 'dec', st)
    Mr = {r: 1 for r in range(1, 17)}; Mr.update({17: 9, 18: 9, 19: 9, 20: 12})
    for r, st in Mr.items():
        A = {n: sum(comb(n, k)**r * comb(n + k, k)**r for k in range(n + 1)) for n in range(1, min(N, 300) + 1)}
        check_seq('3.10 A^(%d) M=%d' % (r, st), A, None, None, 'dec', st)

def mr_prime(n):
    if n < 2: return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in small:
        if n % p == 0: return n == p
    d = n - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in small:
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True
def E_N0_N(q, label, claimE, claimN0, claimN):
    # q: list q_1..q_M (1-indexed via q[n-1])
    M = len(q)
    l = [None] + [ln(q[n-1]) / n for n in range(1, M + 1)]
    E = [n for n in range(1, M) if not l[n] > l[n + 1]]
    Q = [0]
    for x in q: Q.append(Q[-1] + x)
    bad0 = [n for n in range(1, M) if not 2 * Q[n] < (n - 1) * q[n]]   # q[n] = q_{n+1}
    N0 = (max(bad0) + 1) if bad0 else 1
    lq = [None] + [ln(Q[n]) / n for n in range(1, M + 1)]
    badN = [n for n in range(1, M - 1) if not (lq[n + 2] - 2 * lq[n + 1] + lq[n]) > 0]
    Nv = (max(badN) + 1) if badN else 1
    print(label, 'M=%d' % M, 'E=', E, '(claim %s)' % claimE, 'N0=', N0, '(claim %s)' % claimN0, 'N=', Nv, '(claim %s)' % claimN, flush=True)

if part == 'nt1':
    # 2.8-2.11
    LIMQ = NMAX
    PR = primes_upto(LIMQ * 2 + 10)
    PS = set(PR)
    tw = [p for p in PR if p + 2 in PS and p < LIMQ]
    E_N0_N(tw, '2.8 (x,x+2)', '{}', 4, 9)
    t1 = [p for p in PR if p + 2 in PS and p + 6 in PS and p < LIMQ]
    t2 = [p for p in PR if p + 4 in PS and p + 6 in PS and p < LIMQ]
    E_N0_N(t1, '2.9 (x,x+2,x+6)', '{}', 3, 13)
    E_N0_N(t2, '2.9 (x,x+4,x+6)', '{}', 6, 13)
    sg = [p for p in PR if 2 * p + 1 in PS and p < LIMQ]
    E_N0_N(sg, '2.10 (x,2x+1)', '{3,4}', 3, 13)
    sg2 = [p for p in PR if 2 * p - 1 in PS and p < LIMQ]
    E_N0_N(sg2, '2.10 (x,2x-1)', '{2,3,6}', 3, 9)
    q3 = [p for p in PR if p < LIMQ // 10 and mr_prime(p*p + p + 1)]
    E_N0_N(q3, '2.11 (x,x^2+x+1)', '{3,4,12,14}', 3, 17)
    q4 = [x for x in range(1, LIMQ // 100) if mr_prime(x**4 + 1)]
    E_N0_N(q4, '2.11 (x^4+1)', '{1,2,4}', 4, 10)
if part == 'nt2':
    # 2.6 primes x^2+1 and x^2+x+1 ; 2.5 Proth
    X = NMAX
    a = [x*x + 1 for x in range(1, X) if mr_prime(x*x + 1)]
    l = [None] + [ln(a[n-1]) / n for n in range(1, len(a) + 1)]
    Eset = [n for n in range(1, len(a)) if not l[n] > l[n+1]]
    Q = [0]
    for v in a: Q.append(Q[-1] + v)
    lq = [None] + [ln(Q[n]) / n for n in range(1, len(a) + 1)]
    badN = [n for n in range(1, len(a) - 1) if not (lq[n+2] - 2*lq[n+1] + lq[n]) > 0]
    print('2.6(i) x^2+1: #', len(a), 'failures of q_{n+1}<q_n^{1+1/n}:', Eset, '(claim 1,2,4,351); Q-ratio failures (claim none from 13):', [n for n in badN if n >= 13], 'last below:', [n for n in badN if n < 13][-3:])
    b = sorted(set(x*x + x + 1 for x in range(0, X) if mr_prime(x*x + x + 1)))
    l = [None] + [ln(b[n-1]) / n for n in range(1, len(b) + 1)]
    Eset = [n for n in range(1, len(b)) if not l[n] > l[n+1]]
    Q = [0]
    for v in b: Q.append(Q[-1] + v)
    lq = [None] + [ln(Q[n]) / n for n in range(1, len(b) + 1)]
    badN = [n for n in range(1, len(b) - 1) if not (lq[n+2] - 2*lq[n+1] + lq[n]) > 0]
    print('2.6(ii) x^2+x+1: #', len(b), 'failures:', Eset, '(claim 3,6); Q-ratio failures (claim none from 20):', [n for n in badN if n >= 20], 'last below:', [n for n in badN if n < 20][-3:])
    # Proth primes k*2^m+1, k odd, 0<k<2^m
    prs = set()
    for m in range(1, 40):
        for k in range(1, 2**m, 2):
            v = k * 2**m + 1
            if v > 10**8: break
            if mr_prime(v): prs.add(v)
    pr = sorted(prs)[:4000]
    l = [None] + [ln(pr[n-1]) / n for n in range(1, len(pr) + 1)]
    Eset = [n for n in range(1, len(pr)) if not l[n] > l[n+1]]
    PRs = [0]
    for v in pr: PRs.append(PRs[-1] + v)
    bad3 = [n for n in range(1, len(pr) + 1) if not 3 * PRs[n] < n * pr[n-1]]
    lq = [None] + [ln(PRs[n]) / n for n in range(1, len(pr) + 1)]
    badN = [n for n in range(1, len(pr) - 1) if not (lq[n+2] - 2*lq[n+1] + lq[n]) > 0]
    print('2.5 Proth: #', len(pr), 'first', pr[:8], 'root failures:', Eset, '(claim 2,4,5); PR<nPr/3 failures:', [n for n in bad3 if n > 50], 'last<=50:', [n for n in bad3 if n <= 50][-2:], '; ratio failures from 34:', [n for n in badN if n >= 34], 'last below:', [n for n in badN if n < 34][-2:])
if part == 'nt3':
    # 2.12 squarefree, 2.13 sums of two squares, 2.14 partitions, 2.16 harmonic
    N = NMAX
    sf = [n for n in range(1, 4 * N) if all(n % (p*p) for p in range(2, isqrt(n) + 1))][:N]
    S = [0]
    for v in sf: S.append(S[-1] + v)
    check_seq('2.12 squarefree s_n', {n: sf[n-1] for n in range(1, N + 1)}, 'dec', 7)
    check_seq('2.12 S(n)', {n: S[n] for n in range(1, N + 1)}, None, None, 'inc', 7)
    def s2(n):
        a = 0
        while a * a <= n:
            b = isqrt(n - a*a)
            if a*a + b*b == n: return True
            a += 1
        return False
    ss = [n for n in range(1, 6 * N) if s2(n)][:N]
    A = [0]
    for v in ss: A.append(A[-1] + v)
    check_seq('2.13 sum2sq a_n', {n: ss[n-1] for n in range(1, N + 1)}, 'dec', 6)
    check_seq('2.13 A(n)', {n: A[n] for n in range(1, N + 1)}, None, None, 'inc', 6)
    p = [1] + [0] * N
    for k in range(1, N + 1):
        for s in range(k, N + 1): p[s] += p[s - k]
    q = [1] + [0] * N
    for k in range(1, N + 1):
        for s in range(N, k - 1, -1): q[s] += q[s - k]
    check_seq('2.14 p(n)', {n: p[n] for n in range(1, N + 1)}, 'dec', 6, 'inc', 26)
    check_seq('2.14 p*(n)', {n: q[n] for n in range(1, N + 1)}, 'dec', 9, 'inc', 45)
    for m in range(1, 13):
        H = Fr(0); Hs = {}
        for n in range(1, min(N, 300) + 1):
            H += Fr(1, n**m); Hs[n] = H
        check_seq('2.16 H^(%d)' % m, Hs, 'dec', 2, 'inc', 3)
if part == 'nt4':
    # 2.17 finite fields
    def mob(n):
        r = 1; m = n; d = 2
        while d * d <= m:
            if m % d == 0:
                m //= d
                if m % d == 0: return 0
                r = -r
            d += 1
        if m > 1: r = -r
        return r
    def Nirr(q, d):
        return sum(mob(e) * q**(d // e) for e in range(1, d + 1) if d % e == 0) // d
    pps = [q for q in range(2, 260) if len([p for p in range(2, q + 1) if q % p == 0 and is_prime(p)]) == 1]
    bad = []
    NN = NMAX
    for q in pps:
        Mv = [0]
        for n in range(1, NN + 3): Mv.append(Mv[-1] + Nirr(q, n))
        for n in range(1, NN):
            # (i) M(n+1)/M(n) < M(n+2)/M(n+1)  <=> M(n+1)^2 < M(n) M(n+2)
            if not Mv[n+1]**2 < Mv[n] * Mv[n+2]:
                if not (q < 5 and n in (2, 4, 6, 8, 10, 12)): bad.append(('i', q, n))
            # (ii) n>2: M(n)^{1/n} < M(n+1)^{1/(n+1)} <=> M(n)^{n+1} < M(n+1)^n
            if n > 2 and not Mv[n]**(n+1) < Mv[n+1]**n:
                if not (q < 7 and n in (3, 5)): bad.append(('ii', q, n))
            # (iii) n>3: r_n > r_{n+1}
            if n > 3:
                l0, l1, l2 = ln(Mv[n]) / n, ln(Mv[n+1]) / (n+1), ln(Mv[n+2]) / (n+2)
                if not (l2 - 2*l1 + l0) < 0:
                    if not ((q < 8 and n in (5, 7, 9, 11, 13)) or (9 < q < 14 and n == 4)): bad.append(('iii', q, n))
    print('2.17 finite-field failures outside exception lists (q<260 prime powers, n<%d):' % NN, bad)

if part == 'nt5':
    # 2.3: primes in AP, q even <= 45, a odd, gcd 1, 1<=a<=q: (p_n^{1/n}) strictly decreasing for n>=2
    LIMP = NMAX
    PR = primes_upto(LIMP)
    bad = []
    for q in range(2, 46, 2):
        for a in range(1, q + 1, 2):
            if gcd(a, q) != 1: continue
            seq = [p for p in PR if p % q == a % q][:3000]
            l = [None] + [ln(seq[n-1]) / n for n in range(1, len(seq) + 1)]
            for n in range(2, len(seq)):
                if not l[n + 1] < l[n]: bad.append((a, q, n)); break
    print('2.3 failures (a,q,n) for q<=45, n<=3000:', bad)
    # 2.2 Firoozbakht refinement n from 5 to 20000
    bad2 = []
    for n in range(5, 20000):
        lhs = ln(PR[n]) / (n + 1) - ln(PR[n - 1]) / n
        rhs = (1 - D(math.log(math.log(n))) / (2 * n * n)).ln()
        if not lhs < rhs: bad2.append(n)
    print('2.2 failures n in [5,20000):', bad2)
if part == 'nt6':
    # 2.4 primes by least positive / largest negative primitive root
    LIMP = NMAX
    PR = primes_upto(LIMP)
    def fac(n):
        fs = []; d = 2
        while d * d <= n:
            if n % d == 0:
                fs.append(d)
                while n % d == 0: n //= d
            d += 1
        if n > 1: fs.append(n)
        return fs
    pos = {}; neg = {}
    for p in PR:
        if p == 2: continue
        fs = fac(p - 1)
        def isg(g):
            g %= p
            return g != 0 and all(pow(g, (p - 1) // f, p) != 1 for f in fs)
        g = 1
        while not isg(g): g += 1
        pos.setdefault(g, []).append(p)
        g = -1
        while not isg(g): g -= 1
        neg.setdefault(g, []).append(p)
    def perfect_power(a):
        b = abs(a)
        if b == 1: return True
        for m in range(2, 64):
            r = round(b ** (1.0 / m))
            for c in (r - 1, r, r + 1):
                if c > 1 and c ** m == b and (a > 0 or m % 2 == 1): return True
        return False
    res = []
    for a, lst in sorted(list(pos.items()) + list(neg.items())):
        if perfect_power(a) or len(lst) < 4: continue
        l = [None] + [ln(lst[n-1]) / n for n in range(1, len(lst) + 1)]
        fails = [n for n in range(1, len(lst)) if not l[n] > l[n + 1]]
        P = [0]
        for v in lst: P.append(P[-1] + v)
        lq = [None] + [ln(P[n]) / n for n in range(1, len(lst) + 1)]
        f3 = [n for n in range(3, len(lst) - 1) if not (lq[n+2] - 2*lq[n+1] + lq[n]) > 0]
        if fails or f3: res.append((a, len(lst), fails[:5], f3[:5]))
    print('2.4: (a, #primes<%d, root failures, P-ratio failures n>=3); claim only a=-2,n=13:' % LIMP, res)
if part == 'nt7':
    # 3.1 general u_n: A>1, B!=0, A^2>4B, (A>2 or B>=-9): ratio decreasing for n>=4
    bad = []
    for A in range(2, 31):
        for B in range(-60, (A*A + 3) // 4):
            if B == 0 or not A*A > 4*B or not (A > 2 or B >= -9): continue
            u = [0, 1]
            for n in range(1, NMAX + 3): u.append(A*u[n] - B*u[n-1])
            if any(x <= 0 for x in u[1:]): bad.append((A, B, 'nonpos')); continue
            l = [None] + [ln(u[n]) / n for n in range(1, NMAX + 3)]
            for n in range(4, NMAX):
                if not (l[n+2] - 2*l[n+1] + l[n]) < 0: bad.append((A, B, n)); break
    print('3.1 general Lucas failures (A<=30, -60<=B, n<=%d):' % NMAX, bad[:20], len(bad))
