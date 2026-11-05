import sys, math
from math import gcd, isqrt, comb, factorial
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from lib import *
part = sys.argv[1]
LIM = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
PRIMES = primes_upto(max(LIM, 10))
_pf = {}
def pfac(n):
    if n in _pf: return _pf[n]
    m = n; fs = []; d = 2
    while d * d <= m:
        if m % d == 0:
            fs.append(d)
            while m % d == 0: m //= d
        d += 1 if d == 2 else 2
    if m > 1: fs.append(m)
    _pf[n] = fs
    return fs
def is_pr(g, p):
    """g (int or Fraction) primitive root mod prime p; p=2: g odd"""
    if isinstance(g, Fr):
        if g.denominator % p == 0: return False
        g = g.numerator * pow(g.denominator, -1, p)
    g %= p
    if g == 0: return False
    if p == 2: return True
    for q in pfac(p - 1):
        if pow(g, (p - 1) // q, p) == 1: return False
    return True
def is_sq_mod(f, n, sqcache={}):
    if n not in sqcache:
        sqcache.clear(); sqcache[n] = set(x * x % n for x in range(n))
    return f % n in sqcache[n]
fails = []
def F(name, *a):
    fails.append((name,) + a)

if part == 'A':
    # 2.1(i), 2.2(i): composite moduli n
    fib = [0, 1]
    while fib[-1] < 10**7: fib.append(fib[-1] + fib[-2])
    luc = [2, 1]
    while luc[-1] < 10**7: luc.append(luc[-1] + luc[-2])
    for n in range(3, 2001):
        sq = set(x * x % n for x in range(n))
        if n > 4 and not any(f % n not in sq for f in fib if 2 * f < n): F('2.1(i)', n)
        if n > 2 and not any(L % n not in sq for L in luc if L < n): F('2.2(i)', n)
    # 2.1(iii), 2.2(iii)
    for p in PRIMES:
        a, b = 0, 1  # F_0, F_1
        okF = False
        k = 0
        while k <= math.sqrt(p + 2) + 2:
            if k >= 1 and is_pr(a + 1, p): okF = True; break
            a, b = b, a + b; k += 1
        if not okF: F('2.1(iii)', p)
        L0, L1 = 2, 1; k = 0; okL = False
        while k < math.sqrt(p) + 2:
            if k >= 1 and is_pr(L0 + 1, p): okL = True; break
            L0, L1 = L1, L0 + L1; k += 1
        if not okL: F('2.2(iii)', p)
        # 2.3
        if p % 3 == 1:
            k = 1; ok = False
            while 2 * (2**k - 1) < p:
                v = (2**k - 1) % p
                if v and pow(v, (p - 1) // 3, p) != 1: ok = True; break
                k += 1
            if not ok: F('2.3', p)
        # 2.4
        if p > 7 and not any(leg(pow(2, q, p) - 1, p) == 1 for q in PRIMES if q < p): F('2.4a', p)
        if p > 5 and not any(leg(pow(2, q, p) + 1, p) == -1 for q in PRIMES if q < p): F('2.4b', p)
    print('A fails:', fails)

if part == 'B':
    for p in PRIMES:
        # 3.1(i): g = k^2+1 < p
        if not any(is_pr(k*k + 1, p) for k in range(0, isqrt(p) + 1) if k*k + 1 < p): F('3.1(i)', p)
        if p > 3 and not any(is_pr(k*(k+1)//2, p) for k in range(1, 2*isqrt(p) + 3) if 0 < k*(k+1)//2 < p): F('3.1(ii)tri', p)
        if p > 7 and not any(is_pr(k*(k+1), p) for k in range(1, isqrt(p) + 2) if k*(k+1) < p): F('3.1(ii)pron', p)
        # 3.2(i): S_n < p
        if p > 2:
            s = 0; ok = False
            for q in PRIMES:
                s += q
                if s >= p: break
                if is_pr(s, p): ok = True; break
            if not ok: F('3.2(i)', p)
    # 3.2(ii),(iii): indexed by n
    for n in range(2, min(len(PRIMES), 3000)):
        pn = PRIMES[n - 1]
        s = 0; ok = False
        for k in range(1, n + 1):
            s = PRIMES[k - 1] - s
            if is_pr(s, pn): ok = True; break
        if not ok: F('3.2(ii)', n)
        if n > 4:
            pr = 1; ok = False
            for k in range(1, n // 2 + 1):
                pr *= PRIMES[k - 1]
                if is_pr(pr, pn): ok = True; break
            if not ok: F('3.2(iii)', n)
    print('B fails:', fails)

if part == 'C':
    for p in PRIMES:
        if p > 3:
            if not any(is_pr(pow(2, q, p) - 1, p) for q in PRIMES if 2 * q < p): F('3.3(i)', p)
            ok = False
            fact = 1  # (g-1)!
            for g in range(1, p):
                if g > 1: fact = fact * (g - 1) % p
                if is_pr(g, p) and is_pr(pow(2, g, p) - 1, p) and is_pr(fact, p): ok = True; break
            if not ok: F('3.3(ii)', p)
        if p > 2 and not any(is_pr(q, p) and is_pr(pow(2, q, p) - q, p) for q in PRIMES if q < p): F('3.4', p)
        if p > 7:
            ok = False; fact = 1; last = 1
            for q in PRIMES:
                if q >= p: break
                for t in range(last + 1, q + 1): fact = fact * t % p
                last = q
                if is_pr(q, p) and is_pr(fact, p): ok = True; break
            if not ok: F('3.5(i)qq!', p)
        if p > 2:
            ok = False
            for i, q in enumerate(PRIMES):
                if q >= p: break
                if q - 1 < len(PRIMES) and is_pr(q, p) and is_pr(q, PRIMES[q - 1]): ok = True; break
            if not ok: F('3.5(ii)', p)
    print('C fails:', fails)

if part == 'D':
    fib = [0, 1]
    while fib[-1] < 10**8: fib.append(fib[-1] + fib[-2])
    fibp = sorted(set(f for f in fib if f > 0))
    luc = [2, 1]
    while luc[-1] < 10**8: luc.append(luc[-1] + luc[-2])
    lucs = sorted(set(luc))
    cb = [comb(2*k, k) for k in range(40)]
    cat = [comb(2*k, k)//(k+1) for k in range(40)]
    for p in PRIMES:
        prods = set(a*b for a in fibp for b in fibp if a*b < p)
        if not any(is_pr(g, p) for g in prods): F('3.6(i)a', p)
        if not any(is_pr(l*f, p) for l in (1, 2, 3) for f in fibp if l*f < p): F('3.6(i)b', p)
        if not any(is_pr(a*b, p) for a in lucs for b in lucs if a*b < p): F('3.6(ii)', p)
        if p > 7:
            vals = set(5**k + 10**m for k in range(0, 12) for m in range(0, 9) if 5**k + 10**m < p)
            if not any(is_pr(g, p) for g in vals): F('3.7', p)
        if p > 2:
            if not any(is_pr(a + b, p) for a in cb for b in cb if a + b < p): F('3.8(i)', p)
            if not any(is_pr(a + b, p) for a in cat for b in cat if a + b < p): F('3.8(ii)', p)
    print('D fails:', fails)

if part == 'E':
    # Bernoulli/Euler/tangent mod p via exact rationals for small p
    N = LIM
    # Bernoulli numbers B_0..B_N exact (Akiyama-Tanigawa is slow); use recurrence sum_{k<m} C(m+1,k) B_k = -(m+1) B_m
    B = [Fr(1)]
    for m in range(1, N + 1):
        s = sum(comb(m + 1, k) * B[k] for k in range(m))
        B.append(-s / (m + 1))
    # Euler numbers E_n (secant): E_0=1, sum_{k=0}^{n} C(2n,2k) E_{2k} = 0 for n>=1
    E = {0: 1}
    for n in range(1, N // 2 + 2):
        E[2*n] = -sum(comb(2*n, 2*k) * E[2*k] for k in range(n))
    def t(n):  # tangent number t_n = (-1)^{n-1} 2^{2n}(2^{2n}-1) B_{2n}/(2n)
        v = (-1)**(n-1) * 2**(2*n) * (2**(2*n) - 1) * B[2*n] / (2*n)
        assert v.denominator == 1
        return v.numerator
    for p in PRIMES:
        if p > N: break
        if p > 3 and not any(is_pr(B[q - 1], p) for q in PRIMES if q < p): F('BET(i)', p)
        if p > 13 and not any(is_pr(E[q - 1], p) for q in PRIMES if q < p and q > 2): F('BET(ii)', p)
        if p > 2 and not any(is_pr(t(q), p) for q in PRIMES if q < p and 2*q <= N): F('BET(iii)', p)
        if p > 5:
            H = Fr(0); H2 = Fr(0); ok1 = ok2 = False
            PS = set(PRIMES)
            for m in range(1, p):
                # here H = H_{m-1}, H2 = H^{(2)}_{m-1}
                if m in PS:
                    q = m
                    if 2*q <= p + 1 and is_pr(H, p): ok1 = True
                    if 2*q <= p - 1 and is_pr(H2, p): ok2 = True
                H += Fr(1, m); H2 += Fr(1, m*m)
            if not ok1: F('H(i)', p)
            if not ok2: F('H(ii)', p)
    print('E fails:', fails, 'checked p <=', N)
if part == 'F':
    N = LIM
    # Bell numbers via triangle
    bell = [1]; row = [1]
    for i in range(N + 2):
        nr = [row[-1]]
        for x in row: nr.append(nr[-1] + x)
        row = nr; bell.append(row[0])
    for p in PRIMES:
        if p > N: break
        if p > 3:
            if not any(is_pr(comb(2*q, q)//(q+1), p) for q in PRIMES if 2*q < p): F('mix(i)Catalan', p)
            if not any(is_pr(bell[q], p) for q in PRIMES if 2*q < p): F('mix(ii)Bell', p)
            if not any(is_pr(sum(comb(q-1, k)**3 for k in range(q + 1)), p) for q in PRIMES if 2*q < p): F('mix(iii)Franel', p)
            def T(q): return sum(comb(q, 2*k) * comb(2*k, k) for k in range(q//2 + 1))
            if not any(is_pr(T(q), p) for q in PRIMES if q < p): F('mix(iv)T', p)
    print('F fails:', fails, 'checked p <=', N)

if part == 'G':
    # 4.1(i) pairs
    N = LIM
    P = [p for p in PRIMES if p <= N]
    exc = {frozenset(s) for s in [(2,3),(2,11),(2,13),(2,59),(2,131),(2,181),(3,7),(3,31),(3,79),(3,191),(3,199),(5,271),(7,11),(7,13),(7,71)]}
    found_exc = []
    for i, p in enumerate(P):
        for q in P[i:]:
            g = 1
            while not (is_pr(g, p) and is_pr(g, q)): g += 1
            if g * g > 4 * p * q + 1: F('4.1(i)a', p, q, g)
            if g * g >= p * q:
                found_exc.append((p, q, g))
                if frozenset((p, q)) not in exc: F('4.1(i)b', p, q, g)
    print('G fails:', fails, '| pairs with least common g >= sqrt(pq):', found_exc, 'checked primes <=', N)
if part == 'H':
    # 4.2 (i),(ii); right triangles; three consecutive; gap
    N = LIM
    # pythagorean triples with legs/hyp < N
    trip = []
    for m in range(2, isqrt(N) + 2):
        for n_ in range(1, m):
            if (m - n_) % 2 == 1 and gcd(m, n_) == 1:
                a0, b0, c0 = m*m - n_*n_, 2*m*n_, m*m + n_*n_
                k = 1
                while k * c0 < N:
                    trip.append((k*a0, k*b0, k*c0)); k += 1
    exc21 = []; exc_three = []
    for p in PRIMES:
        if p > N: break
        if p > 13:
            ok = False
            prl = [g for g in range(1, p) if is_pr(g, p)]
            prs = set(prl)
            for a in prl:
                for b in prl:
                    if a != b and a + b < p and (a + b) in prs and is_pr(a * b * (a + b), p): ok = True; break
                if ok: break
            if not ok: F('4.2(i)', p)
            if p != 31:
                if not any(is_pr(a*b*c, p) for (a, b, c) in trip if c < p): F('4.2(ii)', p)
        if p > 7:
            if not any((a*b) % 2 == 0 and is_pr(a*b//2, p) for (a, b, c) in trip if c < p): F('tri(i)', p)
        if p > 31:
            if not any(leg(a+b+c, p) == 1 and leg(a*b//2, p) == 1 for (a, b, c) in trip if c < p): F('tri(ii)', p)
        three = any(is_pr(x, p) and is_pr(x+1, p) and is_pr(x+2, p) for x in range(1, p - 2))
        if not three: exc_three.append(p)
        # gap
        if p > 2:
            run = 0; maxrun = 0
            for x in range(0, 2 * p):
                if is_pr(x, p): run = 0
                else:
                    run += 1; maxrun = max(maxrun, run)
            if maxrun > 2 * isqrt(p + 2) + 2: F('gap', p, maxrun)
    print('H fails:', fails, '| primes without 3 consecutive prim roots:', exc_three, 'checked p <=', N)
if part == 'I':
    # S_p(a,b,c) bound and g-quad, p up to LIM
    for p in PRIMES:
        if p > LIM: break
        if p < 5: continue
        prl = [g for g in range(1, p) if is_pr(g, p)]
        L = [0] + [leg(t, p) for t in range(1, p)]
        worst = 0; badq = False
        # normalize: a in {0,1}; (1,b,c) with b^2-4c != 0 ; (0,b,1) b != 0; also (a,b,c) with a not 0 but c... scaling only
        triples = [(1, b, c) for b in range(p) for c in range(p) if (b*b - 4*c) % p] + [(0, b, 1) for b in range(1, p)]
        for (a, b, c) in triples:
            vals = [L[(a*g*g + b*g + c) % p] for g in prl]
            S = sum(vals)
            worst = max(worst, abs(S))
            if p > 13 and p not in (19, 31):
                if 1 not in vals or -1 not in vals: badq = True; badt = (a, b, c)
        if p > 11 and not worst < math.sqrt(p) / 2 * math.log(p): F('Sp-bound', p, worst, math.sqrt(p)/2*math.log(p))
        if badq: F('g-quad', p, badt)
    print('I fails:', fails, 'checked p <=', LIM)
if part == 'J':
    # g +- g^{-1} primitive roots mod p^2 (a>=2 case), p up to LIM
    def is_pr_mod(g, p, a):
        m = p**a
        if g % p == 0: return False
        order = p**(a-1) * (p - 1)
        qs = set(pfac(p - 1)) | ({p} if a > 1 else set())
        return all(pow(g, order // q, m) != 1 for q in qs)
    exc1 = []; exc2 = []
    for p in PRIMES:
        if p > LIM: break
        if p < 3: continue
        m = p * p
        # primitive roots mod p^2: g prim root mod p and g^{p-1} != 1 mod p^2
        ok1 = ok2 = False
        for g in range(1, m):
            if not is_pr_mod(g, p, 2): continue
            gi = pow(g, -1, m)
            if not ok1 and is_pr_mod((g + gi) % m, p, 2): ok1 = True
            if not ok2 and is_pr_mod((g - gi) % m, p, 2): ok2 = True
            if ok1 and ok2: break
        if not ok1: exc1.append(p)
        if not ok2: exc2.append(p)
    for p in exc1:
        if p > 7 and p != 13: F('g+1/g a=2', p)
    for p in exc2:
        if p > 5 and p not in (13, 61): F('g-1/g a=2', p)
    print('J fails:', fails, '| exceptions g+1/g:', exc1, 'g-1/g:', exc2, 'checked p <=', LIM)

if part == 'K':
    N = LIM
    def phi(n):
        r = n; m = n; d = 2
        while d * d <= m:
            if m % d == 0:
                while m % d == 0: m //= d
                r -= r // d
            d += 1
        if m > 1: r -= r // m
        return r
    PH = [0] + [phi(i) for i in range(1, 2 * N + 2)]
    BIGP = primes_upto(20 * N + 1000)
    PS = set(primes_upto(200 * N + 5000))
    def two_pr(p): return p in PS and p > 2 and is_pr(2, p)
    bad_i = []; bad_ii = []
    for n in range(2, N + 1):
        ok = any((BIGP[k-1] + PH[n-k]) in PS and 2*(BIGP[k-1] + PH[n-k]) + 1 in PS and two_pr(BIGP[k-1] + PH[n-k]) for k in range(1, n))
        if not ok: bad_i.append(n)
        ok2 = any(m % 2 == 0 and two_pr(PH[n-m] + PH[m//2] - 1) for m in range(1, n))
        if not ok2: bad_ii.append(n)
    print('4.4(i) n without rep (2..%d):' % N, bad_i, ' (conj: n>37)')
    print('4.4(ii) n without rep (m even):', bad_ii, ' (conj: n>7)')
if part == 'L':
    # primitive prime divisors via gcd stripping
    N = LIM
    def ppd_fail(seq, start, label, lo_index=1):
        # seq: dict n -> integer term (abs), for n from lo_index..N ; report n>=start without primitive prime divisor
        bad = []
        ns = sorted(seq)
        for idx, n in enumerate(ns):
            if n < start: continue
            R = abs(seq[n])
            if R == 0: bad.append(n); continue
            for k in ns[:idx]:
                a = abs(seq[k])
                if a == 0: continue
                g = gcd(R, a)
                while g > 1:
                    R //= g
                    g = gcd(R, g)
                if R == 1: break
            if R == 1: bad.append(n)
        print(label, 'no primitive prime divisor for n in [%d..%d]:' % (start, ns[-1]), bad, flush=True)
    seq = {n: 2**n - n for n in range(1, N + 1)}
    ppd_fail(seq, 2, '5.1 2^n-n (exceptions allowed 5,16)')
    B = [Fr(1)]
    for m in range(1, 2 * N + 1):
        B.append(-sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1))
    ppd_fail({n: B[2*n].numerator for n in range(1, N + 1)}, 5, '5.2 B_2n numerators (n>4)')
    E = {0: 1}
    for n in range(1, N + 1):
        E[2*n] = -sum(comb(2*n, 2*k) * E[2*k] for k in range(n))
    ppd_fail({n: E[2*n] for n in range(1, N + 1)}, 2, '5.2 E_2n (n>=2)')
    t = {n: ((-1)**(n-1) * 2**(2*n) * (2**(2*n) - 1) * B[2*n] / (2*n)).numerator for n in range(1, N + 1)}
    ppd_fail(t, 4, '5.2 tangent t_n (n>=4)')
    H = Fr(0); H2 = Fr(0); Hs = {}; H2s = {}
    for n in range(1, N + 1):
        H += Fr(1, n); H2 += Fr(1, n*n); Hs[n] = H.numerator; H2s[n] = H2.numerator
    ppd_fail(Hs, 2, '5.3 H_n numerators (n>1, n!=7)')
    ppd_fail(H2s, 2, '5.3 H2_n numerators (n>1)')
if part == 'M':
    N = LIM
    bell = [1]; row = [1]
    for i in range(N + 2):
        nr = [row[-1]]
        for x in row: nr.append(nr[-1] + x)
        row = nr; bell.append(row[0])
    exec(open('c1405.py').read().split("if part == 'L':")[1].split("    seq = {n: 2**n - n")[0].replace('\n    ', '\n'))
    ppd_fail({n: bell[n] for n in range(1, N + 1)}, 2, 'Bell(n) n>=2')
    ppd_fail({n: sum(comb(n, k)**3 for k in range(n + 1)) for n in range(1, N + 1)}, 1, 'Franel f_n n>=1')
    ppd_fail({n: sum(comb(n, k)**4 for k in range(n + 1)) for n in range(1, N + 1)}, 1, 'Franel4 n>=1')
    ppd_fail({n: sum(comb(n, 2*k) * comb(2*k, k) for k in range(n//2 + 1)) for n in range(1, N + 1)}, 2, 'T_n n>1')
    M = [1, 1]
    for n in range(2, N + 1): M.append(((2*n + 1) * M[n-1] + (3*n - 3) * M[n-2]) // (n + 2))
    ppd_fail({n: M[n] for n in range(1, N + 1)}, 4, 'Motzkin n>=4')
    ppd_fail({n: sum(comb(n, k) * comb(n + k, k) for k in range(n + 1)) for n in range(1, N + 1)}, 1, 'Delannoy n>=1')
    ppd_fail({n: sum(comb(n, k)**2 * comb(n + k, k)**2 for k in range(n + 1)) for n in range(1, N + 1)}, 1, 'Apery n>=1')
    d = [1, 0]
    for n in range(2, N + 1): d.append((n - 1) * (d[n-1] + d[n-2]))
    ppd_fail({n: d[n] for n in range(1, N + 1)}, 3, 'derangement n>=3')
    ppd_fail({n: sum(comb(n, k)**2 * comb(2*k, k) * comb(2*(n-k), n-k) for k in range(n + 1)) for n in range(1, N + 1)}, 4, 'Domb n>=4')
    # q(n) distinct partitions
    qq = [1] + [0] * (4 * N)
    for k in range(1, 4 * N + 1):
        for s in range(4 * N, k - 1, -1): qq[s] += qq[s - k]
    ppd_fail({n: qq[n] for n in range(1, 4 * N + 1)}, 1, 'q(n) (conj n>203)')
