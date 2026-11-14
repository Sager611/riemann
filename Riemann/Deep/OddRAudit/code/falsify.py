#!/usr/bin/env python3
# falsify.py -- oddraudit hunt for violators of (OR) sigma(n)/n < (e^g/2) log log n among odd n > N:
#  (1) neighbours of N: N p^k / q^j and N p r / q (odd primes p, q, r <= 200, k <= 3, j <= 2), every size;
#  (2) neighbours of every odd CA number Nk >= c0 with x_1 <= XCA: Nk p/q, Nk p, Nk/q, Nk p r/q^j (p, q, r <= 200, j <= 2);
#  (3) the odd superabundant records (sigma(m)/m > all smaller odd m) in (N, c0), from S.txt of enum_s.py.
# Screening in float64 via exact per-prime log increments (|error| < 1e-13); every n with margin > -1e-9 rechecked in mpmath 50 digits.
#   uv run --no-project --with mpmath python3 falsify.py S.txt XCA
import sys, math, time
from mpmath import mp, mpf, log, euler
t0 = time.time(); mp.dps = 50
Sfile = sys.argv[1]; XCA = int(float(sys.argv[2])) if len(sys.argv) > 2 else 3000
def oddprimes(n):
    s = bytearray([1])*(n + 1); s[0] = s[1] = 0
    for i in range(2, math.isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(3, n + 1) if s[i]]
SP = oddprimes(200)
LC = float(euler) - math.log(2)
def sg(p, a): return (p**(a + 1) - 1)//(p - 1)
def lr(p, a): return math.log(sg(p, a)) - a*math.log(p)          # log(sigma(p^a)/p^a)
def mexact(f):
    n = 1; s = 1
    for p, a in f.items(): n *= p**a; s *= sg(p, a)
    return log(mpf(s)/n) - (euler - log(2)) - log(log(log(mpf(n)))), n
def scan(base, lim_lo, tag, out, SP=SP):
    # base: dict exponent vector; all perturbations; returns count, worst
    ls = sum(lr(p, a) for p, a in base.items()); lt = sum(a*math.log(p) for p, a in base.items())
    cnt = 0; worst = (-9, None)
    def test(d):                     # d: dict of exponent changes
        nonlocal cnt, worst
        s_ = ls; t_ = lt
        for p, k in d.items():
            a0 = base.get(p, 0); a1 = a0 + k
            if a1 < 0: return
            s_ += lr(p, a1) - lr(p, a0); t_ += k*math.log(p)
        if t_ <= lim_lo: return
        cnt += 1
        m = s_ - LC - math.log(math.log(t_))
        if m > worst[0]: worst = (m, d)
        if m > -1e-9:
            f = dict(base)
            for p, k in d.items(): f[p] = f.get(p, 0) + k
            me, n = mexact({p: a for p, a in f.items() if a})
            out.append((tag, d, float(me), n))
    for p in SP:
        for k in (1, 2, 3): test({p: k})
        for q in SP:
            if q == p: continue
            for k in (1, 2, 3):
                for j in (1, 2): test({p: k, q: -j})
            for r in SP:
                if r < p or r == q: continue
                d = {p: 1, q: -1}; d[r] = d.get(r, 0) + 1
                test(d)
    for q in SP:
        for j in (1, 2): test({q: -j})
    return cnt, worst
SP100 = [p for p in SP if p <= 100]
fac = {3: 4, 5: 3, 7: 2}
for p in oddprimes(61):
    if p >= 11: fac[p] = 1
N = 277093800961598968947975375; c0 = 67*N
assert math.prod(p**a for p, a in fac.items()) == N
lN = math.log(N)
flag = []
cnt, worst = scan(fac, lN + 1e-9, 'N', flag)
print(f'(1) neighbours of N above N: {cnt} tested, max float margin {worst[0]:.6g} at change {worst[1]}; flagged {len(flag)}')
# odd CA numbers from c0 up to x_1 <= XCA (events sorted in mpmath)
P = oddprimes(XCA)
def F(p, a): return log(1 + mpf(1)/(p*sg(p, a - 1)))/log(p)
ev = []
for p in P:
    a = 1
    while F(p, a) >= F(XCA, 1): ev.append((F(p, a), p, a)); a += 1
ev.sort(key=lambda t: -t[0])
f = {}; tot = 0; wst = (-9, None); started = False; nca = 0
for e, p, a in ev:
    f[p] = a
    if (p, a) == (67, 1): started = True
    if not started: continue
    nca += 1
    c, w = scan(f, lN + 1e-9, f'CA x1={p}', flag, SP100); tot += c
    if w[0] > wst[0]: wst = (w[0], w[1], p, a)
print(f'(2) {nca} odd CA numbers from c0 (events up to x_1 = {XCA}): {tot} neighbours > N tested, max float margin {wst[0]:.6g} ({wst[1:]}); flagged total {len(flag)}')
for t in flag[:20]: print('   flagged:', t[0], t[1], f'margin {t[2]:.6g}', 'n/N =', t[3]/N)
viol = [t for t in flag if t[2] >= 0]
print(f'   exact violators among flagged with n > N: {len(viol)}')
# (3) odd superabundant records in (N, c0)
S = sorted((int(a), int(b)) for a, b in (l.split() for l in open(Sfile)))
best = None; rec = []
for n, s in S:
    if best is None or s*best[0] > best[1]*n:
        best = (n, s); rec.append((n, s))
recg = [(n, s) for n, s in rec if N < n < c0]
mr = [(float(log(mpf(s)/n) - (euler - log(2)) - log(log(log(mpf(n))))), n) for n, s in recg]
print(f'(3) odd SA records among S: {len(rec)} (N is one: {N in [r[0] for r in rec]}); in (N, c0): {len(recg)}, their max margin {max(mr)[0] if mr else None:.6g}; all satisfy (OR): {all(m < 0 for m, _ in mr)}')
for m, n in mr: print(f'   record n = {n} = N*{n//N if n % N == 0 else (str(n/N))}, margin {m:.6g}')
print(f'done {time.time()-t0:.1f}s')
