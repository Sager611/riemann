#!/usr/bin/env python3
# enum_s.py -- independent recount of the gap (oddraudit): best-first search over COMPLETE odd exponent vectors
# (each node fixes the exponent of the next odd prime, 0 allowed), exact integers and an exact rational threshold T0
# (from a hard-coded 60-digit gamma, cross-checked with mpmath), no arb, no eps-cost identity.
# Enumerates S = {odd m < c0 : sigma(m)/m >= T0}, T0 <= (e^g/2) log log N. Any violator of (OR) in (N, c0) lies in S.
# Then certifies the sign of the (OR) log-margin of every member of S with mpmath interval arithmetic (iv, 80 digits).
#   uv run --no-project --with mpmath python3 enum_s.py [out.txt]
import sys, heapq, time, math
from fractions import Fraction
from mpmath import mp, iv, mpf
t0 = time.time()
mp.dps = 80; iv.dps = 80
G60 = '0.577215664901532860606512090082402431042159335939923598805767'   # gamma, 60 digits (OEIS A001620)
assert abs(mp.euler - mpf(G60)) < mpf(10)**-59
gam = iv.mpf([mpf(G60) - mpf(10)**-58, mpf(G60) + mpf(10)**-58])   # contains gamma
def oddprimes(n):
    s = bytearray([1])*(n + 1); s[0] = s[1] = 0
    for i in range(2, math.isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(3, n + 1) if s[i]]
P = oddprimes(300000)
fac = {3: 4, 5: 3, 7: 2}
for p in P:
    if 11 <= p <= 61: fac[p] = 1
def sig(f):
    n = 1; s = 1
    for p, a in f.items(): n *= p**a; s *= (p**(a + 1) - 1)//(p - 1)
    return n, s
N, sN = sig(fac); c0 = 67*N; X = c0 - 1
assert N == 277093800961598968947975375 and c0 == 18565284664427130919514350125
# independent factor check of N by trial division
m = N; ff = {}
for p in [2] + P[:30]:
    while m % p == 0: ff[p] = ff.get(p, 0) + 1; m //= p
assert m == 1 and ff == fac
def margin_iv(n, s):        # log(s/n) - log((e^g/2) log log n), interval
    return iv.log(iv.mpf(s)/iv.mpf(n)) - (gam + iv.log(iv.log(iv.log(iv.mpf(n)))) - iv.log(2))
thr = iv.exp(gam)/2*iv.log(iv.log(iv.mpf(N)))
lo = mpf(thr.a)
Td = 10**40; Tn = int(mp.floor(lo*Td)) - 1           # T0 = Tn/Td < (e^g/2) log log N, exact
assert iv.mpf(Tn)/Td < thr
def G(i, R):                                         # (num, den) >= sigma(k)/k for odd k <= R with primes >= P[i]
    num = 1; den = 1; prod = 1; j = i
    while True:
        p = P[j]
        if prod*p > R: return num, den
        prod *= p; num *= p; den *= p - 1; j += 1
S = {}; nodes = 0; pushes = 0; maxheap = 0
def ok_bound(i, m, s):
    num, den = G(i, X//m)
    return s*num*Td >= m*den*Tn, (s*num)/(m*den)
heap = []; cnt = 0
okb, fb = ok_bound(0, 1, 1)
heapq.heappush(heap, (-fb, cnt, 0, 1, 1))
while heap:
    negb, _, i, m, s = heapq.heappop(heap); nodes += 1
    p = P[i]; R = X//m
    if p > R: continue                               # nothing can be added (m already recorded at creation)
    # zero exponent for p
    okb, fb = ok_bound(i + 1, m, s)
    if okb: cnt += 1; heapq.heappush(heap, (-fb, cnt, i + 1, m, s)); pushes += 1
    pe = p; se = 1 + p
    while m*pe <= X:
        m2 = m*pe; s2 = s*se
        if s2*Td >= m2*Tn: S[m2] = s2              # recorded once: created only here (exponent of p >= 1)
        okb, fb = ok_bound(i + 1, m2, s2)
        if okb: cnt += 1; heapq.heappush(heap, (-fb, cnt, i + 1, m2, s2)); pushes += 1
        pe *= p; se += pe
    maxheap = max(maxheap, len(heap))
print(f'best-first: {nodes} nodes popped, {pushes} pushes, max heap {maxheap}; |S| = {len(S)}; {time.time()-t0:.1f}s', flush=True)
gapS = sorted(m for m in S if N < m < c0); below = sorted(m for m in S if m <= N)
print(f'S in (N, c0): {len(gapS)};  S with m <= N: {len(below)}; N in S: {N in S}')
viol_gap = []; undec = []; worst = None; viol_below = []
for m in sorted(S):
    mg = margin_iv(m, S[m])
    if mg.b < 0: pass
    elif mg.a > 0:
        (viol_gap if m > N else viol_below).append(m)
    else: undec.append(m)
    if N < m < c0 and (worst is None or mg.b > worst[0].b): worst = (mg, m)
print(f'violators in (N, c0): {len(viol_gap)}; undecided: {len(undec)}; violators <= N in S: {len(viol_below)} (largest {max(viol_below) if viol_below else None})')
print(f'largest margin in gap: {mp.nstr(mpf(worst[0].b), 10)} at n = {worst[1]} = N*{Fraction(worst[1], N)}')
if len(sys.argv) > 1:
    with open(sys.argv[1], 'w') as f:
        for m in sorted(S): f.write(f'{m} {S[m]}\n')
print('PASS' if not viol_gap and not undec else 'FAIL', f'{time.time()-t0:.1f}s')
