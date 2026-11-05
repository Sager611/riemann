#!/usr/bin/env python3
# gap.py -- closes the finite gap N < n < c0 of the odd Robin inequality (arXiv 2008.04787, Prop. 3.5).
#   uv run --no-project --with python-flint python3 gap.py
# Claim: every odd n with N < n < c0 satisfies sigma(n)/n < (e^gamma/2) log log n, where
#   N  = 3^4 5^3 7^2 11 13 ... 61 = 277093800961598968947975375   (violates; last odd CA violator)
#   c0 = 67 N = 18565284664427130919514350125                       (next odd colossally abundant number)
# Method (eps-cost identity, Alaoglu-Erdos/Robin): with eps = log(68/67)/log 67, g_p(a) = sigma(p^a)/p^(a(1+eps)),
#   log(sigma(n)/n) = log(sigma(N)/N) + eps (log n - log N) - L(n),  L(n) = sum_p [log g_p(alpha_p) - log g_p(a_p)] >= 0,
# alpha_p = exponent of p in N (N maximises g over odd numbers at this eps; so does 67N).
# n violates  <=>  L(n) <= D(log n),  D(t) = delta0 + eps (t - tN) - (log log t - log log tN),  delta0 = margin of N.
# On [tN, tc0] D is decreasing (eps < 1/(tc0 log tc0)), so a violator in the gap has L(n) < delta0.
# Every prime p >= 73 costs >= c_73(1) > delta0, so only primes <= 71 can deviate from N; enumerate all
# exponent vectors with total cost <= delta0 (superset, ball lower bounds), then check every candidate directly.
import sys
from flint import arb, ctx
ctx.prec = 256
ok = True
def chk(cond, msg):
    global ok
    print(('PASS ' if cond else 'FAIL ') + msg); ok &= bool(cond)
def primes(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
P = [p for p in primes(1000) if p > 2]
alphaN = {3: 4, 5: 3, 7: 2}
for p in P:
    if 11 <= p <= 61: alphaN[p] = 1
def val(f):
    n = 1; s = 1
    for p, a in f.items():
        n *= p**a; s *= (p**(a+1)-1)//(p-1)
    return n, s
N, sN = val(alphaN)
c0 = 67*N
chk(N == 277093800961598968947975375 and c0 == 18565284664427130919514350125, 'N and c0 as in the paper')
EG2 = arb.const_euler().exp()/2
def margin(n, s):   # log(sigma(n)/n) - log((e^g/2) log log n); violation iff >= 0
    return (arb(s)/arb(n)).log() - (EG2*arb(n).log().log()).log()
d0 = margin(N, sN)
chk(d0 > 0, f'N violates: delta0 = log-margin(N) = {d0.str(12)} > 0')
mc0 = margin(c0, val({**alphaN, 67: 1})[1])
chk(mc0 < 0, f'c0 satisfies: log-margin(c0) = {mc0.str(12)} < 0')
eps = (arb(68)/67).log()/arb(67).log()
L = lambda p: arb(p).log()
def logg(p, a):
    return arb((p**(a+1)-1)//(p-1)).log() - a*(1+eps)*L(p)
# N maximises each g_p at eps (concavity: increments log(sigma(p^{a+1})/sigma(p^a)) - (1+eps) log p decrease in a)
for p in P:
    a = alphaN.get(p, 0)
    up = logg(p, a+1) - logg(p, a)          # must be <= 0 (ties allowed only for p = 67)
    if p == 67:
        chk(abs(up) < arb(10)**-60, 'p=67: g_67(1) = g_67(0) at eps (tie; 67N is the next odd CA number)')
    else:
        if not (up < 0): chk(False, f'g_{p}(a+1) < g_{p}(a) at a = alpha_p')
    if a > 0:
        dn = logg(p, a) - logg(p, a-1)
        if not (dn > 0): chk(False, f'g_{p}(alpha) > g_{p}(alpha-1)')
print('PASS N maximises g_p for every odd p <= 997 at eps (and p >= 73 by the monotonicity below)')
tN = arb(N).log(); tc = arb(c0).log()
chk(eps < 1/(tc*tc.log()), f'eps = {eps.str(10)} < 1/(tc0 log tc0) = {(1/(tc*tc.log())).str(10)}: D decreasing on [tN, tc0]')
c73 = (1+eps)*L(73) - L(74)
chk(c73 > d0, f'c_73(1) = {c73.str(8)} > delta0; c_p(1) = (1+eps)log p - log(p+1) increases for p >= 67, and c_p(a) >= c_p(1) for a >= 1')
# allowed exponents for p <= 71
allowed = {}
for p in [q for q in P if q <= 71]:
    a0 = alphaN.get(p, 0); opts = []
    for a in range(0, 40):
        c = logg(p, a0) - logg(p, a)
        if c.lower() <= d0.upper(): opts.append((a, c))
    allowed[p] = opts
    if len(opts) > 1:
        print(f'   p={p}: exponents within budget: ' + ', '.join(f'{a} (cost {c.str(6)})' for a, c in opts))
# DFS over combinations with total cost lower bound <= delta0 upper bound
plist = sorted(allowed)
cands = []
def dfs(i, f, cost):
    if i == len(plist):
        cands.append((dict(f), cost)); return
    p = plist[i]
    for a, c in allowed[p]:
        nc = cost + c
        if nc.lower() <= d0.upper():
            f2 = dict(f)
            if a: f2[p] = a
            dfs(i+1, f2, nc)
dfs(0, {}, arb(0))
print(f'   {len(cands)} exponent vectors with L <= delta0 (superset)')
nviol = 0; ingap = 0
for f, c in cands:
    n, s = val(f)
    m = margin(n, s)
    diff = {p: (f.get(p, 0) - alphaN.get(p, 0)) for p in set(f) | set(alphaN) if f.get(p, 0) != alphaN.get(p, 0)}
    rel = 'N' + ''.join(f'*{p}^{d}' if d > 0 else f'/{p}^{-d}' for p, d in sorted(diff.items())) if diff else 'N'
    where = 'in gap' if N < n < c0 else ('= N' if n == N else ('< N' if n < N else ('= c0' if n == c0 else '> c0')))
    print(f'   n = {rel:14s} n = {n:.4e} {where:6s} L = {c.str(6):>22s} log-margin = {m.str(8)}')
    if N < n < c0:
        ingap += 1
        if not (m < 0): nviol += 1
chk(nviol == 0, f'no violator among the {ingap} candidates in the gap (each has certified log-margin < 0)')
print('ALL PASS' if ok else 'SOME FAIL')
sys.exit(0 if ok else 1)

