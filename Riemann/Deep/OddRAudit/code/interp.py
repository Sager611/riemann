#!/usr/bin/env python3
# interp.py -- oddraudit tests of PROOF.md §2 (WY Lemma 3.4 / Robin Prop. 1) and of the odd CA sequence below c0, built from scratch:
#  (1) events (p,a), p <= PMAX, sorted by F(p,a) (mpmath 60 digits), certified gaps; the odd CA sequence; checks that N and c0 are
#      consecutive, that g = sigma(m)/m^(1+eps) ties at each event eps, and lists which odd CA numbers <= c0 satisfy (OR);
#  (2) for every odd n <= XS (divisor sieve) and every n of S (enum_s output): with N1 <= n <= N2 consecutive odd CA numbers and
#      eps their tie, log(sigma(n)/n) <= log(sigma(N1)/N1) + eps (log n - log N1) (the §2 line), and Lemma 3.4: if n violates (OR),
#      N1 or N2 violates.
#   uv run --no-project --with mpmath --with numpy python3 interp.py XS S.txt
import sys, math, time, bisect
import numpy as np
from mpmath import mp, mpf, log, euler, exp
t0 = time.time(); mp.dps = 60
XS = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**6
Sfile = sys.argv[2] if len(sys.argv) > 2 else None
PMAX = 2000
def oddprimes(n):
    s = bytearray([1])*(n + 1); s[0] = s[1] = 0
    for i in range(2, math.isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(3, n + 1) if s[i]]
P = oddprimes(PMAX)
def sg(p, a): return (p**(a + 1) - 1)//(p - 1)
def F(p, a): return log(1 + mpf(1)/(p*sg(p, a - 1)))/log(p)
ev = []
for p in P:
    a = 1
    while True:
        f = F(p, a)
        if f < F(PMAX, 1): break
        ev.append((f, p, a)); a += 1
ev.sort(key=lambda t: -t[0])
gaps = [(ev[i][0] - ev[i+1][0])/ev[i][0] for i in range(len(ev) - 1)]
print(f'{len(ev)} events with F >= F({PMAX},1); min relative gap between consecutive events {mp.nstr(min(gaps), 5)} (no ties)')
LC = euler - log(2)
def marg(n, s): return log(mpf(s)/n) - LC - log(log(log(mpf(n))))
seq = []; expo = {}; n = 1; s = 1          # odd CA numbers in order: (n, s, eps_entry) eps_entry = event eps at which it starts
for f, p, a in ev:
    n *= p; s = s//sg(p, a - 1)*sg(p, a); expo[p] = a
    seq.append((n, s, f))
N = 277093800961598968947975375; c0 = 67*N
iN = [k for k, t in enumerate(seq) if t[0] == N][0]
print(f'N is odd CA #{iN} (event {ev[iN][1:]}), next: {seq[iN+1][0] == c0} = c0 (event {ev[iN+1][1:]})')
sat = [k for k in range(len(seq)) if seq[k][0] > 3 and seq[k][0] <= c0 and marg(seq[k][0], seq[k][1]) < 0]
print('odd CA numbers <= c0 satisfying (OR):', [(seq[k][0], ev[k][1:]) for k in sat])
viol_after = [k for k in range(iN + 1, len(seq)) if marg(seq[k][0], seq[k][1]) >= 0]
print(f'odd CA numbers after N (x_1 <= {PMAX}) violating (OR): {len(viol_after)}')
# tie check: g(N_k) = g(N_{k+1}) at eps = event of N_{k+1}
bad = 0
for k in range(len(seq) - 1):
    (n1, s1, _), (n2, s2, e) = seq[k], seq[k+1]
    d = log(mpf(s1)/n1) - e*log(n1) - (log(mpf(s2)/n2) - e*log(n2))
    if abs(d) > mpf(10)**-50: bad += 1
print(f'tie identity g(N_k) = g(N_k+1) at the event eps: {len(seq)-1-bad} of {len(seq)-1} pairs exact to 1e-50')
CA = [t[0] for t in seq]
def line_check(n, s, tol=mpf(10)**-40):
    k = bisect.bisect_right(CA, n) - 1             # CA[k] <= n < CA[k+1]
    n1, s1, _ = seq[k]; n2, s2, e = seq[k+1]
    lhs = log(mpf(s)/n); rhs = log(mpf(s1)/n1) + e*(log(mpf(n)) - log(mpf(n1)))
    ok_line = lhs <= rhs + tol
    v = marg(n, s) >= 0
    ok34 = (not v) or marg(n1, s1) >= 0 or marg(n2, s2) >= 0
    return ok_line, ok34, v
# (2a) sieve
sig = np.zeros(XS + 1, dtype=np.int64)
for d in range(1, XS + 1, 2): sig[d::2*d] += d
nl = n34 = nv = 0; worst = None
for n in range(3, XS + 1, 2):
    s = int(sig[n])
    k = bisect.bisect_right(CA, n) - 1
    n1, s1, _ = seq[k]; n2, s2, e = seq[k+1]
    lhs = math.log(s/n); rhs = math.log(s1/n1) + float(e)*(math.log(n) - math.log(n1))
    if lhs > rhs + 1e-12: nl += 1
    th = math.exp(float(euler))/2*math.log(math.log(n))
    if s/n >= th:
        nv += 1
        if not (s1/n1 >= math.exp(float(euler))/2*math.log(math.log(n1)) or s2/n2 >= math.exp(float(euler))/2*math.log(math.log(n2))): n34 += 1
print(f'odd 3 <= n <= {XS}: line violations {nl}; (OR) violators {nv}; Lemma 3.4 failures {n34}; {time.time()-t0:.1f}s')
if Sfile:
    nl = n34 = nv = 0; cnt = 0
    for line in open(Sfile):
        n, s = map(int, line.split()); cnt += 1
        a, b, v = line_check(n, s); nl += (not a); n34 += (not b); nv += v
    print(f'S ({cnt} odd n < c0, sigma/n >= T0): line violations {nl}; (OR) violators {nv}; Lemma 3.4 failures {n34}')
print(f'done {time.time()-t0:.1f}s')
