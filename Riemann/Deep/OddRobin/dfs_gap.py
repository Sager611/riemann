#!/usr/bin/env python3
# dfs_gap.py -- independent check of the gap N < n < c0 by exhaustive branch-and-bound (exact integers; no eps-cost identity).
#   uv run --no-project --with python-flint python3 dfs_gap.py
# Lists every odd n with N < n < c0 and sigma(n)/n >= T0 := a rational lower bound of (e^gamma/2) log log N (any violator in the
# gap has sigma(n)/n >= (e^g/2) log log n > T0), then checks each against its own threshold in arb.
# Pruning: for n = m*k, k built from primes >= p_i, n < c0: sigma(n)/n <= (sigma(m)/m) * prod of the w smallest primes >= p_i of
# p/(p-1), w = the largest number of distinct primes >= p_i whose product is <= (c0-1)/m (k has at most w distinct prime factors).
import sys, time, math
from flint import arb, ctx
ctx.prec = 200
t0 = time.time()
N = 277093800961598968947975375; c0 = 67*N; X = c0 - 1
def primes(n):
    s = bytearray([1])*(n+1); s[0] = s[1] = 0
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(3, n+1) if s[i]]
P = primes(200000)
EG2 = arb.const_euler().exp()/2
thr = EG2*arb(N).log().log()
from decimal import Decimal, getcontext
getcontext().prec = 80
D = 10**40
Tn = int(Decimal(thr.mid().str(60, radius=False))*D) - 10**6
assert arb(Tn)/D < thr and (thr - arb(Tn)/D) < arb(10)**-30
Td = D                                                  # T0 = Tn/Td <= (e^g/2) log log N
def ub(i, R):                                           # (num, den): prod_{j<w} p_{i+j}/(p_{i+j}-1), w max with prod p_{i+j} <= R
    num = 1; den = 1; prod = 1; j = i
    while True:
        if j >= len(P): raise RuntimeError('prime list too short')
        p = P[j]
        if prod*p > R: return num, den
        prod *= p; num *= p; den *= p - 1; j += 1
cands = []; nodes = 0
sys.setrecursionlimit(10000)
def dfs(i, m, s):
    global nodes
    nodes += 1
    if m > N and s*Td >= m*Tn: cands.append((m, s))
    R = X//m
    j = i
    while j < len(P):
        p = P[j]
        if p > R: return
        num, den = ub(j, R)
        if s*num*Td < m*den*Tn: return                    # nothing with primes >= p_j can reach T0 (ub decreases in j)
        pe = p; se = 1 + p
        while pe <= R:
            dfs(j + 1, m*pe, s*se)
            pe *= p; se += pe
        j += 1
    raise RuntimeError('prime list too short')
dfs(0, 1, 1)
print(f'{nodes} search nodes, {len(cands)} odd n in (N, c0) with sigma(n)/n >= T0 = (e^g/2) log log N - <1e-30; {time.time()-t0:.1f}s')
bad = 0; best = None
for m, s in cands:
    marg = (arb(s)/m).log() - (EG2*arb(m).log().log()).log()
    if not (marg < 0): bad += 1; print('VIOLATOR or undecided:', m)
    if best is None or marg.upper() > best[0].upper(): best = (marg, m)
for m, s in sorted(cands)[:12]:
    print(f'   n = {m}  n/N = {m/N:.6g}  sigma/n = {s/m:.8f}')
if best: print(f'largest log-margin among them: {best[0].str(8)} at n = {best[1]} (n/N = {best[1]/N:.6g})')
print('PASS: no odd n with N < n < c0 violates (OR)' if bad == 0 else 'FAIL')
sys.exit(0 if bad == 0 else 1)
