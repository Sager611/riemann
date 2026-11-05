#!/usr/bin/env python3
# oddca.py -- direct check of the odd Robin inequality sigma(N)/N < (e^gamma/2) log log N at every odd colossally abundant
# (odd CA) number N with x_1(eps) <= XD (default 5e8), in arb ball arithmetic (128 bits).
#   uv run --no-project --with python-flint --with numpy python3 oddca.py [XD]
# Odd CA numbers (arXiv 2008.04787, Lemma 2.2): N_eps = prod_{p odd} p^{alpha_p(eps)}. As eps decreases, the exponent of p
# rises from a-1 to a at eps(p,a) = F(p,a) = log(1 + 1/(p + ... + p^a))/log p. Sorting all events (p,a) by decreasing eps
# lists every odd CA number; x_1(eps) is defined by F(x_1,1) = eps, so x_1 <= XD <=> eps >= F(XD,1).
# Events (p,1) come in increasing p; the few events with a >= 2 are merged in at certified positions.
# Below p = 1e5 every number is checked exactly; above, a run of consecutive (p,1) events p_i..p_j is checked as a block:
# every number in the block has S <= S_after and T >= T_before (S = log sigma(N)/N, T = log N both increase), so
# S_after - log(e^gamma/2) - log log T_before < 0 proves them all. Blocks that fail are halved down to single numbers.
import sys, time, math
import numpy as np
from flint import arb, ctx
t0 = time.time()
ctx.prec = 128
XD = int(float(sys.argv[1])) if len(sys.argv) > 1 else 5*10**8
assert XD < 2**32
def primes_upto(X):                     # odd primes <= X, segmented sieve, uint32
    r = math.isqrt(X) + 1
    b = np.ones(r + 1, dtype=bool); b[:2] = False
    for i in range(2, math.isqrt(r) + 1):
        if b[i]: b[i*i::i] = False
    base = np.nonzero(b)[0].tolist()
    out = []; seg = 10**7
    for lo in range(0, X + 1, seg):
        hi = min(lo + seg, X + 1)
        s = np.ones(hi - lo, dtype=bool)
        if lo == 0: s[:3] = False        # 0, 1, 2 (2 excluded: odd primes only)
        for p in base:
            if p*p >= hi: break
            st = max(p*p, ((lo + p - 1)//p)*p)
            s[st - lo::p] = False
        out.append((np.nonzero(s)[0] + lo).astype(np.uint32))
    a = np.concatenate(out)
    return a[a % 2 == 1]
PR = primes_upto(XD)
print(f'{len(PR)} odd primes <= {XD:.3g}; {time.time()-t0:.1f}s', flush=True)
def sig1(p, a): return (p**(a + 1) - 1)//(p - 1)          # sigma(p^a)
def F(p, a): return math.log1p(1/(p*sig1(p, a - 1)))/math.log(p)
def Farb(p, a): return (1/arb(p*sig1(p, a - 1))).log1p()/arb(p).log()
epsD = F(XD, 1); FD = Farb(XD, 1)
# events with a >= 2, eps >= F(XD,1)
hi_ev = []
for p in PR[:20000].tolist():               # F(p,2) < F(XD,1) long before the 20000th prime for XD < 2^32
    if F(p, 2) < epsD*(1 - 1e-9): break
    a = 2
    while F(p, a) >= epsD*(1 - 1e-9):
        if Farb(p, a) >= FD or not (Farb(p, a) < FD): hi_ev.append((F(p, a), p, a))
        a += 1
assert F(int(PR[min(len(PR), 20000) - 1]), 2) < epsD*(1 - 1e-9) or len(PR) < 20000
hi_ev.sort(key=lambda t: -t[0])
for i in range(len(hi_ev) - 1):              # certify the order of near-ties
    if hi_ev[i][0] - hi_ev[i+1][0] < 1e-9*hi_ev[i][0]:
        if not (Farb(*hi_ev[i][1:]) - Farb(*hi_ev[i+1][1:]) > 0):
            print('FAIL: order of', hi_ev[i], hi_ev[i+1]); sys.exit(1)
def x1_of(e):                                # float x with F(x,1) = e (F decreasing)
    lo, hi = 2.0, 1e12
    for _ in range(200):
        mid = math.sqrt(lo*hi)
        if math.log1p(1/mid)/math.log(mid) > e: lo = mid
        else: hi = mid
    return lo
pos = []                                     # index in PR of the first prime q with F(q,1) < eps_e
for e, p, a in hi_ev:
    k = int(np.searchsorted(PR, np.uint32(min(math.ceil(x1_of(e)), 2**32 - 1))))
    for _ in range(5):
        E = Farb(p, a)
        okL = k == 0 or Farb(int(PR[k-1]), 1) > E
        okR = k == len(PR) or Farb(int(PR[k]), 1) < E
        if okL and okR: break
        k += -1 if not okL else 1
    else:
        print('FAIL: position of', (p, a)); sys.exit(1)
    pos.append(k)
assert all(pos[i] <= pos[i+1] for i in range(len(pos) - 1))
print(f'{len(hi_ev)} events with exponent >= 2 merged at certified positions; {time.time()-t0:.1f}s', flush=True)
LC = (arb.const_euler().exp()/2).log()
S = arb(0); T = arb(0); started = False; lastviol = None; nexact = 0; nblock = 0; nnum = 0; kblk = 8; worst = None; X1 = [3]
def exact(p, a):
    global S, T, started, lastviol, nexact, worst
    S += (1/arb(p*sig1(p, a - 1))).log1p(); T += arb(p).log(); expo[p] = a
    if a == 1: X1[0] = p                     # largest prime of N
    if (p, a) == (67, 1): started = True
    if T.lower() <= 1: return
    m = S - LC - T.log().log()
    nexact += 1
    if m > 0:
        if started: print('FAIL: odd CA violator after c0 at event', (p, a)); sys.exit(1)
        lastviol = (p, a, math.prod(q**b for q, b in expo.items()) if p < 1000 else None); return
    if not (m < 0): print('FAIL: undecided margin at event', (p, a)); sys.exit(1)
    if started and X1[0] > 100:
        P = X1[0]; v = float((m*math.sqrt(P)*math.log(P)).upper())
        if worst is None or v > worst[0]: worst = (v, p, a, P)
expo = {}
def run1(i0, i1):                            # (p,1) events for PR[i0:i1]
    global S, T, nblock, nnum, kblk
    if i1 > i0: X1[0] = max(X1[0], int(PR[i1 - 1]))
    i = i0
    while i < i1 and (int(PR[i]) <= 10**5 or not started):
        p = int(PR[i]); exact(p, 1); expo[p] = 1; i += 1
    while i < i1:
        j = min(i + kblk, i1)
        seg = PR[i:j].astype(np.uint64)
        lP = arb(math.prod(seg.tolist())).log()
        S2 = S + arb(math.prod((seg + 1).tolist())).log() - lP; T2 = T + lP
        if S2 - LC - T.log().log() < 0:
            S, T = S2, T2; nblock += 1; nnum += j - i; i = j; kblk = kblk + kblk//4 + 1
        elif j - i == 1:
            if S2 - LC - T2.log().log() < 0: S, T = S2, T2; nblock += 1; nnum += 1; i = j
            else: print('FAIL: single number not certified at p =', int(PR[i])); sys.exit(1)
        else:
            kblk = max(1, (j - i)//2)
cur = 0
for (e, p, a), k in zip(hi_ev, pos):
    run1(cur, k); cur = k
    exact(p, a); expo[p] = a
run1(cur, len(PR))
ok = started and lastviol is not None and lastviol[2] == 277093800961598968947975375 and T.rad() < 1e-20
print(f'last odd CA violator: event {lastviol[:2]}, N = {lastviol[2]} (= 3^4 5^3 7^2 11 ... 61, as in the paper)')
print(f'{"PASS" if ok else "FAIL"}: every odd CA number N >= c0 = 67N with x_1 <= {XD:.3g} satisfies the odd Robin inequality '
      f'({nexact} numbers checked one by one, {nnum} more in {nblock} blocks)')
print(f'   closest normalised margin among exact checks (P = largest prime > 100): m*sqrt(P)*log(P) = {worst[0]:.4f} at event {worst[1:3]}, P = {worst[3]}')
print(f'   log N at the last odd CA number checked: {T.str(15)}; {time.time()-t0:.1f}s')
sys.exit(0 if ok else 1)

