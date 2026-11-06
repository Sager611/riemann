#!/usr/bin/env python3
"""sanity.py -- numerical sanity checks (not part of the proof except item 1 for 2 <= q < 100).
  1. exact least primes p(a,q) for all 2 <= q <= 10^4 and all (a,q) = 1 (sieve), against (q log q)^(1/(1-B)) and,
     for q >= 3, (phi(q) log q)^(1/(1-B));
  2. exact psi(x;q,a) for q <= 100, all a, at 400 geometric points x in [2, 10^8], against 1.5 x^B log^2(qx)
     and the x >= 64 constant 0.174;
  3. exact max_chi n_chi (least N such that {n <= N, (n,q)=1} generates (Z/q)^*) for 3 <= q <= 10^4.
Run: OMP_NUM_THREADS=1 uv run --no-project --with numpy python3 sanity.py   (~1 GB peak at most, one thread)
"""
import math, time, sys
import numpy as np
B = 0.8733819666132076; b = 1 - B; INVB = 1 / b
t0 = time.time()
N = 10 ** 8
sv = np.ones(N + 1, dtype=np.bool_); sv[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if sv[i]: sv[i * i::i] = False
primes = np.nonzero(sv)[0].astype(np.int64)
del sv
print("sieve: %d primes <= 1e8 (%.1f s)" % (len(primes), time.time() - t0), flush=True)

def phi(n):
    r = n; m = n; p = 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0: m //= p
            r -= r // p
        p += 1
    if m > 1: r -= r // m
    return r

# ---- 1. least primes
worst_exp = (0, None); worst_ratio = (0, None); fail = 0; failphi = 0; worst_phi = (0, None)
for q in range(2, 10 ** 4 + 1):
    ph = phi(q); m = int(3 * ph * math.log(ph + 2)) + 200
    while True:
        r = primes[:m] % q
        u, idx = np.unique(r, return_index=True)
        cop = np.array([math.gcd(int(x), q) == 1 for x in u])
        if cop.sum() == ph: break
        m *= 2
        if m > len(primes): raise SystemExit("not all classes found q=%d" % q)
    P = int(primes[idx[cop]].max())
    e = math.log(P) / math.log(q)
    if e > worst_exp[0]: worst_exp = (e, (q, P))
    bound_log = INVB * math.log(q * math.log(q))
    rr = math.log(P) / bound_log
    if rr > worst_ratio[0]: worst_ratio = (rr, (q, P))
    if math.log(P) >= bound_log: fail += 1
    if q >= 3:
        bphi = INVB * math.log(ph * math.log(q))
        if math.log(P) >= bphi: failphi += 1
        if math.log(P) / bphi > worst_phi[0]: worst_phi = (math.log(P) / bphi, (q, P))
print("least primes 2 <= q <= 1e4: max log P/log q = %.4f at (q, P) = %s; max log P / log((q log q)^(1/(1-B))) = %.4f at %s; failures %d (%.1f s)"
      % (worst_exp[0], worst_exp[1], worst_ratio[0], worst_ratio[1], fail, time.time() - t0), flush=True)
print("  phi-version 3 <= q <= 1e4: max log P / log((phi(q) log q)^(1/(1-B))) = %.4f at %s; failures %d" % (worst_phi[0], worst_phi[1], failphi), flush=True)
# certificate part: 2 <= q < 100 (q = 2: (2 log 2)^(1/b) = 13.2 > 3); phi-version for 3 <= q < 100
assert fail == 0 and failphi == 0

# ---- 2. psi(x;q,a)
lp = np.log(primes.astype(np.float64))
# prime powers p^k <= 1e8, k >= 2
pp = []; 
for p in primes[primes <= 10 ** 4]:
    p = int(p); v = p * p
    while v <= N: pp.append((v, math.log(p))); v *= p
ppv = np.array([v for v, _ in pp], dtype=np.int64); ppl = np.array([l for _, l in pp])
allv = np.concatenate([primes, ppv]); alll = np.concatenate([lp, ppl])
o = np.argsort(allv, kind='stable'); allv = allv[o]; alll = alll[o]
xs = np.unique(np.concatenate([np.geomspace(2, N, 400), np.arange(2, 200, 1.0)]))
worst = (0, None); worst64 = (0, None); worstsq = (0, None)
for q in range(1, 101):
    ph = phi(q)
    res = allv % q
    for a in range(q):
        if math.gcd(a, q) != 1: continue
        sel = res == (a % q)
        v = allv[sel]; c = np.cumsum(alll[sel])
        k = np.searchsorted(v, xs, side='right')
        psi = np.where(k > 0, c[np.maximum(k - 1, 0)], 0.0)
        # left limits at x (psi(x-)) cover the interval sup, together with the right values
        E = np.abs(psi - xs / ph)
        den = xs ** B * np.log(q * xs) ** 2
        r = E / den
        i = int(np.argmax(r))
        if r[i] > worst[0]: worst = (float(r[i]), (q, a, float(xs[i])))
        m64 = xs >= 64
        i2 = int(np.argmax(np.where(m64, r, 0)))
        if r[i2] > worst64[0]: worst64 = (float(r[i2]), (q, a, float(xs[i2])))
        rs = E / (np.sqrt(xs) * np.log(q * xs) ** 2)
        i3 = int(np.argmax(np.where(xs >= 100, rs, 0)))
        if rs[i3] > worstsq[0]: worstsq = (float(rs[i3]), (q, a, float(xs[i3])))
print("psi(x;q,a), q <= 100, %d points in [2,1e8]: max |E|/(x^B log^2 qx) = %.4f at %s (bound 1.5); x >= 64: %.4f at %s (bound 0.174);"
      " max |E|/(sqrt(x) log^2 qx), x >= 100: %.4f at %s (%.1f s)"
      % (len(xs), worst[0], worst[1], worst64[0], worst64[1], worstsq[0], worstsq[1], time.time() - t0), flush=True)
assert worst[0] < 1.5 and worst64[0] < 0.174

# ---- 3. character non-residues: max over nonprincipal chi of n_chi = least N with <n <= N coprime> = (Z/q)^*
worstn = (0, None)
for q in range(3, 10 ** 4 + 1):
    ph = phi(q)
    inH = np.zeros(q, dtype=np.bool_); inH[1] = True; size = 1; elems = [1]
    n = 1
    while size < ph:
        n += 1
        if math.gcd(n, q) != 1 or inH[n % q]: continue
        # add n: H <- H * <n>
        new = []
        cur = list(elems)
        g = n % q; mult = g
        while not inH[mult]:
            for h in cur:
                e = (h * mult) % q
                if not inH[e]: inH[e] = True; new.append(e)
            mult = (mult * g) % q
        elems += new; size = len(elems)
    rr = math.log(n) / (INVB * math.log(math.log(q)))
    if n > worstn[0]: worstn = (n, q)
print("non-residues 3 <= q <= 1e4: max_chi n_chi <= %d (at q = %d); (log q)^(1/(1-B)) >= (log 3)^7.898 = 2.10 for every q >= 3 (%.1f s)"
      % (worstn[0], worstn[1], time.time() - t0), flush=True)
