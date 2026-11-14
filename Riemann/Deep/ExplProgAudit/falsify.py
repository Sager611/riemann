#!/usr/bin/env python3
"""falsify.py (explaudit): smallest cases and boundaries of items (1)-(4) by exact computation (floats for the bounds;
margins are large). Not part of the proof."""
import math, numpy as np, mpmath as mp
B = 0.8733819666132076; EXP = 7.8978
N = 10 ** 6
sv = np.ones(N + 1, dtype=bool); sv[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if sv[i]: sv[i * i::i] = False
primes = np.nonzero(sv)[0]
lam = np.zeros(N + 1)
for p in primes:
    v = int(p)
    while v <= N: lam[v] = math.log(p); v *= int(p)
def phi(n): return sum(1 for k in range(1, n + 1) if math.gcd(k, n) == 1)
# (2) least primes, q = 2..300, incl. the sieve/analytic boundary q ~ 100
w = (0, None); wphi = (0, None); near100 = []
for q in range(2, 301):
    ph = phi(q); need = {a for a in range(q) if math.gcd(a, q) == 1}; P = 0
    for p in primes:
        r = int(p) % q
        if r in need: need.discard(r); P = int(p)
        if not need: break
    rr = math.log(P) / (EXP * math.log(q * math.log(q)))
    if rr > w[0]: w = (rr, q, P)
    if q >= 3:
        rp = math.log(P) / (EXP * math.log(ph * math.log(q)))
        if rp > wphi[0]: wphi = (rp, q, P)
    if 95 <= q <= 105: near100.append((q, P, round(rr, 4)))
print("(2) least primes q=2..300: worst log p/log((q log q)^7.8978) = %.4f (q=%d, p=%d); phi-version worst %.4f (q=%d, p=%d)" % (w + wphi))
print("    q=95..105 (q, max_a p(a,q), ratio):", near100)
# (3) non-residues: max over nonprincipal chi of n_chi = least N with <n<=N, (n,q)=1> = (Z/q)^*
wn = (1e9, None)
for q in range(3, 401):
    ph = phi(q); H = {1}; n = 1
    while len(H) < ph:
        n += 1
        if math.gcd(n, q) != 1 or (n % q) in H: continue
        new = set(H); g = n % q; cur = set(H)
        while True:
            cur = {(h * g) % q for h in cur}
            if cur <= new: break
            new |= cur
        H = new
    margin = EXP * math.log(math.log(q)) - math.log(n)
    if margin < wn[0]: wn = (margin, q, n)
print("(3) n_chi, q=3..400: tightest log-margin (log q)^7.8978 vs max_chi n_chi = %.4f at q=%d (n=%d)" % wn)
# (1) exact psi(x;q,a), q <= 10, all a, x at every integer n and n+1- in [max(2,q), 1e6]
xs = np.arange(2, N + 1, dtype=np.float64)
for q in range(1, 11):
    ph = phi(q); worst = 0; w64 = 0; w10 = 0
    for a in range(1, q + 1):
        if math.gcd(a, q) != 1: continue
        m = np.zeros(N + 1); idx = np.arange(N + 1) % q == (a % q); m[idx] = lam[idx]
        c = np.cumsum(m)[2:]
        for xx in (xs, xs + 1 - 1e-9):
            sel = xx >= max(2, q)
            r = np.abs(c - xx / ph) / (xx ** B * np.log(q * xx) ** 2)
            r = np.where(sel & (xx <= N), r, 0)
            worst = max(worst, r.max()); w64 = max(w64, np.where(xx >= 64, r, 0).max()); w10 = max(w10, np.where(xx >= math.exp(10), r, 0).max())
    print("(1) q=%2d: max R over x in [max(2,q),1e6] = %.4f (<=1.5); x>=64: %.5f (<=0.174); x>=e^10: %.6f (<=0.0324)" % (q, worst, w64, w10))
# (4) |L'/L(1,chi)| for primitive chi mod 3, 4, 5, 7, 8 (Hurwitz zeta off s = 1, symmetric differences; pole cancels)
mp.mp.dps = 50
def Lfun(q, chi, s): return sum(chi[a] * mp.zeta(s, mp.mpf(a) / q) for a in chi) / mp.mpf(q) ** s
tests = {3: {1: 1, 2: -1}, 4: {1: 1, 3: -1}, 5: {1: 1, 2: 1j, 4: -1, 3: -1j}, 8: {1: 1, 3: -1, 5: -1, 7: 1}, 7: {1: 1, 2: 1, 4: 1, 3: -1, 5: -1, 6: -1}}
h = mp.mpf(10) ** -15
for q, chi in tests.items():
    fp, fm = Lfun(q, chi, 1 + h), Lfun(q, chi, 1 - h)
    val = ((fp - fm) / (2 * h)) / ((fp + fm) / 2)
    bound = EXP * math.log(math.log(q)) + 10.53
    print("(4) q=%d: |L'/L(1,chi)| = %s <= %.3f" % (q, mp.nstr(abs(val), 6), bound))
