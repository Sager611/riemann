#!/usr/bin/env python3
"""gstar.py (ppraudit) -- exact least prime primitive root g*(p) for all primes 3 <= p <= N, segmented, vectorized.
Independent of PrimRootNum/small.py.  Usage: python3 gstar.py N [blocksize_log2]
Compares with Theorem A (kappa_w from PrimRootNum/logs/kappa_table.txt), B' (0.814, 3/b) and D (1.956, 3.4/b)."""
import sys, math, time
import numpy as np

N = int(float(sys.argv[1])); LB = int(sys.argv[2]) if len(sys.argv) > 2 else 20
DUMP = open(sys.argv[3], "w") if len(sys.argv) > 3 else None
S = 1 << LB
INV_B = 1.0 / (1.0 - 0.8733819666132076)
kap = {}
for line in open("../PrimRootNum/logs/kappa_table.txt"):
    if line.startswith("#"): continue
    f = line.split(); kap[int(f[0])] = float(f[6])
def base_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)
R = int(math.isqrt(N)) + 1
BP = base_primes(max(R, 3000))
SP = BP[BP <= R]
CAND = [int(x) for x in BP[:400]]          # candidate prime roots r (up to 2741)
def modpow(base, e, m):
    res = np.ones_like(m); b = base % m; e = e.copy()
    while True:
        nz = e > 0
        if not nz.any(): break
        odd = (e & 1) == 1
        res = np.where(odd, (res * b) % m, res)
        b = (b * b) % m; e >>= 1
    return res
t0 = time.time()
rec = 0; records = []; nprimes = 0
worstA = (0.0, None); worstB = (0.0, None); worstD = (0.0, None)
per_w = {}                                   # w -> (count, max g*, p at max)
hist = {}
lo = 0
while lo <= N:
    hi = min(lo + S, N + 1)
    m = np.arange(lo, hi, dtype=np.int64)
    isp = np.ones(hi - lo, dtype=bool)
    isp[m < 2] = False
    for l in SP:
        l = int(l)
        if l * l >= hi: break
        st = max(l * l, ((lo + l - 1) // l) * l)
        isp[st - lo::l] = False
    isp[m < 3] = False
    pidx = np.nonzero(isp)[0]
    if len(pidx):
        p = m[pidx]; n = p - 1
        rem = n.copy(); fac = np.zeros((len(p), 10), dtype=np.int64); cnt = np.zeros(len(p), dtype=np.int64)
        for l in SP:
            l = int(l)
            if l > hi: break
            sel = np.nonzero(rem % l == 0)[0] if l < 64 else np.nonzero(n % l == 0)[0]
            if len(sel) == 0: continue
            fac[sel, cnt[sel]] = l; cnt[sel] += 1
            while True:
                dv = sel[rem[sel] % l == 0]
                if len(dv) == 0: break
                rem[dv] //= l; sel = dv
        big = np.nonzero(rem > 1)[0]
        fac[big, cnt[big]] = rem[big]; cnt[big] += 1
        # check full factorization: product of prime powers = n is implied by rem -> 1 after big; verify primes only once
        g = np.zeros(len(p), dtype=np.int64)
        U = np.arange(len(p))
        for r in CAND:
            if len(U) == 0: break
            pu = p[U]; ok = (pu % r) != 0
            for j in range(int(cnt[U].max())):
                has = cnt[U] > j
                lj = np.where(has, fac[U, j], 1)
                e = np.where(has, (pu - 1) // lj, 0)
                pw = modpow(np.full(len(U), r, dtype=np.int64), e, pu)
                ok &= ~(has & (pw == 1))
            g[U[ok]] = r; U = U[~ok]
        if len(U): raise SystemExit("candidate list exhausted at p = %d" % p[U[0]])
        nprimes += len(p)
        L = np.log(p.astype(np.float64)); lg = np.log(g.astype(np.float64))
        for w in np.unique(cnt):
            w = int(w); sw = cnt == w
            rA = lg[sw] / (INV_B * np.log(kap[w] * L[sw]))
            i = int(np.argmax(rA))
            if rA[i] > worstA[0]: worstA = (float(rA[i]), (int(p[sw][i]), int(g[sw][i]), w))
            c0, gm, pm = per_w.get(w, (0, 0, 0)); j = int(np.argmax(g[sw]))
            if g[sw][j] > gm: gm, pm = int(g[sw][j]), int(p[sw][j])
            per_w[w] = (c0 + int(sw.sum()), gm, pm)
        s5 = p >= 5
        if s5.any():
            rB = lg[s5] / (3 * INV_B * np.log(0.814 * L[s5])); i = int(np.argmax(rB))
            if rB[i] > worstB[0]: worstB = (float(rB[i]), (int(p[s5][i]), int(g[s5][i])))
        rD = lg / (3.4 * INV_B * np.log(1.956 * L)); i = int(np.argmax(rD))
        if rD[i] > worstD[0]: worstD = (float(rD[i]), (int(p[i]), int(g[i])))
        if DUMP:
            DUMP.write("".join("%d %d\n" % (int(a), int(c)) for a, c in zip(p, g)))
        for k in range(len(p)):
            if g[k] > rec: rec = int(g[k]); records.append((int(p[k]), rec))
        for v, c in zip(*np.unique(g, return_counts=True)): hist[int(v)] = hist.get(int(v), 0) + int(c)
    lo = hi
print("N = %d: %d primes 3 <= p <= N, max g* = %d, %.1f s" % (N, nprimes, rec, time.time() - t0))
print("records:", ", ".join("%d:%d" % t for t in records))
print("per omega(p-1): count, max g*, at p:", {w: per_w[w] for w in sorted(per_w)})
print("worst log g*/log (kappa_w log p)^(1/b) [Thm A] = %.5f at (p, g*, w) = %s" % worstA)
print("worst log g*/log (0.814 log p)^(3/b) [Thm B', p >= 5] = %.5f at %s" % worstB)
print("worst log g*/log (1.956 log p)^(3.4/b) [Thm D] = %.5f at %s" % worstD)
print("g* histogram (top 12 values):", sorted(hist.items())[-12:])
