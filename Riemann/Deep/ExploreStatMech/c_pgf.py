
# C: Lee-Yang zeros of the squarefree prime-factor gas  P_x(z) = sum_{n<=x} mu(n)^2 z^omega(n)
#    (and A_x(z) = sum_{n<=x} z^omega(n)); test real-rootedness, locate the zero near z=-1 vs Mertens M(x).
import numpy as np, sys, time
X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**7
EXH = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**6
t0 = time.time()
def primes_upto(n):
    s = np.ones(n+1, dtype=bool); s[:2] = False
    for p in range(2, int(n**0.5)+1):
        if s[p]: s[p*p::p] = False
    return np.nonzero(s)[0]
P = primes_upto(X)
om = np.zeros(X+1, dtype=np.int8)
for p in P: om[p::p] += 1
sqf = np.ones(X+1, dtype=bool); sqf[0] = False
for p in P[P*P <= X]: sqf[p*p::p*p] = False
mu = np.where(sqf, np.where(om % 2 == 0, 1, -1), 0).astype(np.int8); mu[0] = 0
print('sieve done', round(time.time()-t0,1), 's; max omega', om.max(), flush=True)
K = int(om.max()) + 1
def rr_check(C, tol=1e-7):
    # C: (m, K) integer coeff arrays (coef of z^k); returns boolean real-rooted for each row (numpy eig) and max rel imag
    out = np.zeros(len(C), dtype=bool); mx = np.zeros(len(C))
    deg = np.array([np.max(np.nonzero(row)[0]) for row in C])
    for d in np.unique(deg):
        idx = np.nonzero(deg == d)[0]
        if d == 0: out[idx] = True; continue
        c = C[idx, :d+1].astype(np.float64)
        # companion matrices for monic poly z^d + (c_{d-1}/c_d) z^{d-1} + ... ; roots negative reals expected
        M = np.zeros((len(idx), d, d))
        M[:, 1:, :-1] = np.eye(d-1)
        M[:, :, -1] = -(c[:, :d] / c[:, d:d+1])
        ev = np.linalg.eigvals(M)
        rel = np.abs(ev.imag) / np.maximum(np.abs(ev), 1e-300)
        mx[idx] = rel.max(axis=1); out[idx] = mx[idx] < tol
    return out, mx
def counts_at(queries, weight_mask):
    # coefficient vectors c_k(x) = #{n<=x: mask(n), omega(n)=k} at sorted queries
    q = np.asarray(queries); res = np.zeros((len(q), K), dtype=np.int64)
    cum = np.zeros(K, dtype=np.int64); start = 1; CH = 2_000_000
    j = 0
    while j < len(q):
        end = min(start + CH, q[-1] + 1)
        om_c = om[start:end].astype(np.int64); mk = weight_mask[start:end]
        cs = np.zeros((end-start, K), dtype=np.int64)
        for k in range(K): cs[:, k] = np.cumsum(mk & (om_c == k))
        while j < len(q) and q[j] < end:
            res[j] = cum + cs[q[j]-start]; j += 1
        cum = cum + cs[-1]; start = end
    return res
# 1) exhaustive: every squarefree x <= EXH (the polynomial changes only there)
qx = np.nonzero(sqf[:EXH+1])[0]; qx = qx[qx >= 1]
C = counts_at(qx, sqf)
ok, mx = rr_check(C)
bad = qx[~ok]
print('P_x exhaustive x<=%d: %d polys, non-real-rooted: %d' % (EXH, len(qx), len(bad)), 'first:', bad[:12].tolist(), flush=True)
# all-n variant A_x
qa = np.arange(1, min(EXH, 2*10**5)+1)
Ca = counts_at(qa, np.ones(X+1, dtype=bool))
oka, mxa = rr_check(Ca)
print('A_x (all n) x<=%d: non-real-rooted: %d' % (qa[-1], (~oka).sum()), 'first:', qa[~oka][:12].tolist(), flush=True)
# 2) sampled up to X: squarefree points on a log grid plus random
rng = np.random.default_rng(1)
samp = np.unique(np.concatenate([np.round(np.logspace(np.log10(EXH), np.log10(X), 4000)).astype(np.int64),
                                rng.integers(EXH, X, 20000)]))
samp = samp[sqf[samp]]
Cs = counts_at(samp, sqf)
oks, mxs = rr_check(Cs)
print('P_x sampled %d squarefree x in (%d,%d]: non-real-rooted %d' % (len(samp), EXH, X, (~oks).sum()), samp[~oks][:10].tolist(), flush=True)
# 3) zero near -1 vs Mertens: for real-rooted P_x, nearest root r to -1; compare (r+1) with -M(x)/P_x'(-1)
Mcum = np.cumsum(mu.astype(np.int64))
for x in [10**k for k in range(3, int(np.log10(X))+1)] + [X]:
    c = counts_at([x], sqf)[0].astype(float); d = np.max(np.nonzero(c)[0])
    r = np.roots(c[:d+1][::-1]); rr = r[np.argmin(np.abs(r+1))]
    k = np.arange(d+1); Pm1 = np.sum(c[:d+1]*(-1.0)**k); dP = np.sum(c[1:d+1]*k[1:]*(-1.0)**(k[1:]-1))
    print('x=%.0e M(x)=%d P_x(-1)=%d  root near -1: %.12g (+1 = %.3e)  -M/P\'(-1)=%.3e  |r+1|*sqrt(x)=%.3f' % (x, Mcum[x], Pm1, rr.real, rr.real+1, -Pm1/dP, abs(rr.real+1)*np.sqrt(x)))
print('roots at x=X:', np.sort(np.roots(Cs[-1][:np.max(np.nonzero(Cs[-1])[0])+1][::-1]).real))
print('time', round(time.time()-t0,1), 's')
