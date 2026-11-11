# Barrett-Jarvis conjecture (LAA 162 (1992); restated in arXiv 2511.13627 l. 524): the "small" eigenvalues of the
# Redheffer matrix R_n (all but the eigenvalue 1 and the two of size ~ +-sqrt n) have negative real part.
# det(xI - R_n) = y^{n-s-1} (y^{s+1} - sum_{k=1}^s F_k(n) y^{s-k}),  y = x - 1, s = floor(log2 n),
# F_k(n) = #{ordered factorizations m_1...m_k <= n, m_i >= 2}.  Exhaustive in n <= N.
import sys, numpy as np
N = int(sys.argv[1])
s_max = N.bit_length() - 1
f = np.zeros((s_max + 1, N + 1), dtype=np.float64)   # f[k][j] = # ordered factorizations of j into k factors >= 2
f[0][1] = 1
for k in range(1, s_max + 1):
    prev = f[k-1]; cur = f[k]
    nz = np.nonzero(prev)[0]
    for d in nz:
        cur[2*d::d] += prev[d]
F = np.cumsum(f, axis=1)       # F[k][n] = sum_{j<=n} f_k(j)
worst = (-1e9, None); recs = []
for n in range(4, N + 1):
    s = n.bit_length() - 1
    coef = [1.0, 0.0] + [-F[k][n] for k in range(1, s + 1)]   # y^{s+1} + 0 y^s - F_1 y^{s-1} - ... - F_s
    r = np.roots(coef)
    r = r[np.argsort(-np.abs(r))]
    small = r[2:]
    if len(small) == 0: continue
    mx = float(np.max(small.real + 1.0))
    if mx > worst[0]: worst = (mx, n); recs.append((n, mx))
print('max over n<=%d of max Re(small eigenvalue) = %.6f at n = %d' % (N, worst[0], worst[1]))
print('record progression (n, max Re):', [(a, round(b, 4)) for a, b in recs[-12:]])
