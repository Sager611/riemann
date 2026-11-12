# explore37 u6: rigorous (Arb, python-flint) floors eta_M(S)^2 for the thinned Baez-Duarte criterion.
# Same construction as u1_exact.py: cylinders I_j = (1/(j+1),1/j], local Gram of {1, log x} (exact moments), upper Cholesky R_j,
# rows R_j (alpha, beta); f(kx) = n(1+log k+log x)+log n!, n = floor(j/k); e = 1 + log x. eta^2 = y.y - r^T G^{-1} r, G = A^T A.
# Every arithmetic step is a ball operation; the printed interval contains the exact floor. usage: python3 u6_arb.py
import sys, time, resource
from flint import arb, arb_mat, ctx
ctx.prec = 256
t0 = time.time()
def chol2(j):
    a = arb(1) / (j + 1); b = arb(1) / j
    F1 = lambda x: x * x.log() - x
    F2 = lambda x: x * x.log() ** 2 - 2 * x * x.log() + 2 * x
    m0 = b - a; m1 = F1(b) - F1(a); m2 = F2(b) - F2(a)  # moments of 1, log x, log^2 x on I_j
    r00 = m0.sqrt(); r01 = m1 / r00; r11 = (m2 - r01 * r01).sqrt()
    return r00, r01, r11
def floor_eta2(S, M):
    ks = [k for k in S if k <= M]; n = len(ks)
    A = arb_mat(2 * M, n); y = arb_mat(2 * M, 1)
    lgam = [arb(0)] + [None] * (M + 1)
    acc = arb(0)
    for i in range(1, M + 2):
        acc = acc + arb(i).log(); lgam[i] = acc  # log i!
    for j in range(1, M + 1):
        r00, r01, r11 = chol2(j); lx = None
        # e = 1 + log x: alpha = 1, beta = 1 in basis {1, log x}
        y[2 * (j - 1), 0] = r00 * 1 + r01 * 1; y[2 * (j - 1) + 1, 0] = r11 * 1
        for c, k in enumerate(ks):
            nn = j // k
            if nn == 0: continue
            alpha = nn * (1 + arb(k).log()) + lgam[nn]; beta = arb(nn)
            A[2 * (j - 1), c] = r00 * alpha + r01 * beta; A[2 * (j - 1) + 1, c] = r11 * beta
    At = A.transpose(); G = At * A; r = At * y; yy = (y.transpose() * y)[0, 0]
    x = G.solve(r)
    return yy - (r.transpose() * x)[0, 0]
primes = [p for p in range(2, 2001) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
jobs = [('{1}+primes', [1] + primes, 1000)]
for nm, S, M in jobs:
    t1 = time.time(); v = floor_eta2(S, M)
    print(f'{nm:12s} M={M:4d}  eta^2 in {v.str(12, radius=True)}   ({time.time()-t1:.1f}s, rss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20:.0f} MB)', flush=True)
print(f'done {time.time()-t0:.1f}s')
