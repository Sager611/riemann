# chaudit2: direct test of Lemma edesmo (l.4328) on exact Q(n), n <= 1e7:
# sup_{0<=t1<t2<=1e7} Q(t2) - Q(t1) - c1 (t2 - t1) = c1 + max_{j<m} (A(m) - A(j)),  A(n) = Q(n) - c1 n.
import numpy as np, sys
from fractions import Fraction as Fr
N = 10**7
mu = np.ones(N + 1, dtype=np.int8); mu[0] = 0
isp = np.ones(N + 1, dtype=bool); isp[:2] = False
for p in range(2, int(N**0.5) + 1):
    if isp[p]: isp[p*p::p] = False
for p in np.nonzero(isp)[0]:
    p = int(p); mu[p::p] *= -1
    if p*p <= N: mu[p*p::p*p] = 0
del isp
Q = np.cumsum(mu != 0, dtype=np.int32); del mu
pairs = [(Fr(1), Fr(1)), (Fr(3, 4), Fr(3, 2)), (Fr(2, 3), Fr(8, 3)), (Fr(16, 25), Fr(114, 25)), (Fr(768, 1225), Fr(9458, 1225)),
         (Fr(18432, 29645), Fr(361192, 29645)), (Fr(442368, 715715), Fr(14328304, 715715))]
for c1, c2 in pairs:
    c1f = float(c1); best = -1e18; runmin = 0.0   # A(0) = 0
    for lo in range(1, N + 1, 10**6):
        n = np.arange(lo, min(N, lo + 10**6 - 1) + 1, dtype=np.float64)
        A = Q[lo:lo + n.size].astype(np.float64) - c1f*n
        rm = np.minimum.accumulate(np.concatenate(([runmin], A)))[:-1]   # min of A over [0, m-1]
        best = max(best, float(np.max(A - rm)))
        runmin = min(runmin, float(A.min()))
    print("pair (%s, %s): sup Q(t2)-Q(t1)-c1(t2-t1) on [0,1e7] = %.6f   c2 = %.6f   ok: %s" % (c1, c2, best + c1f, float(c2), best + c1f <= float(c2) + 1e-9)); sys.stdout.flush()
