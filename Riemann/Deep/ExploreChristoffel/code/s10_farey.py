# explore38 s10: Farey relaxation (all dilations a/b in [1/N,1], b <= Bd) = upper estimate of the continuous NB distance D(1/N)^2,
# for zeta (Vasyunin, double + mpmath check) and DH (digamma). Gram <A(al/t),A(be/t)> = a c g(b c, a d) (homogeneity, as explore13's farey.py).
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 s10_farey.py which N Bd [dps]
import sys, math, time
import numpy as np
from math import gcd
import mpmath as mp
from scipy.special import digamma
t0 = time.time()
which = sys.argv[1]; N = int(sys.argv[2]); Bd = int(sys.argv[3]); dps = int(sys.argv[4]) if len(sys.argv) > 4 else 0
EG = 0.57721566490153286061
cache = {}
def Vsum(h, k):  # sum_{m=1}^{k-1} {m h/k} cot(pi m/k), double
    if k == 1: return 0.0
    m = np.arange(1, k)
    return float(np.sum(((m * h) % k) / k / np.tan(np.pi * m / k)))
def Vsum_mp(h, k):
    return mp.fsum(mp.mpf((m * h) % k) / k * mp.cot(mp.pi * m / k) for m in range(1, k))
def gz(h, k):
    g = gcd(h, k); h //= g; k //= g
    if h > k: h, k = k, h
    if (h, k) not in cache:
        if dps:
            cache[(h, k)] = ((mp.log(2 * mp.pi) - mp.euler) * (h + k) / (2 * h * k) + mp.mpf(k - h) / (2 * h * k) * mp.log(mp.mpf(h) / k) - mp.pi / (2 * h * k) * (Vsum_mp(h, k) + Vsum_mp(k, h)))
        else:
            cache[(h, k)] = ((math.log(2 * math.pi) - EG) * (h + k) / (2 * h * k) + (k - h) / (2 * h * k) * math.log(h / k) - math.pi / (2 * h * k) * (Vsum(h, k) + Vsum(k, h)))
    return cache[(h, k)] / g
kap = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
S = [0.0, 1.0, 1.0 + kap, 1.0, 0.0] if which != 'chi4' else [0.0, 1.0, 1.0, 0.0]
QP = len(S)
def gd(h, k):
    g = gcd(h, k); a_, b_ = h // g, k // g
    if a_ > b_: a_, b_ = b_, a_
    if (a_, b_) not in cache:
        P = QP * a_ * b_
        br = np.union1d(np.arange(0, P + 1, a_), np.arange(0, P + 1, b_)).astype(np.float64)
        a, b = br[:-1], br[1:]
        c = np.array(S)[(a // a_).astype(np.int64) % QP] * np.array(S)[(a // b_).astype(np.int64) % QP]
        mm = c != 0
        cache[(a_, b_)] = float(np.sum(c[mm] / P * (digamma(b[mm] / P) - digamma(a[mm] / P))))
    return cache[(a_, b_)] / g
betaDH = sum(S[r] * (digamma((r + 1) / QP) - digamma(r / QP)) / QP for r in range(1, QP))
def run(dic):
    n = len(dic)
    if which == 'zeta':
        g = gz; tgt = lambda th: th * (math.log(1 / th) + 1 - EG)
    else:
        g = gd; tgt = lambda th: th * betaDH
    if dps:
        mp.mp.dps = dps
        G = mp.matrix(n, n); bb = mp.matrix(n, 1)
        for i, (a, b) in enumerate(dic):
            th = mp.mpf(a) / b
            bb[i] = th * (mp.log(1 / th) + 1 - mp.euler)
            for j in range(i, n):
                c, d = dic[j]
                G[i, j] = G[j, i] = a * c * g(b * c, a * d)
        Lc = mp.cholesky(G); y = mp.lu_solve(Lc, bb)
        return float(1 - mp.fsum(y[i] ** 2 for i in range(n)))
    G = np.zeros((n, n)); bb = np.zeros(n)
    for i, (a, b) in enumerate(dic):
        bb[i] = tgt(a / b)
        for j in range(i, n):
            c, d = dic[j]
            G[i, j] = G[j, i] = a * c * g(b * c, a * d)
    w, U = np.linalg.eigh(G)
    keep = w > w.max() * 1e-15
    y = U[:, keep].T @ bb
    return float(1 - np.sum(y * y / w[keep])), float(w.min() / w.max())
far = sorted({(a // gcd(a, b), b // gcd(a, b)) for b in range(1, Bd + 1) for a in range(1, b + 1) if a * N >= b} | {(1, k) for k in range(1, N + 1)}, key=lambda t: t[0] / t[1])  # Farey(Bd) U {1/k}: D(1/N)^2 <= value <= d_N^2
bd = [(1, k) for k in range(1, N + 1)]
r1 = run(bd); r2 = run(far)
L = math.log(N)
f = lambda r: (r if isinstance(r, float) else r[0])
print(f'{which} N={N} Bd={Bd} dps={dps}: BD n={len(bd)} d2 logN={f(r1)*L:.6f} | Farey n={len(far)} d2 logN={f(r2)*L:.6f} (cond^-1 {r2[1] if not isinstance(r2, float) else 0:.1e}) | {time.time()-t0:.1f}s', flush=True)
