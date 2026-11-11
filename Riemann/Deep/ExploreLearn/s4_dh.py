# explore31 s4 (adversarial): the same designs for the Davenport-Heilbronn function f (functional equation, no Euler product,
# off-line zeros).  Features A(Q/(k t)), A(x) = sum_{n<=x} a_n, a = (1, kap, -kap, -1, 0) mod 5 (Mellin f(s)/s); exact digamma Gram
# copied from Riemann/Deep/ExploreHankel/nbgram.py (explore13).  Gram Q G^f_{jk}; target <chi, A(theta/t)> = theta*beta (theta <= 1).
# Compares natural (Q j), q-design {Q..QN} (-> continuous relaxation) and random+pinned designs with the LC formula.
# usage: uv run --no-project --with numpy --with scipy python3 s4_dh.py
import sys, math, time
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh, cho_factor, cho_solve
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from memcap import check, mem
kap = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
S = np.array([0.0, 1.0, 1.0 + kap, 1.0, 0.0])
CF = 0.19527  # Sum 1/(rho(1-rho)) over DH zeros (explore13 dhconst.py: 0.1952...)
cache = {}
def gdh(j, k):
    d = math.gcd(j, k); h, kk = j // d, k // d
    if h > kk:
        h, kk = kk, h
    if (h, kk) not in cache:
        P = 5 * h * kk
        br = np.union1d(np.arange(0, P + 1, h), np.arange(0, P + 1, kk)).astype(np.float64)
        a, b = br[:-1], br[1:]
        c = S[(a // h).astype(np.int64) % 5] * S[(a // kk).astype(np.int64) % 5]
        m = c != 0
        cache[(h, kk)] = float(np.sum(c[m] / P * (digamma(b[m] / P) - digamma(a[m] / P))))
    return cache[(h, kk)] / d
beta = float(sum(S[r] * (digamma((r + 1) / 5) - digamma(r / 5)) / 5 for r in range(1, 5)))
def gram(idx):
    n = len(idx); Gs = np.empty((n, n))
    for a in range(n):
        for b in range(a, n):
            Gs[a, b] = Gs[b, a] = gdh(int(idx[a]), int(idx[b]))
    return Gs
def dist2(Gs, t, Q):
    cf = cho_factor(Gs, lower=True, check_finite=False)
    return 1 - Q * float(t @ cho_solve(cf, t, check_finite=False))
t0 = time.time()
rng = np.random.default_rng(31)
for N, Qs in [(2, [1, 2, 4, 8, 16, 32, 64, 128]), (5, [1, 2, 4, 8, 16, 32, 64]), (10, [1, 2, 4, 8, 16, 32]), (20, [1, 2, 4, 8, 16])]:
    L = math.log(N)
    for Q in Qs:
        idx = np.arange(Q, Q * N + 1)
        Gs = gram(idx); t = beta / idx
        v = dist2(Gs, t, Q)
        line = f'DH N={N:3d} Q={Q:4d} M={len(idx):4d} d2={v:.8f} d2L={v*L:.6f} L(d2L-Cf)={L*(v*L-CF):+.5f}'
        if Q == Qs[-1] and N >= 5:  # random + pinned theta=1, and LC on the pool
            rs = []
            for r in range(100):
                sel = np.concatenate(([0], 1 + np.sort(rng.choice(len(idx) - 1, size=N - 1, replace=False))))
                rs.append(dist2(Gs[np.ix_(sel, sel)], t[sel], Q))
            Gam = Q * Gs; bet = Q * t
            g00 = Gam[0, 0]; g10 = Gam[1:, 0]
            Gp = Gam[1:, 1:] - np.outer(g10, g10) / g00; bp = bet[1:] - g10 * bet[0] / g00; y2 = 1 - bet[0] ** 2 / g00
            lam, U = eigh(Gp); c = (U.T @ bp) ** 2; lam = np.maximum(lam, 1e-300)
            lo, hi = -60.0, 10.0
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if float(np.sum(lam / (lam + math.exp(mid)))) > N - 1 else (lo, mid)
            ELC = y2 - float(np.sum(c / (lam + math.exp(0.5 * (lo + hi)))))
            line += f' | random+rho1 {np.mean(rs):.6f} (sd {np.std(rs):.6f}) E_LC {ELC:.6f}'
        print(line + f'  t={time.time()-t0:.0f}s', flush=True)
        check(250, f'N={N} Q={Q}')
print('footprint/resident MB', mem())
