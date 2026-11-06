# explore13: Baez-Duarte distances d_N for zeta (Vasyunin Gram) and for the Davenport-Heilbronn function f
# (exact digamma Gram), N up to Nmax.  d_N^2 = inf_c || chi - sum_k c_k F(1/(k t)) ||^2_{L^2(0,inf)},
#   zeta: F(x) = {x} (Mellin -zeta(s)/s);  DH: F(x) = A(x) = sum_{n<=x} a_n, a = (1,k,-k,-1,0) mod 5 (Mellin f(s)/s).
# Usage: uv run --no-project --with mpmath --with numpy --with scipy python3 nbgram.py Nmax_zeta Nmax_dh
#   (reported run: 1280 400, about 90 s single-threaded; results in results.txt)
import sys, math
import numpy as np
from math import gcd
from scipy.special import digamma
import mpmath as mp

EG = 0.5772156649015329
NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 200
ND = int(sys.argv[2]) if len(sys.argv) > 2 else 200

# ---------- zeta: Vasyunin's formula ----------
def V(h, k):  # sum_{m=1}^{k-1} {m h / k} cot(pi m / k)
    if k == 1:
        return 0.0
    m = np.arange(1, k)
    return float(np.sum(((m * h) % k) / k / np.tan(np.pi * m / k)))

def gz_coprime(h, k):
    return ((math.log(2 * math.pi) - EG) * (h + k) / (2 * h * k) + (k - h) / (2 * h * k) * math.log(h / k)
            - math.pi / (2 * h * k) * (V(h, k) + V(k, h)))

def gram_zeta(N):
    G = np.zeros((N, N))
    cache = {}
    for h in range(1, N + 1):
        for k in range(h, N + 1):
            g = gcd(h, k); a, b = h // g, k // g
            if (a, b) not in cache:
                cache[(a, b)] = gz_coprime(a, b)
            G[h - 1, k - 1] = G[k - 1, h - 1] = cache[(a, b)] / g
    b = np.array([(math.log(k) + 1 - EG) / k for k in range(1, N + 1)])
    return G, b

# direct check: int_0^inf F(x/h) F(x/k) dx / x^2, exact on pieces up to R periods, tail by moment expansion
def direct_zeta(h, k, R=400):
    P = h * k // gcd(h, k)
    br = sorted(set([m * h for m in range(R * P // h + 1)] + [m * k for m in range(R * P // k + 1)]))
    tot = mp.mpf(0)
    for a, bb in zip(br[:-1], br[1:]):
        p, q = a // h, a // k
        a_, b_ = mp.mpf(a), mp.mpf(bb)
        if a == 0:  # x^2/(hk)/x^2 on [0, b)
            tot += b_ / (h * k); continue
        tot += (b_ - a_) / (h * k) - (mp.mpf(q) / h + mp.mpf(p) / k) * mp.log(b_ / a_) + p * q * (1 / a_ - 1 / b_)
    # tail: periodic p(y) on [0,P): moments m0, m1, m2
    m0 = m1 = m2 = mp.mpf(0)
    brp = sorted(set([m * h for m in range(P // h + 1)] + [m * k for m in range(P // k + 1)]))
    for a, bb in zip(brp[:-1], brp[1:]):
        p, q = a // h, a // k
        # p(y) = (y/h - p)(y/k - q) = y^2/(hk) - y (q/h + p/k) + p q
        c2, c1, c0 = mp.mpf(1) / (h * k), -(mp.mpf(q) / h + mp.mpf(p) / k), mp.mpf(p * q)
        def mom(n):
            return sum(cc * (mp.mpf(bb) ** (n + j + 1) - mp.mpf(a) ** (n + j + 1)) / (n + j + 1) for j, cc in enumerate((c0, c1, c2)))
        m0 += mom(0); m1 += mom(1); m2 += mom(2)
    Pm = mp.mpf(P)
    tail = m0 / Pm**2 * mp.zeta(2, R) - 2 * m1 / Pm**3 * mp.zeta(3, R) + 3 * m2 / Pm**4 * mp.zeta(4, R)
    return tot + tail

# ---------- DH ----------
kap = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
S = [0.0, 1.0, 1.0 + kap, 1.0, 0.0]  # A(y) = S[floor(y) mod 5]

def gdh_coprime(h, k):
    P = 5 * h * k  # period of A(x/h) A(x/k) (coprime h, k)
    br = np.union1d(np.arange(0, P + 1, h), np.arange(0, P + 1, k)).astype(np.float64)
    a, b = br[:-1], br[1:]
    c = np.array(S)[(a // h).astype(np.int64) % 5] * np.array(S)[(a // k).astype(np.int64) % 5]
    m = c != 0
    return float(np.sum(c[m] / P * (digamma(b[m] / P) - digamma(a[m] / P))))

def gram_dh(N):
    G = np.zeros((N, N)); cache = {}
    for h in range(1, N + 1):
        for k in range(h, N + 1):
            g = gcd(h, k); a, b = h // g, k // g
            if (a, b) not in cache:
                cache[(a, b)] = gdh_coprime(a, b)
            G[h - 1, k - 1] = G[k - 1, h - 1] = cache[(a, b)] / g
    beta = sum(S[r] * (digamma((r + 1) / 5) - digamma(r / 5)) / 5 for r in range(1, 5))
    b = np.array([beta / k for k in range(1, N + 1)])
    return G, b

def direct_dh(h, k, R=400):
    P = 5 * h * k // gcd(h, k)
    br = sorted(set(list(range(0, R * P + 1, h)) + list(range(0, R * P + 1, k))))
    tot = mp.mpf(0)
    for a, bb in zip(br[:-1], br[1:]):
        c = S[(a // h) % 5] * S[(a // k) % 5]
        if c:
            tot += c * (mp.mpf(1) / a - mp.mpf(1) / bb)
    # tail: periodic step function, moments
    m0 = m1 = m2 = mp.mpf(0)
    brp = sorted(set(list(range(0, P + 1, h)) + list(range(0, P + 1, k))))
    for a, bb in zip(brp[:-1], brp[1:]):
        c = S[(a // h) % 5] * S[(a // k) % 5]
        m0 += c * (bb - a); m1 += c * (mp.mpf(bb)**2 - a**2) / 2; m2 += c * (mp.mpf(bb)**3 - a**3) / 3
    Pm = mp.mpf(P)
    return tot + m0 / Pm**2 * mp.zeta(2, R) - 2 * m1 / Pm**3 * mp.zeta(3, R) + 3 * m2 / Pm**4 * mp.zeta(4, R)

def dN2(G, b, N, prec=False):
    if not prec:
        L = np.linalg.cholesky(G[:N, :N])
        y = np.linalg.solve(L, b[:N])
        return 1.0 - float(y @ y)
    mp.mp.dps = 40
    Gm = mp.matrix(G[:N, :N].tolist()); bm = mp.matrix(b[:N].tolist())
    x = mp.lu_solve(Gm, bm)
    return float(1 - (bm.T * x)[0])

if __name__ == "__main__":
    mp.mp.dps = 30
    print("check zeta Gram (Vasyunin vs direct):")
    for (h, k) in [(1, 1), (1, 2), (2, 3), (3, 5), (4, 7)]:
        print(f"  ({h},{k}) vasyunin={gz_coprime(h,k):.12f} direct={float(direct_zeta(h,k,R=60)):.12f}")
    print("check DH Gram (digamma vs direct):")
    for (h, k) in [(1, 1), (1, 2), (2, 3), (3, 5)]:
        print(f"  ({h},{k}) digamma={gdh_coprime(h,k):.12f} direct={float(direct_dh(h,k,R=60)):.12f}")
    Gz, bz = gram_zeta(NZ)
    print("zeta: N, d_N^2, d_N^2 log N   (C = 2 + gamma - log 4pi = 0.0461914179)")
    for N in [1, 2, 5, 10, 20, 30, 40, 50, 80, 100, 160, 200, 320, 640, 1280]:
        if N <= NZ:
            v = dN2(Gz, bz, N)
            extra = f"  [mp40: {dN2(Gz, bz, N, prec=True):.10f}]" if N in (10, 30) else ""
            print(f"  {N:5d} {v:.8f} {v*math.log(N):.5f}{extra}", flush=True)
    np.save("/Users/adriansager/playground/riemann/.scratch/explore13/Gz.npy", Gz)
    Gd, bd = gram_dh(ND)
    print("DH: N, d_N^2, d_N^2 log N")
    for N in [2, 5, 10, 20, 40, 80, 100, 160, 200, 300, 400]:
        if N <= ND:
            v = dN2(Gd, bd, N)
            extra = f"  [mp40: {dN2(Gd, bd, N, prec=True):.10f}]" if N in (10, 30) else ""
            print(f"  {N:5d} {v:.8f} {v*math.log(N):.5f}{extra}", flush=True)
    np.save("/Users/adriansager/playground/riemann/.scratch/explore13/Gd.npy", Gd)
