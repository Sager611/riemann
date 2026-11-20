# explore61 s8: adversarial Davenport-Heilbronn version of s3_rate + s5_zero.
# DH: A(y) = S[floor(y) mod 5], S = (0, 1, 1+kappa, 1, 0); dilation functions A(theta/x); Gram <A(Q/(jx)), A(Q/(kx))> = Q gd(j,k)
# (merged-breakpoint digamma formula copied from explore38 s10_farey.py; natural Gram from ceilaudit6 gram.c 'per' mode).
# Every A(1/(kx)) is constant on each I_m and 0 on (1,inf): V_N^DH in PC (piecewise constants on the harmonic partition), chi in PC.
# omega^2 = dist(A_theta, PC)^2 from the jumps a_n (x = theta/n); kappa_DH(N) = L (pi_N^2 - omega^2) = L dist(P_PC A_theta, V_N)^2.
# Zero formula prediction: sum_{online} |G(rho)|^2  +  L * sum_{offline, beta>1/2} (2 beta - 1) |G(rho)|^2 (H^2 evaluation floor), G = M[P_perp A_theta].
# usage: uv run --no-project --with numpy --with scipy python3 s8_dh.py NMAX MMAX
import sys, os, math, time, json
import numpy as np
from math import gcd
from scipy.special import digamma
from scipy.linalg import cholesky, solve_triangular
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memcap import check, mem
t0 = time.time()
NMAX, MMAX = int(sys.argv[1]), int(sys.argv[2])
kap = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
S = [0.0, 1.0, 1.0 + kap, 1.0, 0.0]; QP = 5
an = [S[n % 5] - S[(n - 1) % 5] for n in range(5)]  # jump of A at y = n (n mod 5)
cache = {}
def gd(h, k):
    g = gcd(h, k); a_, b_ = h // g, k // g
    if a_ > b_: a_, b_ = b_, a_
    if (a_, b_) not in cache:
        P = QP * a_ * b_
        br = np.union1d(np.arange(0, P + 1, a_), np.arange(0, P + 1, b_)).astype(np.float64)
        a, b = br[:-1], br[1:]
        c = np.array(S)[(a // a_).astype(np.int64) % QP] * np.array(S)[(a // b_).astype(np.int64) % QP]
        mm = c != 0
        v = float(np.sum(c[mm] / P * (digamma(b[mm] / P) - digamma(a[mm] / P))))
        if a_ * b_ < 20000: cache[(a_, b_)] = v
        return v / g
    return cache[(a_, b_)] / g
betaDH = sum(S[r] * (digamma((r + 1) / QP) - digamma(r / QP)) / QP for r in range(1, QP))
G = np.array(np.memmap('/Users/adriansager/playground/riemann/.scratch/explore61/Gdh1000.bin', dtype=np.float64, mode='r', shape=(1000, 1000))[:NMAX, :NMAX])
print('Gram check G11 =', G[0, 0], ' gd(1,1) =', gd(1, 1), ' G23 =', G[1, 2], ' gd(2,3) =', gd(2, 3), flush=True)
Lc = cholesky(G, lower=True, overwrite_a=True, check_finite=False); del G
ks = np.arange(1, NMAX + 1)
y = solve_triangular(Lc, betaDH / ks, lower=True, check_finite=False)
d2 = 1 - np.cumsum(y * y)
dz = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/dh_zeros.json'))
gon = np.array(dz['online']); off = np.array(dz['offline']); TMAX = dz['TMAX']
son = 0.5 + 1j * gon; soff = off[:, 0] + 1j * off[:, 1]
Ns = [n for n in [10, 40, 160, 640, 1000] if n <= NMAX]
print(f'DH: {len(gon)} online zeros < {TMAX}, offline (beta, gamma): ' + ' '.join(f'({b:.3f},{g:.1f})' for b, g in off), flush=True)
print(f'chi: L d_N^2 at N={Ns}: ' + ' '.join(f'{d2[n-1]*math.log(n):.4f}' for n in Ns) + f'   sum_on 1/|rho|^2 = {2*np.sum(1/np.abs(son)**2):.4f} (<{TMAX})', flush=True)
def Gmel(Q, k, s):
    n = np.arange(1, (MMAX * Q) // k + 2, dtype=np.int64); kn = k * n; m = kn // Q
    ok = (kn % Q != 0) & (m <= MMAX) & (n % 5 != 0); n, m, kn = n[ok], m[ok].astype(np.float64), kn[ok]
    x = Q / kn.astype(np.float64); l, r = 1 / (m + 1), 1 / m; q = (x - l) / (r - l); a = np.array(an)[n % 5]
    om2 = float(np.sum(a * a * (q - q * q) / (m * (m + 1))))
    out = np.zeros(len(s), complex); lx, ll, lr = np.log(x), np.log(l), np.log(r)
    for i0 in range(0, len(s), 16):
        ss = s[i0:i0 + 16, None]
        out[i0:i0 + 16] = np.sum(a * ((np.exp(ss * lx) - np.exp(ss * ll)) - q * (np.exp(ss * lr) - np.exp(ss * ll))), axis=1) / s[i0:i0 + 16]
    return out, om2
TH = [(2, 3), (3, 4), (3, 5), (4, 5), (4, 7), (5, 8), (2, 5), (3, 7), (3, 8), (2, 7), (2, 9)]
for Q, k in TH:
    cv = np.array([Q * gd(Q * j, k) for j in range(1, NMAX + 1)])
    z = solve_triangular(Lc, cv, lower=True, check_finite=False)
    pi2 = Q * gd(k, k) - np.cumsum(z * z)
    Gon, om2 = Gmel(Q, k, son); Gof, om2b = Gmel(Q, k, soff)
    Gon2, om2h = Gmel(Q, k, son[:1]) if False else (None, None)
    con = 2 * np.sum(np.abs(Gon) ** 2)
    A = np.mean(np.abs(Gon[len(gon) // 2:]) ** 2 * gon[len(gon) // 2:] ** 2)
    tail = 2 * A * (math.log(TMAX / (2 * math.pi)) + 1) / (2 * math.pi * TMAX)
    coff = 2 * np.sum((2 * off[:, 0] - 1) * np.abs(Gof) ** 2)
    kap_m = [math.log(n) * (pi2[n - 1] - om2) for n in Ns]
    print(f'{Q}/{k}: om^2={om2:.5f} |A|^2={Q*gd(k,k):.4f}  L(pi^2-om^2) at N={Ns}: ' + ' '.join(f'{v:.5f}' for v in kap_m) +
          f' | pred online {con:.5f}+tail {tail:.5f} = {con+tail:.5f}; offline slope {coff:.2e} (x L: {coff*math.log(Ns[-1]):.5f} at N={Ns[-1]}) | t={time.time()-t0:.0f}s', flush=True)
    check(280, f'{Q}/{k}')
print(f'done {time.time()-t0:.0f}s mem {mem()}')
