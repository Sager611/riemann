
# t2_clusters.py -- are the configurations that saturate Guth-Maynard's energy bound (Prop 11.1) realizable as large-value sets?
# W = m blocks of K consecutive integers at random positions in [0,T] ("non-sticky" at the coarse scale).
# GM-extremal family: K*m^{1/2} = N^{2-2s} (then E(W) ~ |W|^2 K = |W| N^{4-4s}); Montgomery family: m*K = N^{2-2s}.
# lam1 = top eigenvalue of Gram G(t,t') = sum_{N<=n<2N} n^{i(t-t')}. For |b_n|<=1: sum_W |D|^2 <= lam1*N, so
# bestL2 = lam1*N/|W| bounds the mean of |D|^2 on W (if < V^2 = N^{2s}, W is NOT a large-value set: rigorous).
# unimod: b_n = phase of (M_W^* v1) (|b_n| = 1), D evaluated on W. HM = Schur/Halasz-Montgomery entrywise bound.
# trace-loss = (tr G^3)^{1/3}/lam1: what the cubic trace (GM Lemma 4.2, r=3) loses on this W.
# usage: python3 t2_clusters.py N tau
import numpy as np, sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
tau = float(sys.argv[2]) if len(sys.argv) > 2 else 1.2
s = 0.75
T = int(round(N**tau))
n = np.arange(N, 2*N, dtype=float); ln = np.log(n)
g = np.empty(T+1, dtype=complex)
for a in range(0, T+1, 512):
    u = np.arange(a, min(a+512, T+1), dtype=float)
    g[a:a+len(u)] = np.exp(1j*np.outer(u, ln)).sum(axis=1)
def gram(W):
    D = W[:, None] - W[None, :]
    G = g[np.abs(D)]
    return np.where(D >= 0, G, np.conj(G))
def place(m, K, rng):
    while True:
        starts = np.sort(rng.integers(0, T-K, size=m))
        if m == 1 or np.min(np.diff(starts)) >= 2*K: break
    return np.concatenate([np.arange(a, a+K) for a in starts]).astype(np.int64)
V2 = N**(2*s); base = N**(2-2*s)
print('N=%d tau=%.2f T=%d T/(2 pi N)=%.1f s=%.2f N^{2-2s}=%.1f threshold N^s=%.0f' % (N, tau, T, T/(2*np.pi*N), s, base, N**s))
ag = np.abs(g)
for lo in (10, 50, N, int(2*np.pi*N)):
    if lo < T: print('max|g(u)|/N^{1/2} over u in [%d,T] = %.2f   (|g(0)|/N^{1/2} = %.1f)' % (lo, ag[lo:].max()/np.sqrt(N), ag[0]/np.sqrt(N)))
print('family   m    K  |W|  lam1/N  HM/N  bestL2/V^2  unimod: min|D|/V  med|D|/V  frac>=V  trace-loss')
rng = np.random.default_rng(1)
for fam, mlist in (('GMext', [1, 4, 16, 64]), ('Mont', [1, 4, 16])):
    for m in mlist:
        K = max(1, int(round(base/np.sqrt(m)))) if fam == 'GMext' else max(1, int(round(base/m)))
        for rep in range(2 if N <= 2000 else 1):
            W = place(m, K, rng)
            G = gram(W)
            w, v = np.linalg.eigh(G)
            lam1 = w[-1]; v1 = v[:, -1]
            b = np.exp(-1j*np.outer(ln, W)) @ v1
            b = b/np.abs(b)
            a = np.abs(np.exp(1j*np.outer(W, ln)) @ b)
            tr3 = np.sum(w**3)
            blk = np.repeat(np.arange(m), K)
            off = np.abs(G)*(blk[:, None] != blk[None, :])
            K0 = np.linalg.eigvalsh(G[:K, :K])[-1]
            HM = K0 + off.sum(axis=1).max()
            print('%-6s %4d %4d %4d  %7.2f %6.1f  %8.3f     %8.3f  %8.3f  %6.2f   %6.2f' % (
                fam, m, K, len(W), lam1/N, HM/N, lam1*N/len(W)/V2, a.min()/N**s, np.median(a)/N**s,
                np.mean(a >= N**s), tr3**(1/3)/lam1), flush=True)
