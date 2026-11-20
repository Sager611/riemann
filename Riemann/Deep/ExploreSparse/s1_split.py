# explore61 s1: integer-dilation penalty split P = R + Q for the head-hybrid optimizer (zeta).
# W1 = {f in L2(0,inf): f = a_m + b/x on I_m = (1/(m+1), 1/m), f = b/x on (1,inf); one b}  contains every rho_{1/k}
# (on I_m, {1/(kx)} = (1/k)(1/x) - floor(m/k)).  Hence V_N in W1 and, for the hybrid optimizer h = sum x_i rho_{1/i} + sum e_j rho_{th_j}:
#   P_hyb = d_N^2 - D_hyb^2 = e^T S e,  S = Gram of P_{V_N}^perp rho_th  (Schur complement)
#   R_hyb = dist(h, W1)^2 = e^T Om e,   Om = Gram of P_{W1}^perp rho_th   (closed form from the jumps of floor(th/x))
#   S - Om = Gram of P_{W1 - V_N} rho_th >= 0 (PSD) for every N.
# Gram: q-design indices k <-> dilation Q/k, <rho_{Q/j}, rho_{Q/k}> = Q G(j,k), G = Vasyunin (code copied from explore31 s3_hybrid.py).
# usage: uv run --no-project --with numpy --with scipy python3 s1_split.py N Q K0 [Mmax]
import sys, os, math, time
import numpy as np
from scipy.linalg import cholesky, solve_triangular, eigh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memcap import check, mem
EG = 0.57721566490153286061
C = 2 + EG - math.log(4 * math.pi)
LG = math.log(2 * math.pi) - EG
GF = '/Users/adriansager/playground/riemann/.scratch/explore61/Gz2000.bin'
NMAXF = 2000
Vc = {}
def V(a, k):  # sum_{m=1}^{k-1} {a m/k} cot(pi m/k) (folded form), a coprime to k
    if k == 1:
        return 0.0
    a %= k
    key = (a, k)
    if key in Vc:
        return Vc[key]
    h = (k - 1) // 2
    m = np.arange(1, h + 1)
    v = float(((2 * ((a * m) % k) - k) @ (1.0 / np.tan(np.pi * m / k))) / k) if h else 0.0
    if k <= 1500:
        Vc[key] = v
    return v
def vas(j, k):
    d = math.gcd(j, k); h, kk = j // d, k // d
    if h > kk:
        h, kk = kk, h
    return (LG * (h + kk) + (kk - h) * math.log(h / kk) - math.pi * (V(h, kk) + V(kk % h, h))) / (2.0 * h * kk) / d

def w1_perp_gram(head, Q, Mmax):
    """Om_ij = <P_{W1}^perp rho_{Q/k_i}, P_{W1}^perp rho_{Q/k_j}> truncated to intervals m <= Mmax, plus B, C'."""
    nh = len(head)
    ms, qs, js, lt = [], [], [], []
    for j, k in enumerate(head):
        nmax = (Mmax * Q) // k + 1
        n = np.arange(1, nmax + 1, dtype=np.int64)
        kn = k * n
        m = kn // Q
        ok = (kn % Q != 0) & (m <= Mmax) & (m >= 1)
        n, m, kn = n[ok], m[ok], kn[ok]
        x = Q / kn.astype(np.float64)
        l = 1.0 / (m + 1); r = 1.0 / m
        q = (x - l) / (r - l)
        ms.append(m); qs.append(q); js.append(np.full(len(m), j)); lt.append(np.log(x / l) - q * np.log(r / l))
    m = np.concatenate(ms); q = np.concatenate(qs); jj = np.concatenate(js); lt = np.concatenate(lt)
    o = np.lexsort((q, m)); m, q, jj, lt = m[o], q[o], jj[o], lt[o]
    B = -np.bincount(jj, weights=lt, minlength=nh)
    Om = np.zeros((nh, nh))
    cuts = np.flatnonzero(np.diff(m)) + 1
    starts = np.concatenate(([0], cuts)); ends = np.concatenate((cuts, [len(m)]))
    for s, t in zip(starts, ends):
        mm = int(m[s]); I = 1.0 / (mm * (mm + 1))
        qq = q[s:t]; ii = jj[s:t]
        Mn = np.minimum.outer(qq, qq) - np.outer(qq, qq)
        np.add.at(Om, (ii[:, None], ii[None, :]), I * Mn)
    mall = np.arange(1, 100_001, dtype=np.float64)  # ||1/x - avg||^2 on I_m = 1/(12 m^2) + O(m^-3): tail 1/(12 M)
    Cp = 1.0 + float(np.sum(1.0 - mall * (mall + 1) * np.log1p(1.0 / mall) ** 2)) + 1.0 / (12 * 100_000)
    return Om, B, Cp

if __name__ == '__main__':
    t0 = time.time()
    N, Q, K0 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    Mmax = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
    L = math.log(N)
    G = np.array(np.memmap(GF, dtype=np.float64, mode='r', shape=(NMAXF, NMAXF))[:N, :N])
    head = [k for k in range(Q + 1, Q * K0) if k % Q]
    nh = len(head)
    th = Q / np.array(head, dtype=np.float64)
    Cb = np.empty((N, nh)); Hb = np.empty((nh, nh))
    for a, k in enumerate(head):
        for b2 in range(a, nh):
            Hb[a, b2] = Hb[b2, a] = Q * vas(k, head[b2])
        for j in range(1, N + 1):
            Cb[j - 1, a] = Q * vas(Q * j, k)
    check(280, 'blocks')
    ks = np.arange(1, N + 1)
    bn = (np.log(ks) + 1 - EG) / ks
    ch = Q * (np.log(np.array(head) / Q) + 1 - EG) / np.array(head)
    Lc = cholesky(G, lower=True, check_finite=False)
    y = solve_triangular(Lc, bn, lower=True, check_finite=False)
    d2 = 1 - float(y @ y)
    Z = solve_triangular(Lc, Cb, lower=True, check_finite=False)
    S = Hb - Z.T @ Z
    rr = ch - Z.T @ y
    e = np.linalg.solve(S, rr)
    P = float(rr @ e)
    Te2 = float(e @ Hb @ e)
    Om, B, Cp = w1_perp_gram(head, Q, Mmax)
    Om2, B2, _ = w1_perp_gram(head, Q, Mmax // 2)
    OmE = 2 * Om - Om2; BE = 2 * B - B2  # Richardson in 1/Mmax (tail ~ 1/M)
    Omf = OmE - np.outer(BE, BE) / Cp
    R = float(e @ Omf @ e)
    ev = eigh(S - Omf, eigvals_only=True)
    pw = float(np.sum(np.abs(e) * np.sqrt(np.diag(S))) ** 2)
    print(f'N={N} Q={Q} K0={K0} nh={nh}: d2L={d2*L:.6f} L^2P={L*L*P:.5f} L^2R={L*L*R:.5f} R/P={R/P:.4f} '
          f'L^2|Te|^2={L*L*Te2:.5f} L^2 powerbound={L*L*pw:.4f} | eig(S-Om) min {ev.min():+.2e} max {ev.max():.2e} '
          f'| diag S/Om: min {np.min(np.diag(S)/np.diag(Omf)):.4f} med {np.median(np.diag(S)/np.diag(Omf)):.4f} '
          f'| Mmax={Mmax} Om-Richardson {np.max(np.abs(OmE-Om))/np.max(np.abs(Om)):.1e} | t={time.time()-t0:.0f}s mem {mem()}', flush=True)
    np.savez(f'/Users/adriansager/playground/riemann/Riemann/Deep/ExploreSparse/logs/s1_N{N}_Q{Q}_K{K0}.npz',
             head=np.array(head), th=th, e=e, S=S, Om=Omf, B=BE, Cp=Cp, P=P, R=R, d2=d2, Te2=Te2)
