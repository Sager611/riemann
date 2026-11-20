# explore61 s13: s1_split with continuous freedom beyond the head: extra dilations = head grid Q/k in (1/K0, 1) plus tail midpoints
# 2/(2k+1), K0 <= k <= K1 (one non-integer dilation in every gap (1/(k+1), 1/k) down to 1/(K1+1)).  General rationals a/b:
# <rho_{a/b}, rho_{c/d}> = a c G(b c, a d) (G = Vasyunin), jumps of floor((a/b)/x) at x = a/(b n).   usage: s13_tail.py N Q K0 K1 [Mmax]
import sys, os, math, time
import numpy as np
from math import gcd
from scipy.linalg import cholesky, solve_triangular, eigh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memcap import check, mem
from s1_split import vas, GF, NMAXF, EG
def perp_gram(rats, Mmax):
    nh = len(rats); ms, qs, js, lt = [], [], [], []
    for j, (a, b) in enumerate(rats):
        n = np.arange(1, (Mmax * a) // b + 2, dtype=np.int64); bn = b * n; m = bn // a
        ok = (bn % a != 0) & (m <= Mmax) & (m >= 1); m, bn = m[ok], bn[ok]
        x = a / bn.astype(np.float64); l = 1.0 / (m + 1); r = 1.0 / m; q = (x - l) / (r - l)
        ms.append(m); qs.append(q); js.append(np.full(len(m), j)); lt.append(np.log(x / l) - q * np.log(r / l))
    m = np.concatenate(ms); q = np.concatenate(qs); jj = np.concatenate(js); lt = np.concatenate(lt)
    o = np.lexsort((q, m)); m, q, jj, lt = m[o], q[o], jj[o], lt[o]
    B = -np.bincount(jj, weights=lt, minlength=nh); Om = np.zeros((nh, nh))
    cuts = np.flatnonzero(np.diff(m)) + 1; st = np.concatenate(([0], cuts)); en = np.concatenate((cuts, [len(m)]))
    for s0, t0_ in zip(st, en):
        mm = int(m[s0]); I = 1.0 / (mm * (mm + 1)); qq = q[s0:t0_]; ii = jj[s0:t0_]
        np.add.at(Om, (ii[:, None], ii[None, :]), I * (np.minimum.outer(qq, qq) - np.outer(qq, qq)))
    mall = np.arange(1, 100_001, dtype=np.float64)
    Cp = 1.0 + float(np.sum(1.0 - mall * (mall + 1) * np.log1p(1.0 / mall) ** 2)) + 1.0 / (12 * 100_000)
    return Om, B, Cp
if __name__ == '__main__':
    t0 = time.time()
    N, Q, K0, K1 = (int(x) for x in sys.argv[1:5]); Mmax = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
    L = math.log(N)
    G = np.array(np.memmap(GF, dtype=np.float64, mode='r', shape=(NMAXF, NMAXF))[:N, :N])
    rats = [(Q // gcd(Q, k), k // gcd(Q, k)) for k in range(Q + 1, Q * K0) if k % Q] + [(2, 2 * k + 1) for k in range(K0, K1 + 1)]
    nh = len(rats)
    Cb = np.empty((N, nh)); Hb = np.empty((nh, nh))
    for i, (a, b) in enumerate(rats):
        for i2 in range(i, nh):
            c, d = rats[i2]; Hb[i, i2] = Hb[i2, i] = a * c * vas(b * c, a * d)
        for j in range(1, N + 1):
            Cb[j - 1, i] = a * vas(j * a, b)  # <rho_{1/j}, rho_{a/b}> = a G(j a, b)
    check(280, 'blocks')
    ks = np.arange(1, N + 1); bn = (np.log(ks) + 1 - EG) / ks
    th = np.array([a / b for a, b in rats]); ch = th * (np.log(1 / th) + 1 - EG)
    Lc = cholesky(G, lower=True, check_finite=False); del G
    y = solve_triangular(Lc, bn, lower=True, check_finite=False); d2 = 1 - float(y @ y)
    Z = solve_triangular(Lc, Cb, lower=True, check_finite=False); S = Hb - Z.T @ Z; rr = ch - Z.T @ y
    e = np.linalg.solve(S, rr); P = float(rr @ e)
    Om, B, Cp = perp_gram(rats, Mmax); Om2, B2, _ = perp_gram(rats, Mmax // 2)
    OmE = 2 * Om - Om2; BE = 2 * B - B2; Omf = OmE - np.outer(BE, BE) / Cp
    R = float(e @ Omf @ e); Qd = P - R; D2 = d2 - P; ev = eigh(S - Omf, eigvals_only=True)
    tail = th < 1 / K0
    Rh = float(e[~tail] @ Omf[np.ix_(~tail, ~tail)] @ e[~tail])
    print(f'N={N} Q={Q} K0={K0} K1={K1} nh={nh}: L^2P={L*L*P:.5f} L^2R={L*L*R:.5f} R/P={R/P:.4f} L^3Q={L**3*Qd:.4f} c_e=LQ/R={L*Qd/R:.3f} '
          f'cos={R/math.sqrt((D2-R)*Qd):.4f} P/D2={P/D2:.4f}<=Q/R={Qd/R:.4f} | eig(S-Om) min {ev.min():+.1e} | tail share of L|e|_1 '
          f'{np.abs(e[tail]).sum()/np.abs(e).sum():.3f} | t={time.time()-t0:.0f}s mem {mem()}', flush=True)
