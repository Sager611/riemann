"""gapjinf: BM main-term spectral moments C_k = <(S-v0)^k> kappa^k, k <= K, as (d+1)x(d+1) forms in ExploreShape's
amplifier basis (same as ExploreShape joint.central_forms; Gm vectorised). Read-only imports: GapBar gcore, ExploreShape moments.
usage: mom.py kappa K out.npz [d] [v0] [check]"""
import sys, os, time
from math import comb, factorial
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'GapBar'))
sys.path.insert(0, os.path.join(HERE, '..', 'ExploreShape'))
import gcore as G
import moments_c as Mo        # centred tau-moments (CSH); CSH = 0 reproduces ExploreShape/moments.py


def central_forms_fast(r, th, d, kap, K, v0, D=40, csh=0.0):
    Mo.CSH[0] = csh
    Nb = d + 2 * r + 6
    H = [h.real for h in Mo.arrays_K(kap * np.pi, r, K, D)]
    N0 = Nb + D + 2
    B = G.Basis(th, N0, nq=N0 + 2 * r * r + 12 + K)
    q0 = r * r - 1
    Wc = [(B.V * (B.w * B.X ** q0 * (B.X - v0 + csh) ** p)[:, None]).T @ B.V for p in range(K + 1)]
    Ipow = [np.eye(N0)]
    for k in range(D):
        Ipow.append(B.I @ Ipow[-1])
    s0 = 2 * r - 1
    fac = np.array([float(factorial(k + s0)) for k in range(D)])
    Ir = np.linalg.matrix_power(B.I, 2 * r)[:Nb, :d + 1]
    Q, R = np.linalg.qr(Ir)
    IPQ = np.array([P[:, :Nb] @ Q for P in Ipow[:D]])          # D x N0 x (d+1)
    def Gm(p, h):
        hp = fac[:, None] * h * fac[None, :]
        Y = np.tensordot(hp, IPQ, axes=([1], [0]))             # D x N0 x (d+1)
        WY = np.matmul(Wc[p][None, :, :], Y)
        out = np.tensordot(IPQ, WY, axes=([0, 1], [0, 1]))
        return 0.5 * (out + out.T)
    C = [sum(comb(k, i) * Gm(i, H[k - i]) for i in range(k + 1)) * kap ** k for k in range(K + 1)]
    return C, R, B


if __name__ == '__main__':
    kap = float(sys.argv[1]); K = int(sys.argv[2]); out = sys.argv[3]
    d = int(sys.argv[4]) if len(sys.argv) > 4 else 80
    v0 = float(sys.argv[5]) if len(sys.argv) > 5 else 1.3617
    csh = float(os.environ.get('CSH', '0'))
    t0 = time.time()
    C, R, B = central_forms_fast(2, 0.2499, d, kap, K, v0, csh=csh)
    print('centre shift CSH = %.3f' % csh)
    print('d=%d kappa=%.4f K=%d v0=%.4f: %.1f s' % (d, kap, K, v0, time.time() - t0), flush=True)
    np.savez(out, C=np.array(C), R=R, norm=B.norm[:d + 1], kap=kap, v0=v0, d=d, K=K)
    if len(sys.argv) > 6:
        sys.path.insert(0, os.path.join(HERE, '..', 'ExploreShape'))
        import joint as Jt
        C2 = Jt.central_forms(2, 0.2499, d, kap, K, v0)
        print('max rel diff vs joint.central_forms:', max(np.abs(a - b).max() / np.abs(b).max() for a, b in zip(C, C2)))
