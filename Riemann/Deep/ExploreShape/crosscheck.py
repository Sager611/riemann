"""explore53 cross-check of the higher moments <S^k>, k <= 4, against an INDEPENDENT implementation of BM's printed
9-fold integrand (GapLambda search/bm/bmfloat.py, extended here from k <= 2 to k <= K; unshifted mollifier, real P).
Compares, at the Wirtinger-optimal P for fixed v, the ratios <(v-S)^k>/<1> with moments.py (twist switched off).
usage: crosscheck.py r theta kappa deg v n"""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'GapLambda', 'search', 'bm'))
sys.path.insert(0, HERE)
import bmfloat as BF


def matsK(th, r, kap, D, n, K):
    z, wz = BF.gl01(n)
    S1, E1, S2, E2 = np.meshgrid(z, z, z, z, indexing='ij')
    W4 = np.einsum('i,j,k,l->ijkl', wz, wz, wz, wz)
    S1, E1, S2, E2, W4 = S1.ravel(), E1.ravel(), S2.ravel(), E2.ravel(), W4.ravel()
    kp = kap * np.pi
    M = [np.zeros((D, D), complex) for _ in range(K + 1)]
    for ix, x in enumerate(z):
        R = 1 - x
        x1 = R * S1; x2 = R * (1 - S1) * E1; x3 = R * S2; x4 = R * (1 - S2) * E2
        jac = (R ** 2 * (1 - S1)) * (R ** 2 * (1 - S2)) * W4 * wz[ix]
        a = 1 - th * (x1 + x3); b = 1 - th * (x2 + x4)
        wt = jac * a * b * x ** (r * r - 1) * (x1 * x2 * x3 * x4) ** (r - 1) * np.exp(1j * th * kp * (x2 - x4))
        F1 = BF.leg01(D, 1 - x - x1 - x2); F2 = BF.leg01(D, 1 - x - x3 - x4)
        base0 = th * (x1 + x2 + x3 + x4 + x)
        for it1, t1 in enumerate(z):
            for it2, t2 in enumerate(z):
                s = a * t1 - b * t2
                g = wt * wz[it1] * wz[it2] * (th * (x1 - x2) + s) * (th * (x3 - x4) + s)
                g = g * BF.E(-th * kp * (x3 - x4) - kp * s) * BF.E(th * kp * (x1 - x2) + kp * s)
                Sv = base0 + a * t1 + b * t2
                for k in range(K + 1):
                    M[k] += (F1 * g[:, None]).T @ F2
                    g = g * Sv
    return [0.5 * (A.real + A.real.T) for A in M]


def ratios(Ms, v, K):
    from math import comb
    C = [sum(comb(k, i) * (-v) ** (k - i) * Ms[i] for i in range(k + 1)) for k in range(K + 1)]
    L = np.linalg.cholesky(C[2]); Li = np.linalg.inv(L)
    ev, U = np.linalg.eigh(Li @ C[0] @ Li.T)
    y = Li.T @ U[:, -1]
    n0 = y @ C[0] @ y
    return ev[-1], [float(y @ C[k] @ y / n0) for k in range(K + 1)]


if __name__ == '__main__':
    r = int(sys.argv[1]); th = float(sys.argv[2]); kap = float(sys.argv[3]); deg = int(sys.argv[4]); v = float(sys.argv[5])
    n = int(sys.argv[6]); K = 4
    MB = matsK(th, r, kap, deg + 1, n, K)
    lamB, mB = ratios(MB, v, K)
    import moments as Mo
    orig = Mo.rho_K
    def rho_untw(X1, X2, X3, X4, c, K):
        tw = np.exp(1j * c * (X3 + X4 - X1 - X2) / 2)
        return [o / tw for o in orig(X1, X2, X3, X4, c, K)]
    Mo.rho_K = rho_untw
    MsG, _ = Mo.setup(r, th, deg, kap, K)
    lamG, mG = ratios(MsG, v, K)
    print('r=%d th=%.4f kappa=%.3f deg=%d v=%.3f n=%d' % (r, th, kap, deg, v, n))
    print('  BM printed integrand (bmfloat ext.): h1 = %.8f  <(v-S)^k>/<1> =' % (lamB / kap ** 2), ' '.join('%.8g' % x for x in mB))
    print('  GapBar kernel (moments.py, no twist): h1 = %.8f  <(v-S)^k>/<1> =' % (lamG / kap ** 2), ' '.join('%.8g' % x for x in mG))
