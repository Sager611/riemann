# t1_mayer.py -- explore3, transfer 1 (N12/N13/N7): Mayer's transfer operator for PSL(2,Z)
# L_b f(z) = sum_{m>=1} (z+m)^{-2b} f(1/(z+m)) on holomorphic f on the disc |z-1| < 3/2.
# Matrix in the basis (z-1)^k, entries via Hurwitz zeta (Mayer's continuation to all b != (1-n)/2):
#   M[n,k] = (-1)^(n+k) sum_{j=0}^{k} C(k,j) (-1)^j (2b+j)_n / n! * zeta(2b+j+n, 2).
# Selberg zeta of PSL(2,Z): Z(b) = det(1 - L_b) det(1 + L_b)  (Mayer 1991, Efrat 1993).
# Run: uv run --no-project --with mpmath --with numpy python t1_mayer.py <mode>
import sys, time, math
import mpmath as mp
import numpy as np

def build_M(b, K, dps=40):
    mp.mp.dps = dps
    b = mp.mpc(b)
    Z = [mp.zeta(2*b + m, 2) for m in range(2*K - 1)]
    P = [[None]*K for _ in range(K)]
    for j in range(K):
        a = 2*b + j
        p = mp.mpf(1)
        for n in range(K):
            P[j][n] = p
            p = p * (a + n) / (n + 1)
    M = mp.matrix(K, K)
    for k in range(K):
        for n in range(K):
            s = mp.mpc(0)
            for j in range(k + 1):
                c = math.comb(k, j)
                term = c * P[j][n] * Z[j + n]
                s = s + term if j % 2 == 0 else s - term
            M[n, k] = s if (n + k) % 2 == 0 else -s
    return M

def to_np(M):
    K = M.rows
    return np.array([[complex(M[i, j]) for j in range(K)] for i in range(K)])

def eigs(b, K, dps=40):
    ev = np.linalg.eigvals(to_np(build_M(b, K, dps)))
    return ev[np.argsort(-abs(ev))]

def dets(b, K, dps=40):
    M = build_M(b, K, dps)
    I = mp.eye(K)
    return mp.det(I - M), mp.det(I + M)

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'sanity'
    t0 = time.time()
    if mode == 'sanity':
        # Gauss-Kuzmin-Wirsing operator b = 1: eigenvalues 1, -0.3036630028987326586...
        for K in (20, 30, 40):
            ev = eigs(1, K)
            print('b=1 K=%d top eigenvalues:' % K, ', '.join('%.15f%+.1ei' % (e.real, e.imag) for e in ev[:3]))
        print('time %.1fs' % (time.time() - t0))
    elif mode == 'zeros':
        g = [14.134725141734693790, 21.022039638771554993, 25.010857580145688763,
             30.424876125859513210, 32.935061587739189691]
        even1, odd1 = 13.779751351891, 9.533695261353
        for K in (30, 40, 50):
            print('--- K =', K)
            for gg in g[:3]:
                b = mp.mpc(0.25, gg / 2)
                dm, dp = dets(b, K)
                bn = mp.mpc(0.25, gg / 2 + 0.05)
                dmn, dpn = dets(bn, K)
                print('b = rho/2, gamma=%.6f: |det(1-L)|=%.3e (at +0.05i: %.3e)  |det(1+L)|=%.3e (at +0.05i: %.3e)'
                      % (gg, abs(dm), abs(dmn), abs(dp), abs(dpn)))
            for name, r in (('even1', even1), ('odd1', odd1)):
                b = mp.mpc(0.5, r)
                dm, dp = dets(b, K)
                bn = mp.mpc(0.5, r + 0.05)
                dmn, dpn = dets(bn, K)
                print('cusp %s r=%.6f: |det(1-L)|=%.3e (+0.05i: %.3e) |det(1+L)|=%.3e (+0.05i: %.3e)'
                      % (name, r, abs(dm), abs(dmn), abs(dp), abs(dpn)))
        print('time %.1fs' % (time.time() - t0))

    elif mode == 'conv':
        for (sg, tt) in ((1.0, 20.0), (0.5, 20.0), (0.25, 20.0), (1.0, 35.0), (0.5, 35.0), (0.25, 35.0)):
            row = []
            for K in (40, 50, 60):
                t1 = time.time()
                ev = eigs(mp.mpc(sg, tt), K, dps=40)
                row.append('K=%d r=%.10f |l2|=%.6f (%.1fs)' % (K, abs(ev[0]), abs(ev[1]), time.time() - t1))
            print('b=%.2f+%.1fi: ' % (sg, tt) + ' | '.join(row))
