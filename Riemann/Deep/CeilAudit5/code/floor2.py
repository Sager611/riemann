#!/usr/bin/env python3
# ceilaudit5 floor2.py (N1): least eigenvalue of Weil's zero-side form q[g] = sum_rho fhat(rho), f = g*g~, on L^2[-a/2, a/2]
# with the TRUE zeros 1/2 +- i gamma_k (gamma_k <= GMAX) plus ONE hypothetical off-line quadruple 1/2 +- delta +- i gamma0.
# on-line zero: |G(i g)|^2 + |G(-i g)|^2 ; off-line pair at +-gamma0: 2Re[G(d+i g0) conj G(-d+i g0)] + same at -g0;  G(z) = int g(x) e^{zx} dx.
import sys, math, json
import numpy as np, mpmath as mp
GMAX = float(sys.argv[1]); g0 = float(sys.argv[2])
zs = []; k = 1
while True:
    z = float(mp.zetazero(k).imag)
    if z > GMAX: break
    zs.append(z); k += 1
print('zeros used: %d (gamma <= %.1f); off-line pair at gamma0 = %.3f' % (len(zs), GMAX, g0), flush=True)
out = []
for d in (0.05, 0.1):
    for a in (5.0, 10.0, 20.0, 30.0, 40.0):
        N = int(max(400, a/0.03)); x = (np.arange(N) + 0.5)*a/N - a/2; w = a/N
        def e(z): return np.exp(z*x)*w           # G(z) = e(z) . g
        H = np.zeros((N, N), dtype=complex)
        for gm in zs:
            for s in (1, -1):
                v = e(1j*s*gm); H += np.outer(v.conj(), v)      # g^* (conj(v) v^T) g = |v.g|^2
        for s in (1, -1):
            u = e(d + 1j*s*g0); v = e(-d + 1j*s*g0)
            # 2 Re[(u.g) conj(v.g)] = g^* (conj(v) u^T + conj(u) v^T) g
            H += np.outer(v.conj(), u) + np.outer(u.conj(), v)
        # generalized problem w.r.t. L^2 norm: ||g||^2 = w * g^* g
        lam = float(np.linalg.eigvalsh(H).min())/w
        pair_only = -(math.sinh(d*a)/d - a)
        out.append(dict(delta=d, a=a, N=N, lam=lam, pair_only=pair_only))
        print('delta %.2f a %4.1f N %4d: least eigenvalue with true zeros %.4e | pair alone %.4e | ratio %.3f' % (d, a, N, lam, pair_only, lam/pair_only), flush=True)
        del H
json.dump(out, open('.scratch/ceilaudit5/floor2.json', 'w'), indent=1)
