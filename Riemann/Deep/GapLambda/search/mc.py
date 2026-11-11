"""Independent Monte Carlo check of c, c1(u), c2(u): raw 16-dim integrand exactly as printed in the TeX
(Thm 1.1 l.75-100 at alpha=beta=0; c1 l.617-640; c2 l.648-663), r = 1, no rewriting."""
import sys, time
sys.path.insert(0, '.')
import numpy as np
from flint import fmpq
import core

th = 1249 / 10000; u = 0.6
Pc = [1.0, -3.0, 1.5]          # P(x) = 1 - 3x + 1.5x^2
P = lambda w: Pc[0] + Pc[1] * w + Pc[2] * w * w
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 12345
rng = np.random.default_rng(seed)
N = 150000; reps = int(sys.argv[2]) if len(sys.argv) > 2 else 50
print('seed', seed, 'samples', N * reps)
acc = {'c': [], 'c1': [], 'c2': []}
t0 = time.time()
for rep in range(reps):
    t = rng.random((N * 7, 4))
    ok = (t[:, 0] + t[:, 1] <= 1) & (t[:, 2] + t[:, 3] <= 1) & (t[:, 0] + t[:, 2] <= 1) & (t[:, 1] + t[:, 3] <= 1)
    t = t[ok][:N]
    t1, t2, t3, t4 = t.T
    R = {1: 1 - t1 - t2, 2: 1 - t3 - t4, 3: 1 - t1 - t3, 4: 1 - t2 - t4}
    def tri(Rk):
        a, b = rng.random(N), rng.random(N)
        f = a + b > 1
        a[f], b[f] = 1 - a[f], 1 - b[f]
        return Rk * a, Rk * b
    x3, z3 = tri(R[1]); x4, z4 = tri(R[2]); x1, z1 = tri(R[3]); x2, z2 = tri(R[4])
    wgt = (R[1] ** 2 / 2) * (R[2] ** 2 / 2) * (R[3] ** 2 / 2) * (R[4] ** 2 / 2) / 6.0
    v1, v2, v3, v4 = rng.random((4, N))
    Sx = x1 + x2 + x3 + x4; Sz = z1 + z2 + z3 + z4
    E12 = v1 - v2 + th * (x1 + x2 - z1 - z2 - v1 * Sx + v2 * Sz)
    E34 = v1 - v2 + th * (x3 + x4 - z3 - z4 - v1 * Sx + v2 * Sz)
    base = (1 - th * Sx) * (1 - th * Sz) * E12 * E34
    PP = P(1 - t1 - t2 - x3 - z3) * P(1 - t3 - t4 - x4 - z4) * P(1 - t1 - t3 - x1 - z1) * P(1 - t2 - t4 - x2 - z2)
    G1 = u - th * (t1 + t2 + x1 + x2 + x3 + z3) - v1 * (1 - th * Sx) + v3 * E12
    G2 = u - th * (t3 + t4 + x4 + z1 + z2 + z4) - v2 * (1 - th * Sz) - v3 * E12
    G3 = u - th * (t1 + t3 + x1 + x3 + x4 + z1) - v1 * (1 - th * Sx) + v4 * E34
    G4 = u - th * (t2 + t4 + x2 + z2 + z3 + z4) - v2 * (1 - th * Sz) - v4 * E34
    acc['c'].append(np.mean(wgt * base * PP))
    acc['c1'].append(np.mean(wgt * base * PP * G1 * G2 * G3 * G4))
    acc['c2'].append(np.mean(wgt * base * PP * G1 * G3))
parts = core.poly_parts(fmpq(1249, 10000))
terms = {k: core.split_terms(p) for k, p in parts.items()}
ex = core.evaluate(terms, 1, [fmpq(1), fmpq(-3), fmpq(3, 2)], fmpq(3, 5))
for k in ('c', 'c1', 'c2'):
    m = np.mean(acc[k]); se = np.std(acc[k]) / np.sqrt(reps)
    print('%s  MC %.6e +- %.1e   exact %.6e   z = %.2f' % (k, m, se, float(ex[k]), (m - float(ex[k])) / se))
print('time %.0f s' % (time.time() - t0))
