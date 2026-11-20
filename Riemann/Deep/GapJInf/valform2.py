"""gapjinf: independent validity check of a saved form on explore53's Legendre bubble basis (validate53.basis), with the
degree-adapted scaling g^(min(deg, J) - 1/2) of the basis (well conditioned as g -> 0; same eigenvalues in exact arithmetic).
min over gaps g in [1e-6 pi, pi] of the smallest eigenvalue of the form relative to the H^J Gram matrix.
usage: valform2.py form.json [N]"""
import sys, json
import numpy as np
from validate53 import basis
F = json.load(open(sys.argv[1])); N = int(sys.argv[2]) if len(sys.argv) > 2 else 30
J = F['J']; a = F['a']; phi = {tuple(int(t) for t in k.split(',')): v for k, v in F['phi'].items()}
m = F['moments']
h = -a[0] * m[0] / sum(a[j] * m[j] for j in range(1, J + 1))
G, B = basis(N, J)
worst = (np.inf, None)
for g in np.pi * np.concatenate([np.geomspace(1e-6, 1, 700), [1.0]]):
    sc = g ** (np.minimum(np.arange(N) + 2, J) - 0.5); S2 = np.outer(sc, sc)
    M = sum(a[j] * g ** (1 - 2 * j) * G[j] * S2 for j in range(J + 1))
    for pq, Bm in B.items():
        M = M + phi.get(pq, 0.0) * g ** (-pq[0] - pq[1]) * Bm * S2
    Nrm = sum(g ** (1 - 2 * j) * G[j] * S2 for j in range(J + 1))
    Li = np.linalg.inv(np.linalg.cholesky(Nrm))
    ev = np.linalg.eigvalsh(Li @ M @ Li.T)[0]
    if ev < worst[0]:
        worst = (ev, g)
print('%s: J=%d h=%.7f (LP %.7f); scaled polynomial-basis check N=%d, g in [1e-6, 1] pi: min rel eig %.3e at g/pi=%.3e' % (sys.argv[1], J, h, F['h'], N, worst[0], worst[1] / np.pi))
