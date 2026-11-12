"""Extend the J=16 covering certificate to J harmonics (same configurations, same lambda): float min eigenvalue of M_J.
Escape-route test (lesson 191): do more cosine harmonics escape the TEN-class ceiling?"""
import sys, json
import numpy as np

J = int(sys.argv[3]) if len(sys.argv) > 3 else 32
S, P = int(sys.argv[1]), int(sys.argv[2])
d = json.load(open(f"cover_{S}_{P}.json"))
OM = np.array([np.sqrt(2.0)] + [2 * np.pi * j for j in range(1, J + 1)])
TH0 = 1 / np.sqrt(2); I1SQ = 2 * np.sin(TH0) ** 2
WJ = np.array([0.0] + [(0.5 - 1 / (4 * np.pi ** 2 * j ** 2)) / I1SQ for j in range(1, J + 1)])
H_MT = 1.5 - 1 / np.tan(TH0) / np.sqrt(2)
xg, wg = np.polynomial.legendre.leggauss(300)


def T(al):
    lo, hi = -0.5, 0.5 - al
    s = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo); w = 0.5 * (hi - lo) * wg
    A = np.cos(np.outer(s, OM)); B = np.cos(np.outer(s + al, OM))
    M = (A * w[:, None]).T @ B
    return 0.5 * (M + M.T)


Ts = np.array([T(m / S) for m in range(S)])
Qb = np.zeros((J + 1, J + 1))
for lam, x in zip(d["lam"], d["x"]):
    if lam <= 1e-14:
        continue
    x = np.array(x)
    mu = np.exp(2j * np.pi * np.outer(np.arange(S), x) / S).sum(axis=1)
    wts = 2 * np.abs(mu) ** 2; wts[0] = P * P
    Qb += lam * np.einsum('m,mij->ij', wts, Ts) / (I1SQ * S * S)
nu0 = P / S
M = -Qb.copy()
M[0, 0] += 2 * nu0 - H_MT
M[1:, 1:] += np.diag(WJ[1:])
ev = np.linalg.eigvalsh(M)
print(f"J={J}: min eig M_J = {ev[0]:.3e} (J=16 block min eig {np.linalg.eigvalsh(M[:17,:17])[0]:.3e});  diag W-Q for j>16: min {np.min(np.diag(M)[17:]) if J>16 else 0:.3e}")
