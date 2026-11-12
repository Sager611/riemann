"""Covering certificate for the TEN class over ALL 16-harmonic windows c' in R^16 (c_0 = 1 fixed).
For an S-periodic configuration X of density nu0 = p/S:  e(X;c) = c^T Q_X c - nu0,  H(c) = H_MT - c'^T W c',
slack_X(c) = nu0 - H(c) - e(X;c) = 2 nu0 - H_MT + c'^T W c' - c^T Q_X c.
If lambda in the simplex makes M(lambda) = [[2nu0 - H_MT - Qb00, -qb^T], [-qb, W - Qb']] PSD (Qb = sum lambda_i Q_i),
then for every window some X_i has slack >= 0, so no TEN-type certificate on any window proves kappa > nu0.
Cutting plane: add configurations relaxed at the SDP's worst window until the SDP value t >= 0."""
import sys, json
import numpy as np
from scipy.optimize import minimize
import cvxpy as cp
from winopt import Ttensor, WJ, H_MT, I1SQ, Hc
from tenwin import TEN_C, MT_C

S, P = int(sys.argv[1]), int(sys.argv[2])
NU0 = P / S
rng = np.random.default_rng(7)
MS = np.arange(1, S)
TT = np.array([Ttensor(m / S) for m in MS])
TT = 0.5 * (TT + np.transpose(TT, (0, 2, 1)))
T0 = Ttensor(0.0); T0 = 0.5 * (T0 + T0.T)
BLK = {'A': '.xx', 'B': '.xx.xx.xxx', 'C': '.xx.xx.xx.xxx', 'D': '..xxx.xx.xxx', 'E': '.xx.xxx'}


def seeds():
    out = []
    L = {k: (len(v), v.count('x')) for k, v in BLK.items()}
    for b in range(0, 9):
        for cc in range(0, 7):
            for d in range(0, 8):
                for e in range(0, 13):
                    rest_len = S - 10 * b - 13 * cc - 12 * d - 7 * e
                    rest_p = P - 7 * b - 9 * cc - 8 * d - 5 * e
                    if rest_len < 0 or rest_len % 3 or rest_p != 2 * (rest_len // 3):
                        continue
                    a = rest_len // 3
                    blocks = ['A'] * a + ['B'] * b + ['C'] * cc + ['D'] * d + ['E'] * e
                    for _ in range(3):
                        rng.shuffle(blocks)
                        s = ''.join(BLK[k] for k in blocks)
                        out.append(np.array([i for i, ch in enumerate(s) if ch == 'x'], float))
    return out


SEEDS = seeds()


def Qmat(x):
    mu = np.exp(2j * np.pi * np.outer(MS, x) / S).sum(axis=1)
    return (T0 * P * P + 2 * np.einsum('m,mij->ij', np.abs(mu) ** 2, TT)) / (I1SQ * S * S)


def relax_at(c, k=4):
    kh = np.einsum('i,mij,j->m', c, TT, c) / I1SQ
    k0 = c @ T0 @ c / I1SQ

    def f(x):
        ph = np.exp(2j * np.pi * np.outer(MS, x) / S)
        mu = ph.sum(axis=1)
        e = (k0 * P * P + 2 * np.sum(kh * np.abs(mu) ** 2)) / S ** 2 - P / S
        g = (4 * np.real(np.conj(mu)[:, None] * (2j * np.pi * MS[:, None] / S) * ph) * kh[:, None]).sum(axis=0) / S ** 2
        return e, g
    res = []
    for x0 in SEEDS:
        r = minimize(f, x0, jac=True, method="L-BFGS-B", options={"maxiter": 3000, "gtol": 1e-14})
        res.append((r.fun, np.sort(np.mod(r.x, S))))
    res.sort(key=lambda t: t[0])
    return res[:k]


def slack(c, Q):
    return 2 * NU0 - H_MT + c[1:] @ (WJ[1:] * c[1:]) - c @ Q @ c


def sdp(Qs):
    lam = cp.Variable(len(Qs), nonneg=True)
    t = cp.Variable()
    Qb = sum(lam[i] * Qs[i] for i in range(len(Qs)))
    W = np.diag(WJ)
    M = cp.bmat([[cp.reshape(2 * NU0 - H_MT - Qb[0, 0], (1, 1)), -cp.reshape(Qb[0, 1:], (1, 16))],
                 [-cp.reshape(Qb[1:, 0], (16, 1)), W[1:, 1:] - Qb[1:, 1:]]])
    Ms = 0.5 * (M + M.T)
    cons = [cp.sum(lam) == 1, Ms - t * np.eye(17) >> 0]
    prob = cp.Problem(cp.Maximize(t), cons)
    prob.solve(solver=cp.CLARABEL)
    Y = cons[1].dual_value
    return t.value, lam.value, Y


if __name__ == "__main__":
    print("S", S, "P", P, "nu0", NU0, "seeds", len(SEEDS), flush=True)
    cs = [np.array(json.load(open("winopt_best.json"))["c"]), TEN_C.copy(), MT_C.copy()]
    lib = []
    for c in cs:
        for e, x in relax_at(c):
            lib.append(x)
        print("c-slack best at start window:", max(slack(c, Qmat(x)) for x in lib), flush=True)
    for it in range(25):
        Qs = [Qmat(x) for x in lib]
        t, lam, Y = sdp(Qs)
        print(f"iter {it}: lib {len(lib)} SDP t = {t:.3e}", flush=True)
        if t > 1e-9:
            break
        Y = 0.5 * (Y + Y.T)
        w, V = np.linalg.eigh(Y)
        v = V[:, -1]
        cw = v / v[0]
        print("   worst window H =", Hc(cw), " max_i slack =", max(slack(cw, Q) for Q in Qs), flush=True)
        for e, x in relax_at(cw):
            lib.append(x)
    json.dump({"S": S, "P": P, "lam": [float(v) for v in lam], "x": [x.tolist() for x in lib], "t": float(t)},
              open(f"cover_{S}_{P}.json", "w"))
    print("done t =", t)
