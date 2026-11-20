"""gapjinf: dictionary (genuine-function) upper bound for the order-J moment method (Lemma C of REPORT.md).
h_dict = max { -a0*m0 : sum_{j>=1} a_j m_j = 1, P_a(x) = sum a_j x^j satisfies N_b^* diag(P_a(xi^2)) N_b >= 0 for every
dictionary block b (Bloch fibre of a periodic zero pattern with gaps <= pi) and P_a >= 0 on a grid of [1, Xbig] }.
Every valid order-J form (any per-gap form, any inequality valid for all real F with zero gaps <= pi) has h_J <= h_dict;
h_dict <= 1 means no such inequality of order J succeeds for these moments."""
import numpy as np
import cvxpy as cp


def pattern_blocks(gaps, nth, Xi, tol=1e-9):
    gaps = np.asarray(gaps, float)
    z = np.concatenate([[0.0], np.cumsum(gaps)[:-1]]); P = gaps.sum(); base = 2 * np.pi / P
    out = []
    for t in range(nth):
        th = (t + 0.5) / nth * base
        n = np.arange(int(np.floor((-Xi - th) / base)) - 1, int(np.ceil((Xi - th) / base)) + 2)
        xi = th + base * n
        k = np.abs(xi) <= Xi
        n, xi = n[k], xi[k]
        E = np.exp(2j * np.pi * np.outer(z, n) / P)
        u, s, vh = np.linalg.svd(E)
        rk = int((s > tol * max(1, s.max())).sum())
        N = vh[rk:].conj().T
        if N.shape[1] == 0:
            continue
        out.append((xi, N))
    return out


def real_rep(H):
    A, B = H.real, H.imag
    return np.block([[A, -B], [B, A]])


def sos_halfline(a, J, c):
    """P(x) = s0(x) + (x - c) s1(x), s0, s1 SOS (Markov-Lukacs): P >= 0 on [c, inf)."""
    n0 = J // 2 + 1; n1 = (J - 1) // 2 + 1
    G0 = cp.Variable((n0, n0), PSD=True); G1 = cp.Variable((n1, n1), PSD=True)
    cons = []
    for k in range(J + 1):
        e0 = sum(G0[i, k - i] for i in range(max(0, k - n0 + 1), min(k, n0 - 1) + 1)) if k <= 2 * (n0 - 1) else 0
        e1 = 0
        if 1 <= k <= 2 * (n1 - 1) + 1:
            e1 = e1 + sum(G1[i, k - 1 - i] for i in range(max(0, k - 1 - n1 + 1), min(k - 1, n1 - 1) + 1))
        if k <= 2 * (n1 - 1):
            e1 = e1 - c * sum(G1[i, k - i] for i in range(max(0, k - n1 + 1), min(k, n1 - 1) + 1))
        cons.append(a[k] == e0 + e1)
    return cons


def h_dict(m, blocks, J, X0=6.25, xbig=400.0, ngrid=400, solver='CLARABEL', ret=False):
    """m: moments m_0..m_J (m_j = int xi^{2j} dnu). Scaled variable x' = x/X0."""
    m = np.asarray(m[:J + 1], float); ms = m / X0 ** np.arange(J + 1)
    a = cp.Variable(J + 1)
    cons = [ms[1:] @ a[1:] == 1]
    for xi, N in blocks:
        x = xi ** 2 / X0
        Ms = [real_rep(N.conj().T @ (x[:, None] ** j * N)) for j in range(J + 1)]
        S = sum(a[j] * Ms[j] for j in range(J + 1))
        cons.append(0.5 * (S + S.T) >> 0)
    cons += sos_halfline(a, J, 1.0 / X0)          # P >= 0 on [1, inf): sin(w x), w >= 1, are genuine atoms
    pr = cp.Problem(cp.Maximize(-a[0] * ms[0]), cons)
    pr.solve(solver=solver)
    if ret:
        return pr.value, a.value * X0 ** -np.arange(J + 1.0), pr.status
    return pr.value, pr.status
