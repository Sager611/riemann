"""explore53: J-th order local criterion, SDP with cutting planes (validated forms).
Base: sdpJ LMIs (poly basis N0, s-grid); cuts: for each violation found on a fine s-grid with a larger basis Nv
(validate.basis), add the linear constraint v^T M(s) v >= tol (v = violating eigenvector).  Repeat until the smallest
eigenvalue (relative to the H^J Gram matrix at scale g) is >= -1e-9 on the fine grid.
usage: from sdpcut import sdp_cut"""
import numpy as np
import cvxpy as cp
import sdpJ
import validate as V


def lmi(a, ph, G, B, s, J):
    M = sum(a[j] * s ** (J - j) * G[j] for j in range(J + 1))
    for (p_, q), Bm in B.items():
        M = M + ph[(p_, q)] * s ** ((2 * J - 1 - p_ - q) // 2) * Bm
    return M


def scaled_form(av, phv, Gv, Bv, g, J):
    M = sum(av[j] * g ** (1 - 2 * j) * Gv[j] for j in range(J + 1))
    for (p_, q), Bm in Bv.items():
        M = M + phv[(p_, q)] * g ** (-p_ - q) * Bm
    Nrm = sum(g ** (1 - 2 * j) * Gv[j] for j in range(J + 1))
    return M, Nrm


def sdp_cut(m, J, N0=14, Nv=28, nsig=40, rounds=8, smin=1e-8, nfine=300, verbose=False, margin=0.0):
    G, B = sdpJ.basis(N0, J)
    Gv, Bv = V.basis(Nv, J)
    a = cp.Variable(J + 1); ph = {k: cp.Variable() for k in B}
    cons = [a[1] == 1, a[J] >= 0] + [cp.abs(v) <= 1e3 for v in ph.values()] + [cp.abs(a) <= 1e3]
    for s in np.concatenate([[0.0], np.pi ** 2 * np.geomspace(smin, 1, nsig)]):
        M = lmi(a, ph, G, B, s, J) - margin * sum(s ** (J - j) * G[j] for j in range(J + 1))
        cons.append(0.5 * (M + M.T) >> 0)
    # exact small-gap conditions (Lemmas 3, 4): C >= 0 and the limit x R(x) -> L_inf >= 0 (Schur form for J = 6)
    am = 1e-8
    if J == 5:
        C5 = cp.bmat([[a[3] + ph[(3, 2)], ph[(3, 4)] / 2], [ph[(3, 4)] / 2, a[5]]])
        L5 = cp.bmat([[a[2] + ph[(1, 2)], (ph[(3, 2)] + ph[(1, 4)]) / 2], [(ph[(3, 2)] + ph[(1, 4)]) / 2, a[4] + ph[(3, 4)]]])
        cons += [C5 >> am * np.eye(2), L5 >> am * np.eye(2)]
    if J == 6:
        C6 = cp.bmat([[a[3] + ph[(3, 2)], (ph[(3, 4)] + ph[(5, 2)]) / 2], [(ph[(3, 4)] + ph[(5, 2)]) / 2, a[5] + ph[(5, 4)]]])
        L6 = cp.bmat([[a[2] + ph[(1, 2)], (ph[(1, 4)] + ph[(3, 2)]) / 2, ph[(5, 2)] / 2],
                      [(ph[(1, 4)] + ph[(3, 2)]) / 2, a[4] + ph[(3, 4)], ph[(5, 4)] / 2],
                      [ph[(5, 2)] / 2, ph[(5, 4)] / 2, a[6]]])
        cons += [C6 >> am * np.eye(2), L6 >> am * np.eye(3)]
    fine = np.concatenate([np.pi ** 2 * np.geomspace(smin, 1, nfine), [np.pi ** 2]])
    cuts = []
    for rd in range(rounds):
        pr = cp.Problem(cp.Minimize(sum(a[j] * m[j] for j in range(J + 1))), cons + cuts)
        pr.solve(solver=cp.CLARABEL)
        av = a.value; phv = {k: float(v.value) for k, v in ph.items()}
        worst = []
        for s in fine:
            g = np.sqrt(s)
            M, Nrm = scaled_form(av, phv, Gv, Bv, g, J)
            L = np.linalg.cholesky(Nrm); Li = np.linalg.inv(L)
            ev, U = np.linalg.eigh(Li @ M @ Li.T)
            if ev[0] < margin / 2 - 1e-9:
                worst.append((ev[0], g, Li.T @ U[:, 0]))
        mn = min([w[0] for w in worst], default=0.0)
        if verbose:
            print('    round %d: value %.6g, violations %d, min %.3g' % (rd, pr.value, len(worst), mn), flush=True)
        if not worst:
            return pr.value, av, phv, mn, rd
        worst.sort(key=lambda w: w[0])
        for ev0, g, vv in worst[:6]:
            # constraint on the unscaled-s form: divide the g-scaled form by g^{1-2J} (positive) -> same sign
            Mexpr = sum(a[j] * g ** (1 - 2 * j) * float(vv @ Gv[j] @ vv) for j in range(J + 1))
            for (p_, q), Bm in Bv.items():
                Mexpr = Mexpr + ph[(p_, q)] * g ** (-p_ - q) * float(vv @ Bm @ vv)
            nrm = sum(g ** (1 - 2 * j) * float(vv @ Gv[j] @ vv) for j in range(J + 1))
            cuts.append(Mexpr >= (1e-9 + margin) * nrm)
    return pr.value, av, phv, mn, rounds
