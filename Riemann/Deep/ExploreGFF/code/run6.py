# run6.py (explore17): adversarial tests for the Phi-ladder deflation.
# (a) perturbed prime weight w_2 -> (1+eps) w_2: first-order shift of the ladder eigenvalues and the deflated complement;
# (b) "moved zeros" form Q_test = Q_zeta + 4 Re F(z)^2 (even block; an off-line quartet +-z, +-conj z added; functional-equation
#     symmetric, no Euler product), its own radical Phi_test = P(iD) Phi with P(t) = (t^2 - z^2)(t^2 - conj z^2).
import os, sys, time, math, numpy as np, mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'WeilSimple'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, arb_mat, ctx, fmpq
from weil_matrix import WeilMat
from run5_lib import phi_derivs, basis_vals
t0 = time.time()
def mpmat(T, n, dps):
    M = mp.matrix(n, n)
    for a in range(n):
        for b in range(n): M[a, b] = mp.mpf(T[a, b].mid().str(dps + 10, radius=False))
    return M
def Fz_even(L, N, z):   # F_j(z) = int_W phi_j(x) e^{izx} dx for the even basis, mpmath complex
    out = []
    for j in range(N + 1):
        w = 2*mp.pi*j/L; Nj = 1/mp.sqrt(L) if j == 0 else mp.sqrt(2/L)
        if j == 0: out.append(Nj*2*mp.sin(z*L/2)/z)
        else: out.append(Nj*(mp.sin((z - w)*L/2)/(z - w) + mp.sin((z + w)*L/2)/(z + w)))
    return out
Jmax = 29
for cq, NN, dps in ((5, 40, 60), (10, 40, 90)):
    L = math.log(cq); ctx.prec = int(dps*3.4) + 64; mp.mp.dps = dps
    Wm = WeilMat(fmpq(cq)); T, idx = Wm.block(NN, 'even'); n = len(idx)
    Q = mpmat(T, n, dps)
    # prime-2 matrix P2 = -2 w_2 C(log 2)
    ln2 = arb(2).log(); w2 = ln2/arb(2).sqrt()
    P2 = mp.matrix(n, n)
    for a in range(n):
        for b in range(n): P2[a, b] = mp.mpf((-2*w2*Wm.C(idx[a], idx[b], ln2, 'even')).mid().str(dps + 10, radius=False))
    E, V = mp.eigsy(Q); order = sorted(range(n), key=lambda k: E[k])
    print('c=%d even ladder:' % cq, ' '.join(mp.nstr(E[k], 3) for k in order[:6]), flush=True)
    # (a) first-order coefficient d lambda_k / d eps = <psi_k, P2 psi_k>
    for k in order[:3]:
        v = V[:, k]; d = (v.T*P2*v)[0]
        print('   (a) k=%d lambda=%s  d lambda/d eps = %s  -> fails at eps ~ %s' % (order.index(k), mp.nstr(E[k], 3), mp.nstr(d, 4), mp.nstr(-E[k]/d, 3)), flush=True)
    # deflated complement for eps = +-0.1
    xg, wg = np.polynomial.legendre.leggauss(500); x = xg*L/2; wq = wg*L/2
    mp.mp.dps = 50
    D = np.array([[float(vv) for vv in phi_derivs(xx, Jmax)] for xx in x])
    mp.mp.dps = dps
    Bv = basis_vals(L, 'even', NN, x)
    Qf = np.array([[float(Q[a, b]) for b in range(n)] for a in range(n)]); P2f = np.array([[float(P2[a, b]) for b in range(n)] for a in range(n)])
    Kd = 4 if cq == 5 else 9
    def compl_min(Qm, cols, K):
        U, s, _ = np.linalg.svd(cols[:, :K]/np.linalg.norm(cols[:, :K], axis=0), full_matrices=True)
        r = int(np.sum(s > 1e-10*s[0])); Pc = U[:, r:]; return np.linalg.eigvalsh(Pc.T @ Qm @ Pc)[0]
    colsZ = np.array([Bv @ (wq*D[:, j]) for j in range(0, Jmax + 1, 2)]).T
    for eps in (0.1, -0.1, 0.01):
        print('   (a) eps=%+g: complement of V_%d: lambda_min = %.3f' % (eps, Kd, compl_min(Qf + eps*P2f, colsZ, Kd)), flush=True)
    # (b) moved zeros
    for z in (mp.mpc(10, 0.3), mp.mpc(30, 0.3), mp.mpc(60, 0.3)):
        Fv = Fz_even(L, NN, z)
        R = mp.matrix(n, n)
        for a in range(n):
            for b in range(n): R[a, b] = 4*mp.re(Fv[a]*Fv[b])
        Et = mp.eigsy(Q + R, eigvals_only=True); lmin = min(Et)
        # own radical: Phi_test^{(2j)} = Phi^{(2j+4)} + 2 Re(z^2) Phi^{(2j+2)} + |z|^4 Phi^{(2j)}
        a2 = float(2*mp.re(z*z)); a0 = float(abs(z)**4)
        colsT = np.array([Bv @ (wq*(D[:, 2*j + 4] + a2*D[:, 2*j + 2] + a0*D[:, 2*j])) for j in range(0, (Jmax - 4)//2 + 1)]).T
        Rf = np.array([[float(R[a, b]) for b in range(n)] for a in range(n)])
        cm = [compl_min(Qf + Rf, colsT, K) for K in (Kd - 2, Kd, Kd + 2)]
        cmz = [compl_min(Qf + Rf, colsZ, K) for K in (Kd, Kd + 2)]
        print('   (b) z=%s: lambda_min(Q_test) = %s | complement of own ladder K=%d,%d,%d: %s | complement of zeta ladder K=%d,%d: %s' % (
            mp.nstr(z, 4), mp.nstr(lmin, 3), Kd - 2, Kd, Kd + 2, ' '.join('%.2e' % v for v in cm), Kd, Kd + 2, ' '.join('%.2e' % v for v in cmz)), flush=True)
    print('   time', round(time.time() - t0, 1), flush=True)
