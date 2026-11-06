# run5.py (explore17): deflation test. V_K = span{Phi^{(j)} 1_W} (even j for the even block, odd j for the odd block),
# j exact via the recursion (2 z d/dz) on z^a e^{-b z}, z = e^{2u}. Report lambda_min of Q restricted to the L2-complement of V_K
# inside E_N (Q from WeilMat Arb closed forms, rounded to float64 after assembly).
import os, sys, time, math, numpy as np, mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'WeilSimple'))
from flint import arb, arb_mat, ctx, fmpq
from weil_matrix import WeilMat
t0 = time.time(); mp.mp.dps = 50
NT = 10   # terms n of the Phi series
def phi_derivs(u, J):   # [Phi^{(j)}(u) for j = 0..J]
    z = mp.e**(2*mp.mpf(u)); out = [mp.mpf(0)]*(J + 1)
    for n in range(1, NT + 1):
        b = mp.pi*n*n
        for (cf, a) in ((2*mp.pi**2*n**4, mp.mpf(9)/4), (-3*mp.pi*n*n, mp.mpf(5)/4)):
            poly = [cf]   # coefficients of z^{a+m} e^{-bz}
            for j in range(J + 1):
                out[j] += mp.e**(-b*z)*z**a*mp.fsum(poly[m]*z**m for m in range(len(poly)))
                new = [mp.mpf(0)]*(len(poly) + 1)
                for m, cm in enumerate(poly):
                    new[m] += 2*(a + m)*cm; new[m + 1] += -2*b*cm
                poly = new
    return out
def basis_vals(L, sector, N, x):
    idx = np.arange(0, N + 1) if sector == 'even' else np.arange(1, N + 1)
    w = 2*math.pi*idx/L
    if sector == 'even':
        Nj = np.where(idx == 0, 1/math.sqrt(L), math.sqrt(2/L)); return np.array([Nj[a]*np.cos(w[a]*x) for a in range(len(idx))])
    return np.array([math.sqrt(2/L)*np.sin(w[a]*x) for a in range(len(idx))])
Jmax = 29
for cq, NN in ((5, 40), (10, 40), (20, 50)):
    L = math.log(cq); ctx.prec = 400
    xg, wg = np.polynomial.legendre.leggauss(500); x = xg*L/2; wq = wg*L/2
    D = np.array([[float(v) for v in phi_derivs(xx, Jmax)] for xx in x])   # (nodes, j)
    Wm = WeilMat(fmpq(cq))
    for sec in ('even', 'odd'):
        T, idx = Wm.block(NN, sec); n = len(idx)
        Q = np.array([[float(T[a, b].mid()) for b in range(n)] for a in range(n)])
        Bv = basis_vals(L, sec, NN, x)
        js = list(range(0 if sec == 'even' else 1, Jmax + 1, 2))
        cols = np.array([Bv @ (wq*D[:, j]) for j in js]).T   # coefficients of Phi^{(j)} 1_W
        ev = np.linalg.eigvalsh(Q)
        line = []
        for K in range(0, min(len(js), 15) + 1):
            if K == 0: lm = ev[0]
            else:
                # orthonormal basis of span via SVD (drop numerically dependent directions)
                U, s, _ = np.linalg.svd(cols[:, :K] / np.linalg.norm(cols[:, :K], axis=0), full_matrices=True)
                r = int(np.sum(s > 1e-10*s[0])); Pc = U[:, r:]
                lm = np.linalg.eigvalsh(Pc.T @ Q @ Pc)[0]
            line.append('%d:%.1e' % (K, lm))
        print('c=%d N=%d %s  bottom eigs: %s' % (cq, NN, sec, ' '.join('%.1e' % e for e in ev[:10])), flush=True)
        print('   lambda_min on complement of V_K:', ' '.join(line), flush=True)
    print('   time', round(time.time() - t0, 1), flush=True)
