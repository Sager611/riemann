# run4.py (explore17): is the near-null ladder of the window Weil form spanned by the Riemann-Polya kernel Phi and its derivatives?
# Phi(u) = sum_n (2 pi^2 n^4 e^{9u/2} - 3 pi n^2 e^{5u/2}) exp(-pi n^2 e^{2u}), FT = Xi(t); Phi^{(j)} has FT (-it)^j Xi(t), zero at every zeta zero.
import os, sys, time, math, numpy as np, mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'WeilSimple'))
from flint import arb, arb_mat, ctx, fmpq
from weil_matrix import WeilMat
t0 = time.time()
mp.mp.dps = 40
def Phi(u, d=0):
    u = mp.mpf(u)
    f = lambda v: mp.fsum((2*mp.pi**2*n**4*mp.e**(4.5*v) - 3*mp.pi*n**2*mp.e**(2.5*v))*mp.e**(-mp.pi*n*n*mp.e**(2*v)) for n in range(1, 12))
    return f(u) if d == 0 else mp.diff(f, u, d)
def basis_vals(L, sector, N, x):
    idx = np.arange(0, N + 1) if sector == 'even' else np.arange(1, N + 1)
    w = 2*math.pi*idx/L
    if sector == 'even':
        Nj = np.where(idx == 0, 1/math.sqrt(L), math.sqrt(2/L)); return np.array([Nj[a]*np.cos(w[a]*x) for a in range(len(idx))])
    return np.array([math.sqrt(2/L)*np.sin(w[a]*x) for a in range(len(idx))])
def proj(L, sector, N, d, nq=400):   # coefficients of Phi^{(d)} 1_W in the orthonormal basis, and the captured norm fraction
    xg, wg = np.polynomial.legendre.leggauss(nq); x = xg*L/2; wq = wg*L/2
    fv = np.array([float(Phi(xx, d)) for xx in x]); Bv = basis_vals(L, sector, N, x)
    co = Bv @ (wq*fv); nrm2 = float(np.sum(wq*fv*fv))
    return co, math.sqrt(float(co @ co)/nrm2)
def principal_cos(A, B):   # cosines of principal angles between column spans
    qa, _ = np.linalg.qr(A); qb, _ = np.linalg.qr(B); return np.linalg.svd(qa.T @ qb, compute_uv=False)
for cq, NN, dps in ((5, 40, 60), (10, 40, 90)):
    L = math.log(cq); ctx.prec = int(dps*3.4) + 64; mpd = mp.mp.dps; mp.mp.dps = dps
    Wm = WeilMat(fmpq(cq)); res = {}
    for sec in ('even', 'odd'):
        T, idx = Wm.block(NN, sec); n = len(idx); M = mp.matrix(n, n)
        for a in range(n):
            for b in range(n): M[a, b] = mp.mpf(T[a, b].mid().str(dps + 10, radius=False))
        E, Vm = mp.eigsy(M); order = sorted(range(n), key=lambda k: E[k])
        ev = [float(E[k]) for k in order]; vecs = np.array([[float(Vm[a, k]) for k in order] for a in range(n)])
        res[sec] = (ev, vecs)
    mp.mp.dps = 40
    print('c=%d N=%d even eigenvalues:' % (cq, NN), ' '.join('%.1e' % e for e in res['even'][0][:12]), flush=True)
    print('c=%d N=%d odd  eigenvalues:' % (cq, NN), ' '.join('%.1e' % e for e in res['odd'][0][:12]), flush=True)
    for sec, ds in (('even', (0, 2, 4, 6)), ('odd', (1, 3, 5, 7))):
        cols = []; caps = []
        for d in ds:
            co, cap = proj(L, sec, NN, d); cols.append(co); caps.append(cap)
        A = np.array(cols).T; ev, vecs = res[sec]
        for k in range(1, len(ds) + 1):
            pc = principal_cos(A[:, :k], vecs[:, :k])
            print('  %s k=%d: cos(principal angles) span{Phi^(d), d in %s} vs bottom-%d eigvecs: %s  (captured norm of Phi^(d) in E_N: %s)' % (
                sec, k, ds[:k], k, ' '.join('%.6f' % v for v in pc), ' '.join('%.4f' % v for v in caps[:k])), flush=True)
        # single overlap of ground state with Phi
        g = vecs[:, 0]; ph = A[:, 0]/np.linalg.norm(A[:, 0])
        print('  %s |<psi_0, Phi^(%d)>| = %.8f' % (sec, ds[0], abs(g @ ph)), flush=True)
    print('  time', round(time.time() - t0, 1), flush=True)
