# run7.py (explore17): deflation test at c = 20 (N = 60) with a high-precision orthonormalised Phi-ladder basis.
import sys, os, time, math, numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'WeilSimple'))
from flint import arb, arb_mat, ctx, fmpq
from weil_matrix import WeilMat
from run5_lib import phi_derivs
t0 = time.time()
cq = int(sys.argv[1]) if len(sys.argv) > 1 else 20; NN = int(sys.argv[2]) if len(sys.argv) > 2 else 60; Jmax = int(sys.argv[3]) if len(sys.argv) > 3 else 47
mp.mp.dps = 60; L = mp.log(cq); Lf = float(L)
nq = int(sys.argv[4]) if len(sys.argv) > 4 else 600
xg, wg = np.polynomial.legendre.leggauss(nq)
X = [mp.mpf(float(v))*L/2 for v in xg]; Wq = [mp.mpf(float(v))*L/2 for v in wg]
D = [phi_derivs(xx, Jmax) for xx in X]
print('derivs done', round(time.time() - t0, 1), flush=True)
ctx.prec = 400; Wm = WeilMat(fmpq(cq))
for sec in ('even', 'odd'):
    T, idx = Wm.block(NN, sec); n = len(idx)
    Q = np.array([[float(T[a, b].mid()) for b in range(n)] for a in range(n)])
    js = list(range(0 if sec == 'even' else 1, Jmax + 1, 2))
    # basis values (mpmath) and projections
    def bval(a, x):
        j = idx[a]; w = 2*mp.pi*j/L
        if sec == 'even': return (1/mp.sqrt(L) if j == 0 else mp.sqrt(2/L))*mp.cos(w*x)
        return mp.sqrt(2/L)*mp.sin(w*x)
    Bv = [[bval(a, x) for x in X] for a in range(n)]
    C = mp.matrix(n, len(js))
    for c_, j in enumerate(js):
        for a in range(n): C[a, c_] = mp.fsum(Wq[q]*Bv[a][q]*D[q][j] for q in range(nq))
    # Gram-Schmidt in mpmath, column by column (modified GS twice)
    Qs = []
    for c_ in range(len(js)):
        v = C[:, c_]; v = v/mp.norm(v)
        for _ in range(2):
            for u in Qs: v = v - (u.T*v)[0]*u
        nv = mp.norm(v)
        Qs.append(v/nv)
    U = np.array([[float(Qs[k][a]) for k in range(len(Qs))] for a in range(n)])
    ev = np.linalg.eigvalsh(Q)
    line = []
    for K in range(0, len(js) + 1, 2):
        if K == 0: lm = ev[0]
        else:
            Pfull, _ = np.linalg.qr(np.hstack([U[:, :K], np.eye(n)]))
            Pc = Pfull[:, K:n]
            lm = np.linalg.eigvalsh(Pc.T @ Q @ Pc)[0]
        line.append('%d:%.2e' % (K, lm))
    print('c=%d N=%d %s: lambda_min on complement of V_K: %s' % (cq, NN, sec, ' '.join(line)), flush=True)
    print('  time', round(time.time() - t0, 1), flush=True)
