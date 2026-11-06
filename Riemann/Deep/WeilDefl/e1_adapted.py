# e1.py (weildefl): feasibility of boundary-adapted ladders. NOT rigorous (float64 Galerkin, mpmath/arb jets).
# usage: e1.py c sector kp p N
import sys, os, time, math, numpy as np
from flint import arb, arb_poly, arb_mat, ctx, fmpq
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/WeilSimple')
from weil_matrix import WeilMat
t0 = time.time()
cq = int(sys.argv[1]); sec = sys.argv[2]; kp = int(sys.argv[3]); p = int(sys.argv[4]); NN = int(sys.argv[5])
ctx.prec = int(sys.argv[6]) if len(sys.argv) > 6 else 800
NT = 6
par = 0 if sec == 'even' else 1
M = 2*(kp + p) + par + p + 4          # max derivative order needed
pi = arb.pi(); L = arb(cq).log(); Lf = float(L.mid())
# P[n][m]: Phi_n^{(m)}(u) = sum over (cf,a): e^{-bz} z^a * poly_m(z), z = e^{2u}; d/du = 2z d/dz on z^{a+i} e^{-bz}
def polys(n):
    b = pi*n*n; out = []
    for (cf, a) in ((2*pi**2*n**4, arb(9)/4), (-3*pi*n*n, arb(5)/4)):
        cur = [cf]; lst = []
        for m in range(M + 1):
            lst.append(arb_poly(cur))
            new = [arb(0)]*(len(cur) + 1)
            for i, ci in enumerate(cur):
                new[i] += 2*(a + i)*ci; new[i + 1] += -2*b*ci
            cur = new
        out.append((b, a, lst))
    return out
PP = [polys(n) for n in range(1, NT + 1)]
def derivs(u, ms):   # Phi^{(m)}(u) for m in ms, u >= 0 arb
    z = (2*u).exp(); res = {m: arb(0) for m in ms}
    for terms in PP:
        for (b, a, lst) in terms:
            f = (-b*z).exp()*z**a
            for m in ms: res[m] += f*lst[m](z)
    return res
print('polys', round(time.time() - t0, 1), flush=True)
# boundary jets at L/2: J[m][j] = Phi^{(2j+par+m)}(L/2)
nj = kp + p
d = derivs(L/2, list(range(0, 2*(nj - 1) + par + p + 3)))
Jm = [[d[2*j + par + m] for j in range(nj)] for m in range(p + 2)]
# null space of the p x nj jet matrix (high precision): basis via exact elimination on columns
import mpmath as mp
mp.mp.dps = int(ctx.prec*0.30)
A = mp.matrix(p, nj)
for m in range(p):
    for j in range(nj): A[m, j] = mp.mpf(Jm[m][j].mid().str(int(ctx.prec*0.3), radius=False))
# scale columns by 1/|jet0| to balance
def nullspace(A):
    rows, cols = A.rows, A.cols
    if rows == 0: return [ [mp.mpf(1) if i == k else mp.mpf(0) for i in range(cols)] for k in range(cols)]
    # Gaussian elimination with partial pivoting -> RREF
    R = A.copy(); piv = []; r = 0
    for cidx in range(cols):
        if r >= rows: break
        best = max(range(r, rows), key=lambda i: abs(R[i, cidx]))
        if abs(R[best, cidx]) == 0: continue
        R.swap_row(r, best) if hasattr(R, 'swap_row') else None
        if best != r:
            for jj in range(cols): R[r, jj], R[best, jj] = R[best, jj], R[r, jj]
        pv = R[r, cidx]
        for jj in range(cols): R[r, jj] /= pv
        for i in range(rows):
            if i != r and R[i, cidx] != 0:
                f = R[i, cidx]
                for jj in range(cols): R[i, jj] -= f*R[r, jj]
        piv.append(cidx); r += 1
    free = [cidx for cidx in range(cols) if cidx not in piv]
    basis = []
    for fcol in free:
        v = [mp.mpf(0)]*cols; v[fcol] = mp.mpf(1)
        for i, pc in enumerate(piv): v[pc] = -R[i, fcol]
        basis.append(v)
    return basis
NS = nullspace(A)
print('nullspace dim', len(NS), round(time.time() - t0, 1), flush=True)
# quadrature on [0, L/2] (even/odd functions): Gauss-Legendre
nq = int(sys.argv[7]) if len(sys.argv) > 7 else 400
xg, wg = np.polynomial.legendre.leggauss(nq)
X = [(arb(float(x)) + 1)*L/4 for x in xg]; Wq = [arb(float(w))*L/4 for w in wg]
ms = [2*j + par for j in range(nj)]
DV = [derivs(x, ms) for x in X]
print('node derivs', round(time.time() - t0, 1), flush=True)
def gvals(coef):   # g(x_q), q nodes on [0, L/2]
    return [sum((arb(str(coef[j]))*DV[q][ms[j]] for j in range(nj)), arb(0)) for q in range(nq)]
idx = list(range(0, NN + 1)) if sec == 'even' else list(range(1, NN + 1))
def bvec(x):
    out = []
    for jj in idx:
        w = 2*pi*jj/L
        if sec == 'even': out.append(((1/L).sqrt() if jj == 0 else (2/L).sqrt())*(w*x).cos())
        else: out.append((2/L).sqrt()*(w*x).sin())
    return out
BV = [bvec(x) for x in X]
rows = []; info = []; GV = []; JP = []; J0 = []
for v in NS:
    gv = gvals(v)
    nrm2 = 2*sum((Wq[q]*gv[q]**2 for q in range(nq)), arb(0))
    co = [2*sum((Wq[q]*gv[q]*BV[q][a] for q in range(nq)), arb(0)) for a in range(len(idx))]
    s = 1/nrm2.sqrt()
    # first nonvanishing boundary derivative: order p (of g), g^{(m)}(L/2) = sum_j v_j Phi^{(2j+par+m)}(L/2)
    gp = sum((arb(str(v[j]))*Jm[p][j] for j in range(nj)), arb(0))*s
    g0 = sum((arb(str(v[j]))*Jm[0][j] for j in range(nj)), arb(0))*s
    rows.append([float((cc*s).mid()) for cc in co]); info.append((float(gp.mid()), float(g0.mid()), float((1 - sum(((cc*s)**2 for cc in co), arb(0))).mid())))
    GV.append([x*s for x in gv]); JP.append(gp); J0.append(g0)
# span measures: Gram of normalized g's, max over unit vectors of |g^{(p)}(L/2)| and |g(L/2)|
kk = len(GV); Gm = arb_mat(kk, kk)
for i in range(kk):
    for j in range(i, kk):
        v_ = 2*sum((Wq[q]*GV[i][q]*GV[j][q] for q in range(nq)), arb(0)); Gm[i, j] = v_; Gm[j, i] = v_
Gi = Gm.inv()
bp = arb_mat([[x] for x in JP]); b0 = arb_mat([[x] for x in J0])
spanp = float(((bp.transpose()*Gi*bp)[0, 0]).mid())**0.5; span0 = float(((b0.transpose()*Gi*b0)[0, 0]).mid())**0.5
Gf = np.array([[float(Gm[i, j].mid()) for j in range(kk)] for i in range(kk)])
print('span: max |g^(p)(L/2)| over unit g in V = %.3e, max |g(L/2)| = %.3e, cond(Gram) = %.2e' % (spanp, span0, np.linalg.cond(Gf)))
print('projections', round(time.time() - t0, 1), flush=True)
C = np.array(rows).T     # (basis, kp)
W = WeilMat(fmpq(cq)); ctxp = ctx.prec; ctx.prec = 128
T, _ = W.block(NN, sec); ctx.prec = ctxp
n = len(idx); Q = np.array([[float(T[a, b].mid()) for b in range(n)] for a in range(n)])
ev = np.linalg.eigvalsh(Q)
U, sv, _ = np.linalg.svd(C, full_matrices=True)
print('sv of projected ladder', ' '.join('%.1e' % x for x in sv[:3]), '...', ' '.join('%.1e' % x for x in sv[-3:]))
line = []
for K in range(0, C.shape[1] + 1, max(1, C.shape[1]//8)):
    if K == 0: lm = ev[0]
    else:
        Ub, sb, _ = np.linalg.svd(C[:, :K], full_matrices=True); r = int(np.sum(sb > 1e-12*sb[0])); Pc = Ub[:, r:]
        lm = np.linalg.eigvalsh(Pc.T @ Q @ Pc)[0]
    line.append('%d:%.2e' % (K, lm))
Ub, sb, _ = np.linalg.svd(C, full_matrices=True); r = int(np.sum(sb > 1e-12*sb[0])); Pc = Ub[:, r:]
lmall = np.linalg.eigvalsh(Pc.T @ Q @ Pc)[0]
print('c=%d %s kp=%d p=%d N=%d: lambda_min complement: %s  all(r=%d): %.3e' % (cq, sec, kp, p, NN, ' '.join(line), r, lmall))
# tail estimates per vector (|g^{(p)}(L/2)|, residual norm outside E_N), R_est at K2 = 427 and 2000
for K2 in (427, 2000):
    rK = math.pi*K2/Lf
    Rest = [(30/math.pi)*spanp**2*rK**(-(2*p + 1))/(2*p + 1)]
    print('K2=%d  r=%.0f  span R_est=%.2e  (|g^(p)(L/2)| max %.2e, |g(L/2)| max %.2e, max ||g - P_N g||^2 %.1e)' % (K2, rK, max(Rest), max(abs(i[0]) for i in info), max(abs(i[1]) for i in info), max(i[2] for i in info)))
print('time', round(time.time() - t0, 1))
