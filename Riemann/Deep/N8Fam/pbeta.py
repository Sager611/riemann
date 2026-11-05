# p(beta;F_1) of Johns (Zenodo 23070760, eq. (pbeta)): 2 - inf Q_F(f), f>=0 even, supp f in [-beta/2,beta/2], int f=1,
# Q_F(f) = int f^2 + iint f(x)f(y)F(|x-y|), F = F_1 = min(alpha,1).  Step-function QP (M cells on [0,beta/2]).
# Also: the paper's (2610.08965) second window dilated to support lam <= beta, with the scaled certificate (x, lam*b).
# Run: OMP_NUM_THREADS=1 uv run --no-project --with numpy --with scipy python3 pbeta.py
import numpy as np
from scipy.optimize import minimize
import sys
sys.path.insert(0, '.')
from cfam import C_model, cos_window, A2, X_CERT, B_CERT

def G2(x, C=1.0):
    # G2(x) = int_0^|x| (|x|-z) F(z) dz with F(z) = min(z, 1) for z<=1, C for z>1  (F_C of Johns)
    x = np.abs(x)
    out = np.where(x <= 1, x ** 3 / 6, 0.0)
    xm = x[x > 1]
    out[x > 1] = (xm / 2 - 1 / 3) + C * (xm - 1) ** 2 / 2
    return out

def quad_matrix(beta, M, C=1.0):
    # cells on [0, beta/2]; f even step function with value f_i on +/- cell i; basis b_i = 1 on both cells
    e = np.linspace(0, beta / 2, M + 1)
    h = e[1] - e[0]
    a, b = e[:-1], e[1:]
    def box(a1, b1, a2, b2):  # iint_{[a1,b1]x[a2,b2]} F(|s-t|) ds dt
        return G2(b1[:, None] - a2[None, :], C) - G2(a1[:, None] - a2[None, :], C) - G2(b1[:, None] - b2[None, :], C) + G2(a1[:, None] - b2[None, :], C)
    # each basis function = cell on positive side + mirrored cell on negative side
    P = box(a, b, a, b) + box(a, b, -b, -a) + box(-b, -a, a, b) + box(-b, -a, -b, -a)
    D = 2 * h * np.eye(M)          # int b_i b_j
    m = 2 * h * np.ones(M)        # int b_i
    return D + P, m, e

def pbeta(beta, M=400, C=1.0):
    A, m, e = quad_matrix(beta, M, C)
    # minimize f^T A f s.t. m.f = 1, f >= 0 ; start from unconstrained solution
    try:
        z = np.linalg.solve(A, m); f0 = z / (m @ z)
    except np.linalg.LinAlgError:
        f0 = np.ones(len(m)) / m.sum()
    if f0.min() >= 0:
        f = f0
    else:
        res = minimize(lambda f: f @ A @ f, np.maximum(f0, 0) / (m @ np.maximum(f0, 0)), jac=lambda f: 2 * A @ f,
                       constraints=[{'type': 'eq', 'fun': lambda f: m @ f - 1, 'jac': lambda f: m}],
                       bounds=[(0, None)] * len(m), method='SLSQP', options={'maxiter': 500, 'ftol': 1e-15})
        f = res.x
    return 2 - f @ A @ f, f, e

def cert_paper(beta):
    # best over lam <= beta of the dilated paper window with scaled certificate; also plain
    f2 = cos_window(A2)
    lams = np.linspace(max(0.9, beta - 0.3), beta, 61)
    vals = [(2 - C_model(f2, l, 'gue') - l * B_CERT) / (1 - X_CERT) for l in lams]
    i = int(np.argmax(vals))
    return vals[i], lams[i], 2 - C_model(f2, beta, 'gue')

if __name__ == '__main__':
    johns = {0: 0.932282, 0.25: 0.912574, 0.5: 0.894985, 1: 0.865673, 2: 0.824355, 3: 0.797213, 5: 0.764149, 10: 0.727484}
    print('kappa  beta     p(beta;F1) ours  Johns     paper-window plain(lam=beta)  dilated paper cert (best lam<=beta)  gain vs p')
    for k in [0, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 50, 100, 1000]:
        beta = (2 + k) / (1 + k)
        p, f, e = pbeta(beta)
        c, lam, plain = cert_paper(beta)
        print('%-6g %.5f  %.6f        %s  %.6f                      %.6f (lam %.4f)                 %+.6f'
              % (k, beta, p, ('%.6f' % johns[k]) if k in johns else '   -    ', plain, c, lam, c - p), flush=True)
    print('beta=1 (kappa=inf): p = %.6f (MT 0.6725007), paper cert 0.6735382' % pbeta(1.0)[0])
