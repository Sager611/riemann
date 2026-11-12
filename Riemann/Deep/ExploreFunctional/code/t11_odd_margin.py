# t11_odd_margin.py (explore39): odd sector beyond c = 2. Ground-state representation with the odd ground state psi > 0 on (0,a]:
# q(f) - lam0 |f|^2 = E+_psi(g) - E-_psi(g), E- = reflected prime atoms (negative jump rates w_n on x + y = ln n).
# Relative margin eta = inf_f [Q(f) - lam0 |f|^2]/Den(f), Den(f) = sum_n w_n int psi(x) psi(h-x) (f(x)/psi(x) - f(h-x)/psi(h-x))^2 dx.
import sys, numpy as np, math, json
sys.path.insert(0, '.scratch/explore39'); from xform import Form
for c in [2.05, 2.1, 2.2, 2.5, 3.0, 3.5, 4.0]:
    F = Form(c); n = 20; a = F.a
    Qm, p = F.Q(n, 'odd'); e, V = np.linalg.eigh(Qm); v = V[:, 0]
    psi = lambda x: F.B(np.atleast_1d(x), n, 'odd') @ v
    xs = np.linspace(1e-3*a, a, 400)
    if np.median(psi(xs)) < 0: v = -v
    pos = bool(np.all(psi(xs) > 0))
    Den = np.zeros((n, n))
    for h, w in F.atoms:
        lo, hi = max(0.0, h - a), min(a, h)
        if hi <= lo: continue
        x, wq = F.gauss(lo, hi, 200); y = h - x
        Px = F.B(x, n, 'odd'); Py = F.B(y, n, 'odd'); r = psi(y)/psi(x)
        Den += w*((Px*(wq*r)[:, None]).T @ Px - (Px*wq[:, None]).T @ Py - (Py*wq[:, None]).T @ Px + (Py*(wq/r)[:, None]).T @ Py)
    Den = (Den + Den.T)/2
    U = np.linalg.svd(np.eye(n) - np.outer(v, v))[0][:, :n-1]
    A = U.T @ (Qm - e[0]*np.eye(n)) @ U; A = (A + A.T)/2; Dd = U.T @ Den @ U
    L = np.linalg.cholesky(A); Li = np.linalg.inv(L); M = Li @ Dd @ Li.T
    mmax = float(np.linalg.eigvalsh((M + M.T)/2)[-1])
    print(json.dumps({'c': c, 'lam0_odd': float('%.6g' % e[0]), 'lam1-lam0': float('%.6g' % (e[1] - e[0])), 'psi>0': pos,
                      'eta=inf (Q-lam0)/Den': float('%.4g' % (1/mmax))}), flush=True)
