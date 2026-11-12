# t8_poincare.py (explore39): even-sector Weil positivity as a Poincare inequality (ground-state representation).
# Q^ev = M + 2 C C^T, M = Q - pole Markov (jump kernel kappa(|x-y|) + prime atoms >= 0, every c), C = cosh(x/2).
# For psi > 0, f = psi g:  Q(f) - int (Q psi/psi) f^2 = E_psi(g) - R Var_nu(g),  E_psi(g) = 1/2 iint j_M psi psi (g(x)-g(y))^2,
#   R = 2 <C,psi>^2, nu = C psi dx/<C,psi>.  With psi = ground state: Poincare constant P(E_psi, nu) = R + mu_W, where
#   mu_W = min over f not ~ psi of [Q(f) - lam0 |f|^2]/W(f),  W(f) = (1/<C,psi>)[int (C/psi) f^2 - <C,f>^2/<C,psi>]  (= Var_nu(g)).
import sys, numpy as np, math, json
from scipy.linalg import eigh
sys.path.insert(0, '.scratch/explore39'); from xform import Form
rho = 1.324717957244746
for c in [1.2, rho, 1.5, rho**2, 2.0, 2.5, 3.0, 3.5, 4.0]:
    F = Form(c); n = 20; a = F.a
    Qm, p = F.Q(n, 'even'); e, V = np.linalg.eigh(Qm); v = V[:, 0]
    xs, ws = F.gauss(-a, a, 400); Bx = F.B(xs, n, 'even'); psi = Bx @ v
    if np.median(psi) < 0: v = -v; psi = -psi
    Cv = p['em']; Cpsi = float(Cv @ v); R = 2*Cpsi**2
    Wm = ((Bx*(ws*np.cosh(xs/2)/psi)[:, None]).T @ Bx - np.outer(Cv, Cv)/Cpsi)/Cpsi
    # restrict to psi-perp
    Pm = np.eye(n) - np.outer(v, v); U = np.linalg.svd(Pm)[0][:, :n-1]
    A = U.T @ (Qm - e[0]*np.eye(n)) @ U; Bw = U.T @ Wm @ U; Bw = (Bw + Bw.T)/2; A = (A + A.T)/2
    mw = eigh(A, Bw, eigvals_only=True)
    phi1 = V[:, 1]; W1 = float(phi1 @ Wm @ phi1)
    print(json.dumps({'c': round(c, 6), 'lam0': float('%.6g' % e[0]), 'lam1-lam0': float('%.6g' % (e[1]-e[0])), 'psi>0': bool(np.all(psi > 0)),
                      'R=2<C,psi>^2': round(R, 6), 'mu_W=P-R': float('%.5g' % mw[0]), '(P-R)/R': float('%.4g' % (mw[0]/R)),
                      'W(phi1)': float('%.5g' % W1), 'bound (lam1-lam0)/W(phi1)': float('%.5g' % ((e[1]-e[0])/W1))}), flush=True)
