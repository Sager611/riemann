# ceilaudit8 B3 (explore39): does another ground-state transform psi escape the Poincare margin?  Identity C (explore39):
# Q(f) - int V f^2 = E_psi(g) - R Var_nu(g), V = Q psi/psi, g = f/psi, R = 2<C,psi>^2.  Certificate: V >= 0 pointwise and P >= R.
# Lemma (here): if V >= 0 then for every f not prop. to psi, P - R <= Q(f)/W_psi(f), W_psi(f) = Var_nu(f/psi); so
# P - R <= min(lam0/W_psi(psi0), lam1/W_psi(phi1)): the margin is capped by Q's own eigenvalues for every psi.
import sys, numpy as np, math, json
from scipy.linalg import eigh
sys.path.insert(0, 'Riemann/Deep/ExploreFunctional/code'); from xform import Form
def run(c, n=20):
    F = Form(c); a = F.a; Qm, p = F.Q(n, 'even'); Qm = (Qm+Qm.T)/2; e, V = np.linalg.eigh(Qm); Cv = p['em']
    Mm = Qm - 2*np.outer(Cv, Cv); mu, U = np.linalg.eigh(Mm)
    xs, ws = F.gauss(-a, a, 400); Bx = F.B(xs, n, 'even'); ch = np.cosh(xs/2)
    def W(v, f):  # Var_nu(f/psi) as quadratic form, nu = C psi dx/<C,psi>
        psi = Bx @ v; Cpsi = Cv @ v
        Wm = ((Bx*(ws*ch/psi)[:, None]).T @ Bx - np.outer(Cv, Cv)/Cpsi)/Cpsi
        return Wm, float(f @ Wm @ f)
    def margin(v, name):
        v = v/np.linalg.norm(v); psi = Bx @ v
        if np.median(psi) < 0: v = -v; psi = -psi
        if not np.all(psi > 0): return {'psi': name, 'psi>0': False}
        Vx = (Bx @ (Qm @ v))/psi; Vm = (Bx*(ws*Vx)[:, None]).T @ Bx; R = 2*float(Cv @ v)**2
        Wm, _ = W(v, v); Pm = np.eye(n) - np.outer(v, v); Uq = np.linalg.svd(Pm)[0][:, :n-1]
        A = Uq.T @ (Qm - Vm) @ Uq; Bw = Uq.T @ Wm @ Uq; A = (A+A.T)/2; Bw = (Bw+Bw.T)/2
        PR = float(eigh(A, Bw, eigvals_only=True)[0])
        b0 = W(v, V[:, 0])[1]; b1 = W(v, V[:, 1])[1]
        cap = min(e[0]/b0 if b0 > 1e-14 else np.inf, e[1]/b1)
        return {'psi': name, 'infV': float('%.4g' % Vx.min()), 'supV': float('%.4g' % Vx.max()), 'Barta ok (0<=infV<=lam0)': bool(0 <= Vx.min() <= e[0]*(1+1e-6)),
                '(P-R)/R': float('%.4g' % (PR/R)), 'cap/R': float('%.4g' % (cap/R)), 'PR<=cap': bool(PR <= cap*(1+1e-6)+1e-15)}
    out = {'c': c, 'lam0': float('%.4g' % e[0]), 'lam1': float('%.4g' % e[1]), 'mu0(M)': float('%.5g' % mu[0]), 'mu1(M)': float('%.4g' % mu[1]),
           'M-Perron margin mu1/|mu0|': float('%.4g' % (mu[1]/abs(mu[0])))}
    print(json.dumps(out), flush=True)
    cand = [(V[:, 0], 'psi0(Q)'), (U[:, 0], 'Perron(M)'), (Cv, 'pole C=cosh(x/2)')]
    for s in (1e-6, 1e-4, 1e-2): cand.append((V[:, 0] + s*V[:, 1], f'psi0+{s:g}phi1'))
    for s in (1e-3, 1e-1): cand.append((V[:, 0] + s*V[:, 2], f'psi0+{s:g}phi2'))
    for v, nm in cand: print('   ', json.dumps(margin(v, nm)), flush=True)
for c in (2.5, 3.0, 3.5): run(c)
