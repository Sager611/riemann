# e1arb: optimise the Levinson-Conrey detector in E1's class (float). For fixed (R, P) c is a quadratic form in Q's free
# coefficients x (Q = (1-y) + sum x_j (T_{2j+1}(1-2y) - (1-2y))): one linear solve; for fixed (R, Q) the best P is
# sinh(lam x)/sinh(lam), lam = sqrt(b/a), value sqrt(ab) coth(lam). Alternate to convergence, then maximise kappa in R.
# usage: opt.py nu m1,m2,...   -> opt_<nu>.json
import numpy as np, json, sys
from numpy.polynomial import chebyshev as C
from scipy.optimize import minimize_scalar
xg, wg = np.polynomial.legendre.leggauss(400); xg = (xg+1)/2; wg = wg/2
def basis(m):
    u = 1-2*xg; F = [1-xg]; dF = [-np.ones_like(xg)]
    for j in range(1, m+1):
        c = np.zeros(2*j+2); c[2*j+1] = 1
        F.append(C.chebval(u, c)-u); dF.append(-2*(C.chebval(u, C.chebder(c))-1))
    return np.array(F), np.array(dF)
def grams(R, F, dF):
    e = wg*np.exp(2*R*xg); G = R*F+dF
    return (F*e) @ F.T, (G*e) @ G.T
def PB(l):  # C_P = int P'^2, B_P = int P^2 for P = sinh(l x)/sinh(l)
    s = np.sinh(l); s2 = np.sinh(2*l)/(2*l)
    return l*l*(s2+1)/(2*s*s), (s2-1)/(2*s*s)
def cstar(R, nu, F, dF):
    A, B = grams(R, F, dF); l = 1.0; cp = np.inf
    for it in range(500):
        CP, BP = PB(l); M = (CP/nu)*A+nu*BP*B
        x = np.linalg.solve(M[1:, 1:], -M[1:, 0]); z = np.concatenate([[1.0], x])
        a = (z@A@z)/nu; b = nu*(z@B@z); l = np.sqrt(b/a); c = 0.5+np.sqrt(a*b)/np.tanh(l)
        if cp-c < 1e-16*c: break
        cp = c
    return c, x, l, it
def best(nu, m, lo=1.0, hi=16.0):
    F, dF = basis(m); f = lambda R: -(1-np.log(cstar(R, nu, F, dF)[0])/R)
    Rs = np.linspace(lo, hi, 61); v = [f(R) for R in Rs]; i = int(np.argmin(v))
    r = minimize_scalar(f, bounds=(Rs[max(i-1, 0)], Rs[min(i+1, 60)]), method='bounded', options={'xatol': 1e-10})
    c, x, l, it = cstar(r.x, nu, F, dF)
    return {'nu': nu, 'm': m, 'kappa': -r.fun, 'R': r.x, 'lam': l, 'c': c, 'x': list(x), 'iters': it}
if __name__ == "__main__":
    nu = float(sys.argv[1]); out = []
    for m in [int(s) for s in sys.argv[2].split(',')]:
        d = best(nu, m); out.append(d)
        print(f"nu={nu} m={m} degQ={2*m+1} kappa={d['kappa']:.10f} R={d['R']:.6f} lam={d['lam']:.6f} c={d['c']:.6f} it={d['iters']} x[:4]={np.round(d['x'][:4],6).tolist()} x[-1]={d['x'][-1]:.2e}", flush=True)
    json.dump(out, open(f"opt_{sys.argv[1]}.json", 'w'))
