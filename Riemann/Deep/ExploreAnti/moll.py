# explore55: 1D (complex, CGG) vs 2D (real) second-moment small-ball bounds at zeta zeros. Single-threaded numpy.
import numpy as np, sys, json
from rs import zprime, theta
def mobius(n):
    mu = np.ones(n+1, int); mu[0] = 0; isp = np.ones(n+1, bool); isp[:2] = False
    for p in range(2, n+1):
        if isp[p]:
            isp[2*p::p] = False; mu[p::p] *= -1; mu[p*p::p*p] = 0
    return mu
def feats(g, zp, y, J, chunk=1000):
    """returns D = zeta'(rho) (complex), Mp[k,j] = M_j(rho_k), Mm[k,j] = M_j(1-rho_k) = conj(M_j(rho_k))"""
    mu = mobius(int(y)); n = np.nonzero(mu[1:int(y)+1])[0] + 1
    x = np.log(y/n)/np.log(y); w = mu[n]*n**-0.5
    Wt = np.stack([w*x**j for j in range(1, J+1)], 1)          # (n, J)
    ln = np.log(n.astype(float)); Mp = np.empty((len(g), J), complex)
    for a in range(0, len(g), chunk):
        E = np.exp(-1j*np.outer(g[a:a+chunk], ln)); Mp[a:a+chunk] = E @ Wt
    D = -1j*np.exp(-1j*theta(g))*zp
    return D, Mp, np.conj(Mp)
def bound_complex(Y):            # max_c |sum c.Y|^2/(N sum |c.Y|^2), Y: (N, J) complex
    s = Y.sum(0); G = Y.T @ np.conj(Y)                       # G_ij = sum Y_i conj(Y_j)
    return float(np.real(np.conj(s) @ np.linalg.solve(G, s)))/len(Y)   # = s^H Gbar^{-1}... see NOTES
def bound_real(Y):               # max over real-linear functionals of (Re Y, Im Y)
    W = np.concatenate([Y.real, Y.imag], 1); m = W.sum(0); G = W.T @ W
    return float(m @ np.linalg.solve(G, m))/len(W)
def run(g, zp, th, J):
    T = g[-1]; L = np.log(T/(2*np.pi)); y = T**th
    D, Mp, Mm = feats(g, zp, y, J); N = len(g)
    YA = D[:, None]*Mp; YC = D[:, None]*Mm
    r = dict(N=N, T=round(float(T), 3), L=round(float(L), 4), theta=th, y=round(float(y), 1), J=J)
    r['CGG_1D'] = bound_complex(YA); r['2D_real'] = bound_real(YA)
    r['A+C_1D'] = bound_complex(np.concatenate([YA, YC], 1)); r['A+C_2D'] = bound_real(np.concatenate([YA, YC], 1))
    y1 = YA[:, 0]; mu1 = y1.mean(); ph = mu1/abs(mu1); z = y1/ph          # rotate so mean is real positive
    r['P=x: mu/L'] = float(abs(mu1)/L); r['A/L^2'] = float(np.mean(abs(z)**2)/L**2)
    r['B/L^2 (rot)'] = complex(np.mean(z**2)/L**2)
    r['var_par/var'] = float(np.var(z.real)/(np.var(z.real)+np.var(z.imag)))
    U = 1/D; Jm1 = float(np.sum(abs(U)**2))                                   # Milinovich-Ng control
    def proj(V):
        b = np.conj(V).T @ U; G = np.conj(V).T @ V; return float(np.real(np.conj(b) @ np.linalg.solve(G, b)))
    r['Jm1_emp'] = Jm1; r['MN_1D_frac'] = proj(Mm)/Jm1; r['MN_+Mp_frac'] = proj(np.concatenate([Mm, Mp], 1))/Jm1
    return r
if __name__ == '__main__':
    g_all = np.loadtxt(sys.argv[1]); out = []
    for (a, b) in [(0, 100), (1000, 2000), (10000, 20000), (50000, 100000)]:
        g = g_all[a:b]; zp = zprime(g)
        for th in [0.25, 0.5]:
            for J in [1, 3]:
                r = run(g, zp, th, J); r['window'] = [a+1, b]; out.append(r)
                print(json.dumps({k: (round(v, 5) if isinstance(v, float) else (str(np.round(v, 4)) if isinstance(v, complex) else v)) for k, v in r.items()}), flush=True)
