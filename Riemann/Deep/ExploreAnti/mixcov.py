# explore55 follow-up: antisymmetric mixed covariance Im sum_gamma Y_P conj(Y_Q) and whether Im-parts help the real 2D bound.
import numpy as np, sys, json, re
from moll import mobius, bound_real
from rs import zprime, theta as rstheta
TWOPI = 2*np.pi
def polys(x, th):
    return {'x': x, 'x2': x**2, 'x3': x**3, 'BHB': -th*x**2 + (1+th)*x, 'P2D': ((0.5+th)*x - th*x**2/2)/(0.5+th/2)}
def analyse(D, M, L, th, tag):
    Y = D[:, None]*M; ph = Y[:, 0].mean()/abs(Y[:, 0].mean()); Z = Y/ph; N = len(Z)
    H = (Z.T @ np.conj(Z))/N/L**2; B = (Z.T @ Z)/N/L**2; U = Z.mean(0)/L
    C = np.real(Z).T @ np.imag(Z)/N/L**2 - np.outer(U.real, U.imag)          # Cov(Re Y_P, Im Y_Q)
    k = len(U); off = [(i, j) for i in range(k) for j in range(k) if i < j]
    Zi = Z[:, :3]                                             # x, x^2, x^3 independent (BHB, P2D are combinations of x, x^2)
    q = lambda W: float(W.sum(0) @ np.linalg.lstsq(W.T @ W, W.sum(0), rcond=None)[0])/N
    re_only = q(np.real(Zi)); re_im = q(np.concatenate([np.real(Zi), np.imag(Zi)], 1))
    out = dict(tag=tag, L=round(float(L), 3), theta=th, N=N,
               maxImH_offdiag=round(float(max(abs(H[i, j].imag) for i, j in off)), 5),
               maxReH_offdiag=round(float(max(abs(H[i, j].real) for i, j in off)), 4),
               maxImB=round(float(np.max(abs(B.imag))), 5), maxImU=round(float(np.max(abs(U.imag))), 5),
               max_abs_Cov_ReIm=round(float(np.max(abs(C))), 5),
               bound_Re_only_J3=round(re_only, 6), bound_Re_and_Im_J3=round(re_im, 6))
    print(json.dumps(out), flush=True)
if __name__ == '__main__':
    src = sys.argv[1]
    if 'zeros1' in src:
        g = np.loadtxt(src)[50000:100000]; T = g.max(); L = np.log(T/TWOPI); D = -1j*np.exp(-1j*rstheta(g))*zprime(g)
        for th in [0.35, 0.45]:
            y = T**th; mu = mobius(int(y)); n = np.nonzero(mu[1:int(y)+1])[0] + 1; x = np.log(y/n)/np.log(y)
            Wt = np.stack([mu[n]*n**-0.5*v for v in polys(x, th).values()], 1); ln = np.log(n.astype(float))
            M = np.concatenate([np.exp(-1j*np.outer(g[a:a+500], ln)) @ Wt for a in range(0, len(g), 500)])
            analyse(D, M, L, th, 'zeros1[50000:100000]')
    else:
        import mpmath as mp; mp.mp.dps = 40; t0 = 267653395647
        d = np.array([float(l) for l in open(src) if re.fullmatch(r'\s*[0-9]+\.[0-9]+\s*', l)]); t = t0 + d; T = t.max(); L = np.log(T/TWOPI)
        zp, _ = np.load(src.replace('.txt', '_zp10000.npy'))
        thm = np.array([float(mp.fmod(mp.siegeltheta(mp.mpf(t0)+mp.mpf(str(v))), 2*mp.pi)) for v in d]); D = -1j*np.exp(-1j*thm)*zp
        for th in [0.35, 0.45]:
            y = T**th; mu = mobius(int(y)); n = np.nonzero(mu[1:int(y)+1])[0] + 1; x = np.log(y/n)/np.log(y)
            r = np.array([float(mp.fmod(t0*mp.log(int(k)), 2*mp.pi)) for k in n]); ln = np.log(n.astype(float))
            Wt = np.stack([mu[n]*n**-0.5*v for v in polys(x, th).values()], 1)
            M = np.concatenate([np.exp(-1j*(r[None, :] + d[a:a+8, None]*ln[None, :])) @ Wt for a in range(0, len(t), 8)])
            analyse(D, M, L, th, 'zeros3 (1e12)')
