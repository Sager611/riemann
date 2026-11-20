# hmean: modulus-localization test of H.  At zeros, zeta'(rho) ~ Et_tau(rho), Et_K(s) = sum_{n<=K} (ell - log n) n^{-s}
# (ell = log(gamma/2pi), tau = gamma/2pi), so S_B ~ sum zeta'(rho) Et_tau(rho) B(rho)^2: a first moment whose twist has length
# tau*y^2 >= T.  F(kappa) = mean zeta'(rho) Et_K(rho) B(rho)^2 / ell^2 with K = tau^kappa keeps only the moduli n <= tau^kappa.
# Recipe (NOTES.md sec. 3): F(kappa) -> 1/6 + mu + mu^2 + mu (1-kappa)^2 + (1-kappa)^3/3, mu = (log y/ell) * int P  (P = x: int P = 1/2).
# Single-threaded numpy; sizes: chunk x Nmax complex arrays (<= 30 MB each).
import numpy as np, sys, json, re
TWOPI = 2*np.pi

def mobius(n):
    mu = np.ones(n+1, int); mu[0] = 0; isp = np.ones(n+1, bool); isp[:2] = False
    for p in range(2, n+1):
        if isp[p]:
            isp[2*p::p] = False; mu[p::p] *= -1; mu[p*p::p*p] = 0
    return mu

def pred(kap, mu):
    return 1/6 + mu + mu**2 + mu*(1-kap)**2 + (1-kap)**3/3

def mollcols(nmax, y, mu_arr):
    """column vector c_k = mu(k) P(log(y/k)/log y) k^{-1/2} on k <= y (P = x), zero elsewhere, length nmax"""
    c = np.zeros(nmax); k = np.arange(1, int(y)+1); x = np.log(y/k)/np.log(y)
    c[:len(k)] = mu_arr[k]*x*k**-0.5; return c


def split3(D, Et5, G5, chi, Bs, ell):
    """zeta'(rho) = Et_K(rho) + chi(rho) G_M(1-rho), K = M = sqrt(tau): S_B = P0 (chi^0) + 2 X1 (chi^1) + X2 (chi^2)"""
    out = {}
    for nm, (Bv, mu) in Bs.items():
        B2 = Bv**2; w2 = ell**2; c = lambda z: [round(z.real, 5), round(z.imag, 5)]
        P0 = (Et5**2*B2/w2).mean(); X1 = (chi*Et5*G5*B2/w2).mean(); X2 = (chi**2*G5**2*B2/w2).mean(); ex = (D*D*B2/w2).mean()
        out[nm] = dict(P0=c(P0), X1=c(X1), X2=c(X2), sum=c(P0 + 2*X1 + X2), exact=c(ex))
    return out

def summarize(D, Et, Bs, ell, kappas, info):
    out = dict(info); w2 = ell**2; res = []
    for name, (Bv, mu) in Bs.items():
        B2 = Bv**2; F = [(D*Et[:, j]*B2/w2).mean() for j in range(len(kappas))]
        exact = (D*D*B2/w2).mean(); u = (D*Bv/ell).mean()
        res.append(dict(moll=name, mu=round(mu, 5), u=[round(u.real, 5), round(u.imag, 5)],
                        exact_S_B=[round(exact.real, 5), round(exact.imag, 5)],
                        F=[[k, round(f.real, 5), round(f.imag, 5), round(pred(k, mu), 5)] for k, f in zip(kappas, F)]))
    out['res'] = res; print(json.dumps(out), flush=True)

def run_low(path, a, b, thetas, kappas, chunk=int(__import__('os').environ.get('HMEAN_CHUNK', '40'))):
    sys.path.insert(0, __file__.rsplit('/', 1)[0]); from rs import zprime, theta
    g = np.loadtxt(path)[a:b]; D = -1j*np.exp(-1j*theta(g))*zprime(g)
    tau = g/TWOPI; ell = np.log(tau); Nmax = int(np.floor(tau.max())); taumax = tau.max(); ellmax = np.log(taumax)
    n = np.arange(1, Nmax+1, dtype=float); ln = np.log(n); w = n**-0.5
    mu_arr = mobius(Nmax+1); Bs = {}
    for th in thetas:
        if th == 0: Bs['theta=0'] = (np.ones(len(g), complex), 0.0); continue
        y = taumax**th; c = mollcols(Nmax, y, mu_arr)[:int(y)]; Bv = np.empty(len(g), complex)
        for s in range(0, len(g), 500):
            Bv[s:s+500] = np.exp(-1j*np.outer(g[s:s+500], ln[:int(y)])) @ c
        Bs['theta=%g' % th] = (Bv, th*0.5*np.log(y)/(th*np.sqrt(np.mean(ell**2))) if th else 0.0)
    Kidx = np.floor(np.exp(np.outer(ell, kappas)) + 1e-9).astype(int)
    Et = np.empty((len(g), len(kappas)), complex); G5 = np.empty(len(g), complex); M5 = np.floor(np.sqrt(tau) + 1e-9).astype(int)
    for s in range(0, len(g), chunk):
        E = np.exp(-1j*np.outer(g[s:s+chunk], ln))*w
        C = np.cumsum(E*(ell[s:s+chunk, None] - ln[None, :]), axis=1)
        Et[s:s+chunk] = np.take_along_axis(C, np.clip(Kidx[s:s+chunk]-1, 0, Nmax-1), axis=1)
        Cg = np.cumsum(E[:, :M5.max()]*ln[None, :M5.max()], axis=1); G5[s:s+chunk] = np.conj(Cg[np.arange(len(Cg)), M5[s:s+chunk]-1])
    sp3 = split3(D, Et[:, list(kappas).index(0.5)], G5, np.exp(-2j*theta(g)), Bs, ell) if 0.5 in kappas else None
    summarize(D, Et, Bs, ell, kappas, dict(window=[a+1, b], T=float(g.max()), ell_rms=float(np.sqrt(np.mean(ell**2))), Nmax=Nmax, split3=sp3))

def run_high(path, t0, thetas, kappas, nz=None, chunk=8):
    import mpmath as mp; mp.mp.dps = 40
    d = np.loadtxt(path)[:nz]; t = t0 + d; tau = t/TWOPI; ell = np.log(tau)      # ell constant to 1e-8 here
    N = int(np.floor(np.sqrt(tau.min()))); assert N == int(np.floor(np.sqrt(tau.max())))
    n = np.arange(1, N+1, dtype=float); ln = np.log(n); w = n**-0.5
    r = np.array([float(mp.fmod(t0*mp.log(k), 2*mp.pi)) for k in range(1, N+1)])
    thm = np.array([float(mp.fmod(mp.siegeltheta(mp.mpf(t0)+mp.mpf(repr(float(x)))), 2*mp.pi)) for x in d])
    dth = 0.5*ell; e0 = float(ell.mean()); mu_arr = mobius(N+1)
    cols = [w*(dth.mean() - ln)]                       # Z'(t) main sum weights (dth constant to 1e-8)
    K = [int(np.floor(np.exp(k*e0) + 1e-9)) for k in kappas]; assert max(K) <= N
    for Kj in K:
        c = np.zeros(N); c[:Kj] = (e0 - ln[:Kj])*w[:Kj]; cols.append(c)
    names = []
    for th in thetas:
        if th == 0: continue
        y = float(np.exp(th*e0)); cols.append(mollcols(N, y, mu_arr)); names.append(('theta=%g' % th, 0.5*np.log(y)/e0))
    tk = [0.5, 0.6, 0.7, 0.8, 0.9]; Mj = [int(np.floor(np.exp((1-k)*e0) + 1e-9)) for k in tk]; ng = len(names)
    for M in Mj:
        c = np.zeros(N); c[:M] = ln[:M]*w[:M]; cols.append(c)          # G_M(1-rho) = conj(sum_{n<=M} log n n^{-1/2} e^{-i gamma log n})
    Tm = float(t.max())
    Cm = np.stack(cols, 1); M = np.empty((len(t), Cm.shape[1]), complex)
    for s in range(0, len(t), chunk):
        ph = r[None, :] + d[s:s+chunk, None]*ln[None, :]
        M[s:s+chunk] = np.exp(-1j*ph) @ Cm
    zp = -2*np.imag(np.exp(1j*thm)*M[:, 0]); D = -1j*np.exp(-1j*thm)*zp
    Et = M[:, 1:1+len(K)]; Bs = {'theta=0': (np.ones(len(t), complex), 0.0)} if 0 in thetas else {}
    for j, (nm, mu) in enumerate(names):
        Bs[nm] = (M[:, 1+len(K)+j], mu)
    chi = np.exp(-2j*thm); G = np.conj(M[:, 1+len(K)+ng:1+len(K)+ng+len(Mj)]); tails = {}
    for nm, (Bv, mu) in Bs.items():   # S_B - F(tau^k) ~ sum chi(rho) zeta'(rho) B^2 G_{tau^(1-k)}(1-rho)  (AFE of ell*zeta+zeta' at a zero)
        tails[nm] = [[k, round(v.real, 5), round(v.imag, 5), round(mu*(1-k)**2 + (1-k)**3/3, 5)] for k, v in zip(tk, (chi[:, None]*D[:, None]*(Bv**2)[:, None]*G/ell[:, None]**2).mean(0))]
    sp3 = split3(D, Et[:, list(kappas).index(0.5)], G[:, 0], chi, Bs, ell) if 0.5 in kappas else None
    summarize(D, Et, Bs, ell, kappas, dict(window='zeros3 (10^12+1..)', nz=len(t), T=Tm, ell=e0, N=N, K=K, M=Mj, medZp=float(np.median(abs(zp))), tails_chi2=tails, split3=sp3))

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'low':
        run_low(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), [float(x) for x in sys.argv[5].split(',')], [float(x) for x in sys.argv[6].split(',')])
    else:
        run_high(sys.argv[2], int(sys.argv[3]), [float(x) for x in sys.argv[4].split(',')], [float(x) for x in sys.argv[5].split(',')], int(sys.argv[6]) if len(sys.argv) > 6 else None)
